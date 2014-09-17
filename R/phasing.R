# phasing.R

whichmaxNA <- function(x) {
  # which.max that returns NA if all data is NA
  ifelse(any(is.na(x)), NA, which.max(x))
}

pruneWhichIndividuals <- function(geno, na_thresh) {
  # Return which individuals (columns) do not meet NA threshold
  na_rates  <- colSums(!is.na(geno))/nrow(geno)
  which(na_rates > na_thresh)
}

clusterMAP <- function(x) {
  apply(x, 2, which.max)
}

haploDiff <- function(x) {
  diff <- x[,1] == x[,2]
  sum(diff, na.rm=TRUE)/nrow(x)
}

#
#logSumExpPosterior <- function(liks, pi) {
#  # log sum exp trick for list of matrices of loci x individuals, each for a
#  # particular haplotype
#  bc <- t(mapply(function(x, p) colSums(log(x)) + log(p), liks, pi))
#
#  # denominator (same for both haplotype probs), uses TLOP
#  max_bc <- apply(bc, 2, max)
#  scaled_bc  <- sweep(bc, 2, max_bc, FUN='-') # subtract B
#  denom <- log(colSums(exp(scaled_bc))) + max_bc
#  lse <- sweep(bc, 2, denom, FUN='-')
#  browser()
#  return(lse)
#}

logSumExpPosterior <- function(liks, pi) {
  # log sum exp trick for list of matrices of loci x individuals, each for a
  # particular haplotype
  # note: cluster membership should be biased based on missingess since the same loci
  # are missingin both ll calss
  bc <- t(mapply(function(x, p) colSums(log(x), na.rm=TRUE) + log(p), liks, pi))

  # denominator (same for both haplotype probs), uses TLOP
  max_bc <- apply(bc, 2, max)
  scaled_bc  <- sweep(bc, 2, max_bc, FUN='-') # subtract B
  denom <- log(colSums(exp(scaled_bc))) + max_bc
  sweep(bc, 2, denom, FUN='-')
}

isConverged <- function(ll, ll_last, eps, debug=FALSE) {
  conv <- abs(na.exclude(ll_last - ll)) < eps
  if (debug) message(sprintf("%d of %d individuals converged", sum(conv), length(conv)))
  if (all(conv)) return(TRUE)
  return(FALSE)
}

reshapeIndLL <- function(lls_inds) {
  # reshape the individual likelihoods from each iteration so that they're a
  # dataframe
  ll <- do.call(rbind, mapply(function(x, hap) {
              d <- as.data.frame(do.call(rbind, x))
              colnames(d) <- paste('ind', seq_len(ncol(d)), sep="_")
              d$iter <- seq_len(nrow(d))
              d$hap <- hap
              d
  }, lls_inds, 1L:2L, SIMPLIFY=FALSE))
  setNames(melt(ll, id.vars=c('iter', 'hap')), c("iter", "hap", "ind", "ll"))
}

reshapeHapLL <- function(lls_haps) {
  # reshape the loci likelihoods, from each iteration
  do.call(rbind, mapply(function(x, i) {
          d <- as.data.frame(x)
          d$loci <- seq_len(nrow(x))
          d$iter <- i
          setNames(d, c("h1", "h2", "loci", "iter"))
  }, lls_haps, seq_along(lls_haps), SIMPLIFY=FALSE))
}


#' Use EM Haplotype Phasing by Imputation Algorithm
#'
#' @param pgeno progeny genotypes
#' @param fgeno father (or other parent) genotypes
#' @param freqs allele frequencies
#' @param other_parents indices to corresponding father genotype columns in fgeno
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
#' @param free_pi whether to allow pi to vary (can cause identifiability issues
#' @param na_thresh how much missingness to tolerate before pruning individual
#' @param max_iter maximum number of EM iterations before quiting
#' @param eps epsilon in log likelihood before calling convergence
emHaplotypeImpute <- function(pgeno, fgeno, other_parents, freqs, ehet, ehom,
                              free_pi=FALSE, na_thresh=0.8, init=NULL,
                              max_iter=60L, eps=1e-9) {
  nloci <- length(pgeno)
  nind <- ncol(pgeno)
  pi <- c(0.5, 0.5)
  resp <- matrix(0, nrow=nind, ncol=2)

  # prune some individuals with high missingness, so we can get an initial
  # estimate of cluster responsibility
  remove_inds <- pruneWhichIndividuals(pgeno, na_thresh)

  # initialize responsibilities
  if (is.null(init))
    a <- runif(nind)
  else
    a <- init
  init_resp <- resp <- rbind(a, 1-a)

  # calculate likelihoods for each shared parent allele being (0, 1) ONCE for
  # all individuals, progeny
  liks <- geno_ll(pgeno, fgeno, other_parents, freqs, ehet, ehom)

  thetas <- list()
  lls_inds <- list(list(), list())
  lls_haps <- list()
  #clusters <- list()
  converged <- FALSE

  for (i in seq_len(max_iter)) {
    if (i > 1) {
      # get the total likehood given last MLE for all individuals' MLE alleles, for k in (1, 2)
      # TODO: there is a VERY nasty bug in this the C++ func called below.
      #mlm1 <- max_ll_matrix(liks[[1]], liks[[2]], theta_mle[, 1])
      #mlm2 <- max_ll_matrix(liks[[2]], liks[[2]], theta_mle[, 2])
      mlm1 <- do.call(rbind, lapply(1:nrow(theta_mle), function(i) liks[[1L+theta_mle[i, 1]]][i,]))
      mlm2 <- do.call(rbind, lapply(1:nrow(theta_mle), function(i) liks[[1L+theta_mle[i, 2]]][i,]))
      #browser()
      # individual log likelihoods
      ll1 <- colSums(log(mlm1), na.rm=TRUE)
      ll2 <- colSums(log(mlm2), na.rm=TRUE)
      # append these log likeloods to list
      lls_inds[[1]] <- c(lls_inds[[1]], list(ll1))
      lls_inds[[2]] <- c(lls_inds[[2]], list(ll2))
      #browser()

      # check for convergence, ll1(t), ll2(t), ll1(t-1), ll2(t-1)
      # use t - 2 here, since we just appened these LLs to list
      if (i > 2) {
        converged <- isConverged(ll1, lls_inds[[1]][[i-2]], eps) && isConverged(ll2, lls_inds[[2]][[i-2]], eps)
        if (converged)
          break
      }

     # E step
     ll_resp <- logSumExpPosterior(list(mlm1, mlm2), pi)
     resp <- exp(logSumExpPosterior(list(mlm1, mlm2), pi))
     #clusters <- c(clusters, list(resp))
     #browser()
    }

    # M step
    # compute pi weights
    if (free_pi) # allow free-varying pi
      pi <- rowMeans(resp)

    # compute likelihoods for each (loci, individual). Returns a list of two matrices:
    # likelihood for 0 and 1 alleles of theta_k (k = haplotype)
    # TODO name, not log lik
    ll_weighted <- lapply(1L:2L, function(k) {
                          lapply(liks, function(lik) {
                                 # loops over likelihoods for an allele in (0, 1)
                                 sweep(log(lik), MARGIN=2, STATS=resp[k,], FUN=`*`)
                          })
    })
    # log likelihoods of each loci's allele, for both haplotypes k in (1, 2)
    # TODO: missingness handled acceptably here?
    lls <- lapply(ll_weighted, function(hap) do.call(cbind, lapply(hap, function(x) rowSums(x, na.rm=TRUE))))
    #browser()
    lls_haps <- c(lls_haps, lls)

    # MLE alleles for each haplotype k in (1, 2)
    theta_mle <- do.call(cbind, lapply(lls, function(x) apply(x, 1, whichmaxNA)-1L))
    thetas <- c(thetas, list(theta_mle))
  }
  iter_ind_lls <- reshapeIndLL(lls_inds)
  iter_hap_lls <- reshapeHapLL(lls_haps)
  return(list(haplos=theta_mle, cluster=resp, pi=pi, niter=i-1, converged=converged,
         ll=rbind(ll1, ll2), iter_theta=thetas, iter_ind_lls=iter_ind_lls, iter_hap_lls=iter_hap_lls,
         init=a))
}

matchLabels <- function(x, y) {
  # return version of x, with most consistent labels compared to
  # x.
  tmp <- na.exclude(cbind(x, y))
  a <- tmp[, 1] == tmp[, 2]
  b <- tmp[, 1] == ifelse(tmp[, 2] == 1, 2L, 1L)
  if (sum(a) > sum(b))
    return(y)
  return(ifelse(y == 1L, 2L, 1L))
}

simpleLigation <- function(res, buffer=0.2) {
  tmp <- lapply(res, function(x) {
                # if not clearly near 0 or 1, make NA
                membership <- apply(x$cluster, 2, which.max)
                probs <- apply(x$cluster, 2, max)
                membership <- ifelse(abs(probs - 0.5) < buffer, NA, membership)
                list(probs, membership)
    })
  probs <- do.call(rbind, lapply(tmp, '[[', 1))
  membership <- do.call(rbind, lapply(tmp, '[[', 2))
  matched_labels <- list(membership[1, ])
  for (i in seq(2, nrow(membership))) {
    ml <- matchLabels(matched_labels[[i-1]], membership[i,])
    matched_labels <- c(matched_labels, list(ml))
  }
  do.call(rbind, matched_labels)
}

physicalTiles <- function(x, chrom, tilewidth) {
  tiles <- unlist(tileGenome(seqlengths(x@ranges)[chrom], tilewidth=tilewidth))
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  # convert these tiles to a list of overlapping row loci indices
  ranges_chrom <- x@ranges[this_chrom]
  lapply(split(tiles, start(tiles)), function(tile_i) {
         hits <- findOverlaps(ranges_chrom, tile_i)
         loci <- queryHits(hits)
         loci
  })
}

snpTiles <- function(x, chrom, nsnps) {
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  ranges_chrom <- x@ranges[this_chrom]
  nloci <- length(ranges_chrom)
  m <- floor(nloci/nsnps)
  group <- rep(seq_len(m), each=nsnps)
  # append small end group
  group <- c(group, rep(m+1, nloci - length(group)))
  stopifnot(length(group) == length(ranges_chrom))
  split(seq_len(nloci), group)
}

snpTiles2Physical <- function(phasing, x, chrom, tiles) {
  # check that we have big enough enough SNP tiles to trust
  # the length of first tile as nsnps
  stopifnot(length(tiles[[1]]) == length(tiles[[2]]))
  nsnps <- length(tiles[[1]])
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  ranges_chrom <- x@ranges[this_chrom]
  tmp <- do.call(rbind, lapply(tiles, function(x) {
    start <- start(ranges_chrom[min(x)])
    end <- end(ranges_chrom[max(x)])
    c(start, end)
  }))
  data.frame(loci=1:length(tiles), start=tmp[, 1], end=tmp[, 2])
}

plotPhasingSnpTiles <- function(x, phasing, chrom, tiles, buffer=1e4, reverse_colors=FALSE) {

  if (is.na(buffer)) # choose a good one
    buffer <- seqlengths(x@ranges)[chrom] * 0.0003318355 # this looked good
  d <- ligation2dataframe(simpleLigation(phasing))
  pos <- snpTiles2Physical(phasing, x, chrom, tiles)
  d <- cbind(d, pos[match(d$loci, pos$loci), ])
  if (buffer > 0) {
    buf <- floor(buffer/2)
    d$start <- d$start + buf
    d$end <- d$end - buf
  }
  if (reverse_colors)
    d$hap <- ifelse(d$hap == "1", "2", "1")
  ggplot(d) + geom_segment(aes(x=start, xend=end, y=ind, yend=ind, color=as.factor(hap)), size=3)
}

accuracy <- function(x, y) {
  tbl <- table(x, y)
  sum(tbl[1,1] + tbl[2,2])/sum(tbl)
}

reorderByClosest <- function(x, y) {
  # given two loci x 2 matrices of haplotypes (LL or alleles) create an index
  # that reorders y so that it is the closest corresponding haplotype
  # someCompareFun(x, y[, reoderByClosest(x, y)])
  ac1 <- accuracy(x[, 1], y[, 1]) # correct orderr
  ac2 <- accuracy(x[, 1], y[, 2])
  if (ac2 > ac1)
    return(2L:1L)  # swamp columns
  return(1L:2L)
}

#' Agreement Rate (across both haplotypes) for all windows
#'
#' @param ph1 phasing results from haplotype imputation 1
#' @param ph2 phasing results from haplotype imputation 2
hapAgreementStats <- function(ph1, ph2, x=NULL, tiles=NULL, chrom=NULL) {
  # ff(): set levels so homozygous regions can be compared
  ff <- function(x) factor(x, levels=0:1)
  if (!is.null(tiles) && !is.null(chrom))
    pos <- snpTiles2Physical(ph1, x, chrom, tiles)
  accs <- lapply(1:length(ph1), function(i) {
                 # haplos
                 x <- ph1[[i]]$haplos
                 y <- ph2[[i]]$haplos
                 progeny_na <- ph1[[i]]$progeny_na # same for both; same underlying geno data
                 father_na <- ph1[[i]]$father_na
                 # reorder columns to match haplotype labels so they're consistent
                 swap <- reorderByClosest(x, y)
                 y <- y[, swap]
                 is_swapped <- all(swap == 2L:1L)
                 data.frame(window=i, h1=accuracy(x[, 1], y[, 1]), h2=accuracy(x[, 2], y[, 2]),
                            swapped=is_swapped, progeny_na, father_na)
  })
  d <- do.call(rbind, accs)
  if (!is.null(tiles))
    d$window <- rowMeans(pos[, 2:3])
  d
}

getHapLikelihoods <- function(x) {
  # convenience function to get haplotype log likelihoods from last iteration
  d <- x$iter_hap_lls
  d[x$niter == d$iter, ]
}

hapAgreementPerTile <- function(ph1, ph2, x=NULL, tiles=NULL, chrom=NULL) {
  # version of aggrement that returns raw matrix of
  # when imputed haplotypes agree
  ff <- function(x) factor(x, levels=0:1)
  if (!is.null(tiles) && !is.null(chrom))
    pos <- snpTiles2Physical(ph1, x, chrom, tiles)
  ff <- function(x) factor(x, levels=0:1)
  accs <- lapply(1:length(ph), function(i) {
                 # haplos
                 x <- ph1[[i]]$haplos
                 y <- ph2[[i]]$haplos
                 progeny_na <- ph1[[i]]$progeny_na # same for both; same underlying geno data
                 father_na <- ph1[[i]]$father_na
                 # reorder columns to match haplotype labels so they're consistent
                 swap <- reorderByClosest(x, y)
                 y <- y[, swap]
                 is_swapped <- all(swap == 2L:1L)
                 matches <- x == y
                 # extract last iteration's max likelihoods, then swap accordingly
                 lls1 <- getHapLikelihoods(ph1[[i]])[, c("h1", "h2")]
                 lls2 <- getHapLikelihoods(ph2[[i]])[, c("h1", "h2")][, swap]
                 data.frame(window=i, hm1=matches[, 1], hm2=matches[, 2],
                            ph1_h1=lls1[, 1], ph1_h2=lls1[, 2],
                            ph2_h1=lls2[, 1], ph2_h2=lls2[, 2])
  })
  do.call(rbind, accs)
}



#' Phase Sibling Family by Haplotype Imputation

#'
#' @param x a ProgenyArray object
#' @param parent parent for which to impute haplotypes
#' @param chrom which chromosome to use
#' @param tiles a list of loci indices for each phasing block
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
#' @param which_parent which parent (column in parents(x)) to use
#' @param na_thresh missingness threshold to remove individual from clustering
phaseSibFamily <- function(x, parent, chrom, tiles,
                           ehet=0.8, ehom=0.1,
                           which_parent="parent_1",
                           na_thresh=0.8) {
  if (!(chrom %in% as.character(seqnames(x@ranges))))
    stop(sprintf("chromosome '%s' not in loci", as.character(chrom)))

  # get progeny for the supplied parent (according to which parent it is) --
  pars <- parents(x)
  d <- pars[pars[, as.character(which_parent)] == parent, ]
  other_parent <- setdiff(c("parent_1", "parent_2"), which_parent)

  # get genotype data for this chromosome
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  pgeno <- progenyGenotypes(x)[this_chrom, , drop=FALSE]
  fgeno <- parentGenotypes(x)[this_chrom, , drop=FALSE]
  stopifnot(nrow(pgeno) == nrow(fgeno))
  freqs <- freqs(x)[this_chrom]

  other_parent <- setdiff(c("parent_1", "parent_2"), which_parent)
  # check that all other parents are within bounds parent genotype matrix
  stopifnot(d[, other_parent] %in% 1:ncol(fgeno))

  other_parents <- d[, other_parent] - 1L # since C++ is 0 indexed

  # do we trust these tiles? largest index in last tile should be
  # total length of rows
  stopifnot(max(tiles[[length(tiles)]]) == nrow(pgeno))
  # impute each window
  out <- lapply(tiles, function(loci) {
                message(sprintf("on loci number %d", min(loci)))
                res <- emHaplotypeImpute(pgeno[loci, d$progeny],
                                         fgeno[loci, ],
                                         other_parents,
                                         freqs[loci], ehet, ehom, na_thresh=na_thresh)
                res$nloci <- length(loci)
                res$progeny_na <- sum(is.na(pgeno[loci, d$progeny]))/length(pgeno[loci, d$progeny])
                res$progeny_na_per_loci <- rowMeans(is.na(pgeno[loci, d$progeny]))
                res$father_na <- sum(is.na(fgeno[loci, ]))/length(fgeno[loci, ])
                res$completeness <- sum(apply(!is.na(res$haplos), 1, all))/nrow(res$haplos)
                return(res)
  })
  out
}

#' A diagnostic plot of haplotype imputation of simulated data
#'
#' @param res phasing results
#' @param sim simulation data
emPlot <- function(res, sim) {
  real_haps <- sim$haplos[[1]]
  r1 <- confusionMatrix(res$haplos[, 1], real_haps[, 1])$overall['Accuracy']
  r2 <- confusionMatrix(res$haplos[, 1], real_haps[, 2])$overall['Accuracy']
  hap_order <- 1:2
  if (r2 > r1)
    hap_order <- 2:1
  #ani.record(reset = TRUE)
  var_sites <- !apply(real_haps == 0, 1, all) # all zero sites
  #orig <- haps2dataframe(data.frame(rep(FALSE, nrow(real_haps)), rep(FALSE, nrow(real_haps))))
  #levelplot(z ~ x + y, orig, col.regions=c("red","white"))
  iters <- as.integer(as.character(res$all_lls$iter))
  saveGIF({
    mapply(function(tt, i) {
         #print(table(tt[var_sites, ] == real_haps[var_sites, hap_order]))
         #p <- levelplot(z ~ x + y, haps2dataframe((tt[var_sites, ] == real_haps[var_sites, hap_order])), col.regions=c("red","white"))
         d <- haps2dataframe(1L + (tt[var_sites, ] == real_haps[var_sites, hap_order]))
         d$z <- factor(d$z, levels=1:2)
         #p <- ggplot(d) + geom_tile(aes(x=x,y=y, fill=z))
         p1 <- levelplot(z ~ x + y, d, col.regions=c("red","grey"), colorkey=FALSE, xlab="locus", ylab="haplotype\n(incorrect alleles are red)",
                         scales=list(y=list(at=1:2, labels=c("h1", "h2"))))
         p2 <- xyplot(ll ~ iter | ind, res$all_lls[iters <= i , ], group=hap, type='l', xlim=c(0, max(iters)),
                      ylab="log likelihood", xlab="iteration", ylim=c(1.1*min(res$all_lls$ll), 0.9*max(res$all_lls$ll)))
         grid.arrange(p1, p2, ncol=2)
         #ani.record()
    }, res$all_thetas, 1:max(iters))
    #oopts = ani.options(interval = 0.5)
  }, ani.width = 1000, ani.height=800)
}



#
#
#
#haploImputeWindow <- function(geno, parent_geno, fathers, freqs,
#                              na_thresh=0.6, ehet=0.8, ehom=0.1) {
#  # haplotype imputation, first by IBS clustering into two shared
#  # haplotype groups, and then ML imputation of these group's
#  # alleles
#  predict <- FALSE # feature not added yet
#  stopifnot(nrow(geno) == nrow(parent_geno))
#
#  # prune individuals with high missingness
#  remove_inds <- pruneWhichIndividuals(geno, na_thresh)
#  geno_na <- na.omit(geno[, -remove_inds])
#  if (nrow(geno_na) < 1) {
#    warning("cannot run k-means on this window: all loci contian NA values")
#    browser()
#    return(NULL)
#  }
#  nremoved <- length(na.action(geno_na))
#  if (ncol(geno_na) < 1) {
#    # no individuals available to cluster (all NA); return NA likelihood
#    empty <- list(matrix(NA, nrow=nrow(geno), ncol=2),
#                  matrix(NA, nrow=nrow(geno), ncol=2))
#    return(list(ll=empty, haps=NA, pruned=remove_inds))
#  }
#  # IBS cluter individuals into two groups, by shared haplotype
#  km <- cclust(t(geno_na), 2)
#  clusters <- data.frame(ind=1:ncol(geno), hap=NA)
#  clusters$hap[-remove_inds] <- km$cluster # for complete individuals
#  if (predict) {
#    # Note: this feature not stable yet
#    # predict haplotype of missing individuals
#    # for some reason, the number of columns cclust uses changes.
#    # so we extract these from km
#    unclustered_genos <- na.exclude(t(geno)[remove_inds, colnames(km$centers)])
#    if (nrow(unclustered_genos) > 0) {
#      fit_clust <- predict(km, unclustered_genos)
#      clusters$hap[remove_inds] <- fit_clust$cluster
#    }
#  }
#  # impute alleles for both (1, 2) shared haplotypes
#  hap_geno <- list(geno[, clusters$hap == 1], geno[, clusters$hap == 2])
#  out <- lapply(1:2, function(which_hap) {
#                # iterate haplotype imputation for each cluster
#                hap_geno <- geno[, which(clusters$hap == which_hap), drop=FALSE]
#                hap_fathers <- fathers[which(clusters$hap == which_hap)]
#                ll <- lapply(1:nrow(hap_geno), function(i) {
#                             allele_ll(hap_geno[i, ], ehet, ehom, freqs[i],
#                                       parent_geno[cbind(i, hap_fathers)])
#               })
#               return(do.call(rbind, ll))
#  })
#  return(list(ll=out, haps=clusters$hap, pruned=remove_inds))
#}
#
#whichmaxNA <- function(x) {
#  # which.max that returns NA if all data is NA
#  ifelse(any(is.na(x)), NA, which.max(x))
#}
#
#reshapeResults <- function(x) {
#  # Munging function, used to clean/restructure data from phaseSibFamily
#  ll <- lapply(x, function(res_per_tile) {
#         llscores <- res_per_tile[[1]]
#         hap1 <- apply(llscores[[1]], 1, whichmaxNA)-1L
#         hap2 <- apply(llscores[[2]], 1, whichmaxNA)-1L
#         data.frame(hap1=hap1, hap2=hap2,
#                    ll1_0=llscores[[1]][, 1],
#                    ll1_1=llscores[[1]][, 2],
#                    ll2_0=llscores[[2]][, 2],
#                    ll2_1=llscores[[2]][, 2])
#  })
#  haps <- lapply(x, function(res_per_tile) {
#                 res_per_tile[[2]]
#  })
#  removed <- sapply(x, function(res_per_tile) {
#                 res_per_tile[[3]]
#  })
#  list(ll=ll, haps=haps, removed_inds=removed)
#}
#
