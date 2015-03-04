# phasing.R

#
# TODO: na_thresh no longer used
#
# Refactor ideas
# - would be good to rework nested lapply()s into func

whichmaxNA <- function(x) {
  # which.max that returns NA if any data is NA
  ifelse(any(is.na(x)), NA, which.max(x))
}

pruneWhichIndividuals <- function(geno, na_thresh) {
  # Return which individuals (columns) do not meet NA threshold
  # NOTE: this is no longer used
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

rowSumsNA <- function(x) {
  apply(x, 1, function(x) {
        if (all(is.na(x)))
          return(NA)
        sum(x, na.rm=TRUE)
  })
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
          d[, 1] <- ifelse(d[, 1] == 0, NA, d[, 1])
          d[, 2] <- ifelse(d[, 2] == 0, NA, d[, 2])
          setNames(d, c("h1", "h2", "loci", "iter"))
  }, lls_haps, seq_along(lls_haps), SIMPLIFY=FALSE))
}


#' Use EM Haplotype Phasing by Imputation Algorithm (internal function)
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
#' @param init optional initial reponsibilities (advanced use)
#' @param extra include optional debugging output (advanced use)
#' @param eps epsilon in log likelihood before calling convergence
emHaplotypeImpute <- function(pgeno, fgeno, other_parents, freqs, ehet, ehom,
                              free_pi=FALSE, na_thresh=0.8, init=NULL, extra=FALSE,
                              max_iter=60L, eps=1e-9) {
  # initialize EM pi and responsibility
  nloci <- length(pgeno)
  nind <- ncol(pgeno)
  pi <- c(0.5, 0.5)
  resp <- matrix(0, nrow=nind, ncol=2)

  # prune some individuals with high missingness, so we can get an initial
  # estimate of cluster responsibility
  # REMOVED; not needed anymore
  #remove_inds <- pruneWhichIndividuals(pgeno, na_thresh)

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
      # earlier version useding C++ func; 
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
    lls <- lapply(ll_weighted, function(hap) do.call(cbind, lapply(hap, function(x) rowSumsNA(x))))
    #browser()
    lls_haps <- c(lls_haps, lls)

    # MLE alleles for each haplotype k in (1, 2)
    #browser()
    theta_mle <- do.call(cbind, lapply(lls, function(x) apply(x, 1, whichmaxNA)-1L))
    thetas <- c(thetas, list(theta_mle))
  }
  iter_ind_lls <- reshapeIndLL(lls_inds)
  iter_hap_lls <- reshapeHapLL(lls_haps)
  stopifnot(ncol(resp) == ncol(pgeno)) # this should *always* be the case
  colnames(resp) <- colnames(pgeno)
  out <- list(haplos=theta_mle, cluster=resp, pi=pi, niter=i-1, converged=converged,
              ll=rbind(ll1, ll2), init=a, haplos_lls=lls)
  if (extra)
      out <- c(out, list(iter_theta=thetas, iter_ind_lls=iter_ind_lls,
                         iter_hap_lls=iter_hap_lls))
  out
}

#' Phase all chromosomes
#'
#'
#phase <- function(x, 

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
    #ml <- relabel(rbind(matched_labels[[i-1]], membership[i,]))$cls
    matched_labels <- c(matched_labels, list(ml))
  }
  do.call(rbind, matched_labels)
}

snpTiles2Physical <- function(x, chrom, tiles) {
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
  data.frame(tile=1:length(tiles), start=tmp[, 1], end=tmp[, 2],
             midpoint=rowMeans(tmp[, 1:2]))
}

plotPhasingSnpTiles <- function(x, phasing, chrom, tiles, buffer=1e4, reverse_colors=FALSE) {

  if (is.na(buffer)) # choose a good one
    buffer <- seqlengths(x@ranges)[chrom] * 0.0003318355 # this looked good
  d <- ligation2dataframe(simpleLigation(phasing))
  pos <- snpTiles2Physical(x, chrom, tiles)
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
    pos <- snpTiles2Physical(x, chrom, tiles)
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
    pos <- snpTiles2Physical(x, chrom, tiles)
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
                 ls2 <- getHapLikelihoods(ph2[[i]])[, c("h1", "h2")][, swap]
                 data.frame(window=i, hm1=matches[, 1], hm2=matches[, 2],
                            ph1_h1=lls1[, 1], ph1_h2=lls1[, 2],
                            ph2_h1=lls2[, 1], ph2_h2=lls2[, 2])
  })
  do.call(rbind, accs)
}


estimateTileGenotypeErrors <- function(x, ph, parent, tiles, chrom) {
  pgeno <- progenyGenotypes(x, seqname=chrom)[, parent]
  lapply(seq_along(tiles), function(i) {
         tile <- tiles[[i]]
         raw_geno <- pgeno[tile]
         imputed <- rowSums(ph[[i]]$haplos)
         table(imputed, raw_geno)
  })
}

errorRates <- function(tile_errors) {
  #TODO
}



#' Phase Sibling Family by Haplotype Imputation

#'
#' @param x a ProgenyArray object
#' @param parent parent for which to impute haplotypes
#' @param chrom which chromosome to use
#' @param tiles a PhasingTiles object
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
#' @param which_parent which parent (column in parents(x)) to use
#' @param na_thresh missingness threshold to remove individual from clustering
phaseSibFamily <- function(x, parent, chrom, tiles,
                           ehet=0.8, ehom=0.1,
                           which_parent="parent_1",
                           na_thresh=0.8) {
  if (!(chrom %in% seqlevels(x@ranges)))
    stop(sprintf("chromosome '%s' not in loci", as.character(chrom)))

  # get progeny for the supplied parent (according to which parent it is) --
  pars <- parents(x)
  # Grab the correct column (according to which_parent):
  d <- pars[pars[, as.character(which_parent)] == parent, ]
  other_parent <- setdiff(c("parent_1", "parent_2"), which_parent)

  # get genotype data for this chromosome
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  pgeno <- progenyGenotypes(x, seqname=chrom)
  fgeno <- parentGenotypes(x, seqname=chrom)
  stopifnot(nrow(pgeno) == nrow(fgeno))
  freqs <- freqs(x)[this_chrom]

  other_parent <- setdiff(c("parent_1", "parent_2"), which_parent)
  # check that all other parents are within bounds parent genotype matrix
  stopifnot(d[, other_parent] %in% 1:ncol(fgeno))

  other_parents <- d[, other_parent] - 1L # since C++ is 0 indexed
  
  # impute each window using the indices in a tile
  out <- lapply(tiles@tiles[[chrom]], function(loci) {
                #message(sprintf("on loci number %d", min(loci)))
                res <- emHaplotypeImpute(pgeno[loci, d$progeny],
                                         fgeno[loci, ],
                                         other_parents,
                                         freqs[loci], ehet, ehom, na_thresh=na_thresh)
                # Store diagnostic information
                res$nloci <- length(loci)
                res$progeny_na <- sum(is.na(pgeno[loci, d$progeny]))/length(pgeno[loci, d$progeny])
                res$progeny_na_per_loci <- rowMeans(is.na(pgeno[loci, d$progeny]))
                res$father_na <- sum(is.na(fgeno[loci, ]))/length(fgeno[loci, ])
                res$completeness <- sum(apply(!is.na(res$haplos), 1, all))/nrow(res$haplos)
                return(res)
  })
  out
}

phasingMetadata <- function(tiles, ehet, ehom, na_thresh) {
    list(tile_type=tiles@type, tile_width=tiles@width,
         ehet=ehet, ehom=ehom, na_thresh=na_thresh)
}

#' Phase all Parents in a ProgenyArray object
#'
#' @param x a ProgenyArray object
#' @param tiles a PhasingTiles object
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
#' @param na_thresh how much missingness to tolerate before pruning individual
#' @param parallel whether to run across multiple cores (use option mc.cores to set number of cores)
#' @param verbose report status messages
setMethod("phaseParents", c(x="ProgenyArray"),
          function(x, tiles, ehet=0.8, ehom=0.1, na_thresh=0.8,
                   parallel=FALSE, verbose=TRUE) {
              ncores <- getOption("mc.cores")
              lpfun <- if (!is.null(ncores) && ncores > 1) mclapply else lapply
              x@tiles <- tiles # add tiles to ProgenyArray Object
              # first, note that some mothers may be inconsistent; that is
              # the inferred parent may not be a mother included. We use NA for
              # these.
              pars <- seq_len(ncol(parentGenotypes(x)))
              x@phased_parents <- lpfun(pars, function(par) {
                  chroms <- names(x@tiles@tiles)# uses tile chromosome not slot @ranges!
                  setNames(lpfun(chroms, function(chr) {
                      if (verbose) message(sprintf("phasing parent %d, chrom %s", par, chr))
                      phaseSibFamily(x, par, chr, tiles, ehet=ehet,
                                     ehom=ehom, na_thresh=na_thresh)
                  }), chroms)
              })
              names(x@phased_parents) <- colnames(parentGenotypes(x))
              x@phased_parents_metadata <- phasingMetadata(tiles, ehet, ehom, na_thresh)
              return(x)
          })

#' A function that gets the reference alleles used in the tiles.  This is an
#' important difference from ref() which returns all reference alleles.
tileRef <- function(x) {
    # TODO; this is messy; a lot of this should be cleaned up with higher level
    # funcs. In other words, this is fragile hacky code that bears the signature
    # of Vance Bison, not Vince Buffalo.
    # TODO: this could be moved to the tiling code;
    stopifnot(length(x@ranges) == length(x@ref))
    all_chroms <- names(x@tiles@tiles)
    # ditch chromosomes not in tiles:
    # have to convert to character beacuse stupid
    # S4 Rle vec won't work in indexing later
    seqnames_char <- as.character(seqnames(x@ranges))  
    chrs_tiles <- seqnames_char %in% all_chroms
    chrs <- seqnames_char[chrs_tiles]
    refs <- x@ref[chrs_tiles]
    stopifnot(length(x@tiles@tiles) == length(all_chroms))
    unlist(Map(function(ref_chr, tile_chr) {
        ref_chr[unlist(tile_chr)]
    }, split(refs, chrs), x@tiles@tiles))
}

tileAlt <- function(x) {
    # some as above for alt; TODO refactor
    stopifnot(length(x@ranges) == length(x@alt))
    all_chroms <- names(x@tiles@tiles)
    # ditch chromosomes not in tiles:
    seqnames_char <- as.character(seqnames(x@ranges))  
    chrs_tiles <- seqnames_char %in% all_chroms
    chrs <- seqnames_char[chrs_tiles]
    alts <- x@alt[chrs_tiles]
    stopifnot(length(x@tiles@tiles) == length(all_chroms))
    unlist(Map(function(alt_chr, tile_chr) {
        alt_chr[unlist(tile_chr)]
    }, split(alts, chrs), x@tiles@tiles))
}


bindProgenyHaplotypes <- function(x, parent, progeny) {
    # create a haplotype for a progeny given the parent's phased haplotypes
    do.call(rbind, lapply(x@phased_parents[[parent]], function(chrom) {
        # loop through all chromosomes of one parent
        do.call(rbind, lapply(chrom, function(tile) {
            # loop through all tiles of one chromosome
            # TODO this is quite inefficient
            clts <- apply(tile$cluster, 2, which.max)[progeny]
            out <- tile$haplos[, clts, drop=FALSE]
            if (all(is.na(out))) browser()
            out
        }))
    }))
}

#' Extract all progeny haplotypes from all parents
extractProgenyHaplotypes <- function(x) {
    progeny_names <- colnames(progenyGenotypes(pa))
    parent_names <- colnames(parentGenotypes(pa))
    parent_1 <- parent_names[x@parents$parent_1]
    parent_2 <- parent_names[x@parents$parent_2]
    progeny <- progeny_names[x@parents$progeny]
    Map(function(prog, p1, p2) {
        par1 <- bindProgenyHaplotypes(pa, p1, prog)
        par2 <- bindProgenyHaplotypes(pa, p2, prog)
        browser()
        stopifnot(dim(par1) == dim(par2))
        if (length(x@ref) > 0 && length(x@alt) > 0) {
            ref <- tileRef(x)
            alt <- tileAlt(x)
            stopifnot(length(ref) == nrow(par1) && length(alt) == nrow(par2))
        } else {
            ref <- NULL; alt <- NULL
        }
        do.call(cbind, lapply(seq_len(ncol(par1)), function(i) {
            encodeHaplotype(par1[, i], par2[, i], ref, alt)
        }))
    }, progeny, parent_1, parent_2)
}



# Take two haplotype vectors and merge them into strings like "0|1"
encodeHaplotype <- function(hap1, hap2, ref=NULL, alt=NULL) {
    stopifnot(length(hap1) == length(hap2))
    stopifnot(is.null(ref) == is.null(alt)) # both should be same
    if (!is.null(ref) && !is.null(alt)) {
        alleles <- cbind(ref, alt)
        rrows <- 1:nrow(alleles)
        # TODO check this
        hap1 <- alleles[cbind(rrows, hap1 + 1L)]
        hap2 <- alleles[cbind(rrows, hap2 + 1L)]
    }
    ifelse(is.na(hap1) | is.na(hap2), NA, paste(hap1, hap2, sep ="|"))
}

#' A function that takes a ProgenyArray object and a progeny ID, and creates
#' creates this progeny's haplotypes by looking at the cluters inferred by EM.
#' With these, this function grabs the corresponding haplotype from the parents
#' to build the progeny's haplotype.
#'
#' Note that we use the names in the cluster matrix, since this is safer than
#' using indexes.
getProgenyHaplotypesFromClusters <- function(x, progeny) {
    prog <- x@parents[x@parents$progeny == progeny, ]
    par1 <- extractProgenyHaplotypes(x, prog$parent_1, prog$progeny)
    par2 <- extractProgenyHaplotypes(x, prog$parent_2, prog$progeny)
    stopifnot(dim(par1) == dim(par2))
    if (length(x@ref) > 0 && length(x@alt) > 0) {
        ref <- tileRef(x)
        alt <- tileAlt(x)
        stopifnot(length(ref) == nrow(par1) && length(alt) == nrow(par2))
    } else {
        ref <- NULL; alt <- NULL
    }
    do.call(cbind, lapply(seq_len(ncol(par1)), function(i) {
        encodeHaplotype(par1[, i], par2[, i], ref, alt)
    }))
}


llratio <- function(x) max(x) - min(x) # TODO think about this

#' Merge phasing data
#'
#' @param x a ProgenyArray object with phased parents
#' @param include_ll a logical indicating whether to include the haplotype log-likelihoods
setMethod("phases", c(x="ProgenyArray"),
          function(x, include_ll=FALSE) {
              if (length(x@ref) == 0 || length(x@alt) == 0)
                  warning("ref/alt not set; reporting haplotypes in 0|1 format.")
              # loads of data reshaping to do. Slot phased_parents contains
              # nested lists:
              #
              # Parents
              #  Chromosomes
              #    Tiles
              #      Phasing info
              # This is a bit crufty; lots of cruft to deal with indices
              # to grab names which should have passed earlier TODO
              parent_names <- colnames(parentGenotypes(x))
              all_chroms <- names(x@tiles@tiles)
              pars <- lapply(seq_along(x@phased_parents), function(parent_i) {
                  parent <- x@phased_parents[[parent_i]]
                  do.call(rbind, lapply(seq_along(parent), function(chrom_i) {
                      chrom_name <- all_chroms[chrom_i]
                      this_chrom <- parent[[chrom_i]]
                      do.call(rbind, lapply(seq_along(this_chrom), function(tile_i) {
                          tile <- x@tiles@tiles[[chrom_i]][[tile_i]]
                          phases <- this_chrom[[tile_i]] # for this tile
                          par_haps <- encodeHaplotype(phases$haplos[, 1],
                                                      phases$haplos[, 2],
                                                      x@ref[tile], x@alt[tile])
                          if (include_ll) {
                              # TODO
                              ll0 <- apply(phases$haplos_ll[[1]], 1, llratio)
                              ll1 <- apply(phases$haplos_ll[[2]], 1, llratio)
                          }
                          out <- cbind(par_haps)
                          colnames(out) <- parent_names[parent_i]
                          out
                     }))
                  }))
              })

              if (include_ll) {
                  ll <- lapply(seq_along(x@phased_parents), function(parent_i) {
                     parent <- x@phased_parents[[parent_i]]
                     do.call(rbind, lapply(seq_along(parent), function(chrom_i) {
                         chrom_name <- all_chroms[chrom_i]
                         this_chrom <- parent[[chrom_i]]
                         do.call(rbind, lapply(seq_along(this_chrom), function(tile_i) {
                             tile <- x@tiles@tiles[[chrom_i]][[tile_i]]
                             phases <- this_chrom[[tile_i]] # for this tile
                             ll0 <- apply(phases$haplos_ll[[1]], 1, llratio)
                             ll1 <- apply(phases$haplos_ll[[2]], 1, llratio)
                             pn <- paste(c("ll0", "ll1"),
                                         parent_names[parent_i], sep="_")
                             out <- cbind(ll0, ll1)
                             colnames(out) <- pn
                             out
                         }))
                     }))
                 })
              }

              # now, we need to reconstruct each kid's phase
              prog <- extractProgenyHaplotypes(pa)
              # get positions from tiles, bind everything together.
              pos <- x@tiles@info$smoothed_genetic_map[, c("seqnames", "position")]
              pos <- setNames(pos, c("chr", "position"))
              prog_binded <- do.call(cbind, prog)
              colnames(prog_binded) <- colnames(progenyGenotypes(x))
              cbind(pos, do.call(cbind, pars), prog_binded)
          })

writePhases <- function(x, file) {
    if (regexpr(".*\\.gz", file) != -1)
        file <- gzfile(file)
    write.table(x, file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
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

extractFromPhases <- function(x, item) {
    # extract each parents phase from a list, reshaope to dataframe
    # parents, chromosomes, tiles
    unname(mapply(function(par, nm) {
        tmp <- lapply(par, function(chr) {
            do.call(rbind, lapply(chr, function(t) {
                as.data.frame(t[[item]])
            }))
        })
        d <- do.call(rbind, tmp)
        colnames(d) <- paste(nm, c("1", "2"), sep="_")
        d
    }, x, names(x), SIMPLIFY=FALSE))
}

hapLR <- function(x) {
    # create a likelihood ratio from most likely allele / least likely
    apply(lls, 1, function(x) {
        i <- which.max(x)
        x[i]/x[setdiff(1:2, i)]
    })
}

reshapeParentPhases <- function(x, tiles) {
    haps <- do.call(cbind, extractFromPhases(x@phased_parents, 'haplos'))
    lls <- extractFromPhases(x@phased_parents, 'haplos_ll')
    lrt1 <- hapLR(lls[[1]])
    lrt2 <- hapLR(lls[[2]])
    colnames(ll) <- paste("ll", colnames(ll), sep="_")
    cbind(haps, ll, lrt1, lrt2)
}

