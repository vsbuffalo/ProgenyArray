# phasing-utilities.R -- These are older tools I used during testing phasing


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
                 data.frame(window=i, h1=accuracy(x[, 1], y[, 1]),
                            h2=accuracy(x[, 2], y[, 2]),
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


