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
  # TODO notes on output object.
  
  # initialize EM pi and responsibility
  nloci <- length(pgeno)
  nind <- ncol(pgeno)
  if (nloci == 0 || nind == 0)
    stop(sprintf("emHaplotypeImpute(): nind or nloci = 0"))
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
  # NOTE: we decrement other_parents since C++ is 0 indexed.
  liks <- geno_ll(pgeno, fgeno, other_parents-1L, freqs, ehet, ehom)

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
      # append these log likeloods to list; yes this is inefficient.... TODO
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
                                 sweep(log(lik), MARGIN=2, STATS=resp[k, , drop=FALSE], FUN=`*`)
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
  stopifnot(ncol(resp) == length(other_parents)) # this should *always* be the case
  colnames(resp) <- colnames(pgeno)
  out <- list(haplos=theta_mle, cluster=resp, pi=pi, niter=i-1, converged=converged,
              ll=rbind(ll1, ll2), init=a, haplos_lls=lls)
  if (extra)
      out <- c(out, list(iter_theta=thetas, iter_ind_lls=iter_ind_lls,
                         iter_hap_lls=iter_hap_lls))
  out
}

#' Create a sibling family for a single parent.
#'
#' @param x a ProgenyArray object
#' @param parent the parent to bring the sib family
#' @param ignore_selfs should we ignore selfed individuals?
#'
#' @export
sibFamily <- function(x, parent, ignore_selfs) {
  pars <- parents(x)
  # This builds a full-sib family
  tmp <- pars[pars$parent_1 == parent | pars$parent_2 == parent,
              c("progeny", "parent_1", "parent_2")]
  if (nrow(tmp) == 0) {
    warning(sprintf("sibFamily(): parent '%d' has no offspring", parent))
    return(NULL)
  }
  prog_i <- tmp$progeny
  other_parents <- apply(tmp[, -1], 1, function(x) {
      if (x[1] == x[2]) return(parent) # selfed ind; return any parent
      setdiff(x, parent)
    })
  out <- data.frame(focal_parent=parent, other_parent=other_parents, progeny=prog_i)
  if (ignore_selfs) {
    return(out[out$focal_parent != out$other_parent, ])
  }
  out
}

#' Get all sibling families for a ProgenyArray object
allSibFamilies <- function(x, ignore_selfs=FALSE) {
  pars <- seq_len(ncol(parentGenotypes(x)))
  parnames <- colnames(parentGenotypes(x))
  setNames(lapply(pars, function(p) sibFamily(x, p, ignore_selfs)), parnames)
}

#' Filter sibling families by minimum number of offspring
filterSibFamilies <- function(x, sibfamilies, min_child,
                              remove_full_sibs=TRUE) {
  # parent_names is used for friendlier messages
  msg <- "filterSibFamilies(): parent %d ('%s') has fewer than %d offspring (%d); not including in phasing/imputation."
    filtered_sibfams <- Map(function(fam, n) {
      if (is.null(fam)) return(NULL)
      if (nrow(fam) < min_child) {
        warning(sprintf(msg, fam$focal_parent[1], n, min_child, nrow(fam)))
        return(NULL)
      }

      if (remove_full_sibs) {
        # get the individual with the least missingness
        tmp <- lapply(split(fam, list(fam$focal_parent, fam$other_parent)),
                 function(d) {
                   miss <- colSums(is.na(progenyGenotypes(x)[, d$progeny, drop=FALSE]))
                   d[which.min(miss), , drop=FALSE]
               })
        nofull <- do.call(rbind, tmp)
        return(nofull)
      }
      
      return(fam)
    }, sibfamilies, names(sibfamilies))

  
  
  return(filtered_sibfams)
}

#' Phase Sibling Family by Haplotype Imputation

#'
#' @param x a ProgenyArray object
#' @param sibfam a dataframe from \code{sibFamily()}
#' @param chrom which chromosome to use
#' @param tiles a PhasingTiles object
#' @param ehet heterozygous error rate
#' @param ehom homozygous error rate
#' @param na_thresh missingness threshold to remove individual from clustering
phaseSibFamily <- function(x, sibfam, chrom, tiles,
                           ehet=0.8, ehom=0.1,
                           na_thresh=0.8) {

  if (!(chrom %in% seqlevels(x@ranges)))
    stop(sprintf("chromosome '%s' not in loci", as.character(chrom)))

  # get genotype data for this chromosome
  this_chrom <- as.logical(seqnames(x@ranges) == chrom)
  pgeno <- progenyGenotypes(x, seqname=chrom)
  fgeno <- parentGenotypes(x, seqname=chrom)
  stopifnot(nrow(pgeno) == nrow(fgeno))
  freqs <- freqs(x)[this_chrom]

  # impute each window using the indices in a tile
  out <- lapply(tiles@tiles[[chrom]], function(loci) {
                #message(sprintf("on loci number %d", min(loci)))
                res <- emHaplotypeImpute(pgeno[loci, sibfam$progeny, drop=FALSE],
                                         fgeno[loci, , drop=FALSE],
                                         sibfam$other_parent,
                                         freqs[loci], ehet, ehom, na_thresh=na_thresh)
                stopifnot(ncol(res$cluster) == length(sibfam$other_parent))
                # Store diagnostic information
                res$nloci <- length(loci)
                res$progeny_na <- sum(is.na(pgeno[loci, sibfam$progeny]))/length(pgeno[loci, sibfam$progeny])
                res$progeny_na_per_loci <- rowMeans(is.na(pgeno[loci, sibfam$progeny]))
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
#' @param min_child minimum number of children to include a parent for phasing
#' @param ignore_selfs whether to ignore selfed individuals to avoid bias (should be true, unless debugging)
#' @param verbose report status messages
#'
#' @export
setMethod("phaseParents", c(x="ProgenyArray"),
          function(x, tiles, ehet=0.8, ehom=0.1, na_thresh=0.8, min_child=8,
                   ignore_selfs=TRUE, 
                   verbose=TRUE) {
  no_full_sibs=FALSE ## for debugging
  ncores <- getOption("mc.cores")
  lpfun <- if (!is.null(ncores) && ncores > 1) mcMap else Map
  x@tiles <- tiles # add tiles to ProgenyArray Object
  # first, note that some mothers may be inconsistent; that is
  # the inferred parent may not be a mother included. We use NA for
  # these.
  
  # create sibling families, for all parents; filter out parents
  # that don't have sufficient children for phasing/imputation
  sibfams <- filterSibFamilies(x, allSibFamilies(x, ignore_selfs=ignore_selfs),
                                              remove_full_sibs=no_full_sibs,
                                              min_child)
  #browser()
  x@sibfams <- sibfams

  phaseFun <- function(sibfam, parname, i, n) {
    if (is.null(sibfam)) {
      warning(sprintf("phaseParents(): skipping parent '%s'; no sib family data", parname))
      return(NULL) # nothing to phase, either no progeny or too few
    }
    chroms <- names(x@tiles@tiles)# uses tile chromosome not slot @ranges!
    setNames(lpfun(function(chr) {
                     msg <- sprintf("phasing parent '%s' (%d of %d), chrom %s", parname, i, n, chr)
                     if (verbose) message(msg)
                     phaseSibFamily(x, sibfam, chr, tiles, ehet=ehet,
                                    ehom=ehom, na_thresh=na_thresh)
                   }, chroms), chroms)
  }
  
  x@phased_parents <- lpfun(phaseFun, sibfams, names(sibfams),
                            seq_along(sibfams), length(sibfams))

  names(x@phased_parents) <- colnames(parentGenotypes(x))
  x@phased_parents_metadata <- phasingMetadata(tiles, ehet, ehom, na_thresh)
  return(x)
})
