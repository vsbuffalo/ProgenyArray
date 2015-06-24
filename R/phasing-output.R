#' A function that gets the reference alleles used in the tiles.  This is an
#' important difference from ref() which returns all reference alleles.
tileRef <- function(x) {
    if (length(x@ref) == 0) return(NULL)
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
    if (length(x@alt) == 0) return(NULL)
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

#' the LL ratio is defined as the most likeley allele's log likelihood minus the
#' least likely.
llratio <- function(x) max(x) - min(x)


getTileAlleles <- function(x, split=FALSE) {
  # get ref/alt alleles optionally splitting by chromosome
  all_alleles <- as.data.frame(cbind(tileRef(x), tileAlt(x)), stringsAsFactors=FALSE)
  colnames(all_alleles) <- c("ref", "alt")
  if (!split)
    return(all_alleles)
  grp <- unlist(Map(function(x, i) rep(i, length(unlist(x))), x@tiles@tiles, seq_along(x@tiles@tiles)))
  split(all_alleles, grp)
}

# Take two haplotype vectors and merge them into strings like "0|1"
encodeHaplotype <- function(hap1, hap2, alleles=NULL) {
    # these are painfully slow to turn on in real implementation, but are good checks:
    #stopifnot(!is.null(alleles)) # not implemented
    #stopifnot(length(hap1) == length(hap2))
    #stopifnot(sum(sapply(hap1, function(x) sum(sapply(x, length)))) == nrow(alleles))
    if (!is.null(alleles)) {
        rrows <- 1:nrow(alleles)
        # TODO check this
        hap1 <- alleles[cbind(rrows, hap1 + 1L)]
        hap2 <- alleles[cbind(rrows, hap2 + 1L)]
    }
    ifelse(is.na(hap1) | is.na(hap2), NA, paste(hap1, hap2, sep ="|"))
}

bindProgenyHaplotypes <- function(x, parent, progeny) {
    # create a haplotype for a progeny given the parent's phased haplotypes
    # unlist, not sapply; simplify2array screw everything up
    all_chroms <- names(x@tiles@tiles)
    getProgenyHaplotypeType <- function(tile) {
            # loop through all tiles of one chromosome
            # grabbing the tile for the appropriate progeny
            # moved outside of lapply() for speed
            clts <- apply(tile$cluster, 2, which.max)[progeny]
            tile$haplos[, clts, drop=FALSE]
    }
    out <- lapply(all_chroms, function(chrom) {
        # loop over all tiles for this parent/chromosome, returning haplotypes unlisted into a single vec
    	unlist(lapply(x@phased_parents[[parent]][[chrom]], getProgenyHaplotypeType))
    }) 
    out
}

createSelfedHaplotypeCombinations <- function(haplos) {
  # for two sets of imputed haplotypes 0, 1, create a set of haplotype
  # combinations that could be inherited in selfed individuals, 00, 01, 11 for
  # ML checking.
  list(haplos[, c(1, 1)], haplos[, c(1, 2)], haplos[, c(2, 2)])
}

bindSelfedProgenyHaplotypes <- function(x, parent, progeny, error_matrix) {
  # infer the correct pair of haplotypes for a selfed individual, using ML
  # approach.
  all_chroms <- names(x@tiles@tiles)
  # create haps for all chroms
  all_haps <- lapply(all_chroms, function(chrom) {
                       lapply(x@phased_parents[[parent]][[chrom]], function(x) {
                                createSelfedHaplotypeCombinations(x$haplos)
                              })
  })

  # log likelihoods of observed progeny genotypes under each of the selfing
  # combinations
  prog_geno <- progenyGenotypes(x)[, progeny]
  haplo_ll <- lapply(all_haps, function(haps) {
                       # across chromosomes
                       lapply(haps, function(hap) {
                         # across tiles
                         # calculate genotype likelihoods across haplotype combinations
                         sapply(hap, function(hapcomb) {
                           # calculate log likelihood across loci
                           sum(log(error_matrix[cbind(rowSums(hapcomb, na.rm=TRUE), prog_geno)]), na.rm=TRUE)
                         })
                       })
                     })
  max_ll_haps <- lapply(haplo_ll, function(chrom) {
                          lapply(chrom, function(tile) {
                                   tile[[which.max(tile)]]
                                 })
                        })
  out <- list(haplo_ll=haplo_ll, haps=max_ll_haps)
  browser()
}

extractProgenyHaplotypes <- function(x, included_parents, verbose=TRUE) {
  vmessage <- function(x) {
    if (verbose)
      message(x, appendLF=FALSE)
  }

  ERROR_MAT <- genotypingErrorMatrix(x@phased_parents_metadata$ehet,
                                     x@phased_parents_metadata$ehom)
 
  pars <- parents(x, use_names=TRUE)
  # only keep those parents with full sib fams
  kp <- pars$parent_1 %in% included_parents & pars$parent_2 %in% included_parents
  pars <- pars[kp, ]
  parent_1 <- pars$parent_1
  parent_2 <- pars$parent_2
  progeny <- pars$progeny
  if (length(x@ref) > 0 && length(x@alt) > 0) {
    # alleles set
    alleles <- getTileAlleles(x, split=TRUE)
  } else {
    stop("not implemented")
    alleles <- NULL
  }

  ncores <- getOption("mc.cores")
  mapfun <- if (!is.null(ncores) && ncores > 1) mcMap else Map
  n <- length(progeny)
  out <- mapfun(function(prog, p1, p2, i) {
                  message(sprintf("extracting progeny haplotypes for '%s' (parents: '%s' + '%s') %d of %d", prog, p1, p2, i, n))
                  ## new routine for selfed progeny -- can't rely on cluster responsibility probs
                  is_selfed <- p1 == p2
                  if (is_selfed) {
                    ## evaluate the three haplotypes using ML to find most
                    ## likely one. This is identical to creating parent
                    ## genotypes from the reconstructed haplotypes and finding
                    ## the most likely parent.
                    par_self <- bindSelfedProgenyHaplotypes(x, p1, prog, ERROR_MAT)
                  } else {
                    ### not selfed
                    # prog, p1, and p2 are names of progeny and both parents
                    par1 <- bindProgenyHaplotypes(x, p1, prog)
                    par2 <- bindProgenyHaplotypes(x, p2, prog)
                    # par1 and par2 are lists (of chromsoomes) of parent haplotypes for each progeny
                    # these are then bound together with encodeHaplotype
                    #stopifnot(all(sapply(par1, length) == sapply(par2, length)))
                    # by chrom:
                  }
                  Map(function(p1, p2, chrom_alleles) encodeHaplotype(p1, p2, chrom_alleles), par1, par2, alleles)
                }, progeny, parent_1, parent_2, seq_along(progeny))

  # now, collate these into a list of chrom by progeny
  all_chroms <- names(x@tiles@tiles)
  collated_out <- vector('list', length(all_chroms))
  vmessage("collating progeny haplotypes into lists by chromosome... ")
  for (chrom_i in seq_along(all_chroms)) {
    vmessage(sprintf("  collating chromosome number %d... ", chrom_i))
    collated_out[[chrom_i]] <- lapply(out, '[[', chrom_i)
  }
  vmessage("done.\n")
  collated_out
}

tileCol <- function(x) {
  # make a vector of tile numbers
  tiles <- x@tiles@tiles
  lapply(tiles, function(x) unlist(lapply(seq_along(x), function(i) rep(i, length(x[[i]])))))
}

#' Merge phasing data
#'
#' @param x a ProgenyArray object with phased parents
#' @param include_ll a logical indicating whether to include the haplotype log-likelihoods
#'
#' @export
setMethod("phases", c(x="ProgenyArray"),
          function(x, include_ll=FALSE, verbose=TRUE) {
            vmessage <- function(x) {
              if (verbose)
                message(x, appendLF=FALSE)
            }
            
            if (length(x@ref) == 0 || length(x@alt) == 0)
              warning("ref/alt not set; reporting haplotypes in 0|1 format.")
            # loads of data reshaping to do. Slot phased_parents contains
            # nested lists:
            #
            # Parents
            #  Chromosomes
            #    Tiles
            #      Phasing info
            #
            # This is a bit crufty; lots of cruft to deal with indices to grab
            # names which should have passed earlier TODO might be good to
            # refactor into separate functions at some point.
            vmessage("binding all alleles together... ")
            # get alleles at tile sites
            alleles <- getTileAlleles(x)
            vmessage("done.\n")
            vmessage("extracting all parent haplotypes... ")
            parent_names <- names(x@sibfams)
            all_chroms <- names(x@tiles@tiles)
            inc_pars <- names(x@sibfams)[!sapply(x@sibfams, is.null)]

            lpfun <- lapply # to support parallel in future; currently this is
            # disabled as shared objects cause memory issues for large
            # ProgenyArray objects

            pars <- lpfun(seq_along(all_chroms), function(chrom_i) {
                       chrom_name <- all_chroms[chrom_i]
                       all_pars_for_chrom <- lapply(seq_along(x@phased_parents), function(parent_i) {
                           parent <- x@phased_parents[[parent_i]]
                           this_chrom <- parent[[chrom_i]] 
                           if (is.null(parent)) return(NULL) # didn't phase this parent for some reason
                           o <- do.call(rbind, lapply(seq_along(this_chrom), function(tile_i) {
                               tile <- x@tiles@tiles[[chrom_i]][[tile_i]]
                               phases <- this_chrom[[tile_i]] # for this tile
                               par_haps <- encodeHaplotype(phases$haplos[, 1],
                                                           phases$haplos[, 2],
                                                           alleles[tile, ])
                               out <- as.data.frame(cbind(par_haps))
                               colnames(out) <- parent_names[parent_i]
                               out
                             }))
                           return(as.data.frame(o))
                         })
                       do.call(cbind, all_pars_for_chrom[!sapply(all_pars_for_chrom, is.null)])
                     })
            vmessage("done.\n")

            ll <- NULL # in case not included; consistent list object out
            if (include_ll) {
              vmessage("extracting all parent likelihoods... ")
              # TODO: could be cleaned up with parent/chrom/tile iterator
              ll <- lpfun(seq_along(all_chroms), function(chrom_i) {
                            chrom_name <- all_chroms[chrom_i]
                            all_pars_for_chrom <- lapply(seq_along(x@phased_parents), function(parent_i) {
                                parent <- x@phased_parents[[parent_i]]
                                this_chrom <- parent[[chrom_i]] 
                                if (is.null(parent)) return(NULL) # didn't phase this parent for some reason
                                
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
                              })
                            do.call(cbind, all_pars_for_chrom[!sapply(all_pars_for_chrom, is.null)])
                          })
              vmessage("done.\n")
            }

            # now, we need to reconstruct each kid's phase
            vmessage("extracting all progeny haplotypes... ")
            inc_pars <- names(x@sibfams)[!sapply(x@sibfams, is.null)]
            prog <- extractProgenyHaplotypes(x, inc_pars, verbose=verbose)
            vmessage("done.\n")
            # get positions from tiles, bind everything together.
            pos <- x@tiles@info$smoothed_genetic_map[, c("seqnames", "position")]
            pos <- setNames(pos, c("chr", "position"))
            pos_by_chroms <- split(pos, pos$chr)
            tiles <- tileCol(x)
            list(pos=pos_by_chroms, tiles=tiles, ll=ll, prog=prog, pars=pars)

          })
