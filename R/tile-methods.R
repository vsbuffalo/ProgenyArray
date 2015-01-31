# tile-methods.R

## TODO: validation methods?

# physical tiles for a single chromosome
.physicalTiles <- function(x, chrom, tilewidth) {
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

# SNP tiles for a single chromosome
.snpTiles <- function(x, chrom, nsnps) {
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

#' Load genetic map of four columns: marker, chrom, position, interval, cumulative
loadGeneticMap <- function(mapfile) {
  cols <- c(marker="factor", seqnames="factor", position="integer",
            interval="numeric", cumulative="numeric")
  read.delim(mapfile, colClasses=cols, col.names=names(cols))
}

#' This uses a monotonic polynomial regression (9th order by default) to smooth
#' over the genentic map marker windows, and then extrapolate the genetic
#' distance between genotyping markers.
approxGeneticMap <- function(genmap, degree=9) {
  # check if multicore available
  ncores <- getOption("mc.cores", 2L)
  lpfun <- lapply
  if (ncores > 1) lpfun <- mclapply

  gm <- genmap
  chrs <- split(gm, gm$seqnames)
  if (!all(!sapply(chrs, function(x) is.unsorted(x$position) & is.unsorted(x$cumulative))))
      stop("Genetic map unsorted; either position or cumulative distance not sorted.")

  approx <- lpfun(chrs, function(d) monpol(cumulative ~ position, data=d,
                                           degree=degree))
  names(approx) <- unique(gm$seqnames)
  list(fits=approx, map=genmap)
}

predictGeneticMap <- function(approx, ranges, as_df=TRUE) {
  # make data dataframe from range object
  data <- data.frame(seqnames=seqnames(ranges), position=start(ranges))
  tmp <- unique(data$seqnames)
  keep_chrs <- tmp %in% unique(approx$map$seqnames)
  if (!all(keep_chrs)) {
      msg <- "%d chromosomes in ProgenyArray object do not matching chromosomes in genetic map: %s"
      warning(sprintf(msg, sum(keep_chrs), paste(tmp[!keep_chrs], sep=", ")))
  }
  # drop chroms not in gen map
  data <- data[data$seqnames %in% unique(approx$map$seqnames), ]
  data$seqnames <- droplevels(data$seqnames)
  pred <- unlist(lapply(split(data, data$seqnames), function(d) {
      fit <- approx$fits[[as.character(d$seqnames[1])]]
      unname(predict(fit, d)[, 1])
  }))
  if (!as_df) return(pred)
  data$pred <- pred
  data
}

plotGeneticMap <- function(genmap, approx) {
  p <- ggplot(approx) 
  p <- p + facet_wrap(~seqnames, scales="free")
  p <- p + geom_point(data=genmap, aes(x=position, y=cumulative), size=1)
  p <- p + geom_line(aes(x=position, y=gendist), color="blue")
  p
}

#' Create a PhasingTiles object for tiles based on physical distances
#' @param x a ProgenyArray object
#' @param width width in bases of each tile (except last)
#'
#' @export
physicalTiles <- function(x, width) {
  nchroms <- nlevels(seqnames(x@ranges))
  tiles <- lapply(seqlevels(x@ranges),
                  function(chrom) .physicalTiles(x, chrom, width))
  names(tiles) <- seqlevels(x@ranges)
  new("PhasingTiles", "physical", as.integer(width), tiles)
}

#' Create a PhasingTiles object for tiles based on fixed number of SNPs per tile
#' @param x a ProgenyArray object
#' @param nsnps number of SNPs per tile
#'
#' @export
snpTiles <- function(x, nsnps) {
  nchroms <- nlevels(seqnames(x@ranges))
  tiles <- lapply(seqlevels(x@ranges),
                  function(chrom) .snpTiles(x, chrom, nsnps))
  names(tiles) <- seqlevels(x@ranges)
  new("PhasingTiles", type="snp", width=as.integer(nsnps), tiles=tiles)
}

#' Create a PhasingTiles object for tiles based on genetic distance
#'
#' @param x a ProgenyArray object
#' @param map a \code{data.frame} of three columns: chrom, position, genetic distance
#' @param width the genetic distance (in centiMorgans) of each tile
#'
#' @export
geneticTiles <- function(x, map, width) {
    stop("Not implemented yet.")
}


#' Pretty-print a PhasingTiles object
#'
#' @param object a PhasingTiles object
#'
#' @export
setMethod("show",
          c(object="PhasingTiles"),
          function(object) {
              cat(sprintf("PhasingTiles (%s)\n", object@type))
              cat(sprintf(" number of chromosomes: %d\n", length(object@tiles)))
              cat(sprintf(" width: %d\n", object@width))
          })
