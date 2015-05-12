## phasing-plots.R

splitJointHaplotype <- function(x) {
    do.call(rbind, lapply(strsplit(x, "|", fixed=TRUE)))
}
