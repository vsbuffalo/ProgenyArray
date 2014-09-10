
sampleHaplotypeFromSFS <- Vectorize(function(freqs) rbinom(1, 1, freqs), vectorize.args='freqs')

randomGametes <- function(freqs, F=0) {
  inbred <- sample(0:1, 1, prob=c(1-F, F))
  h1 <- sampleHaplotypeFromSFS(freqs)
  if (inbred) {
    return(cbind(h1, h1))
  }
  h2 <- sampleHaplotypeFromSFS(freqs)
  cbind(h1, h2)
}

createIndividuals <- function(n, freqs, F=0) {
  replicate(n, randomGametes(freqs, F),
            simplify=FALSE)
}

haplotypes2NoisyGenotypes <- function(x, ehet, ehom, na_rate=0) {
  d <- do.call(cbind, lapply(x, function(y) addGenotypeError(rowSums(y), ehet, ehom)))
  if (na_rate == 0)
    return(d)
  d[as.logical(rbinom(length(d), 1, na_rate))] <- NA_integer_
  d
}

#' Simulate a half-sib family, fixed number of sites
#'
#' @export n size of halfsib family
#' @export nsites number of sites
#' @export mgamete_probs probability of getting particular mother gamete
#' @export rates selfing, halfsib, fullsib proportions
#' @export nfullsib_fathers number of fathers contributing to fullsib families
#' @export father_F proportion of fathers that are inbred
#' @export mom_inbred whether mom is inbred
#' @export ehet heterozygous error rate
#' @export ehom homozygous error rate
#'
#' @return a list full of simulation data
#'
#' \enumerate{
#'   \item \code{progeny} progeny genotypes (with error)
#'   \item \code{mom} mother's haplotypes
#'   \item \code{anno} dataframe containing haplotypes form each parent, fathers
#'   \item \code{freqs} allele frequencies at loci
#'   \item \code{hs_haplo} halfsib father's haplotypes
#'   \item \code{fs_haplo} fullsib father's haplotypes
#'   \item \code{hs_geno} halfsib father's genotypes (with error)
#'   \item \code{fs_geno} fullsib father's genotypes (with error)
#' }
sibFamily <- function(nprogeny, nsites, mgamete_probs=c(0.5, 0.5),
                      rates=c(selfing=0.5, halfsib=0.2, fullsib=0.3),
                      nfullsib_fathers=2, mom_inbred=FALSE,
                      father_F=0.1, ehet=0.8, ehom=0.1,
                      na_rate=0) {

  # n is number of parents, an overestimate. Assuming no selfing or fullsibs,
  # this is at most 2*nprogeny (each kid has separate parents)
  n <- 2*nprogeny
  freqs <- sampleSiteFreqs(n, nsites)
  types <- factor(c("selfed", "halfsib", "fullsib"))
  ind_type <- sample(types, nprogeny, prob=rates, replace=TRUE)

  # generate mother (who is a father in selfs)
  mom <- randomGametes(freqs, F=as.numeric(mom_inbred))

  # generate fathers, fullsib fathers + number of halfsib families
  nhalfsib_families <- sum(ind_type == "halfsib")
  nfullsib_families <- sum(ind_type == "fullsib")
  parents <- createIndividuals(nfullsib_fathers+nhalfsib_families, freqs, father_F)
  parents <- c(list(mom), parents) # mom is potential father (selfing), kept in pos 1

  # assign labels to fathers, used to sample subsets for particular types
  # of fathers
  father_type <- rep(c("self", "halfsib", "fullsib"),
                     c(1, nhalfsib_families, nfullsib_fathers))
  fathers <- integer(nprogeny)
  fathers[ind_type == "selfed"] <- 1L # mother is going to be set to father 1
  fathers[ind_type == "halfsib"] <- sample(which(father_type == "halfsib"),
                                           nhalfsib_families)
  fathers[ind_type == "fullsib"] <- sample(which(father_type == "fullsib"),
                                           nfullsib_families, replace=TRUE)
  stopifnot(!any(is.na(fathers)))

  # generate which father gamete we will receive.
  father_gamete_i <- sample(1:2, nprogeny, replace=TRUE)
  mother_gamete_i <- sample(1:2, nprogeny, prob=mgamete_probs, replace=TRUE)

  # metadata dataframe to keep track of everything
  d <- data.frame(type=ind_type, mother=1L, father=fathers,
                  mgamete=mother_gamete_i, fgamete=father_gamete_i)
  stopifnot(all(fathers[d$type == "selfed"] == 1))
  progeny <- mapply(function(mg, fg, which_father) {
                      father_gamete <- parents[[which_father]][, fg]
                      mother_gamete <- mom[, mg]
                      cbind(mother_gamete, father_gamete)
                    }, mother_gamete_i, father_gamete_i, fathers,
                    SIMPLIFY=FALSE)
  genos <- haplotypes2NoisyGenotypes(progeny, ehet, ehom, na_rate)

  # parent genotypes
  parent_genos <- haplotypes2NoisyGenotypes(parents, ehet, ehom, na_rate)
  list(progeny=genos, metadata=d, freqs=freqs,
       haplos=parents, genos=parent_genos)
}

