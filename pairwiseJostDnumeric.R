# Function to calculate pairwise Jost's D between populations using
# numeric SNP data (0,1,2).  Results are returned per locus.

# Lindsay V. Clark, May 9, 2016
# This function uses equation 12 from Jost
# (2008; http://dx.doi.org/10.1111/j.1365-294X.2008.03887.x).  Because we are
# only looking at pairs of populations with this function, n = 2.

# Polyploid samples are allowed, and all samples need not be the same ploidy.
# To account for this, the parts of the equation that use twice the harmonic
# mean of the number of individuals per subpopulation (2Ñ) instead use the
# harmonic mean of the number of gene copies per subpopulation.  For diploids,
# this is the same number.

# This function has been tested on a real dataset, and produces per-locus values
# identical to those produced by mmod's D_Jost.

### Arguments:
# genmat: a matrix of numeric genotypes with individuals in rows and
# loci in columns.  Missing data should be NA.  For diploids, genotypes are 0, 1,
# or 2; for tetraploids, genotypes are 0, 1, 2, 3, or 4; etc.

# pops: a factor, integer, or character vector, with one value for each
# individual, indicating the population to which that individual belongs.

# ploidy: a single value, or one value per individual.

# freq: optional. A matrix of second allele frequencies, with population names
# as the names of the first dimension.  There should be one column per locus.
# If freq = NULL, allele frequencies will be estimated from genmat.  The freq
# argument is handy if allele frequencies have been estimated by polyfreqs or
# some other method, but otherwise is not necessary.

pairwise.JostD.numeric <- function(genmat, pops, ploidy = 2, freq = NULL){
    if(length(pops) != dim(genmat)[1]){
        stop("Number of individuals does not match between pops and genmat.")
    }
    if(min(genmat, na.rm = TRUE) < 0){
        stop("genmat contains negative values")
    }
    if(any(ploidy < 1)){
        stop("ploidy must be greater than zero")
    }
    if(length(ploidy) == 1){
        if(max(genmat, na.rm = TRUE) > ploidy){
            stop("genmat contains values greater than ploidy")
        }
        # expand ploidy
        ploidy <- rep(ploidy, dim(genmat)[1])
    } else {
        if(length(ploidy) != length(pops)){
            stop("Number of individuals does not match between ploidy and pops.")
        }
        if(any(apply(genmat, 1, max, na.rm = TRUE) > ploidy)){
            stop("genmat contains values greater than ploidy")
        }
    }

    # names of all populations
    pops <- as.character(pops)
    allpops <- unique(pops)
    # number of loci
    nLoc <- dim(genmat)[2]
    # number of pops
    nPops <- length(allpops)

    # deal with pre-calculated allele frequencies
    if(is.null(freq)){
        estfreq <- TRUE
    } else {
        estfreq <- FALSE
        if(dim(freq)[2] != nLoc){
            stop("genmat and freq must have the same number of columns.")
        }
        if(dim(freq)[1] != nPops){
            stop("freq needs one row per population.")
        }
        if(!all(allpops %in% dimnames(freq)[[1]])){
            stop("Population names must match between pops and freq.")
        }
    }

    # items that have to be calculated for each pop
    nonmissing <- Hs <- matrix(NA, nrow = nPops, ncol = nLoc,
                                       dimnames = list(allpops,
                                           dimnames(genmat)[[2]]))
    if(estfreq){
        freq <- Hs
    }
    for(p in allpops){
        theseind <- pops == p # individuals in this population
        # number of non-missing alleles (not individuals) per locus
        nonmissing[p,] <- apply(genmat[theseind,], 2,
                                function(x) sum(ploidy[theseind][!is.na(x)]))
        if(estfreq){
            # second allele frequency per locus
            freq[p,] <- colSums(genmat[theseind,], na.rm=TRUE)/nonmissing[p,]
        }
        # expected heterozygosity per locus
        Hs[p,] <- 1 - (freq[p,]^2 + (1 - freq[p,])^2)
    }

    # pairwise calculations
    result <- list()
    for(i in 1:(nPops - 1)){
        pop1 <- allpops[i]
        for(j in (i+1):nPops){
            pop2 <- allpops[j]
            # non-weighted mean allele frequencies across a pair of pops
            meanfreq <- (freq[pop1,] + freq[pop2,])/2
            # harmonic mean of subpopulation sample sizes (in terms of alleles)
            meanNonmissing <- 2/(1/nonmissing[i,] + 1/nonmissing[j,])
            # estimated expected het per locus averaged across two pops
            Hsest <- meanNonmissing/(meanNonmissing - 1) *
                (Hs[i,] + Hs[j,])/2
            # expected het per locus when two pop are combined
            Ht <- 1 - (meanfreq^2 + (1 - meanfreq)^2)
            # estimated expeted het per locus when two pop combined
            Htest <- Ht + Hsest/(2*meanNonmissing)

            # estimate D
            D <- (Htest - Hsest)/(1 - Hsest) * 2

            # add to results
            result[[length(result) + 1]] <- D
            names(result)[length(result)] <- paste(pop1, pop2, sep="_")
        }
    }

    return(result)
}
