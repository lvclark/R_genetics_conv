# Lindsay V. Clark, 18 March 2013
# R function to import data from TASSEL's HapMap format to adegenet's
# genlight format.
# Also available at: http://dx.doi.org/10.13012/C5CC0XMJ
# On my computer (64 bit OS, 8Gb RAM, 2.7 GHz processor) the function takes
# approximately 8e-06 seconds per individual per locus to run.
# (Big datasets could take a few minutes.)
# At this time it does not import chromosome position information, but it could
# be modified to do so fairly easily.

# Example use:
# mydata <- hapMap2genlight("HapMap.hmp.txt")

hapMap2genlight <- function(file){
  require(adegenet)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,-(2:10)]
  samples <- scan(file, what = character(), nlines = 1)[-(1:11)]
  loci <- row.names(hapmap)
  
  # set up conversion table
  s <- as.integer(c(0,1,2,NA))
  ac <- s
  ag <- s
  at <- s
  cg <- s
  ct <- s
  gt <- s
  names(ac) <- c("A","M","C","N")
  names(ag) <- c("A","R","G","N")
  names(at) <- c("A","W","T","N")
  names(cg) <- c("C","S","G","N")
  names(ct) <- c("C","Y","T","N")
  names(gt) <- c("G","K","T","N")
  conv <- list(ac,ac,ag,ag,at,at,cg,cg,ct,ct,gt,gt)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G")
  
  # Pull out and convert genotypes
  S <- length(samples)
  SBlist <- vector(mode="list",S)   # set up list of SNPbin objects
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=TRUE, USE.NAMES=FALSE)
    # create SNPbin object for this individual
    SBlist[[i]] <- new("SNPbin", mygen)
  }
  
  # make genlight object
  x <- new("genlight", SBlist)
  locNames(x) <- loci
  indNames(x) <- samples
  
  return(x)
}
