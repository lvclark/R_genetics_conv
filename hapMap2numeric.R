# Lindsay V. Clark, 3 August 2015
# R function to import data from TASSEL's HapMap format to numeric (0,1,2) format.
# Example use:
# mydata <- hapMap2numeric("HapMap.hmp.txt")
# "shortnames" indicates whether individual names should be shortened to remove everything
# after the first underscore.

hapMap2numeric <- function(file, shortnames=TRUE){
  hapmap <- as.matrix(read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,-(2:10)])
  samples <- dimnames(hapmap)[[2]][-1]
  loci <- dimnames(hapmap)[[1]]
  
  # shorten sample names if desired
  if(shortnames){
    samples <- sapply(samples, function(x) strsplit(x, split="_")[[1]][1],
                      USE.NAMES=FALSE)
  }
  
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
  
  # matrix to hold output
  x <- matrix(NA, nrow=length(samples), ncol=length(loci),
              dimnames=list(samples, loci))
  # convert genotypes
  for(L in 1:length(loci)){
    thisconv <- conv[[hapmap[L, 1]]]
    x[,L] <- thisconv[hapmap[L, -1]]
  }
  
  return(x)
}
