# Lindsay V. Clark, 26 July 2015
# R function to import data from TASSEL's HapMap format to adegenet's
# genind format.
# Previous version available at: http://dx.doi.org/10.13012/C5CC0XMJ
# Example use:
# mydata <- hapMap2genind("HapMap.hmp.txt")
hapMap2genind <- function(file){
  require(adegenet)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE, comment.char = "")
  samples <- scan(file, what = character(), nlines = 1)[-(1:11)]
  conv <- c("AA", "CC", "TT", "GG", "AC", "AT", "AG", "CT", "CG", "TG", NA)
  names(conv) <- c("A", "C", "T", "G", "M" , "W", "R", "Y", "S",  "K", "N")
  mydf <- matrix(NA, nrow=dim(hapmap)[1], ncol=length(samples),
                 dimnames=list(row.names(hapmap), samples))
  for(i in 1:length(samples)){
    mydf[,i] <- conv[hapmap[[i+10]]]
  }
  mydf <- as.data.frame(t(mydf))
  x <- df2genind(mydf, type="codom", ncode=1)
  return(x)
}