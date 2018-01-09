# Function to export to STRUCTURE format from genind object.
# genind objects are created in the R package adegenet.  The function below is an R function.
# Lindsay V. Clark, 26 July 2015

# obj: genind object
# file: file name to write
# pops: whether to include population info in the file
# Function is flexible with regards to ploidy, although genotypes are
# considered to be unambiguous.
# Missing data must be recorded as NA in obj@tab.

# example use: 
# data(nancycats)
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)

genind2structure <- function(obj, file="", pops=FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:nPop(obj)
    names(popnums) <- as.character(unique(pop(obj)))
    popcol <- rep(popnums[as.character(pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=nLoc(obj), dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), dimnames(obj@tab)[[2]])] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}
