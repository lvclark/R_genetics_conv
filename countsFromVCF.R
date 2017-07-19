# Function for getting read counts (allelic depth), as well as genotypes and some marker metadata, from VCF file.
# This function has been tested with VCF files from the TASSEL GBSv2 pipeline.
# Unlike the vcfR package, this function does not read the whole VCF file at once, resulting in considerable
# memory savings and allowing much larger files to be read.

# Lindsay V. Clark, July 19, 2017

# filename = VCF file name.  Can be uncompressed or GZIPed.  If zipped, should end in ".gz".
# nread = number of lines to read at once
# indToKeep and indToExclude = optional, character or integer vector indicating which
#     individuals to keep OR which individuals to exclude
# min.ind.minor.allele = minimum number of individuals with minor allele, if marker is to be retained
# max.proportion.missing = maximum proportion of individuals that can have missing data, if marker is to be retained
# additionalFilter = function that indicates whether to keep a SNP, using two vectors of depths that have
#                    already been subsetted using indToExclude or indToKeep (if provided).
#                    SNPs will be excluded if they don't pass min.ind.minor.allele, max.proportion.missing, OR
#                    additionalFilter.
# snpsPreallocated = the approximate number of SNPs that are expected to pass the filters.  Does not need to be 
#                    accurate, but will save computational time and memory if it is.
countsFromVCF <- function(filename, nread = 1000L, indToKeep = NULL, indToExclude = NULL,
                          min.ind.minor.allele = 2, max.proportion.missing = 1, verbose = TRUE,
                          additionalFilter = function(dep1, dep2){ return(TRUE) },
                          snpsPreallocated = 1e5L){
  if(!is.null(indToKeep) && !is.null(indToExclude)){
    stop("Cannot provide both indToKeep and indToExclude.")
  }
  if(substring(filename, nchar(filename) - 2, nchar(filename)) == ".gz"){
    vcfFile <- gzfile(filename, "rt") # GZIP compressed VCF
  } else {
    vcfFile <- file(filename, "rt") # uncompressed VCF
  }
   
  # get past header
  theselines <- readLines(vcfFile, 1)
  while(substring(theselines[1], 1, 2) == "##"){
    theselines <- readLines(vcfFile, 1)
  }
  stopifnot(substring(theselines[1], 1, 6) == "#CHROM") # expected column header line
  headerrow <- strsplit(theselines[1], split = "\t")[[1]]
  stopifnot(identical(headerrow[1:9], c("#CHROM", 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                                        'FILTER', 'INFO', 'FORMAT')))
  indNames <- headerrow[-(1:9)]
  # convert indToKeep or indToExclude to a numeric vector indicating which columns
  # have individuals to keep.
  if(is.numeric(indToKeep)){
    indToKeep <- indToKeep + 9
    if(max(indToKeep) > length(headerrow)){
      stop("indToKeep includes numbers beyond the number of individuals in the file.")
    }
  }
  if(is.numeric(indToExclude)){
    if(max(indToExclude) + 9 > length(headerrow)){
      stop("indToExclude includes numbers beyond the number of individuals in the file.")
    }
    indToKeep <- (10:length(headerrow))[-indToExclude]
  }
  if(is.character(indToKeep)){
    indToKeep <- match(indToKeep, indNames) + 9
    if(any(is.na(indToKeep))){
      stop("Some individual names from indToKeep not found in file.")
    }
  }
  if(is.character(indToExclude)){
    indToKeep <- (10:length(headerrow))[!indNames %in% indToExclude]
  }
  if(is.null(indToKeep)){
    indToKeep <- 10:length(headerrow)
  }
  # subset indivdual names for output
  if(!is.null(indToKeep)){
    indNames <- indNames[indToKeep - 9]
  }
  nind <- length(indNames) # number of individuals
  maxMissingInd <- floor(nind * max.proportion.missing) # max number ind. missing
  
  nLinesRead <- nread # for checking if we are at end of file
  
  AllName <- character(snpsPreallocated)   # marker names
  AllChr <- character(snpsPreallocated)    # Chromosome names
  AllPos <- integer(snpsPreallocated)      # positions
  AllDep1 <- matrix(NA, nrow = snpsPreallocated, ncol = nind) # depth for allele 1
  AllDep2 <- matrix(NA, nrow = snpsPreallocated, ncol = nind) # depth for allele 2
  AllGen <- matrix(NA, nrow = snpsPreallocated, ncol = nind)
  AllREF <- character(snpsPreallocated)    # reference allele
  AllALT <- character(snpsPreallocated)    # alternative allele
  testcnt <- 0 # for tracking progress through file
  currentSNP <- 0 # the number of the last SNP to be entered
  
  # start reading genotypes
  while(nLinesRead == nread){
    theselines <- strsplit(readLines(vcfFile, nread), split = "\t")
    nLinesRead <- length(theselines) # for stopping loop if at end of file
    TheseName <- character(0)
    TheseChr <- character(0)
    ThesePos <- integer(0)
    TheseDep1 <- matrix(NA, nrow = 0, ncol = nind)
    TheseDep2 <- matrix(NA, nrow = 0, ncol = nind)
    TheseGen <- matrix(NA, nrow = 0, ncol = nind)
    TheseREF <- character(0)
    TheseALT <- character(0)
    for(line in theselines){
      # get position of allele depth and genotype in format
      thisFormat <- strsplit(line[9], split = ":")[[1]]
      ADpos <- match('AD', thisFormat) # position for allelic depth in format
      GTpos <- match('GT', thisFormat) # position for genotype in format
      # extract read depth
      thisdata <- strsplit(line[indToKeep], split = ":") # list, by genotype, with all the data for each genotype
      thisad <- strsplit(sapply(thisdata, function(x) x[ADpos]), ",")
      thisdepth1 <- sapply(thisad, function(x) as.integer(x[1]))
      thisdepth2 <- sapply(thisad, function(x) as.integer(x[2]))
      # skip if too many missing, or minor allele too infrequent, or just one allele, or multiple alternative alleles
      if(!any(is.na(thisdepth2)) && sum(thisdepth1 + thisdepth2 == 0) <= maxMissingInd && 
         sum(thisdepth1 > 0) >= min.ind.minor.allele && sum(thisdepth2 > 0) >= min.ind.minor.allele &&
         additionalFilter(thisdepth1, thisdepth2) && nchar(line[5]) == 1){
        TheseChr <- c(TheseChr, line[1])
        ThesePos <- c(ThesePos, thispos <- as.integer(line[2]))
        TheseName <- c(TheseName, line[3])
        TheseDep1 <- rbind(TheseDep1, thisdepth1)
        TheseDep2 <- rbind(TheseDep2, thisdepth2)
        TheseREF <- c(TheseREF, line[4])
        TheseALT <- c(TheseALT, line[5])
        GTstrings <- sapply(thisdata, function(x) x[GTpos])
        GTnum <- rep(as.integer(NA), nind)
        GTnum[GTstrings == "0/0"] <- 0L
        GTnum[GTstrings %in% c("0/1","1/0")] <- 1L
        GTnum[GTstrings == "1/1"] <- 2L
        TheseGen <- rbind(TheseGen, GTnum)
      }
    }
    nNewSnps <- length(TheseName)
    AllName[currentSNP+(1:nNewSnps)] <- TheseName
    AllChr[currentSNP+(1:nNewSnps)] <- TheseChr
    AllPos[currentSNP+(1:nNewSnps)] <- ThesePos
    AllDep1[currentSNP+(1:nNewSnps),] <- unname(TheseDep1)
    AllDep2[currentSNP+(1:nNewSnps),] <- unname(TheseDep2)
    AllREF[currentSNP+(1:nNewSnps)] <- TheseREF
    AllALT[currentSNP+(1:nNewSnps)] <- TheseALT
    AllGen[currentSNP+(1:nNewSnps),] <- TheseGen
    testcnt <- testcnt + 1
    currentSNP <- currentSNP + nNewSnps
    if(verbose){
      cat(c(paste(testcnt * nread, "lines read"),
            paste(currentSNP, "SNPs retained")), sep = "\n")
    } 
#    if(testcnt == 10) break # for quick tests and debugging of the function
  }
  if(currentSNP < snpsPreallocated){
    AllName <- AllName[1:currentSNP]
    AllChr <- AllChr[1:currentSNP]
    AllPos <- AllPos[1:currentSNP]
    AllDep1 <- AllDep1[1:currentSNP,]
    AllDep2 <- AllDep2[1:currentSNP,]
    AllGen <- AllGen[1:currentSNP,]
    AllREF <- AllREF[1:currentSNP]
    AllALT <- AllALT[1:currentSNP]
  }
  dimnames(AllDep1) <- dimnames(AllDep2) <- dimnames(AllGen) <- list(AllName, indNames)
  close(vcfFile)
  
  return(list(Names = AllName, Chr = AllChr, Pos = AllPos, Ref = AllREF, Alt = AllALT,
              Counts1 = AllDep1, Counts2 = AllDep2, Genotypes = AllGen))
}
