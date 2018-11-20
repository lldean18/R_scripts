################################################################################
### Modification of the origional CalcD function to save the p_value as well ###
############## Modified from CalcD in the evobiR package #######################
################################################################################

# load dependent packages
library(seqinr)

# write the modified function so that it saves the P-value as well as the d value
CalcD_P_value <- function(alignment = "alignment.fasta", 
                  sig.test="N",                                                    # options are "N", "B", "J"
                  ambig="D", #options are D R I
                  block.size = 1000,                                                # size of blocks to drop in jacknife
                  replicate=1000,
                  align.format='fasta'){
  # this function is used regardless of approach
  d.calc <- function(alignment){
    abba <- 0                                                                         #  set up my variables
    baba <- 0                                                                         #  set up my variables
    for(i in 1:ncol(alignment)){                                               #  run through all sites
      if(length(unique(alignment[, i])) == 2){                                 #  unique(c(p1,p2,p3,o))==2 aka biallelic
        if(alignment[1, i] != alignment[2, i]){                         #  p1 != p2   aka different resolutions in p1 and p2
          if(alignment[4, i] != alignment[3, i]){                       #  o != p3    durand says "less likely pattern due to seq. errors
            if(alignment[3, i] == alignment[1, i]) {baba <- baba + 1}   #  add to the count of baba sites
            if(alignment[2, i] == alignment[3, i]) {abba <- abba + 1}   #  add to the count of abba sites
          } 
        }
      }
    }
    d <- (abba - baba) / (abba + baba)   #what its all about
    results <- list()
    results[[1]] <- d
    results[[2]] <- abba
    results[[3]] <- baba
    return(results)
  }
  
  #### Test of empirical data
  alignment <- read.alignment(alignment, format = align.format, forceToLower=T)                         #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }
  #### This section is being added to deal with reccurent 
  #### Requests to deal with ambiguity in sequence data
  # R A or G
  # Y C or T
  # S G or C
  # W A or T
  # K G or T
  # M A or C
  
  ## First we deal with the situation where the user
  ## wishes to simply drop ambig sites
  if(ambig == "D"){
    target <- c("a","c","g","t")
    keep <- vector()
    for(i in 1:ncol(alignment.matrix)){
      keep[i] <- all(alignment.matrix[,i] %in% target)
    }
    alignment.matrix <- alignment.matrix[,keep]
  }
  
  ## Next we deal with the situation where users want to
  ## randomly resolve ambigous sites
  if(ambig == "R"){
    # I still want to limit sites so we first drop 
    # those sites that look like 3 or 4 possibilities
    target <- c("a", "c", "g", "t", "r",
                "y", "s", "w", "k", "m")
    keep <- vector()
    for(i in 1:ncol(alignment.matrix)){
      keep[i] <- all(alignment.matrix[,i] %in% target)
    }
    alignment.matrix <- alignment.matrix[,keep]
    
    # this function will be applied to each site in our data
    # it resolves ambiguities randomly
    resolver <- function(x){
      if(x=="r") z <- sample(c("a", "g"), 1)
      if(x=="y") z <- sample(c("c", "t"), 1)
      if(x=="s") z <- sample(c("g", "c"), 1)
      if(x=="w") z <- sample(c("a", "t"), 1)
      if(x=="k") z <- sample(c("g", "t"), 1)
      if(x=="m") z <- sample(c("a", "c"), 1)
      if(x %in% c("a", "c", "g", "t")) z <- x
      return(z)
    }
    alignment.matrix <- apply(alignment.matrix, c(1,2), resolver)
  }
  
  
  results <- d.calc(alignment.matrix)
  d <- results[[1]]
  if(is.nan(d)) d <- 0
  abba <- results[[2]]
  baba <- results[[3]]
  
  ## THIS SECTION WILL CALCULATE THE P-VAL BASED ON BOOTSTRAPPING
  ## SITES ARE SAMPLED WITH REPLACEMENT TO MAKE A NEW DATASET OF
  ## OF EQUAL SIZE TO THE ORIGINAL DATASET THIS ALLOWS US TO CALCULATE
  ## THE STANDARD DEVIATION AND THUS A Z SCORE.
  if(sig.test=="B"){
    sim.d<-vector()
    foo <- ncol(alignment.matrix)
    sim.matrix<-matrix(,4,foo)
    cat("\nperforming bootstrap")
    for(k in 1:replicate){
      if(k/(replicate/100) == round(k/(replicate/100))) cat(".")
      sim.matrix[1:4,1:foo] <- alignment.matrix[1:4, sample(1:foo, replace=T)]
      results <- d.calc(sim.matrix)
      sim.d[k] <- results[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(d/sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))
    ## NOW WE MAKE THE OUTPUTS  
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic / Z-score = ", d, " / ", z)
    cat("\n\nResults from ", replicate, "bootstraps")
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ",new.pval,"\n\n") #after Eaton and Ree 2013 
  }
  
  ## THIS SECTION WILL CALCULATE THE P-VAL BASED ON JACKKNIFING
  ## THE DATA IS REANALYZED WHILE DROPPING A PORTION OF THE DATA
  ## THE SIZE OF THE DROPPED PORTION IS DETERMINED BY THE BLOCK SIZE
  ## ARGUMENT THIS PROCEDURE ALLOWS US TO CALCULATE
  ## THE STANDARD DEVIATION AND THUS A Z SCORE. THIS APPROACH IS PARTICULARLY 
  ## IMPORTANT WHEN WE WILL BE USING DATA WHERE THE SNPs MAY BE IN LINKAGE WITH
  ## ONE ANOTHER
  if(sig.test=="J"){
    
    #first lets test whether we are limited in the number of reps by block size
    max.rep <- ncol(alignment.matrix) - block.size
    if(block.size >= (ncol(alignment.matrix)/2)){
      stop(call. = F, paste("\nWith a block size of", block.size, 
                            "and an alignment of", ncol(alignment.matrix), 
                            "sites \nsome sites would never be included in the \nanalysis",
                            "\n\nThe maximum block size is 1/2 the alignment length"))
    }
    if(max.rep < replicate){
      stop(call. = F, paste("\nWith a block size of", block.size, 
                            "and an alignment of", ncol(alignment.matrix), 
                            "sites", replicate, "replicates\nare not possible"))
    }
    
    
    if(max.rep >= replicate){
      drop.pos <- seq.int(from=1, to=(max.rep-1), length.out=replicate)
      replicate2 <- replicate
    }
    sim.d<-vector()
    foo <- ncol(alignment.matrix)
    sim.matrix<-matrix(,4,foo-block.size)
    cat("\nperforming jackknife")
    for(k in 1:replicate2){  
      if(k/2 == round(k/2)) cat(".")
      sim.matrix[1:4,1:(foo-block.size-1)] <-alignment.matrix[1:4, -drop.pos[k]:-(drop.pos[k]+block.size)]
      results <- d.calc(sim.matrix)
      sim.d[k] <- results[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(d/sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))
    ## NOW WE MAKE THE OUTPUTS  
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\nD raw statistic", d)
    cat("\nZ-score = ", z)
    cat("\n\nResults from", replicate2, "jackknifes with block size of", block.size)
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ",new.pval,"\n\n") #after Eaton and Ree 2013 
  }
  if(sig.test=="N"){
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic = ", d,"\n\n")
  }
  return(list(D=d, p.value=new.pval))
}

# test the new function
#TEST <- CalcD_P_value(alignment = "C:/Users/mbzlld/Google Drive/Post Doc/RAD data analysis/ABBA BABA analysis/populations_PGDSpider_OBSE_SCAD_OBSM_BEPA/batch_4_ordered.fasta",
#              sig.test = "J", ambig = "R", block.size = 100, replicate = 100)
#TEST
