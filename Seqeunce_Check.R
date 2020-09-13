#Given a list of Uniprot IDs that need to be checked, a file of modifications for each Uniprot ID,
#and FASTA information for each of these protiens, check to see that all modifications are correct. 
library(Biostrings)
library(BiocGenerics)
library(readr)
path <- ""
setwd(path)
#import FASTA file 
fastaName <- ""
unreviewed1 <- readAAStringSet(fastaName)
#import list of unverified Uniprots, one column with Uniprots
textName <- ""
unverified_uniprotID <- read_csv(textName,col_names = FALSE)
#import a txt file (tab separated) with Uniprot in second column, site (EG: A123) in third column
siteInfoName <- ""
ADP_Uniprot_Site_Separated <- read_delim(siteInfoName, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
#create list with names and seqeunces of proteins, information from FASTA file
UNRev1 <- vector(length = length(unreviewed1))
SEQ1 <- vector(length = length(unreviewed1))
names1 <- names(unreviewed1)
sequence1 <- paste(unreviewed1)
for (i in 1:length(sequence1)) {
  SEQ1[i] = sequence1[i]
  UNRev1[i] = names1[i]
  UNRev1[i] = substring(UNRev1[i],4)
}

#trim uncessary characters from protein name
for (i in 1:length(UNRev1)){
  string = toString(UNRev1[i])
  splitString = strsplit(string,split='\\|')
  UNRev1[i] = splitString[[1]][1]
}
#joining column to create a dataframe with columns Uniprot ID and AA Sequence 
DFUnreviewed <- data.frame(UNRev1, SEQ1)
names(DFUnreviewed)[1] <- "Uniprot"
names(DFUnreviewed)[2] <- "Sequence"
#compare the elements of our FASTA and seqeunce data to the list of Uniprots we hope to confirm
intersection <- intersect(DFUnreviewed$Uniprot,unverified_uniprotID$X1)
#Any differneces mean more FASTA data will need to be downloaded for these proteins
difference <- setdiff(unverified_uniprotID$X1,DFUnreviewed$Uniprot)


#given the list of Uniprots we hope to confirm, get the letter and number of the modification from ADP_Uniprot_Site_Separated and 
#compare this letter and position to the acutal sequence data from DFUnreviewed. Report results as true or false.
correctVector <- vector(length = length(intersection))
for (i in 1:length(intersection)) {
  location1 <- which(grepl(intersection[i], DFUnreviewed$Uniprot))
  location2 <- which(grepl(intersection[i], ADP_Uniprot_Site_Separated$X2))
  sites <- unlist(strsplit(ADP_Uniprot_Site_Separated$X3[location2], ","))
  occurances <- length(sites)
  locationOfAA <- substring(sites,2)
  AminoAcid <- substring(sites,0,1)
  verification <- substring(DFUnreviewed$Sequence[location1],locationOfAA,locationOfAA)
  results <- AminoAcid == verification
  if('FALSE' %in% results){
    correctVector[i] <- FALSE
  } else {
    correctVector[i] <- TRUE
  }
}
#Create a list the returns False if the site infromation is confirmed, if the information was incorrect, print the UNIPROT ID 
#of the protein that needs to be reviewed
errors <- vector(length = length(intersection))
for (i in 1:length(intersection)) {
  if(correctVector[i] == FALSE){
    errors[i] <- intersection[i]
    print(intersection[i])
  }
}



