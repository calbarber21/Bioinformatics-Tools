#This code takes input from a .txt file downloaded from Uniprot. For each protein in the file domain information is extracted for
# the follwing domain features: METAL,ZN_FING,COILED,ACT_SITE,REGION,MOTIF,DOMAIN, and length of the protien. 
library(foreign)
library(readr)
library(tidyverse)
library(plyr)
#define path
path <- ""
setwd(path)
#define and import Uniprot .txt file 
fileName <- ""
X1_69_alternative <- read_delim(fileName, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

#Only certain rows are relevant, take only those rows
X1_69_alternative <- X1_69_alternative %>% filter(substr(X1,1,2) == "ID"| substr(X1,1,2) == "AC"|substr(X1,1,2) == "FT")
#Define and empty data frame that we will fill as we loop through the text file
currentID <- ""
METAL <- ""
ZN_FING <- ""
COILED <- ""
ACT_SITE <- ""
REGION <- ""
MOTIF <- ""
DOMAIN <- ""
length <- ""
df <- data.frame(currentID,METAL,ZN_FING,COILED,ACT_SITE,REGION,MOTIF,DOMAIN,length)
#Loop through the text file and extract relevatn features from each row
for (i in 1:nrow(X1_69_alternative)) {
  if(substr(X1_69_alternative$X1[i],1,3) == "ID "){
    if(i != 1){
      #After each loop we want to add this information to the dataframe
      METAL <- substring(METAL,3,)
      ZN_FING <- substring(ZN_FING,3,)
      COILED <- substring(COILED,3,)
      ACT_SITE <- substring(ACT_SITE,3,)
      REGION <- substring(REGION,3,)
      MOTIF <- substring(MOTIF,3,)
      DOMAIN <- substring(DOMAIN,3,)
      df_i <- data.frame(currentID,METAL,ZN_FING,COILED,ACT_SITE,REGION,MOTIF,DOMAIN,length)
      df <- rbind(df,df_i) 
      #now clear the varaibles
      currentID <- ""
      METAL <- ""
      ZN_FING <- ""
      COILED <- ""
      ACT_SITE <- ""
      REGION <- ""
      MOTIF <- ""
      DOMAIN <- ""
      length <- ""
    }
    length_i <- paste(unlist(t(X1_69_alternative$X1[i])), collapse=",")
    length_i <- substring(length_i,str_locate_all(pattern = " \\d", length_i)[[1]][2],)
    length_i <- substring(length_i,1,nchar(length_i) - 4)
    length <- paste(length,length_i)
  }
  #The major problem to be solved is that one row labled ID, may not have just one ID, there may be multiple IDs for the same protien
  #We need to detect when this is the case (count number of semi colons which seperates the IDs)
  if(substr(X1_69_alternative$X1[i],1,3) == "AC "){
  numIDs <- str_count(X1_69_alternative$X1[i],";")
  currentID <- vector(length = numIDs)
  if(numIDs == 1){
    currentID <- substr(X1_69_alternative$X1[i],6,str_locate(pattern =';',X1_69_alternative$X1[i])[1] - 1)
  } else {
    currentID[1] <- substr(X1_69_alternative$X1[i],6,str_locate(pattern =';',X1_69_alternative$X1[i])[1] - 1)
  for (j in 2:numIDs) {
    currentID[j] <- substr(X1_69_alternative$X1[i],str_locate_all(pattern = ";",X1_69_alternative$X1[i])[[1]][j-1] + 2,str_locate_all(pattern = ";",X1_69_alternative$X1[i])[[1]][j] - 1)
  }
  }
  }
  #Below there exists a block of code much like this for each type of domain we want to extract. Each domain feature has a unique 
  #pattern and the requisite regex code exists in each block below. Generally we want to extract the type of domain, starting 
  #and endingAA position, and potentially other information for the features DOMAIN, MOTIF, and REGION 
  #Common features of data extraction will be explained in this first block of code for METAL.
  if(substr(X1_69_alternative$X1[i],6,10) == "METAL"){
    #If a protein has a METAL domain information will be in the row that included "METAL" and the next row
    METAL_i <- paste(unlist(t(X1_69_alternative$X1[i:(i+2)])), collapse=",")
    #We dont want the information that matched this pattern
    METAL_i <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", METAL_i, perl=TRUE)
    #The first 4 characters will be labeled "FT  " we do not want these characters
    METAL_i <- substring(METAL_i,4,)
    METAL_i <- paste(substr(METAL_i,1,str_locate(pattern = "FT",METAL_i)[1] - 2),substring(METAL_i,str_locate(pattern = "note=",METAL_i)[2] + 2,str_locate_all(pattern = "FT",METAL_i)[[1]][4] - 4))
    #Here we extract the starting and ending AA positions in the case that this is in fact a range (say "123..130")
    if(!grepl("\\d \\..\\d",METAL_i)){
      digits <- nrow(str_locate_all(pattern = "\\d", METAL_i)[[1]])
      METAL_i <- paste(substr(METAL_i,1,str_locate_all(pattern = "\\d", METAL_i)[[1]][digits*2])," ",substr(METAL_i,str_locate_all(pattern = "\\d", METAL_i)[[1]][digits+1],str_locate_all(pattern = "\\d", METAL_i)[[1]][digits*2])," ",substring(METAL_i,str_locate_all(pattern = "\\d", METAL_i)[[1]][digits*2] + 2,),".", sep = "")
      #Extract location information in the case that the domain exists at a single amino acid.
      } else {
      METAL_i <- gsub("\\.."," ",METAL_i)
      METAL_i <- paste(METAL_i,".",sep = "")
      }
    #Add this information to the variable METAL outsie this loop, if there are multiple metal features then they will be 
    #separated by a semicolon
    METAL <- paste(METAL,METAL_i,sep = "; ")
  }
  #Zinc Finger
  if(substr(X1_69_alternative$X1[i],6,12) == "ZN_FING"){
    ZN_FING_i <- paste(unlist(t(X1_69_alternative$X1[i:(i+2)])), collapse=",")
    ZN_FING_i <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", ZN_FING_i, perl=TRUE)
    ZN_FING_i <- substring(ZN_FING_i,4,)
    ZN_FING_i <- paste(substr(ZN_FING_i,1,str_locate(pattern = "FT",ZN_FING_i)[1] - 2),substring(ZN_FING_i,str_locate(pattern = "note=",ZN_FING_i)[2] + 2,str_locate_all(pattern = "FT",ZN_FING_i)[[1]][4] - 4))
    if(!grepl("\\d\\..\\d",ZN_FING_i)){
      digits <- nrow(str_locate_all(pattern = "\\d", ZN_FING_i)[[1]])
      ZN_FING_i <- paste(substr(ZN_FING_i,1,str_locate_all(pattern = "\\d", ZN_FING_i)[[1]][digits*2])," ",substr(ZN_FING_i,str_locate_all(pattern = "\\d", ZN_FING_i)[[1]][digits+1],str_locate_all(pattern = "\\d", ZN_FING_i)[[1]][digits*2])," ",substring(ZN_FING_i,str_locate_all(pattern = "\\d", ZN_FING_i)[[1]][digits*2] + 2,),".", sep = "")
    } else {
      ZN_FING_i <- gsub("\\.."," ",ZN_FING_i)
      ZN_FING_i <- paste(ZN_FING_i,".",sep = "")
    }
    ZN_FING <- paste(ZN_FING,ZN_FING_i,sep = "; ")
  }
  #Coiled
  if(substr(X1_69_alternative$X1[i],6,11) == "COILED"){
    COILED_i <- paste(unlist(t(X1_69_alternative$X1[i:(i+1)])), collapse=",")
    COILED_i <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", COILED_i, perl=TRUE)
    COILED_i <- substring(COILED_i,4,)
    COILED_i <- substring(COILED_i,1,str_locate(pattern = "FT",COILED_i)[1] - 2)
    if(!grepl("\\d\\..\\d",COILED_i)){
      digits <- nrow(str_locate_all(pattern = "\\d", COILED_i)[[1]])
      COILED_i <- paste(substr(COILED_i,1,str_locate_all(pattern = "\\d", COILED_i)[[1]][digits*2])," ",substr(COILED_i,str_locate_all(pattern = "\\d", COILED_i)[[1]][digits+1],str_locate_all(pattern = "\\d", COILED_i)[[1]][digits*2]),".", sep = "")
    } else {
      COILED_i <- gsub("\\.."," ",COILED_i)
      COILED_i <- paste(COILED_i,".",sep = "")
    }
    COILED <- paste(COILED,COILED_i,sep = "; ")
  }
  #Active Site
  if(substr(X1_69_alternative$X1[i],6,13) == "ACT_SITE"){
    ACT_SITE_i <- paste(unlist(t(X1_69_alternative$X1[i:(i+1)])), collapse=",")
    ACT_SITE_i <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", ACT_SITE_i, perl=TRUE)
    ACT_SITE_i <- substring(ACT_SITE_i,4,)
    ACT_SITE_i <- paste(substr(ACT_SITE_i,1,str_locate(pattern = "FT",ACT_SITE_i)[1] - 2),substring(ACT_SITE_i,str_locate(pattern = "note=",ACT_SITE_i)[2] + 2,str_locate_all(pattern = "\"",ACT_SITE_i)[[1]][4] - 1))
    if(!grepl("\\d\\..\\d",ACT_SITE_i)){
      digits <- nrow(str_locate_all(pattern = "\\d", ACT_SITE_i)[[1]])
      ACT_SITE_i <- paste(substr(ACT_SITE_i,1,str_locate_all(pattern = "\\d", ACT_SITE_i)[[1]][digits*2])," ",substr(ACT_SITE_i,str_locate_all(pattern = "\\d", ACT_SITE_i)[[1]][digits+1],str_locate_all(pattern = "\\d", ACT_SITE_i)[[1]][digits*2])," ",substring(ACT_SITE_i,str_locate_all(pattern = "\\d", ACT_SITE_i)[[1]][digits*2] + 2,),".", sep = "")
    } else {
      ACT_SITE_i <- gsub("\\.."," ",ACT_SITE_i)
      ACT_SITE_i <- paste(ACT_SITE_i,".",sep = "")
    }
    ACT_SITE <- paste(ACT_SITE,ACT_SITE_i,sep = "; ")
  }
  #Region
  if(substr(X1_69_alternative$X1[i],6,11) == "REGION"){
    REGION_i <- paste(unlist(t(X1_69_alternative$X1[i:(i+1)])), collapse=",")
    REGION_i <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", REGION_i, perl=TRUE)
    REGION_i <- substring(REGION_i,4,)
    REGION_i <- paste(substr(REGION_i,1,str_locate(pattern = "FT",REGION_i)[1] - 2),substring(REGION_i,str_locate(pattern = "note=",REGION_i)[2] + 2,str_locate_all(pattern = "\"",REGION_i)[[1]][4] - 1))
    if(!grepl("\\d\\..\\d",REGION_i)){
      digits <- nrow(str_locate_all(pattern = "\\d", REGION_i)[[1]])
      REGION_i <- paste(substr(REGION_i,1,str_locate_all(pattern = "\\d", REGION_i)[[1]][digits*2])," ",substr(REGION_i,str_locate_all(pattern = "\\d", REGION_i)[[1]][digits+1],str_locate_all(pattern = "\\d", REGION_i)[[1]][digits*2])," ",substring(ACT_SITE_i,str_locate_all(pattern = "\\d", ACT_SITE_i)[[1]][digits*2] + 2,),".", sep = "")
    } else {
      REGION_i <- gsub("\\.."," ",REGION_i)
      REGION_i <- paste(REGION_i,".",sep = "")
    }
    REGION <- paste(REGION,REGION_i,sep = "; ")
  }
  #Motif
  if(substr(X1_69_alternative$X1[i],6,10) == "MOTIF"){
    MOTIF_i <- paste(unlist(t(X1_69_alternative$X1[i:(i+1)])), collapse=",")
    MOTIF_i <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", MOTIF_i, perl=TRUE)
    MOTIF_i <- substring(MOTIF_i,4,)
    MOTIF_i <- paste(substr(MOTIF_i,1,str_locate(pattern = "FT",MOTIF_i)[1] - 2),substring(MOTIF_i,str_locate(pattern = "note=",MOTIF_i)[2] + 2,str_locate_all(pattern = "\"",MOTIF_i)[[1]][4] - 1))
    if(!grepl("\\d\\..\\d",MOTIF_i)){
      digits <- nrow(str_locate_all(pattern = "\\d", MOTIF_i)[[1]])
      MOTIF_i <- paste(substr(MOTIF_i,1,str_locate_all(pattern = "\\d", MOTIF_i)[[1]][digits*2])," ",substr(MOTIF_i,str_locate_all(pattern = "\\d", MOTIF_i)[[1]][digits+1],str_locate_all(pattern = "\\d", MOTIF_i)[[1]][digits*2])," ",substring(MOTIF_i,str_locate_all(pattern = "\\d", MOTIF_i)[[1]][digits*2] + 2,),".", sep = "")
    } else {
      MOTIF_i <- gsub("\\.."," ",MOTIF_i)
      MOTIF_i <- paste(MOTIF_i,".",sep = "")
    }
    MOTIF <- paste(MOTIF,MOTIF_i,sep = "; ")
    
  }
  #Domain
  if(substr(X1_69_alternative$X1[i],6,11) == "DOMAIN"){
    DOMAIN_i <- paste(unlist(t(X1_69_alternative$X1[i:(i+1)])), collapse=",")
    DOMAIN_i <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", DOMAIN_i, perl=TRUE)
    DOMAIN_i <- substring(DOMAIN_i,4,)
    DOMAIN_i <- paste(substr(DOMAIN_i,1,str_locate(pattern = "FT",DOMAIN_i)[1] - 2),substring(DOMAIN_i,str_locate(pattern = "note=",DOMAIN_i)[2] + 2,str_locate_all(pattern = "\"",DOMAIN_i)[[1]][4] - 1))
    if(!grepl("\\d\\..\\d",DOMAIN_i)){
      digits <- nrow(str_locate_all(pattern = "\\d", DOMAIN_i)[[1]])
      DOMAIN_i <- paste(substr(DOMAIN_i,1,str_locate_all(pattern = "\\d", DOMAIN_i)[[1]][digits*2])," ",substr(DOMAIN_i,str_locate_all(pattern = "\\d", DOMAIN_i)[[1]][digits+1],str_locate_all(pattern = "\\d", DOMAIN_i)[[1]][digits*2])," ",substring(DOMAIN_i,str_locate_all(pattern = "\\d", DOMAIN_i)[[1]][digits*2] + 2,),".", sep = "")
    } else {
      DOMAIN_i <- gsub("\\.."," ",DOMAIN_i)
      DOMAIN_i <- paste(DOMAIN_i,".",sep = "")
    }
    DOMAIN <- paste(DOMAIN,DOMAIN_i,sep = "; ")
  }
}
#Relabel columns and remove the first row which will always be empty
df <- df %>% rename("Metal.binding" = METAL, "Zinc.finger" = ZN_FING,"Coiled.coil" = COILED,"Active.site" = ACT_SITE,"Motif" = MOTIF,"Region" = REGION,"Domain..FT." = DOMAIN,"Length" = length,"Entry" = currentID) 
df <- df[-1,]
