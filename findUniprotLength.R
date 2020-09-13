#Given a standard .txt file downloaed from the Uniprot website, create a dataframe with the lengths and labels for each Uniprot ID
library(foreign)
library(readr)
library(tidyverse)
library(plyr)
path <- ""
setwd(path)
#import text file
textFileName <- ""
textFile <- read_delim(textFileName, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
#only certain rows from this file will be relevant, these are labeled ID,AC, and FT and can always be found in the first 
#two characters of the line
textFile <- textFile %>% filter(substr(X1,1,2) == "ID"| substr(X1,1,2) == "AC"|substr(X1,1,2) == "FT")

#prepare empty dataframe to be filled
currentID <- ""
length <- ""
df <- data.frame(currentID,length)
#The Uniprot ID is always found in rows labeled ID
for (i in 1:nrow(textFile)) {
  if(substr(textFile$X1[i],1,3) == "ID "){
    if(i != 1){
      df_i <- data.frame(currentID,length)
      df <- rbind(df,df_i) 
      currentID <- ""
      length <- ""
    }
    length_i <- paste(unlist(t(textFile$X1[i])), collapse=",")
    length_i <- substring(length_i,str_locate_all(pattern = " \\d", length_i)[[1]][2],)
    length_i <- substring(length_i,1,nchar(length_i) - 4)
    length <- paste(length,length_i)
  }
  #The major problem to be solved is that one row labled ID, may not have just one ID, there may be multiple IDs for the same protien
  #We need to detect when this is the case (count number of semi colons which seperates the IDs)
  if(substr(textFile$X1[i],1,3) == "AC "){
    numIDs <- str_count(textFile$X1[i],";")
    currentID <- vector(length = numIDs)
    if(numIDs == 1){
      currentID <- substr(textFile$X1[i],6,str_locate(pattern =';',textFile$X1[i])[1] - 1)
    } else {
      currentID[1] <- substr(textFile$X1[i],6,str_locate(pattern =';',textFile$X1[i])[1] - 1)
      for (j in 2:numIDs) {
        currentID[j] <- substr(textFile$X1[i],str_locate_all(pattern = ";",textFile$X1[i])[[1]][j-1] + 2,str_locate_all(pattern = ";",textFile$X1[i])[[1]][j] - 1)
      }
    }
  }
}
df <- df[-1,]
#this final dataframe will have two columns, Uniprot_ID with labels and Length with length information.
df <- df %>% rename("currentID" = Uniprot_ID,"length" = Length)
