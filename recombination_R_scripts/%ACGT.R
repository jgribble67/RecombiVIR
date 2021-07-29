library(dplyr)
library(Biostrings)
##Load in data, save quantification of rows as variable, and slice start and stop sequences
dat <- read.table("/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/NT_frequency/WTP3A_DVGs_sequences.txt", header = FALSE)
n = nrow(dat)
dat_start <- select(dat, V9)
dat_stop <- select(dat, V10)
##generate matrix of sequences (Start sequences)
new <- matrix(nrow = 41, ncol = n)
for(i in 1:41){
  for(j in 1:n){
    new [i,j] <- substring(dat_start[j,], i, i)
  }
}
##generate matrix of sequences (Stop sequences)
new_stop <- matrix(nrow = 41, ncol = n)
for(x in 1:41){
  for(y in 1:n){
    new_stop [x,y] <- substring(dat_stop[y,], x, x)
  }
}
##Count matrix (Start sequences)
countTable_start <- matrix(nrow = 41, ncol = 4)
for(i in 1:41){
  columnSeq_start <- DNAStringSet(paste0(new[i,], collapse = ""))
  columnCounts_start <- letterFrequency(columnSeq_start, letters = "ACGT", OR = 0)
  countTable_start[i,] <- columnCounts_start
}
##Count matrix (Stop sequences)
countTable_stop <- matrix(nrow = 41, ncol = 4)
for(x in 1:41){
  columnSeq_stop <- DNAStringSet(paste0(new_stop[x,], collapse = ""))
  columnCounts_stop <- letterFrequency(columnSeq_stop, letters = "ACGT", OR = 0)
  countTable_stop[x,] <- columnCounts_stop
}
##Rename columns, calculate frequency, and save for start sequences
colnames(countTable_start) <- c("A", "C", "G", "U")
freqTable_start <- countTable_start/n
df1<- round(t(freqTable_start), digit = 4)
df1 <- df1 * 100
df1 <- as.data.frame(t(df1))
##Rename columns, calculate frequency, and save for stop sequences
colnames(countTable_stop) <- c("A", "C", "G", "U")
freqTable_stop <- countTable_stop/n
df2 <- round(t(freqTable_stop), digit = 4)
df2 <- df2 * 100
df2 <- as.data.frame(t(df2))
##Add position lables. +1 indicates junction-participating nucleotide. Positive positions are upstream of the site, negative positions are downstream of the site.
vec_start <- c("+21", "+20", "+19", "+18", "+17", "+16", "+15", "+14", "+13", "+12", "+11", "+10", "+9", "+8", "+7", "+6", "+5", "+4", "+3", "+2", "+1", "-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15", "-16", "-17", "-18", "-19", "-20")
df1$Position <- vec_start
vec_stop <- c("-20", "-19", "-18", "-17", "-16", "-15", "-14", "-13", "-12", "-11", "-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15", "+16", "+17", "+18", "+19", "+20", "+21")
df2$Position <- vec_stop
df1 <- df1[c(5,1,2,3,4)]
df2 <- df2[c(5,1,2,3,4)]
write.table(df1, file = "/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/NT_frequency/WTP3A_start_%ACGU.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df2, file = "/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/NT_frequency/WTP3A_stop_%ACGU.txt", sep = "\t", quote = FALSE, row.names = FALSE)
