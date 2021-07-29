library(Biostrings)
library(dplyr)
SARS2_TRS = "ACGAAC"
SARS2_fasta <- readDNAStringSet("SARS-CoV-2_new.fa", "fasta")
TRS_positions <- vmatchPattern(SARS2_TRS, SARS2_fasta)
TRS_start <- startIndex(TRS_positions)
TRS_stop <- endIndex(TRS_positions)
TRS_start_window <- sapply(TRS_start, function(x) x - 31) #set window for TRS start positions
TRS_stop_window <- sapply(TRS_stop, function(x) x + 31) #set window for TRS stop positions
colnames(TRS_start_window) <- c("Start")
colnames(TRS_stop_window) <- c("Stop")
TRS_matrix <- merge(TRS_start_window, TRS_stop_window, by="row.names")
TRS_matrix <- as.matrix(TRS_matrix[-1])
if(nrow(TRS_matrix) > 9) {
  stop("Too many TRS positions detected in the genome. Virus may be evolving new TRS locations. Script aborting.")
}
if(nrow(TRS_matrix) < 9) {
  stop("Not enough TRS sequences detected in the genome. Virus may be using a mutated version of the consensus sequence. Script aborting.")
}
df <- read.table("test_nosgmRNA2_forward_junctions.txt", header = TRUE)
df_TRSL <- filter(df, between(df$Start, as.numeric(TRS_matrix[1,1]), as.numeric(TRS_matrix[1,2])))
df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[2,1]), as.numeric(TRS_matrix[2,2])))
df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[3,1]), as.numeric(TRS_matrix[3,2])))
df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[4,1]), as.numeric(TRS_matrix[4,2])))
df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[5,1]), as.numeric(TRS_matrix[5,2])))
df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[6,1]), as.numeric(TRS_matrix[6,2])))
df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[7,1]), as.numeric(TRS_matrix[7,2])))
df_sgmRNA8 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[8,1]), as.numeric(TRS_matrix[8,2])))
df_sgmRNA9 <- filter(df_TRSL, between(df_TRSL$Stop, as.numeric(TRS_matrix[9,1]), as.numeric(TRS_matrix[9,2])))
#Create null dataframes if no observations
if(nrow(df_sgmRNA2) == 0) {
  df_sgmRNA2 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
if(nrow(df_sgmRNA2) == 0) {
  df_sgmRNA3 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
if(nrow(df_sgmRNA4) == 0) {
  df_sgmRNA4 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
if(nrow(df_sgmRNA5) == 0) {
  df_sgmRNA5 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
if(nrow(df_sgmRNA6) == 0) {
  df_sgmRNA6 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
if(nrow(df_sgmRNA7) == 0) {
  df_sgmRNA7 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
if(nrow(df_sgmRNA8) == 0) {
  df_sgmRNA8 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
if(nrow(df_sgmRNA9) == 0) {
  df_sgmRNA9 <- data.frame("Genome" = "MT020881.1", "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
}
#Add column identifying sgmRNA species
df_sgmRNA2$Type <- "sgmRNA2"
df_sgmRNA3$Type <- "sgmRNA3"
df_sgmRNA4$Type <- "sgmRNA4"
df_sgmRNA5$Type <- "sgmRNA5"
df_sgmRNA6$Type <- "sgmRNA6"
df_sgmRNA7$Type <- "sgmRNA7"
df_sgmRNA8$Type <- "sgmRNA8"
df_sgmRNA9$Type <- "sgmRNA9"
#Slice canonical sgmRNA species
sgmRNA2_canonical <- df_sgmRNA2 %>% arrange(desc(df_sgmRNA2$Depth)) %>% slice(1)
sgmRNA3_canonical <- df_sgmRNA3 %>% arrange(desc(df_sgmRNA3$Depth)) %>% slice(1)
sgmRNA4_canonical <- df_sgmRNA4 %>% arrange(desc(df_sgmRNA4$Depth)) %>% slice(1)
sgmRNA5_canonical <- df_sgmRNA5 %>% arrange(desc(df_sgmRNA5$Depth)) %>% slice(1)
sgmRNA6_canonical <- df_sgmRNA6 %>% arrange(desc(df_sgmRNA6$Depth)) %>% slice(1)
sgmRNA7_canonical <- df_sgmRNA7 %>% arrange(desc(df_sgmRNA7$Depth)) %>% slice(1)
sgmRNA8_canonical <- df_sgmRNA8 %>% arrange(desc(df_sgmRNA8$Depth)) %>% slice(1)
sgmRNA9_canonical <- df_sgmRNA9 %>% arrange(desc(df_sgmRNA9$Depth)) %>% slice(1)
#Create concatenated dataframe of canonical sgmRNAs
df_canonical <- rbind(sgmRNA2_canonical, sgmRNA3_canonical, sgmRNA4_canonical, sgmRNA5_canonical, sgmRNA6_canonical, sgmRNA7_canonical, sgmRNA8_canonical, sgmRNA9_canonical)
df_canonical$Total <- sum(df_canonical$Depth)
df_sgmRNA <- rbind(df_sgmRNA2, df_sgmRNA3, df_sgmRNA4, df_sgmRNA5, df_sgmRNA6, df_sgmRNA7, df_sgmRNA8, df_sgmRNA9)
df_sgmRNA$Total <- sum(df_sgmRNA$Depth)
#Print list of alternative sgmRNAs
df_alternative <- anti_join(df_sgmRNA, df_canonical, by = "Depth")
df_alternative$Total <- sum(df_alternative$Depth)
df_alt_summary <- df_alternative %>% group_by(Type) %>% summarise(Sum = sum(Depth))
#Slice DVGs and turn format into BED
df_DVG <- anti_join(df, df_sgmRNA, by = c("Start", "Stop"))
df_DVG$Duplication <- "Duplication"
df_DVG$Strand <- "+"
df_DVG$Start2 <- df_DVG$Start
df_DVG$Stop2 <- df_DVG$Stop
df_DVG <- df_DVG[c(1,2,3,8,4,9,10,11)]
#Write tables
write.table(df_canonical, file = "sample_canonical_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_alternative, file = "sample_alternative_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_sgmRNA, file = "sample_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_alt_summary, file = "sample_alt_sgmRNA_summary.txt", sep = "\t", row.names = FALSE)
write.table(df_DVG, file = "sample_DVGs.bed.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)