#!/usr/bin/env Rscript
print("Usage: Recombination_Analysis.R sample_base_name wd_path virus reference_accession reference_file")
library(dplyr)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
wd <- args[2]
if (dir.exists(paste(wd, "Junction_Plots", "/", sep = "/", collapse = "/"))) {
  cat("Junction_Plots exists in current working directory and is a directory", sep="\n")
} else {
  cat("Junction_Plots does not exist in current working directory - creating", sep="\n")
  dir.create(file.path(wd, "Junction_Plots"))
}
if (dir.exists(paste(wd, "Junction_Files", "/", sep = "/", collapse = "/"))) {
  cat("Junction_Files exists in current working directory and is a directory", sep="\n")
} else {
  cat("Junction_Files does not exist in current working directory - creating", sep="\n")
  dir.create(file.path(wd, "Junction_Files"))
}
data <- read.table(paste(wd, "BED_Files/", name, "_virema_Virus_Recombination_Results.bed", sep = ""), header = FALSE, skip = 1)
data <- data[,c(1,2,3,5)]
data <- data %>% rename(Genome = V1, Start = V2, Stop = V3, Depth = V5)
data <- data[order(data$Depth), ]
#Calculate total recombination depths for Jfreq
junction_depth <- sum(data$Depth)
cat(paste0("The total junction depth is: ", junction_depth, "\n"))
coverage <- read.table(paste(wd, name, "_virema_coverage.txt", sep=""), header = FALSE)
coverage <- coverage %>% rename(Genome = V1, Position = V2, Coverage = V3)
total_depth <- sum(coverage$Coverage)
cat(paste0("The total mapped depth is: ", total_depth, "\n"))
mean_depth <- mean(coverage$Coverage)
cat(paste0("The mean depth is: ", mean_depth, "\n"))
data_forward <- data[which(data$Start < data$Stop), ]
data_forward$Total = sum(data_forward$Depth)
data_forward$Frequency = data_forward$Depth / data_forward$Total
data_forward$logFreq = log10(data_forward$Frequency)
data_graph <- ggplot(data_forward, aes(Stop, Start, colour = logFreq, alpha = logFreq)) + geom_point(size = 2) + theme_linedraw(base_size = 20) + ylim(0,31500) + xlim(0,31500) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = c(0.35,0.80)) + labs(x = "Stop Position (nt)", y = "Start Position (nt)") + guides(alpha = FALSE)
print(data_graph + scale_color_gradientn(guide = guide_colorbar(direction = "horizontal", barwidth = 10), colours = rainbow(6)))
ggsave(filename = paste(name, ".tiff", sep=""), plot = last_plot(), device = "tiff", path = paste(wd, "Junction_Plots/", sep=""), scale = 1, width = 5, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
write.table(data_forward, file=paste(wd, "Junction_Files/", name, "_forward_junctions.txt", sep = ""), sep = "\t", row.names = FALSE)
#Calculate recombination frequency at each genomic position. (Positional Recombination Frequency)
library(tidyr)
options(dplyr.summarise.inform = FALSE)
df_single <- pivot_longer(data_forward, c(2,3), values_to = "Position") 
df_single <- df_single[,c(1,7,2)]
df_agg <- df_single %>% group_by(Position) %>% summarise(Depth = sum(Depth))
df_agg <- df_agg[order(df_agg$Position), ]
df_PRF <- right_join(df_agg, coverage, by = "Position")
df_PRF[is.na(df_PRF)] <- 0
df_PRF$Total = df_PRF$Depth + df_PRF$Coverage
df_PRF$Frequency = df_PRF$Depth / df_PRF$Total
df_PRF <- df_PRF[c(3,1,2,4,5,6)]
write.table(df_PRF, file = paste(wd, "/Junction_Files/", name, "_PRF.txt", sep=""), sep = "\t", row.names = FALSE)
#sgmRNA filtering
if (dir.exists(paste(wd, "sgmRNAs_DVGs", "/", sep = "/", collapse = "/"))) {
  cat("sgmRNAs_DVGs exists in current working directory and is a directory", sep="\n")
} else {
  cat("sgmRNAs_DVGs does not exist in current working directory - creating", sep="\n")
  dir.create(file.path(wd, "sgmRNAs_DVGs"))
}
##Okay so I was thinking I will add another argument delineating which virus is in it and pull the accession number 
if (args[3] == "MHV") {
  df_TRSL <- filter(data_forward, between(data_forward$Start, 31, 103))
  df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$Stop, 21713, 21785))
  df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$Stop, 23888, 23960))
  df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$Stop, 27901, 27973))
  df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$Stop, 28284, 28356))
  df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$Stop, 28924, 28996))
  df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$Stop, 29621, 29693))
  #Create null dataframes if no observations
  if(nrow(df_sgmRNA2) == 0) {
    df_sgmRNA2 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA2) == 0) {
    df_sgmRNA3 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA4) == 0) {
    df_sgmRNA4 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA5) == 0) {
    df_sgmRNA5 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA6) == 0) {
    df_sgmRNA6 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA7) == 0) {
    df_sgmRNA7 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  #Add column identifying sgmRNA species
  df_sgmRNA2$Type <- "sgmRNA2"
  df_sgmRNA3$Type <- "sgmRNA3"
  df_sgmRNA4$Type <- "sgmRNA4"
  df_sgmRNA5$Type <- "sgmRNA5"
  df_sgmRNA6$Type <- "sgmRNA6"
  df_sgmRNA7$Type <- "sgmRNA7"
  #Slice canonical sgmRNA species
  sgmRNA2_canonical <- df_sgmRNA2 %>% arrange(desc(df_sgmRNA2$Depth)) %>% slice(1)
  sgmRNA3_canonical <- df_sgmRNA3 %>% arrange(desc(df_sgmRNA3$Depth)) %>% slice(1)
  sgmRNA4_canonical <- df_sgmRNA4 %>% arrange(desc(df_sgmRNA4$Depth)) %>% slice(1)
  sgmRNA5_canonical <- df_sgmRNA5 %>% arrange(desc(df_sgmRNA5$Depth)) %>% slice(1)
  sgmRNA6_canonical <- df_sgmRNA6 %>% arrange(desc(df_sgmRNA6$Depth)) %>% slice(1)
  sgmRNA7_canonical <- df_sgmRNA7 %>% arrange(desc(df_sgmRNA7$Depth)) %>% slice(1)
  #Create concatenated dataframe of canonical sgmRNAs
  df_canonical <- rbind(sgmRNA2_canonical, sgmRNA3_canonical, sgmRNA4_canonical, sgmRNA5_canonical, sgmRNA6_canonical, sgmRNA7_canonical)
  df_canonical$Total <- sum(df_canonical$Depth)
  df_sgmRNA <- rbind(df_sgmRNA2, df_sgmRNA3, df_sgmRNA4, df_sgmRNA5, df_sgmRNA6, df_sgmRNA7)
  df_sgmRNA$Total <- sum(df_sgmRNA$Depth)
  #Print list of alternative sgmRNAs
  df_alternative <- anti_join(df_sgmRNA, df_canonical, by = "Depth")
  df_alternative$Total <- sum(df_alternative$Depth)
  df_alt_summary <- df_alternative %>% group_by(Type) %>% summarise(Sum = sum(Depth))
  #Slice DVGs and turn format into BED
  df_DVG <- anti_join(data_forward, df_sgmRNA, by = c("Start", "Stop"))
  df_DVG$Duplication <- "Duplication"
  df_DVG$Strand <- "+"
  df_DVG$Start2 <- df_DVG$Start
  df_DVG$Stop2 <- df_DVG$Stop
  df_DVG <- df_DVG[c(1,2,3,8,4,9,10,11)]
  #Write tables
  write.table(df_canonical, file = paste(wd, "sgmRNAs_DVGs/", name, "_canonical_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_alternative, file = paste(wd, "sgmRNAs_DVGs/", name, "_alternative_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_sgmRNA, file = paste(wd, "sgmRNAs_DVGs/", name, "_total_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_alt_summary, file = paste(wd, "sgmRNAs_DVGs/", name, "_alt_sgmRNA_summary.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_DVG, file = paste(wd, "sgmRNAs_DVGs/", name, "_DVGs.bed.txt", sep=""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
} else if (args[3] == "MERS") {
  df_TRSL <- filter(data_forward, between(data_forward$Start, 32, 97))
  df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$Stop, 21374, 21439))
  df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$Stop, 25490, 25555))
  df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$Stop, 25812, 25877))
  df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$Stop, 26802, 26867))
  df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$Stop, 27552, 27617))
  df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$Stop, 27807, 27872))
  df_sgmRNA8 <- filter(df_TRSL, between(df_TRSL$Stop, 28514, 28579))
  #Create null dataframes if no observations
  if(nrow(df_sgmRNA2) == 0) {
    df_sgmRNA2 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA2) == 0) {
    df_sgmRNA3 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA4) == 0) {
    df_sgmRNA4 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA5) == 0) {
    df_sgmRNA5 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA6) == 0) {
    df_sgmRNA6 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA7) == 0) {
    df_sgmRNA7 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA8) == 0) {
    df_sgmRNA8 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  #Add column identifying sgmRNA species
  df_sgmRNA2$Type <- "sgmRNA2"
  df_sgmRNA3$Type <- "sgmRNA3"
  df_sgmRNA4$Type <- "sgmRNA4"
  df_sgmRNA5$Type <- "sgmRNA5"
  df_sgmRNA6$Type <- "sgmRNA6"
  df_sgmRNA7$Type <- "sgmRNA7"
  df_sgmRNA8$Type <- "sgmRNA8"
  #Slice canonical sgmRNA species
  sgmRNA2_canonical <- df_sgmRNA2 %>% arrange(desc(df_sgmRNA2$Depth)) %>% slice(1)
  sgmRNA3_canonical <- df_sgmRNA3 %>% arrange(desc(df_sgmRNA3$Depth)) %>% slice(1)
  sgmRNA4_canonical <- df_sgmRNA4 %>% arrange(desc(df_sgmRNA4$Depth)) %>% slice(1)
  sgmRNA5_canonical <- df_sgmRNA5 %>% arrange(desc(df_sgmRNA5$Depth)) %>% slice(1)
  sgmRNA6_canonical <- df_sgmRNA6 %>% arrange(desc(df_sgmRNA6$Depth)) %>% slice(1)
  sgmRNA7_canonical <- df_sgmRNA7 %>% arrange(desc(df_sgmRNA7$Depth)) %>% slice(1)
  sgmRNA8_canonical <- df_sgmRNA8 %>% arrange(desc(df_sgmRNA8$Depth)) %>% slice(1)
  #Create concatenated dataframe of canonical sgmRNAs
  df_canonical <- rbind(sgmRNA2_canonical, sgmRNA3_canonical, sgmRNA4_canonical, sgmRNA5_canonical, sgmRNA6_canonical, sgmRNA7_canonical, sgmRNA8_canonical)
  df_canonical$Total <- sum(df_canonical$Depth)
  df_sgmRNA <- rbind(df_sgmRNA2, df_sgmRNA3, df_sgmRNA4, df_sgmRNA5, df_sgmRNA6, df_sgmRNA7, df_sgmRNA8)
  df_sgmRNA$Total <- sum(df_sgmRNA$Depth)
  #Print list of alternative sgmRNAs
  df_alternative <- anti_join(df_sgmRNA, df_canonical, by = "Depth")
  df_alternative$Total <- sum(df_alternative$Depth)
  df_alt_summary <- df_alternative %>% group_by(Type) %>% summarise(Sum = sum(Depth))
  #Slice DVGs and turn format into BED
  df_DVG <- anti_join(data_forward, df_sgmRNA, by = c("Start", "Stop"))
  df_DVG$Duplication <- "Duplication"
  df_DVG$Strand <- "+"
  df_DVG$Start2 <- df_DVG$Start
  df_DVG$Stop2 <- df_DVG$Stop
  df_DVG <- df_DVG[c(1,2,3,8,4,9,10,11)]
  #Write tables
  write.table(df_canonical, file = paste(wd, "sgmRNAs_DVGs/", name, "_canonical_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_alternative, file = paste(wd, "sgmRNAs_DVGs/", name, "_alternative_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_sgmRNA, file = paste(wd, "sgmRNAs_DVGs/", name, "_total_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_alt_summary, file = paste(wd, "sgmRNAs_DVGs/", name, "_alt_sgmRNA_summary.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_DVG, file = paste(wd, "sgmRNAs_DVGs/", name, "_DVGs.bed.txt", sep=""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
} else if (args[3] == "SARS2") {
  library(Biostrings)
  SARS2_TRS = "ACGAAC"
  SARS2_fasta <- readDNAStringSet(args[5], "fasta")
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
  df_TRSL <- filter(data_forward, between(data_forward$Start, as.numeric(TRS_matrix[1,1]), as.numeric(TRS_matrix[1,2])))
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
    df_sgmRNA2 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA2) == 0) {
    df_sgmRNA3 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA4) == 0) {
    df_sgmRNA4 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA5) == 0) {
    df_sgmRNA5 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA6) == 0) {
    df_sgmRNA6 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA7) == 0) {
    df_sgmRNA7 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA8) == 0) {
    df_sgmRNA8 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
  }
  if(nrow(df_sgmRNA9) == 0) {
    df_sgmRNA9 <- data.frame("Genome" = args[4], "Start" = 0, "Stop" = 0, "Depth" = 0, "Total" = 0, "Frequency" = 0, "logFreq" = 0)
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
  df_DVG <- anti_join(data_forward, df_sgmRNA, by = c("Start", "Stop"))
  df_DVG$Duplication <- "Duplication"
  df_DVG$Strand <- "+"
  df_DVG$Start2 <- df_DVG$Start
  df_DVG$Stop2 <- df_DVG$Stop
  df_DVG <- df_DVG[c(1,2,3,8,4,9,10,11)]
  #Write tables
  write.table(df_canonical, file = paste(wd, "sgmRNAs_DVGs/", name, "_canonical_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_alternative, file = paste(wd, "sgmRNAs_DVGs/", name, "_alternative_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_sgmRNA, file = paste(wd, "sgmRNAs_DVGs/", name, "_total_sgmRNAs.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_alt_summary, file = paste(wd, "sgmRNAs_DVGs/", name, "_alt_sgmRNA_summary.txt", sep=""), sep = "\t", row.names = FALSE)
  write.table(df_DVG, file = paste(wd, "sgmRNAs_DVGs/", name, "_DVGs.bed.txt", sep=""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
} else {
  stop("Unknown genome! Check reference accession number and change to AY910861.1 (MHV), JX869059.2 (MERS), or SARS-CoV-2.")
}