#!/usr/bin/env Rscript
print("Usage: sgmRNA_filtering_MERS.R input_file output_canonical_sgmRNAs.txt output_alternative_sgmRNAs.txt output_total_sgmRNAs.txt, output_alt_sgmRNAs_summary.txt, output_DVGs.bed.txt")
print("input_file = file path for forward junctions")
print("output_canonical_sgmRNAs.txt = file path for tab-delineated file containing labeled canonical sgmRNA junctions")
print("output_alternative_sgmRNAs.txt = file path for tab-delineated file containing all alternative sgmRNA junctions")
print("output_total_sgmRNAs.txt = file path for tab-delineated file containing all sgmRNA junctions")
print("output_alt_sgmRNA_summary.txt = file path for tab-delineated file containing labeled summary of all alternative sgmRNA species")
print("output_DVGs.bed.txt = file path for BED format file of DVG junctions")
args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
#args[1] is file path and name of forward junctions tab-delineated file
df <- read.table(args[1], header = TRUE)
df_TRSL <- filter(df, between(df$Start, 32, 97))
df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$Stop, 21374, 21439))
df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$Stop, 25490, 25555))
df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$Stop, 25812, 25877))
df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$Stop, 26802, 26867))
df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$Stop, 27552, 27617))
df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$Stop, 27807, 27872))
df_sgmRNA8 <- filter(df_TRSL, between(df_TRSL$Stop, 28514, 28579))
#Add column identifying sgmRNA species
if(nrow(df_sgmRNA2) > 0) {
  df_sgmRNA2$Type <- "sgmRNA2"
}
if(nrow(df_sgmRNA3) > 0) {
  df_sgmRNA3$Type <- "sgmRNA3"
}
if(nrow(df_sgmRNA4) > 0) {
  df_sgmRNA4$Type <- "sgmRNA4"
}
if(nrow(df_sgmRNA5) > 0) {
  df_sgmRNA5$Type <- "sgmRNA5"
}
if(nrow(df_sgmRNA6) > 0) {
  df_sgmRNA6$Type <- "sgmRNA6"
}
if(nrow(df_sgmRNA7) > 0) {
  df_sgmRNA7$Type <- "sgmRNA7"
}
if(nrow(df_sgmRNA8) > 0) {
  df_sgmRNA8$Type <- "sgmRNA8"
}
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
df_DVG <- anti_join(df, df_sgmRNA, by = c("Start", "Stop"))
df_DVG$Duplication <- "Duplication"
df_DVG$Strand <- "+"
df_DVG$Start2 <- df_DVG$Start
df_DVG$Stop2 <- df_DVG$Stop
df_DVG <- df_DVG[c(1,2,3,8,4,9,10,11)]
#Write tables.
write.table(df_canonical, file = args[2], sep = "\t", row.names = FALSE)
write.table(df_alternative, file = args[3], sep = "\t", row.names = FALSE)
write.table(df_sgmRNA, file = args[4], sep = "\t", row.names = FALSE)
write.table(df_alt_summary, file = args[5], sep = "\t", row.names = FALSE)
write.table(df_DVG, file = args[6], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)