library(dplyr)
data_junctions <- read.table("/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/Junction_Files/SARS2-0125-2_forward_junctions_single.txt", header = TRUE)
data_coverage <- read.table("/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/Junction_Files/SARS2-0125-2_virema.coverage.txt", header = FALSE)
data_coverage <- data_coverage %>% rename(Genome = V1, Position = V2, Coverage = V3)
data_agg <- data_junctions %>% group_by(Position) %>% summarise(Depth = sum(Depth))
data_agg <- data_agg[order(data_agg$Position), ]
data_PRF <- right_join(data_agg, data_coverage, by = "Position")
data_PRF[is.na(data_PRF)] <- 0
data_PRF$Total = data_PRF$Depth + data_PRF$Coverage
data_PRF$Frequency = data_PRF$Depth / data_PRF$Total
data_PRF <- data_PRF[c(3,1,2,4,5,6)]
write.table(data_PRF, file = "/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/PRF/SARS2-0125-2_PRF.txt", sep = "\t", row.names = FALSE)
