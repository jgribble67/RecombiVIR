#!/usr/bin/env Rscript
print("Usage: Junction_Pattern_Plots.R sample_base_name")
args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(ggplot2)
name <- args[1]
data <- read.table(paste("./BED_Files/", name, "_virema_Virus_Recombination_Results.bed", sep = ""), header = FALSE, skip = 1)
data <- data[,c(1,2,3,5)]
data <- data %>% rename(Genome = V1, Start = V2, Stop = V3, Depth = V5)
data <- data[order(data$Depth), ]
#Calculate total recombination depths for Jfreq
junction_depth <- sum(data$Depth)
paste0("The total junction depth is: ", junction_depth)
coverage <- read.table(paste(name, "_virema_coverage.txt", sep=""), header = FALSE)
coverage <- coverage %>% rename(Genome = V1, Position = V2, Depth = V3)
total_depth <- sum(coverage$Depth)
paste0("The total mapped depth is: ", total_depth)
data_forward <- data[which(data$Start < data$Stop), ]
data_forward$Total = sum(data_forward$Depth)
data_forward$Frequency = data_forward$Depth / data_forward$Total
data_forward$logFreq = log10(data_forward$Frequency)
data_graph <- ggplot(data_forward, aes(Stop, Start, colour = logFreq, alpha = logFreq)) + geom_point(size = 2) + theme_linedraw(base_size = 20) + ylim(0,31500) + xlim(0,31500) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = c(0.35,0.80)) + labs(x = "Stop Position (nt)", y = "Start Position (nt)") + guides(alpha = FALSE)
print(data_graph + scale_color_gradientn(guide = guide_colorbar(direction = "horizontal", barwidth = 10), colours = rainbow(6)))
ggsave(filename = paste(name, ".tiff", sep=""), plot = last_plot(), device = "tiff", path = "Junction_Plots/", scale = 1, width = 5, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
write.table(data_forward, file=paste(name, "_forward_junctions.txt", sep = ""), sep = "\t", row.names = FALSE)