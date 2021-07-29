library(dplyr)
df <- read.table("/Users/gribblj/Dropbox/NHC_recombination/MERS_HAE/Recombination/Junction_Files/UNTA_forward_junctions.txt", header = TRUE)
df_TRSL <- filter(df, between(df$start, 32, 97))
df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$stop, 21374, 21439))
df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$stop, 25490, 25555))
df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$stop, 25812, 25877))
df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$stop, 26802, 26867))
df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$stop, 27552, 27617))
df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$stop, 27807, 27872))
df_sgmRNA8 <- filter(df_TRSL, between(df_TRSL$stop, 28514, 28579))
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
sgmRNA2_canonical <- df_sgmRNA2 %>% arrange(desc(df_sgmRNA2$depth)) %>% slice(1)
sgmRNA3_canonical <- df_sgmRNA3 %>% arrange(desc(df_sgmRNA3$depth)) %>% slice(1)
sgmRNA4_canonical <- df_sgmRNA4 %>% arrange(desc(df_sgmRNA4$depth)) %>% slice(1)
sgmRNA5_canonical <- df_sgmRNA5 %>% arrange(desc(df_sgmRNA5$depth)) %>% slice(1)
sgmRNA6_canonical <- df_sgmRNA6 %>% arrange(desc(df_sgmRNA6$depth)) %>% slice(1)
sgmRNA7_canonical <- df_sgmRNA7 %>% arrange(desc(df_sgmRNA7$depth)) %>% slice(1)
sgmRNA8_canonical <- df_sgmRNA8 %>% arrange(desc(df_sgmRNA8$depth)) %>% slice(1)
#Create concatenated dataframe of canonical sgmRNAs
df_canonical <- rbind(sgmRNA2_canonical, sgmRNA3_canonical, sgmRNA4_canonical, sgmRNA5_canonical, sgmRNA6_canonical, sgmRNA7_canonical, sgmRNA8_canonical)
df_canonical$Total <- sum(df_canonical$depth)
major_sgmRNA_depth = sum(df_canonical$depth)
df_sgmRNA <- rbind(df_sgmRNA2, df_sgmRNA3, df_sgmRNA4, df_sgmRNA5, df_sgmRNA6, df_sgmRNA7, df_sgmRNA8)
df_sgmRNA$Total <- sum(df_sgmRNA$depth)
total_sgmRNA_depth = sum(df_sgmRNA$depth)
#Print list of alternative sgmRNAs
df_alternative <- anti_join(df_sgmRNA, df_canonical, by = (c("start", "stop")))
df_alternative$Total <- sum(df_alternative$D=depth)
df_alt_summary <- df_alternative %>% group_by(Type) %>% summarise(Sum = sum(depth))
minor_sgmRNA_depth = sum(df_alt_summary$Sum)
#Slice DVGs and turn format into BED
df_DVG <- anti_join(df, df_sgmRNA, by = c("start", "stop"))
DVG_depth = sum(df_DVG$depth)
#Write tables.
#write.table(df_canonical, file = args[2], sep = "\t", row.names = FALSE)
#write.table(df_alternative, file = args[3], sep = "\t", row.names = FALSE)
#write.table(df_sgmRNA, file = args[4], sep = "\t", row.names = FALSE)
#write.table(df_alt_summary, file = args[5], sep = "\t", row.names = FALSE)
#write.table(df_DVG, file = args[6], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)