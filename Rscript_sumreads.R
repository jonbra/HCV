# NB! Slight deviation in the number of mapped reads if I just look in the weeSAM stats file. Why?
#.libPaths("/home/ngs/miniconda3/lib/R/library")
.libPaths("/home/jonbra/miniconda3/envs/HCV/lib/R/library")

args = commandArgs(trailingOnly=TRUE)

df<-read.csv(file=args[1], header=TRUE, sep="\t", dec=".")

# df<-read.delim(file="Q-3a-1x-2-R_S37_ICTV_tanoti_stats.txt", sep="\t", dec=".")

dat1 <- data.frame(do.call(rbind, strsplit(as.vector(df$Ref_Name), split = "_")))
colnames(dat1)<-c("Genotype", "ID")

df_full<-cbind(df, dat1)

sum_reads<-aggregate(Mapped_Reads ~ Genotype, df_full, sum)
sum_reads$Mapped_Reads<-as.numeric(sum_reads$Mapped_Reads)
sum_reads<-sum_reads[order(sum_reads$Mapped_Reads, decreasing=TRUE),]

sum_reads$Percent<-sum_reads$Mapped_Reads / sum(sum_reads$Mapped_Reads)

write.table(sum_reads, file=args[2], row.names=FALSE)
