#' ---
#' title: "Mapping statistics"
#' author: "Benesova"
#' ---

#+ setup1, include=FALSE
library(magrittr)
library(ggplot2)
library(gridExtra)
library(SummarizedExperiment)
library(reshape2)
library(cowplot)
library(grid)
library("readxl")

#+ load annotations, include=FALSE
load("temp/annot_miRXplore.RData")
annot_miRXplore <-annot_miRXplore[(annot_miRXplore$Kit != "NEXTflex_UMI")&
                                    (annot_miRXplore$Kit != "QIAseq_UMI"),]
annot_miRXplore <-droplevels(annot_miRXplore)

load("temp/annot_plasma.RData")
annot_plasma <-annot_plasma[(annot_plasma$Kit != "NEXTflex_UMI")&
                                    (annot_plasma$Kit != "QIAseq_UMI"),]
annot_plasma <-droplevels(annot_plasma)

#+ load mapping miRXplore log df, include=FALSE
path<-"src/log_mapping/"
files<-list.files(path)
idx<-grep("tab",files)
log_files<-files[idx]
log_files<-paste0(path, log_files)
map_names <- c("A","B","C","Ch","D","E","F","G","H","I","J")
miRXplore_mapping <- as.data.frame(c(seq(1,14,by=1)))
counter <- 0
for (i in log_files){
  counter <- counter + 1
  df <- as.data.frame(read.table(file = i, header = TRUE, sep = "\t", fill=TRUE))
  rownames(df) <- df$Sample
  cnames <-rep(map_names[counter], length(colnames(df)))
  colnames(df) <- paste(cnames,colnames(df),sep="_")
  df <- df[rownames(annot_miRXplore),]
  if (!("NA" %in% rownames(df))){
    miRXplore_mapping <- cbind(miRXplore_mapping,df)
    }
}
miRXplore_mapping <- as.data.frame(miRXplore_mapping)

#'**Discarded reads miRXplore**
#+ discarded miRXplore,cache=FALSE, results='show', echo=FALSE
Discarded_reads <- as.data.frame(cbind(miRXplore_mapping$A_Reads.too.short.after.trimming, 
                                       miRXplore_mapping$B_Short_reads,
                                       miRXplore_mapping$B_Long_reads,
                                       miRXplore_mapping$C_Aligned.reads, 
                                       miRXplore_mapping$D_Aligned.reads, 
                                       miRXplore_mapping$A_Input.reads))
rownames(Discarded_reads) <- rownames(miRXplore_mapping)
colnames(Discarded_reads) <- c("Trim","Short","Long","rRNA","Spikes", "Raw_reads")
df <-  c(1,2,3,4,5)
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  means <- colMeans(Discarded_reads[rownames(annot_miRXplore[annot_miRXplore$Kit == i,]),])
  percent_all <- round((means[c(1:(length(means)-1))]/means[length(means)])*100,2)
  percent_discarded <- round((means[c(1:(length(means)-1))]/sum(means[c(1:(length(means)-1))]))*100,2)
  df <- rbind(df, as.numeric(as.character(percent_discarded)))
}
df <- as.data.frame(df)
colnames(df)<- c("Trim","Short","Long","rRNA","Spikes")
df <- df[-1,]
rownames(df) <- levels(annot_miRXplore$Kit)
df$Protocol <- rownames(df)
write.csv(df, file = "outs/miRXplore_discarded.csv")
df <-melt(df)
df$variable <- factor(df$variable, level = c("Trim","Short","Long","rRNA","Spikes"))
df$Protocol <- factor(df$Protocol, level = c("Lexogen","Norgen","QIAseq","NEXTflex",
                                             "RealSeq","SMARTer","EdgeSeq"))
discarded_miRXplore <- ggplot(data=df, aes(x=Protocol, y=value, color = variable, fill = variable))+
  geom_bar(stat="identity")+
  ylab("% of discarded reads")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_brewer(palette="Blues")+
  scale_color_brewer(palette="Blues")+
  ggtitle("Discarded reads miRXplore")
legend <- cowplot::get_legend(discarded_miRXplore)
grid.newpage()
grid.draw(legend)
discarded_miRXplore <- discarded_miRXplore + theme(legend.position = "none")
discarded_miRXplore

#'**Mapping reads miRXplore**
#+ mapping miRXplore,cache=FALSE, results='show', echo=FALSE
miRXplore_mapped <- miRXplore_mapping$E_Unique + miRXplore_mapping$E_Multimapping
false_isomiRs <-as.data.frame(read.table("src/log_mapping/miRXplore_isomiRs.txt", 
                                         header = TRUE, sep = "\t"))
false_isomiRs <-false_isomiRs[,-1]
false_isomiRs <- colSums(false_isomiRs)
top10 <- read.table("outs/top10_reads_miRXplore.csv", header = TRUE, sep =",")
top10 <- top10[top10$Protocol != "QIAseq_UMI",]
top10 <-unlist(lapply(top10$mean,function(x){rep(x,2)}))
mapping <- as.data.frame(cbind(miRXplore_mapped, 
                               top10,
                               false_isomiRs,
                               rowSums(Discarded_reads[,c(1:5)]),
                               miRXplore_mapping$A_Input.reads))
rownames(mapping) <- rownames(miRXplore_mapping)
colnames(mapping) <- c("Unique","top10","Isomirs","Discarded","Raw_reads")
df <-  c(1,2,3,4,5,6,7,8)
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  means <- colMeans(mapping[rownames(annot_miRXplore[annot_miRXplore$Kit == i,]),])
  means["Mapped"] <- means["Unique"] - means["top10"]
  means["smallRNA"] <- 0
  means["Genome"] <- 0
  means["Unmapped"] <- means["Raw_reads"] - (means["Unique"] + means["Discarded"] + means["Isomirs"]
                                             + means["smallRNA"] + means["Genome"])
  percent_all <- round((means[names(means) != "Raw_reads"]/means["Raw_reads"])*100,2)
  df <- rbind(df, as.numeric(as.character(percent_all)))
}
df <- as.data.frame(df)
colnames(df)<- c("Unique","top10","Isomirs","Discarded","Mapped", "smallRNA", "Genome","Unmapped")
df <- df[-1,]
rownames(df) <- levels(annot_miRXplore$Kit)
df$Protocol <- rownames(df)
write.csv(df, file = "outs/miRXplore_mapping.csv")
df <-melt(df)
df <- df[df$variable %in% c("Discarded","Unmapped","Genome","smallRNA","Isomirs","top10","Mapped"),]
df$variable <- factor(df$variable, level = c("Discarded","Unmapped","Genome","smallRNA","Isomirs","top10","Mapped"))
df$Protocol <- factor(df$Protocol, level = c("Lexogen","Norgen","QIAseq","NEXTflex",
                                             "RealSeq","SMARTer","EdgeSeq"))
mapping_miRXplore <- ggplot(data=df, aes(x=Protocol, y=value, color = variable, fill = variable))+
  geom_bar(stat="identity")+
  ylab("% of raw reads")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_brewer(palette="PRGn")+
  scale_color_brewer(palette="PRGn")+
  ggtitle("Mapping statistics miRXplore")
legend <- cowplot::get_legend(mapping_miRXplore)
grid.newpage()
grid.draw(legend)
mapping_miRXplore <- mapping_miRXplore + theme(legend.position = "none")
mapping_miRXplore

#+ load mapping plasma log df, include=FALSE
path<-"src/log_mapping/"
files<-list.files(path)
idx<-grep("tab",files)
log_files<-files[idx]
log_files<-paste0(path, log_files)
map_names <- c("A","B","C","Ch","D","E","F","G","H","I","J")
plasma_mapping <- as.data.frame(c(seq(1,14,by=1)))
counter <- 0
for (i in log_files){
  counter <- counter + 1
  df <- as.data.frame(read.table(file = i, header = TRUE, sep = "\t", fill=TRUE))
  rownames(df) <- df$Sample
  cnames <-rep(map_names[counter], length(colnames(df)))
  colnames(df) <- paste(cnames,colnames(df),sep="_")
  df <- df[rownames(annot_plasma),]
  if (!("NA" %in% rownames(df))){
    plasma_mapping <- cbind(plasma_mapping,df)
  }
}
plasma_mapping <- as.data.frame(plasma_mapping)

#'**Discarded reads plasma**
#+ discarded plasma,cache=FALSE, results='show', echo=FALSE
Discarded_reads <- as.data.frame(cbind(plasma_mapping$A_Reads.too.short.after.trimming, 
                                       plasma_mapping$B_Short_reads,
                                       plasma_mapping$B_Long_reads,
                                       plasma_mapping$C_Aligned.reads, 
                                       plasma_mapping$D_Aligned.reads, 
                                       plasma_mapping$A_Input.reads))

rownames(Discarded_reads) <- rownames(plasma_mapping)
colnames(Discarded_reads) <- c("Trim","Short","Long","rRNA","Spikes", "Raw_reads")
df <-  c(1,2,3,4,5)
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  print(i)
  means <- colMeans(Discarded_reads[rownames(annot_plasma[annot_plasma$Kit == i,]),])
  percent_all <- round((means[c(1:(length(means)-1))]/means[length(means)])*100,2)
  percent_discarded <- round((means[c(1:(length(means)-1))]/sum(means[c(1:(length(means)-1))]))*100,2)
  df <- rbind(df, as.numeric(as.character(percent_discarded)))
}
df <- as.data.frame(df)
colnames(df)<- c("Trim","Short","Long","rRNA","Spikes")
df <- df[-1,]
rownames(df) <- levels(annot_plasma$Kit)
df$Protocol <- rownames(df)
write.csv(df, file = "outs/plasma_discarded.csv")
df <-melt(df)
df$variable <- factor(df$variable, level = c("Trim","Short","Long","rRNA","Spikes"))
df$Protocol <- factor(df$Protocol, level = c("Lexogen","Norgen","QIAseq","NEXTflex",
                                             "RealSeq","SMARTer","EdgeSeq"))
discarded_plasma <- ggplot(data=df, aes(x=Protocol, y=value, color = variable, fill = variable))+
  geom_bar(stat="identity")+
  ylab("% of discarded reads")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_brewer(palette="Blues")+
  scale_color_brewer(palette="Blues")+
  ggtitle("Discarded reads plasma")
legend <- cowplot::get_legend(discarded_plasma)
grid.newpage()
grid.draw(legend)
discarded_plasma <- discarded_plasma + theme(legend.position = "none")
discarded_plasma

#'**Mapping reads plasma**
#+ mapping plasma,cache=FALSE, results='show', echo=FALSE
smallRNA <- as.data.frame(cbind(plasma_mapping$Ch_Unique + plasma_mapping$Ch_Multimapping,
                                 plasma_mapping$I_Unique + plasma_mapping$I_Multimapping,
                                 plasma_mapping$J_Unique + plasma_mapping$J_Multimapping))
rownames(smallRNA) <-rownames(plasma_mapping)
colnames(smallRNA) <- c("tRNA", "piRNA", "ncRNA")
mature_mapped <-plasma_mapping$G_Unique + plasma_mapping$G_Multimapping
genome_mapped <- plasma_mapping$F_Unique + plasma_mapping$F_Multimapping
isomiRs <-as.data.frame(read.table("src/log_mapping/plasma_isomiRs.txt", 
                                         header = TRUE, sep = "\t"))
isomiRs <-isomiRs[,-1]
isomiRs <- colSums(isomiRs)
top10 <- read.table("outs/top10_reads_plasma.csv", header = TRUE, sep =",")
top10 <- top10[top10$Protocol != "QIAseq_UMI",]
top10 <-unlist(lapply(top10$mean,function(x){rep(x,2)}))
mapping <- as.data.frame(cbind(mature_mapped, 
                               top10,
                               isomiRs,
                               rowSums(Discarded_reads[,c(1:5)]), 
                               rowSums(smallRNA),
                               genome_mapped,
                               plasma_mapping$A_Input.reads))
rownames(mapping) <- rownames(plasma_mapping)
colnames(mapping) <- c("Unique","top10","Isomirs","Discarded","smallRNA","Genome_mapped","Raw_reads")
df <-  c(1,2,3,4,5,6,7,8,9)
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  means <- colMeans(mapping[rownames(annot_plasma[annot_plasma$Kit == i,]),])
  means["Mapped"] <- means["Unique"] - means["top10"]
  means["Genome"] <- means["Genome_mapped"] - (means["Unique"] + means["Isomirs"] + means["smallRNA"])
  means["Unmapped"] <- means["Raw_reads"] - (means["Unique"] + means["Discarded"] + means["Isomirs"]
                                             + means["smallRNA"] + means["Genome"])
  percent_all <- round((means[names(means) != "Raw_reads"]/means["Raw_reads"])*100,2)
  df <- rbind(df, as.numeric(as.character(percent_all)))
}
df <- as.data.frame(df)
colnames(df)<- c("Unique","top10","Isomirs","Discarded", "smallRNA","Genome_mapped","Mapped","Genome","Unmapped")
df <- df[-1,]
rownames(df) <- levels(annot_plasma$Kit)
df$Protocol <- rownames(df)
write.csv(df, file = "outs/plasma_mapping.csv")
df <-melt(df)
df <- df[df$variable %in% c("Discarded","Unmapped","Genome","smallRNA","Isomirs","top10", "Mapped"),]
df$variable <- factor(df$variable, level = c("Discarded","Unmapped","Genome","smallRNA","Isomirs","top10", "Mapped"))
df$Protocol <- factor(df$Protocol, level = c("Lexogen","Norgen","QIAseq","NEXTflex",
                                             "RealSeq","SMARTer","EdgeSeq"))
mapping_plasma <- ggplot(data=df, aes(x=Protocol, y=value, color = variable, fill = variable))+
  geom_bar(stat="identity")+
  ylab("% of raw reads")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_brewer(palette="PRGn")+
  scale_color_brewer(palette="PRGn")+
  ggtitle("Mapping statistics plasma")
legend <- cowplot::get_legend(mapping_plasma)
grid.newpage()
grid.draw(legend)
mapping_plasma <- mapping_plasma + theme(legend.position = "none")
mapping_plasma

#'**SmallRNA plasma**
#+ mapping plasma,cache=FALSE, results='show', echo=FALSE
df <-  c(1,2,3)
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  means <- colMeans(smallRNA[rownames(annot_plasma[annot_plasma$Kit == i,]),])
  percent_smallRNA <- round((means/sum(means))*100,2)
  df <- rbind(df, as.numeric(as.character(percent_smallRNA)))
}
colnames(df)<- c("tRNA","piRNA","ncRNA")
df <- as.data.frame(df[-1,])
rownames(df) <- levels(annot_plasma$Kit)
df$Protocol <- rownames(df)
write.csv(df, file = "outs/plasma_smallRNA.csv")
df<-melt(df)
df$Protocol <- factor(df$Protocol, level = c("Lexogen","Norgen","QIAseq","NEXTflex",
                                             "RealSeq","SMARTer","EdgeSeq"))
smallRNA <- ggplot(data=df, aes(x=Protocol, y=value, color = variable, fill = variable))+
  geom_bar(stat="identity")+
  ylab("% reads mapped to other small RNA")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_brewer(palette="BuGn")+
  scale_color_brewer(palette="BuGn")+
  ggtitle("smallRNA")
legend <- cowplot::get_legend(smallRNA)
grid.newpage()
grid.draw(legend)
smallRNA <- smallRNA + theme(legend.position = "none")
smallRNA

grid.arrange(mapping_miRXplore, mapping_plasma, discarded_miRXplore, discarded_plasma, nrow=2)

#rm(list = setdiff(ls(), lsf.str()))
