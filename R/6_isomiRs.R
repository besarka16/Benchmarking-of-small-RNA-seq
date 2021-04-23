#' ---
#' title: "IsomiRs and false isomiRs"
#' author: "Benesova + Androvic"
#' ---

#+ setup1, include=FALSE
library(GenomicFeatures)
library(GenomicAlignments)
library(GenomicRanges)
library(IRanges)
library(Rsubread)
library(AnnotationDbi)
library(magrittr)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(grid)
library(DESeq2)
library(variancePartition)
library(ggseqlogo)
library(RColorBrewer)
library(dplyr)
library(reshape)
library(maditr)

#'**False isomiRs in miRXplore sample**
#+ load isomiR data, cache=FALSE, results='hide', echo=FALSE
con <- file("src/log_mapping/miRXplore_isomiRs.txt","r")
first_line <- readLines(con,n=1)
first_line <- read.table(textConnection(first_line))
first_line <- first_line[2:length(first_line)]
X <-c(paste0("X.",seq(1,20, by=1)))
X <-t(as.data.frame(c("X",X)))
first_line <- t(as.data.frame((cbind(X, first_line))))
close(con)
x <- readLines("src/log_mapping/miRXplore_isomiRs.txt")
x <- gsub( "Xed", "X_seed", x )
x <- gsub( "_", "\t", x )
isomiRs <-read.table(textConnection(x),fill=TRUE, sep = "\t", skip=1)
colnames(isomiRs) <- first_line
rm(X,x,first_line,con)

#+ 3end distribution, cache=FALSE, results='hide', echo=FALSE
load("temp/annot_miRXplore.RData")
annot <- annot_miRXplore[!(annot_miRXplore$Kit %in% c("QIAseq_UMI", "NEXTflex_UMI")),]
df <- c()
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df<- cbind(df,rowMeans(isomiRs[,rownames(annot[annot$Kit == i,])]))
}
df <- as.data.frame(df)
colnames(df) <- unique(annot[order(annot$Kit),"Kit"])
df$end_5_start <- as.factor(substring(isomiRs$X.3,1,1))
df$add_5 <- as.factor(isomiRs$X.1)
df$length_5 <- as.factor(isomiRs$X.2)
df$end_3_end <- as.factor(unlist(lapply(isomiRs$X.7, function(x) substring(x, length(x)-1, length(x)))))
df$add_3 <- as.factor(isomiRs$X.5)
df$length_3 <- as.factor(isomiRs$X.6)

df_3 <- df[df$add_3 %in% c("Add3","Can3"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_3 <- cbind(df_3,(df_3[,i]/sum(df_3[,i]))*100)
}
df_3 <- df_3[,c(8:dim(df_3)[2])]
colnames(df_3) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_3)
colourCount = length(unique(melted$end_3_end))
mycolors <- colorRampPalette(brewer.pal(11, "PRGn"))(colourCount)
miRXplore_3_end <- ggplot(data=melted, aes(x=variable, y=value, color = end_3_end, fill = end_3_end))+
  geom_bar(stat="identity")+
  ylab("% isomiRs raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("False isomiRs distribution on 3'")
legend <- cowplot::get_legend(miRXplore_3_end)
grid.newpage()
grid.draw(legend)
miRXplore_3_end <- miRXplore_3_end + theme(legend.position = "none")
miRXplore_3_end

#+ 3end distribution export table, cache=FALSE, results='hide', echo=FALSE
#Get the % of first 3? isomir bases for each kit. There are 2 ways to do it:
#1st way
temp1 <- group_by(melted, variable, end_3_end)#use dplyr package to group data by Kit names (variable), and base (end_3_end)
temp2 <- temp1 %>% summarise(sum(value))#obtain sums of reads for each base per each kit
end_3_distribution.1 <- dcast(temp2, end_3_end~variable) #reformat it into wide format for visualization
write.table(end_3_distribution.1, file = "outs/false.isomir_3end_distribution.txt", col.names = T, row.names = F, quote = F, sep = "\t")

#2nd way
end_3_distribution.2 <- data.frame(matrix(NA, ncol=1, nrow=5))[-1]
for (i in unique(annot[order(annot$Kit),"Kit"])){
temp1 <- melted[melted$variable ==i,]
temp2 <- as.data.frame(temp1 %>% group_by(end_3_end) %>% summarise(sum(value)))[2]
colnames(temp2) <- i
end_3_distribution.2 <- cbind(end_3_distribution.2,temp2)
rownames(end_3_distribution.2) <- c("A","C","G","T","X")
}

#+ length distribution on 3end, cache=FALSE, results='hide', echo=FALSE
df_3 <- df[df$add_3 %in% c("Add3","Can3","Trim3"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_3 <- cbind(df_3,(df_3[,i]/sum(df_3[,i]))*100)
}
df_3 <- df_3[,c(8:dim(df_3)[2])]
colnames(df_3) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_3)
colourCount = length(unique(melted$length_3))
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(colourCount)
miRXplore_3_end_length <- ggplot(data=melted, aes(x=variable, y=value, color = length_3, fill = length_3))+
  geom_bar(stat="identity")+ 
  ylab("% isomiRs raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("False isomiRs length distribution on 3'")
legend <- cowplot::get_legend(miRXplore_3_end_length)
grid.newpage()
grid.draw(legend)
miRXplore_3_end_length <- miRXplore_3_end_length + theme(legend.position = "none")
miRXplore_3_end_length

#+ 5end distribution, cache=FALSE, results='hide', echo=FALSE
df_5 <- df[df$add_5 %in% c("Add5","Can5"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_5 <- cbind(df_5,(df_5[,i]/sum(df_5[,i]))*100)
}
df_5 <- df_5[,c(8:dim(df_5)[2])]
colnames(df_5) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_5)
colourCount = length(unique(melted$end_5_start))
mycolors <- colorRampPalette(brewer.pal(11, "PRGn"))(colourCount)
miRXplore_5_end <- ggplot(data=melted, aes(x=variable, y=value, color = end_5_start, fill = end_5_start))+
  geom_bar(stat="identity")+
  ylab("% of isomiRs Raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("False isomiRs distribution on 5'")
legend <- cowplot::get_legend(miRXplore_5_end)
grid.newpage()
grid.draw(legend)
miRXplore_5_end <- miRXplore_5_end + theme(legend.position = "none")
miRXplore_5_end

#+ 5end distribution export table, cache=FALSE, results='hide', echo=FALSE
temp1 <- group_by(melted, variable, end_5_start)#use dplyr package to group data by Kit names (variable), and base (end_3_end)
temp2 <- temp1 %>% summarise(sum(value))#obtain sums of reads for each base per each kit
end_5_distribution <- dcast(temp2, end_5_start~variable) #reformat it into wide format for visualization
write.table(end_5_distribution, file = "outs/false.isomir_5end_distribution.txt", col.names = T, row.names = F, quote = F, sep = "\t")

#+ length distribution on 3end, cache=FALSE, results='hide', echo=FALSE
df_5 <- df[df$add_5 %in% c("Add5","Can5","Trim5"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_5 <- cbind(df_5,(df_5[,i]/sum(df_5[,i]))*100)
}
df_5 <- df_5[,c(8:dim(df_5)[2])]
colnames(df_5) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_5)
colourCount = length(unique(melted$length_5))
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(colourCount)
miRXplore_5_end_length <- ggplot(data=melted, aes(x=variable, y=value, color = length_5, fill = length_5))+
  geom_bar(stat="identity")+
  ylab("% of isomiRs Raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("False isomiRs length distribution on 5'")
legend <- cowplot::get_legend(miRXplore_5_end_length)
grid.newpage()
grid.draw(legend)
miRXplore_5_end_length <- miRXplore_5_end_length + theme(legend.position = "none")
miRXplore_5_end_length

#'**IsomiRs in plasma sample**
#+ load isomiRs plasma data, cache=FALSE, results='hide', echo=FALSE
con <- file("src/log_mapping/plasma_isomiRs.txt","r")
first_line <- readLines(con,n=1)
first_line <- read.table(textConnection(first_line))
first_line <- first_line[2:length(first_line)]
X <-c(paste0("X.",seq(1,20, by=1)))
X <-t(as.data.frame(c("X",X)))
first_line <- t(as.data.frame((cbind(X, first_line))))
close(con)
x <- readLines("src/log_mapping/plasma_isomiRs.txt")
x <- gsub( "Xed", "X_seed", x )
x <- gsub( "_", "\t", x )
isomiRs <-read.table(textConnection(x),fill=TRUE, sep = "\t", skip=1)
colnames(isomiRs) <- first_line
rm(X,x,first_line,con)

#+ plasma load annot, cache=FALSE, results='hide', echo=FALSE
load("temp/annot_plasma.RData")
annot <- annot_plasma[!(annot_plasma$Kit %in% c("QIAseq_UMI", "NEXTflex_UMI")),]
df <- c()
for (i in unique(annot[order(annot$Kit),"Kit"])){
  print(i)
  df<- cbind(df,rowMeans(isomiRs[,rownames(annot[annot$Kit == i,])]))
}
df <- as.data.frame(df)
colnames(df)<-unique(annot[order(annot$Kit),"Kit"])
df$miRNAs <- as.factor(isomiRs$X)

#+ plasma top10 contribution, cache=FALSE, results='hide', echo=FALSE
counts_per_miRNA_all <- c()
for (kit in unique(annot[order(annot$Kit),"Kit"])){
  counts_per_miRNA <-c()
  for (i in unique(df$miRNAs)){
    vec <- lapply(df[df$miRNAs == i,kit], as.numeric)
    count <-sum(unlist(vec))
    counts_per_miRNA <- rbind(counts_per_miRNA,c(i,count))
  }
  counts_per_miRNA <- as.data.frame(counts_per_miRNA)
  counts_per_miRNA_all <- cbind(counts_per_miRNA_all,counts_per_miRNA$V2)
}
counts_per_miRNA_all <- as.data.frame(counts_per_miRNA_all)
colnames(counts_per_miRNA_all) <- unique(annot[order(annot$Kit),"Kit"])
rownames(counts_per_miRNA_all) <- unique(df$miRNAs)
top10 <- read.table("outs/top10_names_plasma.csv", header = TRUE, sep = ",")
top10 <- top10[,-1]

data <-c()
for (kit in unique(annot[order(annot$Kit),"Kit"])){
  counts_top10 <-as.numeric(sum(counts_per_miRNA_all[rownames(counts_per_miRNA_all) %in% top10[,kit],kit]))
  counts_others <-as.numeric(sum(counts_per_miRNA_all[!(rownames(counts_per_miRNA_all) %in% top10[,kit]),kit]))
  total <- sum(counts_top10,counts_others)
  percent_top10 <- (counts_top10/total) *100
  percent_others <- (counts_others/total) *100
  data <- rbind(data, c(counts_top10, counts_others,percent_top10,percent_others, kit))
}
data <- as.data.frame(data)
colnames(data) <- c("top10", "others","p_top10","p_others", "Protocol")
melted <-melt(data, id.vars = "Protocol")
melted$Protocol <- factor(melted$Protocol, level = c("Lexogen","Norgen","QIAseq","NEXTflex",
                                                                "RealSeq","SMARTer","EdgeSeq"))
melted$value <-lapply(melted$value,as.numeric)
melted$variable <- factor(melted$variable, level=c("others","top10","p_top10","p_others"))
melted <- melted[melted$variable %in% c("p_top10","p_others"),]
ggplot(data=melted, aes(x=Protocol, y=value, color = variable, fill = variable))+
  geom_bar(stat="identity")+
  ylab("Raw counts")+
  xlab("")+
  scale_fill_brewer(palette="Paired")+
  scale_color_brewer(palette="Paired")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ggtitle("Representation of isomiRs counts belonging to top10 miRNAs")

#+ plasma 3end distribution, cache=FALSE, results='hide', echo=FALSE
df <- c()
for (i in unique(annot[order(annot$Kit),"Kit"])){
  print(i)
  df<- cbind(df,rowMeans(isomiRs[,rownames(annot[annot$Kit == i,])]))
}
df <- as.data.frame(df)
colnames(df) <- unique(annot[order(annot$Kit),"Kit"])
df$end_5_start <- as.factor(substring(isomiRs$X.3,1,1))
df$add_5 <- as.factor(isomiRs$X.1)
df$length_5 <- as.factor(isomiRs$X.2)
df$end_3_end <- as.factor(unlist(lapply(isomiRs$X.7, function(x) substring(x, length(x)-1, length(x)))))
df$add_3 <- as.factor(isomiRs$X.5)
df$length_3 <- as.factor(isomiRs$X.6)

df_3 <- df[df$add_3 %in% c("Add3","Can3"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_3 <- cbind(df_3,(df_3[,i]/sum(df_3[,i]))*100)
}
df_3 <- df_3[,c(8:dim(df_3)[2])]
colnames(df_3) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_3)
colourCount = length(unique(melted$end_3_end))
mycolors <- colorRampPalette(brewer.pal(11, "PRGn"))(colourCount)
plasma_3_end <- ggplot(data=melted, aes(x=variable, y=value, color = end_3_end, fill = end_3_end))+
  geom_bar(stat="identity")+
  ylab("% isomiRs raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("IsomiRs distribution on 3'")
legend <- cowplot::get_legend(plasma_3_end)
grid.newpage()
grid.draw(legend)
plasma_3_end <- plasma_3_end + theme(legend.position = "none")
plasma_3_end

#+ plasma 3end length distribution, cache=FALSE, results='hide', echo=FALSE
df_3 <- df[df$add_3 %in% c("Add3","Can3","Trim3"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_3 <- cbind(df_3,(df_3[,i]/sum(df_3[,i]))*100)
}
df_3 <- df_3[,c(8:dim(df_3)[2])]
colnames(df_3) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_3)
colourCount = length(unique(melted$length_3))
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(colourCount)
plasma_3_end_length <- ggplot(data=melted, aes(x=variable, y=value, color = length_3, fill = length_3))+
  geom_bar(stat="identity")+
  ylab("% isomiRs raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("IsomiRs length distribution on 3'")
legend <- cowplot::get_legend(plasma_3_end_length)
grid.newpage()
grid.draw(legend)
plasma_3_end_length <- plasma_3_end_length + theme(legend.position = "none")
plasma_3_end_length

#+ plasma 5end distribution, cache=FALSE, results='hide', echo=FALSE
df_5 <- df[df$add_5 %in% c("Add5","Can5"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_5 <- cbind(df_5,(df_5[,i]/sum(df_5[,i]))*100)
}
df_5 <- df_5[,c(8:dim(df_5)[2])]
colnames(df_5) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_5)
colourCount = length(unique(melted$end_5_start))
mycolors <- colorRampPalette(brewer.pal(11, "PRGn"))(colourCount)
plasma_5_end <- ggplot(data=melted, aes(x=variable, y=value, color = end_5_start, fill = end_5_start))+
  geom_bar(stat="identity")+
  ylab("% of isomiRs raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("IsomiRs distribution on 5'")
legend <- cowplot::get_legend(plasma_5_end)
grid.newpage()
grid.draw(legend)
plasma_5_end <- plasma_5_end + theme(legend.position = "none")
plasma_5_end

#+ plasma 5end length distribution, cache=FALSE, results='hide', echo=FALSE
df_5 <- df[df$add_5 %in% c("Add5","Can5","Trim5"),]
for (i in unique(annot[order(annot$Kit),"Kit"])){
  df_5 <- cbind(df_5,(df_5[,i]/sum(df_5[,i]))*100)
}
df_5 <- df_5[,c(8:dim(df_5)[2])]
colnames(df_5) <- c("end_5_start", "add_5", "length_5","end_3_end", "add_3","length_3",
                    paste0(unique(annot[order(annot$Kit),"Kit"])))
melted <-melt(df_5)
colourCount = length(unique(melted$length_5))
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(colourCount)
plasma_5_end_length <- ggplot(data=melted, aes(x=variable, y=value, color = length_5, fill = length_5))+
  geom_bar(stat="identity")+
  ylab("% of isomiRs raw counts")+
  xlab("")+
  labs(fill = "",color = "")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = mycolors)+
  scale_color_manual(values = mycolors)+
  ggtitle("IsomiRs length distribution on 5'")
legend <- cowplot::get_legend(plasma_5_end_length)
grid.newpage()
grid.draw(legend)
plasma_5_end_length <- plasma_5_end_length + theme(legend.position = "none")
plasma_5_end_length



