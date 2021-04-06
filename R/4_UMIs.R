#' ---
#' title: "UMI"
#' author: "Benesova + Androvic"
#' ---

#+ setup1, include=FALSE
library(GGally)
library(ggplot2)
library(ggridges)
library(reshape2)
library(variancePartition)
library(gridExtra)
library(grid)
library(pheatmap)
library(matrixStats)
library(magrittr)
library(ggpubr)

#+ load annotations, include=FALSE
load("src/UMI/miRXplore_UMI.RData")
load("src/UMI/plasma_UMI.RData")
annot <- read.table("src/UMI/annotation_UMI.txt", sep = "\t", header = TRUE)
rownames(annot) <- annot$Sample

miRXplore_mean <- c()
for (i in unique(annot$Kit)){
  miRXplore_mean <- cbind(miRXplore_mean,rowMeans(miRXplore_UMI[,rownames(annot[annot$Kit==i,])]))
  }
miRXplore_mean <- as.data.frame(miRXplore_mean)
colnames(miRXplore_mean)<-unique(annot$Kit)


miRXplore_mean$QIAseq_8bp_ratio <- miRXplore_mean$QIAseq/miRXplore_mean$QIAseq_8
miRXplore_mean$QIAseq_10bp_ratio <- miRXplore_mean$QIAseq/miRXplore_mean$QIAseq_10
miRXplore_mean$QIAseq_12bp_ratio <- miRXplore_mean$QIAseq/miRXplore_mean$QIAseq_12
miRXplore_mean$NEXTflex_8bp_ratio <- miRXplore_mean$NEXTflex/miRXplore_mean$NEXTflex_8


plasma_mean <- c()
for (i in unique(annot$Kit)){
  plasma_mean <- cbind(plasma_mean,rowMeans(plasma_UMI[,rownames(annot[annot$Kit==i,])]))
}
plasma_mean <- as.data.frame(plasma_mean)
colnames(plasma_mean)<-unique(annot$Kit)


plasma_mean$QIAseq_8bp_ratio <- plasma_mean$QIAseq/plasma_mean$QIAseq_8
plasma_mean$QIAseq_10bp_ratio <- plasma_mean$QIAseq/plasma_mean$QIAseq_10
plasma_mean$QIAseq_12bp_ratio <- plasma_mean$QIAseq/plasma_mean$QIAseq_12
plasma_mean$NEXTflex_8bp_ratio <- plasma_mean$NEXTflex/plasma_mean$NEXTflex_8

theme_set(theme_bw())
p1<-ggplot(data=miRXplore_mean, aes(x=log2(QIAseq),y=QIAseq_8bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("miRXplore QIAseq 8 bp UMI")+
  ylim(0,17)
p1


p2<-ggplot(data=miRXplore_mean, aes(x=log2(QIAseq),y=QIAseq_10bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("miRXplore QIAseq 10 bp UMI")+
  ylim(0,17)
p2
p3<-ggplot(data=miRXplore_mean, aes(x=log2(QIAseq),y=QIAseq_12bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("miRXplore QIAseq 12 bp UMI")+
  ylim(0,17)
p3

p4<-ggplot(data=miRXplore_mean, aes(x=log2(NEXTflex),y=NEXTflex_8bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("miRXplore NEXTflex 8 bp UMI")+
  ylim(0,17)
p4


p5<-ggplot(data=plasma_mean, aes(x=log2(QIAseq),y=QIAseq_8bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("plasma QIAseq 8 bp UMI")+
  ylim(0,40)
p5
p6<-ggplot(data=plasma_mean, aes(x=log2(QIAseq),y=QIAseq_10bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("plasma QIAseq 10 bp UMI")+
  ylim(0,40)
p6
p7<-ggplot(data=plasma_mean, aes(x=log2(QIAseq),y=QIAseq_12bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("plasma QIAseq 12 bp UMI")+
  ylim(0,40)
p7
p8<-ggplot(data=plasma_mean, aes(x=log2(NEXTflex),y=NEXTflex_8bp_ratio))+
  geom_point(color = "Lightsteelblue", size = 2)+
  geom_smooth(method = "loess", span = 0.4, color = "#8f0021",  alpha = 0.4, size = 0.5)+
  ylab("Correction ratio (counts/UMIs)")+
  xlab("log2(counts)")+
  theme(panel.grid.minor = element_blank())+
  ggtitle("plasma NEXTflex 8 bp UMI")+
  ylim(0,40)
p8

g<-grid.arrange(p1,p2,p3,p5,p6,p7,p4,p8, ncol = 3)
ggsave(file= "outs/UMI_length.pdf",plot = g, device = "pdf", height = 40, width = 40, units = "cm")

