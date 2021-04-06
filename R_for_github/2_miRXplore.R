#' ---
#' title: "miRXplore"
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
library(DESeq2)
library(RColorBrewer)
library(ggpubr)

#+ load data, include=FALSE
load("temp/miRXplore.RData")
load("temp/miRXplore_norm.RData")
load("temp/annot_miRXplore.RData")
targeted_in_miRXplore_EdgeSeq <- as.data.frame(read.table("src/miRXplore_targeted_EdgeSeq.tab", header = TRUE, sep = "\t"))
rownames(targeted_in_miRXplore_EdgeSeq) <- targeted_in_miRXplore_EdgeSeq$miRNA_name
annot_miRXplore <- annot_miRXplore[annot_miRXplore$Kit != "NEXTflex_UMI",]
annot_miRXplore <- droplevels(annot_miRXplore)
miRXplore <- miRXplore[,rownames(annot_miRXplore)]
miRXplore_norm <- miRXplore_norm[,rownames(annot_miRXplore)]

#'**miRXplore total counts per sample.**
#+ count_tables,cache=FALSE, results='show', echo=FALSE
Total_counts <- c("Protocol","Mean_total_counts")
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  Total_counts<-rbind(Total_counts,c(i,mean(colSums(miRXplore[,annot_miRXplore$Kit == i]))))
}
Total_counts
write.csv(miRXplore,file = "outs/miRXplore_raw_counts.csv")
rm(Total_counts)

#'**miRXplore coefficient of variation.**
#'Mean coefficient of variation (standard deviation/mean) for two replicates.
#+ cv,cache=FALSE, results='show', echo=FALSE
cv<-c("Protocol","Mean_CV")
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  if (i =="EdgeSeq"){
    df <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
    df <- df[rownames(df) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
  }
  else {
  df <-miRXplore_norm[,annot_miRXplore$Kit == i]
  }
  CV1 <- sd(df[,1])/mean(df[,1])
  CV2 <- sd(df[,2])/mean(df[,2])
  cv <-rbind(cv,c(i,mean(CV1,CV2)))
}
cv
rm(cv,CV1,CV2)

#'**Between kit correlation**
#+ between kit correlation,cache=FALSE, results='show', echo=FALSE

library(corrplot)
data <- c()
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  mean <- rowMeans(miRXplore_norm[,annot_miRXplore$Kit == i])
  data <-cbind(data,mean)
}
colnames(data) <- unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])
data <- log2(data+1)
data[data==0] <- NA
data <- data[,c("Norgen","Lexogen","QIAseq","QIAseq_UMI","NEXTflex","RealSeq","SMARTer","EdgeSeq")]#Order columns
data_substr <- data[rownames(data) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
corrm_all <- Hmisc::rcorr(as.matrix(data), type = "pearson")
corrm_substr <- Hmisc::rcorr(as.matrix(data_substr), type = "pearson")
corrm <- as.data.frame(corrm_all$r)
corrm_substr <- as.data.frame(corrm_substr$r)
corrm$EdgeSeq <- corrm_substr$EdgeSeq
corrm[8,] <- corrm_substr[8,]
corrm <- as.matrix(corrm)
col1 <- colorRampPalette(brewer.pal(11, "RdBu"))
pdf(height = 7, width = 7, file = "outs/miRXplore_between.kit.correlation_log2counts_nonzeroPairs.pearson.pdf")
corrplot(corrm, type ="upper", method = "color", order = "original", col = rev(col1(40)),
         tl.col = "black", tl.srt = 45, outline = T, diag = F, addCoef.col = T)#Correlations with edgeseq are based only on targeted miRNAs, but correlations between all other kits are calculated based on full set of 962 miRNAs
dev.off()

#'**Replicate correlation**
#+ replicate correlation,cache=FALSE, results='show', echo=FALSE
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  if (i =="EdgeSeq"){
    df <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
    df <- df[rownames(df) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
  }
  else {
    df <-as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
  }
df = log2(df+1)
colnames(df) = c("Replicate_1", "Replicate_2")
pdf(height = 2.5, width = 2.5, file = paste0("outs/",i, "replicate_correlation_miRXplore.pdf"))
print(ggscatter(df, x = "Replicate_1", y="Replicate_2",
          add = "reg.line", conf.int = TRUE, cor.coef = T,
          color = "#2166ac", size = 1,
          add.params = list(color="#8f0021", fill="lightgray", size = 0.5),
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
          title = i, ggtheme = theme_bw()) + rremove("grid"))
dev.off()
}
#'**10 most abundant miRNA**
#+ mtop10_miRNAs,cache=FALSE, results='show', echo=FALSE
ind <- 0
export_top10 <- c("Reads_replicate_1", "Reads_replicate_2", "mean","Protocol")
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  ind <- ind +1
  if (i =="EdgeSeq"){
    df <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
    df <- df[rownames(df) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
  }
  else {
  df <-as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
  }
  df$Mean <- rowMeans(df)
  df<-df[order(df[,3],decreasing=TRUE),]
  df<-df[1:10,]
  miRNA_names <- factor(rownames(df), levels = rownames(df))
  reads_top10 <- c(colSums(miRXplore[rownames(df),annot_miRXplore$Kit == i]),
                   mean(colSums(miRXplore[rownames(df),annot_miRXplore$Kit == i])),i)
  export_top10 <- rbind(export_top10,reads_top10)
  df$miRNA <- miRNA_names
  df<-df[,-c(1,2)]
  colnames(df)<-c("Mean","miRNA")
  df$Percentage_of_normalized_counts <- (as.numeric(as.character(df$Mean/1000000)))*100
  assign(x=paste0("p",ind),value=p<-ggplot(data=df, aes(x=miRNA,y=Percentage_of_normalized_counts))+
           geom_bar(stat="identity", position="dodge",color = "lightsteelblue2", fill = "lightsteelblue2")+
           coord_cartesian(xlim =c(1, 10), ylim = c(0, 10))+
           theme(axis.text.x = element_text(angle = 45, hjust = 1))+
           xlab("")+
           ylab("% normalized counts")+
           ggtitle(i))
}
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow = 3)
colnames(export_top10) <- export_top10[1,]
export_top10 <- as.data.frame(export_top10[-1,])
rownames(export_top10) <- export_top10$Protocol
write.csv(export_top10, file = "outs/top10_reads_miRXplore.csv")

#'**Cumulative frequency of miRXplore**
#+ miRNA-cumulative frequency,cache=FALSE, results='show', echo=FALSE
ind <- 0
CF<-as.data.frame(seq(1,962,1))
colnames(CF)<-"Ranking"
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  df <- as.data.frame((rowMeans(miRXplore_norm[,annot_miRXplore$Kit == i])/1000000))
  colnames(df)<-"Frequency"
  df$miRNA_names <- rownames(df)
  df<-df[order(df[,1],decreasing = FALSE),]
  df$Cumulative_frequency<-cumsum(df$Frequency)
  CF <- cbind(CF,df$Cumulative_frequency)
}
CF <- CF[,-1]
colnames(CF)<-unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])
CF<-melt(CF, value.name = "Cumulative_frequency")
CF1 <- CF[(CF$variable != "EdgeSeq"),]
CF2 <- CF[(CF$variable == "EdgeSeq"),]
CF3 <- CF2[496:962,]
CF <- rbind(CF3,CF1)
rm(CF1,CF2,CF3)
miRNA_rank_seq <-rep(seq(1,962,1),length(unique(annot_miRXplore$Kit))-1)
miRNA_rank_EdgeSeq <-seq(1,length(targeted_in_miRXplore_EdgeSeq$miRNA_name),1)
miRNA_rank <- c(miRNA_rank_seq,miRNA_rank_EdgeSeq)
CF$miRNA_rank <- miRNA_rank
#+ ranking without zero count miRNAs,cache=FALSE, results='show', echo=FALSE
CF <- CF[CF$Cumulative_frequency !=0,]
miRNA_rank <- c()
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  miRNA_rank<-c(miRNA_rank,dim(CF[CF$variable == i, ])[1])
}
miRNA_rank <- lapply(miRNA_rank,function(x) seq(1,x, by = 1))
CF$miRNA_rank <- unlist(miRNA_rank)
#+ percentile miRNAs,cache=FALSE, results='show', echo=FALSE
percent_miRNA <- c()
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  number_of_miRNA <- dim(CF[CF$variable == i,])[1]
  percent <- as.vector(sapply(unlist(CF[CF$variable == i,"miRNA_rank"]), function(s){(s/number_of_miRNA)*100}))
  percent_miRNA <- c(percent_miRNA, percent)
}
CF$percent_miRNA <- percent_miRNA
CF$variable <- factor(CF$variable,level = c("Norgen","Lexogen","QIAseq","QIAseq_UMI","NEXTflex","RealSeq","SMARTer","EdgeSeq"))
theme_set(theme_bw())
plot_CF_miRXplore <-ggplot(CF, aes(x=percent_miRNA, y=Cumulative_frequency, color=variable))+
  geom_line(size=1.5)+
  xlab("miRNA percentile")+
  ylab("Cumulative frequency")+
  ggtitle("Linear cumulative frequency of miRXplore")+
  geom_hline(yintercept = 0.5,linetype = "dashed" )+
  theme(panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#cccccc","#737c82","#b2d4e6","#379ed2","#063a74","#ffc5ad","#e56638", "#8f0021"))
plot_CF_miRXplore
ggsave(file = "outs/cumulative.frequency_miRXplore.pdf", device = "pdf", height = 9, width = 12.5, units = "cm", useDingbats = F)
rm(CF,miRNA_rank_seq,miRNA_rank_EdgeSeq,miRNA_rank)

#'**Number of detected miRNA, dependence on raw counts**
#+ sensitivity of detection (to raw counts),results='show', echo=FALSE
raw_counts <- t(as.data.frame(read.table("src/raw_count_star", header = TRUE, sep = "\t")))
raw_counts <- as.data.frame(raw_counts[grep("miRXplore", rownames(raw_counts), ignore.case = TRUE),])
colnames(raw_counts)<-"raw_counts"
no_miRNA <- colSums((miRXplore > 5)== TRUE)
no_dropouts <- colSums((miRXplore > 5)== FALSE)
sensitivity <- as.data.frame(no_miRNA)
mean_sensitivity <- c()
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  if (i =="EdgeSeq"){
    df <- as.data.frame(miRXplore[,annot_miRXplore$Kit == i])
    df <- df[rownames(df) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
  }
  else {
    df <-miRXplore[,annot_miRXplore$Kit == i]
  }
  df1 <- df[(df[,1])>5,]
  df2 <- df[(df[,2])>5,]
  number_of_detected <- length(intersect(rownames(df1),rownames(df2)))
  mean_sensitivity <-c(mean_sensitivity,as.numeric(number_of_detected))
}
names(mean_sensitivity) <- unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])
mean_sensitivity <- as.data.frame(mean_sensitivity)
sensitivity 
mean_sensitivity

#+ downsampling,  results='hide',include=FALSE
proportions <- 10^seq(-3,0,0.01)
for (proportion in proportions){
  downsampled_matrix <- as.data.frame(apply(miRXplore, 2, function(x) rbinom(nrow(miRXplore), x, proportion)))
  rownames(downsampled_matrix) <- rownames(miRXplore)
  mean_nomiRNA <- c()
  for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
    if (i =="EdgeSeq"){
      df <- as.data.frame(downsampled_matrix[,annot_miRXplore$Kit == i])
      df <- df[rownames(df) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
    }
    else {
      df <-downsampled_matrix[,annot_miRXplore$Kit == i]
    }
    df1 <- df[(df[,1])>5,]
    df2 <- df[(df[,2])>5,]
    number_of_detected <- length(intersect(rownames(df1),rownames(df2)))
    mean_nomiRNA <-c(mean_nomiRNA,as.numeric(number_of_detected))
  }
  mean_nomiRNA <- as.data.frame(mean_nomiRNA)
  colnames(mean_nomiRNA) <- proportion
  mean_sensitivity<-cbind(mean_sensitivity, mean_nomiRNA)
}
mean_sensitivity <- mean_sensitivity[,-1]
raw_counts_downsampled <- c("raw_counts","Protocol")
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  print(i)
  print(as.vector(raw_counts[annot_miRXplore$Kit == i,]))
  mean_of_samples <-mean(as.vector(raw_counts[rownames(annot_miRXplore[annot_miRXplore$Kit == i,]),]))
  raw_counts_down <- as.numeric(as.numeric(mean_of_samples)* as.numeric(proportions))
  raw_counts_down <- cbind(round(raw_counts_down),rep(i,length(raw_counts_down)))
  raw_counts_downsampled <- rbind(raw_counts_downsampled,raw_counts_down)
}
mean_sensitivity$Sample <- factor(rownames(mean_sensitivity), level =c("Lexogen","Norgen","QIAseq",
                                                                       "QIAseq_UMI","NEXTflex","NEXTflex_UMI",
                                                                       "RealSeq","SMARTer","EdgeSeq"))
mean_sensitivity <- melt(mean_sensitivity, value.name = "no_miRNA")
mean_sensitivity <- mean_sensitivity[order(mean_sensitivity$Sample),]
colnames(raw_counts_downsampled) <- raw_counts_downsampled[1,]
raw_counts_downsampled <- as.data.frame(raw_counts_downsampled[-1,])
mean_sensitivity$raw_counts_downsampled <- as.numeric(as.character(raw_counts_downsampled$raw_counts))
mean_sensitivity$Protocol <- raw_counts_downsampled$Protocol
mean_sensitivity <- mean_sensitivity[!(mean_sensitivity$Protocol %in% c("QIAseq_UMI","NEXTflex_UMI")),]
plot_S_miRXplore <- ggplot(data=mean_sensitivity, aes(x=raw_counts_downsampled, y=no_miRNA, colour=Sample)) + 
  geom_line(size=2)+
  xlab("Number of raw reads")+
  ggtitle("miRXplore")+
  ylab("Number of detected miRNAs")+
  geom_hline(yintercept = 962)+
  geom_hline(yintercept = 467,linetype = "dashed" )+
  scale_color_brewer(palette="Set3")
legend <- cowplot::get_legend(plot_S_miRXplore)
  grid.newpage()
  grid.draw(legend)
plot_S_miRXplore <- plot_S_miRXplore + theme(legend.position = "none")
plot_S_miRXplore
rm(raw_counts,no_miRNA,no_dropouts,mean_sensitivity,raw_counts_downsampled,raw_counts_down, mean_nomiRNA,
   mean_of_samples, sensitivity)

#'**Ligation and PCR bias**
#+ miRNA ligation,cache=FALSE, results='show', echo=FALSE
df_bias <- as.data.frame(t(c(1,2,3)))
colnames(df_bias)<-c("log2fc","miRNA","Protocol")
df_bias <- df_bias[-1,]
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  if (i =="EdgeSeq"){
    predicted <- 1000000/dim(targeted_in_miRXplore_EdgeSeq)[1]
    df <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
    df <- df[rownames(df) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
  }
  else {
  predicted <-1000000/962
  df <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
  }
  df$mean <- rowMeans(df)
  df$log2fc <-log2(df$mean/predicted)
  df$miRNA <- unlist(rownames(df), use.names=FALSE)
  df$Protocol <- rep(paste0("",i),dim(df)[1])
  df <- df[,-c(1,2,3)]
  df_bias <- rbind(df_bias,df)
}
df_bias_plot <- as.data.frame(df_bias[df_bias$log2fc != "-Inf",])
df_bias_plot$Protocol <- as.factor(df_bias_plot$Protocol)
df_bias_plot$Protocol <- factor(df_bias_plot$Protocol,level = c("Norgen","Lexogen","QIAseq","QIAseq_UMI","NEXTflex","RealSeq","SMARTer","EdgeSeq"))
level_order <- levels(df_bias_plot$Protocol)
theme_set(theme_bw())
plot_ligation_bias <- ggplot(data=df_bias_plot, aes(x=log2fc,y=Protocol, fill=Protocol))+
  geom_density_ridges(alpha=I(1))+
  ggtitle("Technical bias")+
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank())+
  scale_x_continuous(breaks = c(-10,-5,-2,0,2,5,10))+
  geom_vline(aes(xintercept = -1),linetype = "dashed", color = "292929")+
  geom_vline(aes(xintercept = 1),linetype = "dashed", color = "292929")+
  xlab("Deviation from expected (log2)")+
  scale_fill_manual(values = c("#cccccc","#737c82","#b2d4e6","#379ed2","#063a74","#ffc5ad","#e56638", "#8f0021"))
plot_ligation_bias
ggsave(file= "outs/Technical.bias_miRXplore.pdf", device = "pdf", height = 13, width = 12, units = "cm")
 rm(predicted, df_bias_plot)

#'**Ligation and PCR bias correlation with GC content**
#+ miRNA ligation bias,cache=FALSE, results='show', echo=FALSE
gc_content <-c()
miRXplore_seuqences <- read.table("src/miRXplore_annot.tab", header = TRUE, sep = "\t")
rownames(miRXplore_seuqences) <- miRXplore_seuqences$Description
for (i in rownames(miRXplore_seuqences)){
  seq <- miRXplore_seuqences[i,"Sequence"]
  C <- stringr::str_count(seq, c("C"))
  G <- stringr::str_count(seq, c("G"))
  GC <- C + G
  GC_content <-round((GC /stringr::str_length(seq))*100,2)
  gc_content <- c(gc_content, GC_content)
}
miRXplore_seuqences$GC_content <- gc_content
rm(seq,C,G,GC,GC_content,gc_content)
df_log2fc <- as.data.frame(t(c(1,2,3,4,5,6)))
colnames(df_log2fc) <- c("log2fc", "miRNA", "Protocol", "Description","Sequence", "GC_content")
df_log2fc <- df_log2fc[-1,]
ind <- 0
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  ind <- ind + 1 
  df <- df_bias[df_bias$Protocol == i,]
  df <- cbind(df,miRXplore_seuqences[rownames(miRXplore_seuqences) %in% df$miRNA,])
  df_log2fc <- rbind(df_log2fc,df)
  df <- df[df$log2fc != "-Inf",]
  lm1 <-lm(df$log2fc ~ df$GC_content)
  intercept1 <- lm1$coefficients[1]
  slope1 <- lm1$coefficients[2]
  r2 <- round(summary(lm1)$r.squared,2)
  rsq_label <- paste('R^2 == ', r2)
  assign(x=paste0("p",ind),value=p<-ggplot(data=df, aes(x=GC_content,y=log2fc))+
           geom_point(color="plum4", alpha=0.5)+
           ylab("log2fc")+
           xlab("GC content %")+
           ggtitle(i)+
           xlim(0,100)+
           ylim(-15,10)+
           #geom_text(x=90, y=6, label=paste0("R2=",r2))+
           annotate(geom="text", x=90, y=8, label=rsq_label,color="gray5", parse = TRUE)+
           geom_abline(intercept = intercept1, slope = slope1))
}
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow=2)
df_log2fc <- df_log2fc[,-2]
write.csv(df_log2fc, file = "outs/log2cf_GC_sequence_by_kit.csv")
kits <- unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])
rm(df_log2fc)

#'**Number of unbiased miRNAs**
#+ unbiased miRNAs,cache=FALSE, results='show', echo=FALSE
df_number_of_biased <- as.data.frame(t(c(1,2,3,4,5,6,7,8)))
#df_number_of_biased <- df_number_of_biased[-1,]
for (i in unique(annot_miRXplore[order(annot_miRXplore$Kit),"Kit"])){
  print(i)
  n_unbiased <- dim(df_bias[(df_bias$log2fc < 1) & (df_bias$log2fc > -1) & (df_bias$Protocol == i),])[1]
  n_under <- dim(df_bias[(df_bias$log2fc < -1) & (df_bias$Protocol == i),])[1]
  n_over <- dim(df_bias[(df_bias$log2fc > 1) & (df_bias$Protocol == i),])[1]
  p_unbiased <- (n_unbiased/dim(df_bias[df_bias$Protocol == i,])[1])*100
  p_under <- (n_under/dim(df_bias[df_bias$Protocol == i,])[1])*100
  p_over <- (n_over/dim(df_bias[df_bias$Protocol == i,])[1])*100
  df <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit == i])
  df$mean <- rowMeans(df)
  df1 <- df_bias[(df_bias$Protocol == i) & (df_bias$log2fc < 1) & (df_bias$log2fc > -1),]
  df2 <- df[rownames(df) %in% df1$miRNA,]
  p_r_unbiased <- (sum(df2$mean)/sum(df$mean))*100
  df_number_of_biased <- rbind(df_number_of_biased,c(i,n_unbiased,n_under,n_over,p_unbiased,
                                                     p_under,p_over, p_r_unbiased))
}
colnames(df_number_of_biased) <- c("Protocol", "Unbiased", "Under-estimated","Over-estimated",
                                   "Unbiased%", "Under-estimated%","Over-estimated%", "Unbiased_reads")
df_number_of_biased <- df_number_of_biased[-1,]
write.csv(df_number_of_biased,file = "outs/ligation_bias.csv")
rm(n_unbiased,n_under,n_over, df_number_of_biased,p_unbiased,p_under,p_over, p_r_unbiased)

#'**QIAseq expression vs UMI**
#+ UMI vs expression,cache=FALSE, results='show', echo=FALSE
qiagen_means <- c()
for (i in c("QIAseq","QIAseq_UMI")){qiagen_means <- cbind(qiagen_means,rowMeans(miRXplore_norm[,rownames(annot_miRXplore[annot_miRXplore$Kit == i,])]))}
colnames(qiagen_means) <- c("QIAseq","QIAseq_UMI")
qiagen_means <-as.data.frame(qiagen_means)
ggplot(data = qiagen_means,aes(x=QIAseq_UMI, y=QIAseq))+
  geom_point()

#'**NEXTflex expression vs UMI**
#+ UMI vs expression nextflex,cache=FALSE, results='show', echo=FALSE
load("temp/miRXplore_norm.RData")
load("temp/annot_miRXplore.RData")
nextflex_means <- c()
for (i in c("NEXTflex","NEXTflex_UMI")){nextflex_means <- cbind(nextflex_means,rowMeans(miRXplore_norm[,rownames(annot_miRXplore[annot_miRXplore$Kit == i,])]))}
nextflex_means <- as.data.frame(nextflex_means)
colnames(nextflex_means) <- c("NEXTflex","NEXTflex_UMI")
ggplot(data = nextflex_means,aes(x=NEXTflex_UMI, y=NEXTflex))+
  geom_point()

#'**Heatmap clustering of kits and miRNA**
#+ heatmap,cache=FALSE, results='show', echo=FALSE  
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
scaled <- as.data.frame(t(apply(miRXplore_norm, 1, cal_z_score)))
scaled<-na.omit(scaled)
distance_matrix<-dist(scaled)
clustering<-hclust(distance_matrix)
plot(clustering, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = NULL, ylab = "Height")
top<-scaled[head(clustering$order,15),]
bottom<-scaled[tail(clustering$order,15),]

#'15 most similar miRNA (similar in expression values) 
#+ heatmap_clust,cache=FALSE, results='show', echo=FALSE
pheatmap(top)
pheatmap(bottom)

#'15 miRNA with biggest difference between minimal and maximal expression value 
#+ heatmap_max_diff,cache=FALSE, results='show', echo=FALSE
ordered<-scaled[order(abs(rowMaxs(as.matrix(scaled))-rowMins(as.matrix(scaled)))),]
pheatmap(head(ordered,15))

#'15 miRNA with lowest difference between minimal and maximal expression value 
#+ heatmap_min_diff,cache=FALSE, results='show', echo=FALSE
pheatmap(tail(ordered,15))

#'15 miRNA with lowest padj (deseq2, LRT - tested for affect of kit)
#+ heatmap_deseq,cache=FALSE, results='show', echo=FALSE
load("temp/miRXplore.RData")
miRXplore_cluster <- miRXplore[rownames(targeted_in_miRXplore_EdgeSeq),]
de<-DESeqDataSetFromMatrix(miRXplore_cluster, colData = annot_miRXplore, design = ~ Kit)
res<-results(DESeq(de, test= "LRT", reduced = ~ 1))
exp_val<-as.data.frame(assay(de))
exp_val<-as.data.frame((sweep(exp_val, 2, colSums(exp_val), "/"))*1000000)
scaled <- as.data.frame(t(apply(exp_val, 1, cal_z_score)))
scaled$padj <-res$padj
scaled<-na.omit(scaled)
scaled<-scaled[order(scaled$padj, decreasing = FALSE),]
top <- head(scaled,15)
top <- top[,colnames(top) != "padj"]
miRXplore_heatmap<-pheatmap(top)
miRXplore_seuqences <- read.table("src/miRXplore_annot.tab", header = TRUE, sep = "\t")
rownames(miRXplore_seuqences) <- miRXplore_seuqences[,1]
colnames(miRXplore_seuqences) <- c("miRNA", "Sequence")
miRXplore_top_15 <- c("miRNA", "Sequence", "GC_content")
for (i in rownames(top)){
  seq<-as.character(miRXplore_seuqences[miRXplore_seuqences$miRNA ==i ,2])
  print(seq)
  C <- stringr::str_count(seq, c("C"))
  G <- stringr::str_count(seq, c("G"))
  GC <- C + G
  GC_content <-round((GC /stringr::str_length(seq))*100,2)
  miRXplore_top_15 <- rbind(miRXplore_top_15, c(i, seq, GC_content))
}
write.csv(miRXplore_top_15,"outs/miRXplore_top15_sequence.csv")

#'PCA on normalized expression - variance stabilizing normalization.
#+ pca,cache=FALSE, results='show', echo=FALSE
data <- miRXplore[rownames(miRXplore) %in% targeted_in_miRXplore_EdgeSeq$miRNA_name,]
de<-DESeqDataSetFromMatrix(data, colData = annot_miRXplore, design = ~ Kit)
vsd<-vst(de, nsub = 200, fitType = "local")
plotPCA(vsd, intgroup = "Kit")

#'**Linear model - explained variance by PCR bias and by ligation bias - QIAseq only**
#+ explained variance,cache=FALSE, results='show', echo=FALSE
#'Mixed (with random effects) linear model was used on normalized data.
qiagen <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit %in% c("QIAseq", "QIAseq_UMI")])
qiagen$PredictedQIAseq_miRXplore_1 <- rep(1000000/962,962)
qiagen$PredictedQIAseq_miRXplore_2 <- rep(1000000/962,962)
annot_qiagen <- read.table("src/annot_qiagen_star.tab")
annot_qiagen$Replicate <- as.factor(annot_qiagen[,"Replicate"])
colnames(annot_qiagen) <- c("Protocol", "Sample_source","Replicate", "PCR_bias", "Ligation_bias")
qiagen <- qiagen[(qiagen$QIAseq_miRXplore_1 > 5) & (qiagen$QIAseq_miRXplore_2 > 5),]
form <- ~ (1|PCR_bias) + (1|Ligation_bias) + (1|Replicate)
varPart <- fitExtractVarPartModel( qiagen, form, annot_qiagen )
vp <- sortCols(varPart)
plotPercentBars(vp[1:10,])
plot_var<-plotVarPart(vp, col = c("#b2182b","#2166ac","#b2d4e6","#f4a582")) + ggtitle("Explained variance in miRXplore (QIAseq data)")+
  scale_x_discrete(limits=c("Ligation_bias","Residuals","Replicate","PCR_bias"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plot_var
ggsave(file = "outs/explained.variance_QIAseq.data.pdf", device = "pdf", height = 13, width = 14, units = "cm")
head(varPart)
vpSummaries <- fitVarPartModel( qiagen, form, annot_qiagen, fxn=summary )
vpSummaries[[2]]
form <- ~ PCR_bias + Ligation_bias + Replicate
C = canCorPairs( form, annot_qiagen)
plotCorrMatrix( C )
save(varPart, file = "outs/varPart_qiagen.RData")

#correlation of PC-Rbiased miRNAs with GC content
GC_content <- read.table("outs/log2cf_GC_sequence_by_kit.csv", header = TRUE, sep = ",")
GC_content <- GC_content[GC_content$Protocol == "QIAseq",]
GC_content <- GC_content[GC_content$Description %in% rownames(varPart),]
rownames(GC_content)<-GC_content$Description
GC_content$PCR_bias <- varPart$PCR_bias
lm1 <-lm(GC_content$PCR_bias ~ GC_content$GC_content)
intercept1 <- lm1$coefficients[1]
slope1 <- lm1$coefficients[2]
r2 <- round(summary(lm1)$r.squared,2)
rsq_label <- paste('R^2 == ', r2)
ggplot(data = GC_content, aes(x=GC_content, y=PCR_bias))+
  geom_point()+
  geom_abline(intercept = intercept1, slope = slope1)+
  annotate(geom="text", x=90, y=0.7, label=rsq_label,color="gray5", parse = TRUE)
cor.test(GC_content$PCR_bias,GC_content$GC_content)

#'**Linear model - explained variance by PCR bias and by ligation bias - NEXTflex only**
#+ explained variance 2,cache=FALSE, results='show', echo=FALSE
#'Mixed (with random effects) linear model was used on normalized data.
nextflex <- as.data.frame(miRXplore_norm[,annot_miRXplore$Kit %in% c("NEXTflex", "NEXTflex_UMI")])
nextflex$PredictedNEXTflex_miRXplore_1 <- rep(1000000/962,962)
nextflex$PredictedNEXTflex_miRXplore_2 <- rep(1000000/962,962)
annot_nextflex <- read.table("src/annot_nextflex_star.tab")
annot_nextflex$Replicate <- as.factor(annot_nextflex[,"Replicate"])
colnames(annot_nextflex) <- c("Protocol", "Sample_source","Replicate", "PCR_bias", "Ligation_bias")
nextflex <- nextflex[(nextflex$NEXTflex_miRXplore_1 > 5) & (nextflex$NEXTflex_miRXplore_2 > 5),]
form <- ~ (1|PCR_bias) + (1|Ligation_bias) + (1|Replicate)
varPart <- fitExtractVarPartModel( nextflex, form, annot_nextflex )
vp <- sortCols(varPart)
plotPercentBars(vp[1:10,])
plot_var<-plotVarPart(vp) + ggtitle("Explained variance in miRXplore (NEXTflex data)")
plot_var
head(varPart)
vpSummaries <- fitVarPartModel( nextflex, form, annot_nextflex, fxn=summary )
vpSummaries[[2]]
form <- ~ PCR_bias + Ligation_bias + Replicate
C = canCorPairs( form, annot_nextflex)
plotCorrMatrix( C )
save(varPart, file = "outs/varPart_nextflex.RData")

#correlation of PCRbiased miRNAs with GC content
GC_content <- read.table("outs/log2cf_GC_sequence_by_kit.csv", header = TRUE, sep = ",")
GC_content <- GC_content[GC_content$Protocol == "NEXTflex",]
GC_content <- GC_content[GC_content$Description %in% rownames(varPart),]
rownames(GC_content)<-GC_content$X
GC_content$PCR_bias <- varPart$PCR_bias
lm1 <-lm(GC_content$PCR_bias ~ GC_content$GC_content)
intercept1 <- lm1$coefficients[1]
slope1 <- lm1$coefficients[2]
r2 <- round(summary(lm1)$r.squared,2)
rsq_label <- paste('R^2 == ', r2)
ggplot(data = GC_content, aes(x=GC_content, y=PCR_bias))+
  geom_point()+
  geom_abline(intercept = intercept1, slope = slope1)+
  annotate(geom="text", x=90, y=0.7, label=rsq_label,color="gray5", parse = TRUE)
cor.test(GC_content$PCR_bias,GC_content$GC_content)


sessionInfo()
