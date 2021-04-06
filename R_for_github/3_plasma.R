#' ---
#' title: "Plasma"
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
library(GGally)
library(cowplot)
library(RColorBrewer)
library(ggpubr)

#+ load data, include=FALSE
load("temp/plasma.RData")
load("temp/plasma_norm.RData")
load("temp/annot_plasma.RData")
targeted_in_plasma_EdgeSeq <- as.data.frame(read.table("src/plasma_targeted_EdgeSeq.tab", header = TRUE, sep = "\t"))
rownames(targeted_in_plasma_EdgeSeq) <- targeted_in_plasma_EdgeSeq$miRNA_name
annot_plasma <- annot_plasma[annot_plasma$Kit != "NEXTflex_UMI",]
annot_plasma <- droplevels(annot_plasma)
plasma <- plasma[,rownames(annot_plasma)]
plasma_norm <- plasma_norm[,rownames(annot_plasma)]


#'**Plasma total counts per sample.**
#+ count_tables,cache=FALSE, results='show', echo=FALSE
Total_counts <- c("Protocol","Mean_total_counts")
for (i in unique(annot_plasma$Kit)){
  Total_counts<-rbind(Total_counts,c(i,mean(colSums(plasma[,annot_plasma$Kit == i]))))
}
Total_counts
write.csv(plasma,file = "outs/plasma_raw_counts.csv")
rm(Total_counts)

#'**Plasma coefficient of variation.**
#'Mean coefficient of variation (standard deviation/mean) for two replicates.
#+ cv,cache=FALSE, results='show', echo=FALSE
cv<-c("Protocol","Mean_CV")
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  if (i =="EdgeSeq"){
    df <- as.data.frame(plasma_norm[,annot_plasma$Kit == i])
    df <- df[rownames(df) %in% targeted_in_plasma_EdgeSeq$miRNA_name,]
  }
  else {
    df <-plasma_norm[,annot_plasma$Kit == i]
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
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  mean <- rowMeans(plasma_norm[,annot_plasma$Kit == i])
  data <-cbind(data,mean)
}
colnames(data) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
data <- as.data.frame(log2(data+1))
data[data==0] <- NA
data <- data[c("Norgen","Lexogen","QIAseq","QIAseq_UMI","NEXTflex","RealSeq","SMARTer","EdgeSeq")]#Order columns
data_substr <- data[rownames(data) %in% targeted_in_plasma_EdgeSeq$miRNA_name,]
corrm_all <- Hmisc::rcorr(as.matrix(data), type = "pearson")
corrm_substr <- Hmisc::rcorr(as.matrix(data_substr), type = "pearson")
corrm <- as.data.frame(corrm_all$r)
corrm_substr <- as.data.frame(corrm_substr$r)
corrm$EdgeSeq <- corrm_substr$EdgeSeq
corrm[8,] <- corrm_substr[8,]
corrm <- as.matrix(corrm)
col1 <- colorRampPalette(brewer.pal(11, "RdBu"))
pdf(height = 7, width = 7, file = "outs/plasma_between.kit.correlation_log2counts_nonzeroPairs.pearson.pdf")
corrplot(corrm, type ="upper", method = "color", order = "original", col = rev(col1(40)),
         tl.col = "black", tl.srt = 45, outline = T, diag = F, addCoef.col = T)#Correlations with edgeseq are based only on targeted miRNAs, but correlations between all other kits are calculated based on full set of 962 miRNAs
dev.off()
data$miRNA <- rownames(data)
data <- melt(data)
g <- ggplot(data = data, x = miRNA, y = value, color = variable) + geom_line(aes(x = miRNA, y = value, color = variable, group = variable))
g
#Replicate correlation
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  if (i =="EdgeSeq"){
    df <- as.data.frame(plasma_norm[,annot_plasma$Kit == i])
    df <- df[rownames(df) %in% targeted_in_plasma_EdgeSeq$miRNA_name,]
  }
  else {
    df <-as.data.frame(plasma_norm[,annot_plasma$Kit == i])
  }
  df <- log2(df+1)
  colnames(df) = c("Replicate_1", "Replicate_2")
  pdf(height = 2.5, width = 2.5, file = paste0("outs/",i, "replicate_correlation_plasma.pdf"))
  print(ggscatter(df, x = "Replicate_1", y="Replicate_2",
                  add = "reg.line", conf.int = TRUE, cor.coef = T,
                  color = "#2166ac", size = 1,
                  add.params = list(color="#8f0021", fill="lightgray", size = 0.5),
                  cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                  title = i, ggtheme = theme_bw()) + rremove("grid"))
  dev.off()
}

#'**QIAseq expression vs UMI**
#+ UMI vs expression,cache=FALSE, results='show', echo=FALSE
qiagen_means <- c()
for (i in c("QIAseq","QIAseq_UMI")){qiagen_means <- cbind(qiagen_means,rowMeans(plasma_norm[,rownames(annot_plasma[annot_plasma$Kit == i,])]))}
colnames(qiagen_means) <- c("QIAseq","QIAseq_UMI")
qiagen_means <-as.data.frame(qiagen_means)
qiagen_means <- qiagen_means[(qiagen_means$QIAseq > 0) & (qiagen_means$QIAseq < 200000) ,]
ggplot(data = qiagen_means,aes(x=QIAseq_UMI, y=QIAseq))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()


#'**10 most abundant miRNA and Ligation bias**
#+ mtop10_miRNAs,cache=FALSE, results='show', echo=FALSE
ind <- 0
export_top10_names <- as.data.frame(c(1,2,3,4,5,6,7,8,9,10))
export_top10 <- c("Reads_replicate_1", "Reads_replicate_2", "mean","Protocol")
log2fc <- read.table("outs/log2cf_GC_sequence_by_kit.csv", header = TRUE, sep = ",")
miRXplore_plasma_names <- read.table("src/mirXplore_mature.tab", header = TRUE)
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  ind <- ind +1
  if (i =="EdgeSeq"){
    df <- as.data.frame(plasma_norm[,annot_plasma$Kit == i])
    df <- df[rownames(df) %in% targeted_in_plasma_EdgeSeq$miRNA_name,]
  }
  else {
    df <-as.data.frame(plasma_norm[,annot_plasma$Kit == i])
  }
  df$mean <- rowMeans(df)
  df<-df[order(df[,3],decreasing=TRUE),]
  df<-df[1:10,]
  miRNA_names <- factor(rownames(df), levels = rownames(df))
  reads_top10 <- c(colSums(plasma[rownames(df),annot_plasma$Kit == i]),
                   mean(colSums(plasma[rownames(df),annot_plasma$Kit == i])),i)
  export_top10 <- rbind(export_top10,reads_top10)
  df$miRNA <- miRNA_names
  export_top10_names <-cbind(export_top10_names,df$miRNA)
  df<-df[,-c(1,2)]
  colnames(df)<-c("Mean","miRNA")
  log2fc_vector <- c()
  for (miRNA in df$miRNA){
    if (miRNA %in% miRXplore_plasma_names$name_mature){
      miRXplore_name <- as.character(miRXplore_plasma_names[miRXplore_plasma_names$name_mature == miRNA, "name_mirxplore"])
      log2fc_value <- log2fc[(log2fc$Protocol == i) & (log2fc$Description == miRXplore_name),"log2fc"]
      log2fc_vector <-c(log2fc_vector,log2fc_value)
    }
    else{
      log2fc_vector<-c(log2fc_vector,"NA")
    }
  }
  df$log2fc <- as.numeric(as.character(log2fc_vector))
  df$Percentage_of_normalized_counts <- (df$Mean/1000000)*100
  x1<-ggplot(data=df, aes(x=miRNA,y=Percentage_of_normalized_counts))+
           geom_bar(stat="identity", position="dodge", color = "lightsteelblue2", fill = "lightsteelblue2")+
           coord_cartesian(xlim =c(1, 10), ylim = c(0, 100))+
           theme(axis.text.x=element_blank())+
           ylab("% of normalized counts")+
           xlab("")+
           ggtitle(i)
  x2<-ggplot(data=df, aes(x=miRNA,y=log2fc))+
    geom_point(color = "royalblue4", size = 2)+
    scale_fill_brewer(palette = "Blues")+
    theme_half_open()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    coord_cartesian(xlim =c(1, 10), ylim = c(-6, 6))+
    scale_y_continuous(position = "right", breaks = round(seq(-6, 6, by = 2),0))+
    geom_hline(aes(yintercept = 0),linetype = "solid", size = 0.3)+
    geom_hline(aes(yintercept = 1),linetype = "dashed",size = 0.3)+
    geom_hline(aes(yintercept = -1),linetype = "dashed",size = 0.3)+
    theme(legend.position = "none")+
    ylab("log2fc")+
    xlab("")
  aligned_plots<-align_plots(x1, x2, align="hv", axis="tblr")
  assign(x=paste0("p",ind),value=p<-ggdraw(aligned_plots[[1]])+draw_plot(aligned_plots[[2]]))
}
colnames(export_top10) <- export_top10[1,]
export_top10 <- as.data.frame(export_top10[-1,])
rownames(export_top10) <- export_top10$Protocol
write.csv(export_top10, file = "outs/top10_reads_plasma.csv")
export_top10_names <- export_top10_names[,-1]
colnames(export_top10_names) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
write.csv(export_top10_names, file = "outs/top10_names_plasma.csv")
g<-grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow = 4)
ggsave(file= "outs/Top10_vs_bias.pdf",plot = g, device = "pdf", height = 40, width = 23, units = "cm")
rm(log2fc,log2fc_value,log2fc_vector,miRNA,miRNA_names,miRXplore_name,p1,p2,p3,p4,p5,p6,p7,p8,
   miRXplore_plasma_names, ind,x1,x2, aligned_plots,export_top10,reads_top10)

#'**Cumulative frequency of plasma**
#+ miRNA-cumulative frequency,cache=FALSE, results='show', echo=FALSE
ind <- 0
CF<-as.data.frame(seq(1,2656,1))
colnames(CF)<-"Ranking"
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  df <- as.data.frame((rowMeans(plasma_norm[,annot_plasma$Kit == i])/1000000))
  colnames(df)<-"Frequency"
  df$miRNA_names <- rownames(df)
  df<-df[order(df[,1],decreasing = FALSE),]
  df$Cumulative_frequency<-cumsum(df$Frequency)
  CF <- cbind(CF,df$Cumulative_frequency)
}
CF <- CF[,-1]
colnames(CF)<-unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
CF<-melt(CF, value.name = "Cumulative_frequency")
CF1 <- CF[(CF$variable != "EdgeSeq"),]
CF2 <- CF[(CF$variable == "EdgeSeq"),]
CF3 <- CF2[593:2656,]
CF <- rbind(CF3,CF1)
rm(CF1,CF2,CF3)
miRNA_rank_seq <-rep(seq(1,2656,1),length(unique(annot_plasma$Kit))-1)
miRNA_rank_EdgeSeq <-seq(1,length(targeted_in_plasma_EdgeSeq$miRNA_name),1)
miRNA_rank <- c(miRNA_rank_seq,miRNA_rank_EdgeSeq)
CF$miRNA_rank <- unlist(miRNA_rank)
#+ ranking without zero count miRNAs,cache=FALSE, results='show', echo=FALSE
CF <- CF[CF$Cumulative_frequency !=0,]
miRNA_rank <- c()
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  print(i)
  miRNA_rank<-c(miRNA_rank,dim(CF[CF$variable == i, ])[1])
}
miRNA_rank <- lapply(miRNA_rank,function(x) seq(1,x, by = 1))
CF$miRNA_rank <- unlist(miRNA_rank)
#+ percentile miRNAs,cache=FALSE, results='show', echo=FALSE
percent_miRNA <- c()
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  number_of_miRNA <- dim(CF[CF$variable == i,])[1]
  percent <- as.vector(sapply(unlist(CF[CF$variable == i,"miRNA_rank"]), function(s){(s/number_of_miRNA)*100}))
  percent_miRNA <- c(percent_miRNA, percent)
}
CF$percent_miRNA <- percent_miRNA
CF$variable <- factor(CF$variable,level = c("Norgen","Lexogen","QIAseq","QIAseq_UMI","NEXTflex","RealSeq","SMARTer","EdgeSeq"))
theme_set(theme_bw())
plot_CF_plasma <-ggplot(CF, aes(x=percent_miRNA, y=Cumulative_frequency, color=variable))+
  geom_line(size=1.5)+
  scale_y_log10()+
  xlab("miRNA percentile")+
  ylab("Cumulative frequency")+
  scale_color_manual(values = c("#cccccc","#737c82","#b2d4e6","#379ed2","#063a74","#ffc5ad","#e56638", "#8f0021"))+
  ggtitle("Log10 cumulative frequency of Plasma")+
  geom_hline(yintercept = 0.01,linetype = "dashed" )+
  theme(panel.grid.minor = element_blank())
plot_CF_plasma
ggsave(file = "outs/cumulative.frequency_plasma.pdf", device = "pdf", height = 9, width = 12.5, units = "cm", useDingbats = F)
rm(CF,miRNA_rank_seq,miRNA_rank_EdgeSeq,miRNA_rank)

#'**Number of detected miRNA, dependence on raw counts**
#+ sensitivity of detection (to raw counts),results='show', echo=FALSE
if (!(dir.exists("outs/downsampling_detection_rate"))){
  dir.create("outs/downsampling_detection_rate")
}
annot_plasma <- annot_plasma[annot_plasma$Kit != "QIAseq_UMI",]
annot_plasma <- droplevels(annot_plasma)
plasma <- plasma[,rownames(annot_plasma)]
plasma_norm <- plasma_norm[,rownames(annot_plasma)]
for(number in c(1,5,10)){
  raw_counts <- t(as.data.frame(read.table("src/raw_count_star", header = TRUE, sep = "\t")))
  raw_counts <- as.data.frame(raw_counts[grep("plasma", rownames(raw_counts), ignore.case = TRUE),])
  colnames(raw_counts)<-"raw_counts"
  no_miRNA <- colSums((plasma > number)== TRUE)
  no_dropouts <- colSums((plasma > number)== FALSE)
  sensitivity <- as.data.frame(no_miRNA)
  mean_sensitivity <- c()
  for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
    if (i =="EdgeSeq"){
      df <- as.data.frame(plasma[,annot_plasma$Kit == i])
      df <- df[rownames(df) %in% targeted_in_plasma_EdgeSeq$miRNA_name,]
    }
    else {
      df <-plasma[,annot_plasma$Kit == i]
    }
    df1 <- df[(df[,1])>=number,]
    df2 <- df[(df[,2])>=number,]
    number_of_detected <- length(intersect(rownames(df1),rownames(df2)))
    mean_sensitivity <-c(mean_sensitivity,as.numeric(number_of_detected))
  }
  names(mean_sensitivity) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
  mean_sensitivity <- as.data.frame(mean_sensitivity)
  sensitivity 
  mean_sensitivity
  
  #+ downsampling,  results='hide',include=FALSE
  proportions <- 10^seq(-3,0,0.01)
  for (proportion in proportions){
    downsampled_matrix <- apply(plasma, 2, function(x) rbinom(nrow(plasma), x, proportion))
    rownames(downsampled_matrix) <- rownames(plasma)
    mean_nomiRNA <- c()
    for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
      if (i =="EdgeSeq"){
        df <- as.data.frame(downsampled_matrix[,annot_plasma$Kit == i])
        df <- df[rownames(df) %in% targeted_in_plasma_EdgeSeq$miRNA_name,]
      }
      else {
        df <-downsampled_matrix[,annot_plasma$Kit == i]
      }
      df1 <- df[(df[,1])>=number,]
      df2 <- df[(df[,2])>=number,]
      number_of_detected <- length(intersect(rownames(df1),rownames(df2)))
      mean_nomiRNA <-c(mean_nomiRNA,as.numeric(number_of_detected))
    }
    mean_nomiRNA <- as.data.frame(mean_nomiRNA)
    colnames(mean_nomiRNA) <- proportion
    mean_sensitivity<-cbind(mean_sensitivity, mean_nomiRNA)
  }
  mean_sensitivity <- mean_sensitivity[,-1]
  raw_counts_downsampled <- c("raw_counts","Protocol")
  for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
    print(raw_counts[rownames(annot_plasma[annot_plasma$Kit == i,]),])
    mean_of_samples <-mean(as.vector(raw_counts[rownames(annot_plasma[annot_plasma$Kit == i,]),]))
    raw_counts_down <- as.numeric(as.numeric(mean_of_samples)* as.numeric(proportions))
    raw_counts_down <- cbind(round(raw_counts_down),rep(i,length(raw_counts_down)))
    raw_counts_downsampled <- rbind(raw_counts_downsampled,raw_counts_down)
  }
  mean_sensitivity$Sample <- factor(rownames(mean_sensitivity), level =c("Lexogen","Norgen","QIAseq",
                                                                          "NEXTflex","RealSeq","SMARTer",
                                                                         "EdgeSeq"))
  mean_sensitivity <- melt(mean_sensitivity, value.name = "no_miRNA")
  mean_sensitivity <- mean_sensitivity[order(mean_sensitivity$Sample),]
  colnames(raw_counts_downsampled) <- raw_counts_downsampled[1,]
  raw_counts_downsampled <- as.data.frame(raw_counts_downsampled[-1,])
  raw_counts_downsampled$Protocol <- factor(raw_counts_downsampled$Protocol, level = c("Lexogen","Norgen","QIAseq",
                                                                                       "NEXTflex","RealSeq","SMARTer",
                                                                                       "EdgeSeq"))
  raw_counts_downsampled <- raw_counts_downsampled[order(raw_counts_downsampled$Protocol),]
  mean_sensitivity$raw_counts_downsampled <- as.numeric(as.character(raw_counts_downsampled$raw_counts))
  mean_sensitivity$Protocol <- raw_counts_downsampled$Protocol
  mean_sensitivity <- mean_sensitivity[!(mean_sensitivity$Protocol %in% c("QIAseq_UMI","NEXTflex_UMI")),]
  plot_S_plasma <- ggplot(data=mean_sensitivity, aes(x=raw_counts_downsampled, y=no_miRNA, colour=Sample)) + 
    geom_line(size=2)+
    xlab("Number of raw reads")+
    ggtitle("Plasma")+
    ylab("Number of detected miRNAs")+
    #geom_hline(yintercept = max(mean_nomiRNA))+
    #scale_color_brewer(palette="Set3")
    #theme(legend.position = "none")
    scale_color_manual(values = c("#cccccc","#737c82","#b2d4e6","#063a74","#ffc5ad","#e56638","#8f0021"))
  plot_S_plasma
  path <- "outs/downsampling_detection_rate/"
  ggsave(plot = plot_S_plasma, path = path, filename = paste0("detection_rate_count.thresh_greater.than_",number,".pdf"), device = "pdf", width = 15, height = 12, unit = "cm")
  #rm(raw_counts,no_miRNA,no_dropouts,mean_sensitivity,raw_counts_downsampled,raw_counts_down, mean_nomiRNA,
    # mean_of_samples, sensitivity)
}

#'**Heatmap clustering of kits and miRNA**
#+ heatmap,cache=FALSE, results='show', echo=FALSE 
load("temp/annot_plasma.RData")
load("temp/plasma.RData")
load("temp/plasma_norm.RData")
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
plasma_1 <- plasma_norm[rowSums(plasma_norm) > 32,]
scaled <- as.data.frame(t(apply(plasma_1, 1, cal_z_score)))
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
#+ heatmap_deseq1,cache=FALSE, results='show', echo=FALSE

#deseq
plasma_cluster <- plasma[rownames(plasma) %in% targeted_in_plasma_EdgeSeq$miRNA,]
de<-DESeqDataSetFromMatrix(plasma_cluster, colData = annot_plasma, design = ~ Kit)
res<-results(DESeq(de, test= "LRT", reduced = ~ 1))
exp_val<-as.data.frame(assay(de))
exp_val<-as.data.frame((sweep(exp_val, 2, colSums(exp_val), "/"))*1000000)
scaled <- as.data.frame(t(apply(exp_val, 1, cal_z_score)))
scaled$padj <-res$padj
scaled<-na.omit(scaled)
scaled<-scaled[order(scaled$padj, decreasing = FALSE),]
scaled <- scaled[,colnames(scaled) != "padj"]
top <- head(scaled,15)
heatmap_plasma <-pheatmap(top)

#'15 miRNA with highest variance
#+ heatmap_deseq,cache=FALSE, results='show', echo=FALSE
plasma_1 <- as.data.frame(plasma_norm)
plasma_1$sd <- unlist(apply(plasma_1,1,sd))
ordered <- plasma_1[order(plasma_1$sd, decreasing = TRUE),]
ordered <- ordered[,-17]
pheatmap(head(ordered,15))

#'PCA on normalized expression - variance stabilizing normalization.
#+ pca,cache=FALSE, results='show', echo=FALSE
data <- plasma[rownames(plasma) %in% targeted_in_plasma_EdgeSeq$miRNA,]
de<-DESeqDataSetFromMatrix(data, colData = annot_plasma, design = ~ Kit)
vsd<-vst(de, nsub = 200)
plotPCA(vsd, intgroup = c("Kit"))

#'**Correlation of kits after and before normalization by bias ratio**
#+ heatmap,cache=FALSE, results='show', echo=FALSE  
annot_plasma <- annot_plasma[annot_plasma$Kit != "NEXTflex_UMI",]
annot_plasma <- droplevels(annot_plasma)
plasma <- plasma[,rownames(annot_plasma)]
plasma_norm <- plasma_norm[,rownames(annot_plasma)]
normalized_data <- c()
data <-c()
log2fc <- as.data.frame(read.table("outs/log2cf_GC_sequence_by_kit.csv", header = TRUE, sep = ","))
miRNA_names <- as.data.frame(read.table("src/mirXplore_mature.tab"))
colnames(miRNA_names) <- c("mirXplore", "mature", "sequence")
miRNA_names <- miRNA_names[-1,]
rownames(miRNA_names)<- miRNA_names$mature
miRNA_names <- miRNA_names[rownames(miRNA_names) %in% targeted_in_plasma_EdgeSeq$miRNA_name,]
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  means <- as.data.frame(rowMeans(plasma_norm[rownames(plasma_norm) %in% miRNA_names$mature,annot_plasma$Kit == i]))
  names(means) <- "mean"
  means$miRNA <- rownames(means)
  bias_ratios <- c()
  mirxplore_names <-c()
  for (miRNA in rownames(means)){
    mirxplore_name <- as.character(miRNA_names[miRNA,"mirXplore"])
    mirxplore_names <-c(mirxplore_names,mirxplore_name)
    bias_ratio <- 2^(log2fc[(log2fc$Description == mirxplore_name) & (log2fc$Protocol == i),"log2fc"])
    bias_ratios <-c(bias_ratios,bias_ratio)
  }
  normalized <- means$mean/bias_ratios
  normalized_data <- cbind(normalized_data,normalized)
  data <- cbind(data,means$mean)
}
normalized_data <- as.data.frame(normalized_data)
colnames(normalized_data) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
rownames(normalized_data) <- rownames(means)
normalized_data <- normalized_data[is.finite(rowSums(normalized_data)),]
normalized_data <- normalized_data[c("Norgen","Lexogen","QIAseq","QIAseq_UMI","NEXTflex","RealSeq","SMARTer","EdgeSeq")]
normalized_data <- log2(normalized_data+1)
corrm <- Hmisc::rcorr(as.matrix(normalized_data), type = "pearson")
col1 <- colorRampPalette(brewer.pal(11, "RdBu"))
pdf(height = 7, width = 7, file = "outs/plasma_between.kit.correlation_reduced.set.after.correction_log2counts_pearson.pdf")
corrplot(corrm$r, type ="upper", method = "color", order = "original", col = rev(col1(40)),
         tl.col = "black", tl.srt = 45, outline = T, diag = F, addCoef.col = T)
dev.off()

data <- as.data.frame(data)
colnames(data) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
rownames(data) <- rownames(means)
data <- data[is.finite(rowSums(data)),]
data = data[c("Norgen","Lexogen","QIAseq","QIAseq_UMI","NEXTflex","RealSeq","SMARTer","EdgeSeq")]
data = log2(data+1)
corrm <- Hmisc::rcorr(as.matrix(data), type = "pearson")
col1 <- colorRampPalette(brewer.pal(11, "RdBu"))
pdf(height = 7, width = 7, file = "outs/plasma_between.kit.correlation_reduced.set.before.correction_log2counts_pearson.pdf")
corrplot(corrm$r, type ="upper", method = "color", order = "original", col = rev(col1(40)),
         tl.col = "black", tl.srt = 45, outline = T, diag = F, addCoef.col = T)
dev.off()

write.table(data, file = "outs/before.correction.txt", row.names = T, col.names = NA, quote = F, sep = '\t')
write.table(normalized_data, file = "outs/after.correction.txt", row.names = T, col.names = NA, quote = F, sep = '\t')

sessionInfo()
