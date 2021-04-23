#' ---
#' title: "RT-qPCR"
#' author: "Benesova + Androvic"
#' ---

#+ setup1, include=FALSE
library("readxl")
library(ggplot2)
library(gridExtra)
library(CORREP)
library(effects)
library(readxl)

#+ load data, include=FALSE
load("temp/annot_plasma.RData")
load("temp/plasma_norm.RData")
log2fc <- read.table("outs/log2cf_GC_sequence_by_kit.csv", header = TRUE, sep = ",")
miRXplore_plasma_names <- read.table("src/mirXplore_mature.tab", header = TRUE)
targeted_in_plasma_EdgeSeq <- as.data.frame(read.table("src/plasma_targeted_EdgeSeq.tab", header = TRUE, sep = "\t"))

#+ load qPCR data, include=FALSE
qPCR <- as.data.frame(read_excel("src/qPCR_quantification.xlsx",
                        sheet = "df_for_R2"))
rownames(qPCR) <- qPCR[,1]
qPCR <- qPCR[,-1]
qPCR <-qPCR[,grep("<", qPCR[5,], invert = TRUE)]
qPCR <-qPCR[,grep("let", colnames(qPCR), invert = TRUE)]
#remove miRNAs with unoptimal efficiency
off_eff <- c("hsa-miR-330-3p", "hsa-miR-125b-5p", "hsa-miR-320a-3p", "hsa_miR-222-3p", "hsa_miR-486-5p",
             "hsa_miR-20a-5p", "hsa_miR-29a-3p")
qPCR <- qPCR[,!(colnames(qPCR) %in% off_eff)]
qPCR_for_plot <- qPCR[1:4,]
qPCR_for_plot2<-melt(qPCR_for_plot)
miRNA_order <- names(sort(-colMeans(qPCR_for_plot)))
colnames(qPCR_for_plot2)<-c("miRNA", "Absolute_concentration")
qPCR_for_plot2$Absolute_concentration_log10 <-log10(qPCR_for_plot2$Absolute_concentration)
colorcount <- length(unique(qPCR_for_plot2$miRNA))
getPalette <- colorRampPalette(brewer.pal(11,"RdBu"))
ggplot(data = qPCR_for_plot2, aes(x = miRNA, y= Absolute_concentration_log10, fill = miRNA)) + 
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank())+
  scale_x_discrete(limits = miRNA_order)+
  scale_fill_manual(values = getPalette(19), limits = miRNA_order)
ggsave(file = "outs/absolute.quantities_sorted.pdf",device = "pdf", height = 20, width = 18, units = "cm")
 mir1AV <-colMeans(qPCR[c(1,2),], na.rm = TRUE) 
mir2AV <- colMeans(qPCR[c(3,4),], na.rm = TRUE) 
qPCR <- as.data.frame(t(rbind(mir1AV, mir2AV)))

#+ combine qPCR and RNAseq data - correlation, include=FALSE
RNAseq <- plasma_norm[rownames(qPCR),]
missing_miRNA_EdgeSeq <- rownames(RNAseq[!(rownames(RNAseq) %in% targeted_in_plasma_EdgeSeq$miRNA_name),
                            c("EdgeSeq_plasma_1","EdgeSeq_plasma_2")])
RNAseq[missing_miRNA_EdgeSeq,c("EdgeSeq_plasma_1", "EdgeSeq_plasma_2")] <- NA
combined <- cbind(qPCR, RNAseq)
results_df <- data.frame("Statistic" = c("corr.p_value", "corr.coef", "r2", "sw.test.qPCR", "sw.test.RNAseq",
                                         "sw.test.res"))
ind <- 0
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  ind <- ind + 1
  df <- as.data.frame(log10(rowMeans(combined[,colnames(combined) %in% rownames(annot_plasma[annot_plasma$Kit == i,])]) +1))
  df <- cbind(df,as.data.frame(log10(rowMeans(combined[,c("mir1AV", "mir2AV")]) +1)))
  colnames(df)<-c("RNAseq", "qPCR")
  correlation <-cor.test(df$RNAseq,df$qPCR, method  = "pearson")
  lm1<-lm(RNAseq ~ qPCR, data = df)
  lm_save <- summary(lm1)
  results_df[paste(i)] <-c(correlation$p.value, correlation$estimate, 
                                                                  lm_save$r.squared, shapiro.test(df$qPCR)$p.value,
                                                                  shapiro.test(df$RNAseq)$p.value, 
                                                                  shapiro.test(lm_save$residuals)$p.value)
  intercept1 <- lm1$coefficients[1]
  slope1 <- lm1$coefficients[2]
  r <- round(correlation$estimate,2)
  r2 <- round(lm_save$r.squared,2)
  p <- round(correlation$p.value,4)
  rsq_label <- paste('R^2 == ', r2)
  r_label <- paste('R == ', r)
  p_label <- paste('p-value == ', p)
  assign(x=paste0("p",ind),value=p<-ggplot(data=df, aes(x=qPCR, y=RNAseq))+
           geom_smooth(method = "lm", size = 0.5, color = "#8f0021", fill = "#EDF2F4", alpha = 1)+
           geom_point(color = "#063a74", size = 2)+
           theme(panel.grid = element_blank())+
           coord_cartesian(ylim =c(-1, 6), xlim = c(2.2, 8))+
           annotate(geom="text", x=3, y=5.9, label=rsq_label,color="gray5", parse = TRUE)+
           annotate(geom="text", x=3, y=4.9, label=r_label,color="gray5", parse = TRUE)+
           ylab("RNAseq")+
           ggtitle(i))
}
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3)
ggsave(file = "outs/qPCR_RNASeq_correlation.pdf", device = "pdf", height =15, width = 15, units = "cm")
rownames(results_df) <- results_df[,1]
results_df <- results_df[,-1]
write.csv(results_df, file = "outs/qPCR_RNAseq_corr_mean.csv")
write.csv(combined, file = "outs/correlation_input_data")

#+ correlation-after normalization by miRXplore bias, include=FALSE
annot_plasma <-annot_plasma[(annot_plasma$Kit != "NEXTflex_UMI")&
                              (annot_plasma$Kit != "QIAseq_UMI"),]
annot_plasma <-droplevels(annot_plasma)
combined2 <- combined
df <- c()
ind <- 0
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  print(i)
  log2fc_vector <- c()
  for (miRNA in rownames(combined2)){
    if (miRNA %in% miRXplore_plasma_names$name_mature){
      miRXplore_name <- as.character(miRXplore_plasma_names[miRXplore_plasma_names$name_mature == miRNA, "name_mirxplore"])
      if (miRXplore_name %in% log2fc[(log2fc$Protocol == i),"Description"]){
      log2fc_value <- log2fc[(log2fc$Protocol == i) & (log2fc$Description == miRXplore_name),"log2fc"]
      log2fc_vector <-c(log2fc_vector,log2fc_value)
      }
      else{
        print("missing")
        log2fc_vector<-c(log2fc_vector,NA)
        }
    }
    else{
      log2fc_vector<-c(log2fc_vector,NA)
    }
  }
  print(log2fc_vector)
  df <- cbind(df,log2fc_vector)
}
bias <- as.data.frame(df)
colnames(bias) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
rownames(bias) <- rownames(combined)

bias_log2 <- as.data.frame(bias)
colnames(bias_log2) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
rownames(bias_log2) <- rownames(combined)
bias_ratio <- 2^bias_log2
RNAseq2 <- c()
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){RNAseq2 <- cbind(RNAseq2,rowMeans(RNAseq[,rownames(annot_plasma[annot_plasma$Kit == i,])]))}
RNAseq <- as.data.frame(RNAseq2)
colnames(RNAseq) <- unique(annot_plasma[order(annot_plasma$Kit),"Kit"])
RNAseq_normalized <- RNAseq/bias_ratio


results_df <- data.frame("Statistic" = c("corr.p_value", "corr.coef", "r2", "sw.test.qPCR", "sw.test.RNAseq",
                                         "sw.test.res"))
combined <- cbind(qPCR, RNAseq_normalized)
ind <- 0
for (i in unique(annot_plasma[order(annot_plasma$Kit),"Kit"])){
  ind <- ind + 1
  df <- as.data.frame(log10(combined[,i] +1))
  df <- cbind(df,as.data.frame(log10(rowMeans(combined[,c("mir1AV", "mir2AV")]) +1)))
  colnames(df)<-c("RNAseq", "qPCR")
  correlation <-cor.test(df$RNAseq,df$qPCR, method  = "pearson")
  lm1<-lm(RNAseq ~ qPCR, data = df)
  lm_save <- summary(lm1)
  results_df[paste(i)] <-c(correlation$p.value, correlation$estimate, 
                           lm_save$r.squared, shapiro.test(df$qPCR)$p.value,
                           shapiro.test(df$RNAseq)$p.value, 
                           shapiro.test(lm_save$residuals)$p.value)
  intercept1 <- lm1$coefficients[1]
  slope1 <- lm1$coefficients[2]
  r <- round(correlation$estimate,2)
  r2 <- round(lm_save$r.squared,2)
  p <- round(correlation$p.value,4)
  rsq_label <- paste('R^2 == ', r2)
  r_label <- paste('R == ', r)
  p_label <- paste('p-value == ', p)
  assign(x=paste0("p",ind),value=p<-ggplot(data=df, aes(x=qPCR, y=RNAseq))+
           geom_smooth(method = "lm", size = 0.5, color = "#8f0021", fill = "#EDF2F4", alpha = 1)+
           geom_point(color = "#063a74", size = 2)+
           theme(panel.grid = element_blank())+
           coord_cartesian(ylim =c(-1, 6), xlim = c(2.2, 8))+
           annotate(geom="text", x=3, y=5.9, label=rsq_label,color="gray5", parse = TRUE)+
           annotate(geom="text", x=3, y=4.9, label=r_label,color="gray5", parse = TRUE)+
           ylab("RNAseq")+
           ggtitle(i))
}
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3)

