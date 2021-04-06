#' ---
#' title: "Count tables"
#' author: "Benesova"
#' ---

#+ setup1, include=FALSE
library(GenomicFeatures)
library(Rsubread)

#+ load annotation,cache=FALSE, results='show', echo=FALSE
annot <-read.table("src/annotation.tab",sep = "\t", header = TRUE) 
annot_miRXplore <- annot[annot$Sample_source=="miRXplore",]
rownames(annot_miRXplore) <- annot_miRXplore$Name
save(annot_miRXplore,file="temp/annot_miRXplore.RData")
annot_plasma <- annot[annot$Sample_source=="plasma_gel",]
rownames(annot_plasma) <-  annot_plasma$Name
save(annot_plasma,file="temp/annot_plasma.RData")
annot_ncRNA <- annot_plasma[(!(annot_plasma$Kit %in% c("QIAseq_UMI", "NEXTflex_UMI"))),]
rownames(annot_ncRNA) <- annot_ncRNA$Name  

#'**miRXplore**

#+ get bam files,cache=FALSE, results='show', echo=FALSE
path<-"src/mapped_miRXplore/"
files<-list.files(path)
idx<-grep("miRXplore",files)
files_miRXplore<-files[idx]
miRXplore<-paste0(path, files_miRXplore)
#+ make count table,cache=FALSE, results='show', echo=FALSE
saf <-read.table("src/miRXplore.saf",
                              sep = "\t", header = TRUE) 

f<-featureCounts(miRXplore,annot.ext = saf, nthreads = 8, verbose = TRUE,
                        countMultiMappingReads=FALSE)
miRXplore<-f$counts
colnames(miRXplore)<-annot_miRXplore$Name
save(miRXplore, file = "temp/miRXplore.RData")
#+ normalize count table,cache=FALSE, results='show', echo=FALSE
miRXplore_norm <- (sweep(miRXplore, 2, colSums(miRXplore), "/"))*1000000
save(miRXplore_norm, file = "temp/miRXplore_norm.RData")


#'**mature miRNAs**

#+ get bam files mirna,cache=FALSE, results='show', echo=FALSE
path<-"src/mapped_mature/"
files<-list.files(path)
idx<-grep("plasma",files)
files_idx<-files[idx]
files<-paste0(path, files_idx)
#+ make count table mirna,cache=FALSE, results='show', echo=FALSE
saf <-read.table("src/mature.saf",
                 sep = "\t", header = TRUE) 

f<-featureCounts(files, annot.ext = saf, 
                 nthreads = 8,verbose = TRUE, countMultiMappingReads=FALSE)
plasma<-f$counts
colnames(plasma)<-annot_plasma$Name
save(plasma, file = "temp/plasma.RData")
#+ normalize count table mirna,cache=FALSE, results='show', echo=FALSE
plasma_norm <- (sweep(plasma, 2, colSums(plasma), "/"))*1000000
save(plasma_norm, file = "temp/plasma_norm.RData")


#'**ncRNA**
#+ get bam files ncRNA,cache=FALSE, results='show', echo=FALSE
path<-"src/mapped_ncRNA/"
files<-list.files(path)
idx<-grep("plasma",files)
files_idx<-files[idx]
files<-paste0(path, files_idx)
#+ make count table ncRNA,cache=FALSE, results='show', echo=FALSE
saf <-read.table("src/ncRNA.saf",
                 sep = "\t", header = TRUE) 

f<-featureCounts(files, annot.ext = saf, 
                 nthreads = 8,verbose = TRUE, countMultiMappingReads=FALSE)
ncRNA<-f$counts
colnames(ncRNA)<-annot_ncRNA$Name
save(ncRNA, file = "temp/ncRNA_counts.RData")


sessionInfo()


