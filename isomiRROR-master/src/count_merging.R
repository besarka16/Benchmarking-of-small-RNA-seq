library(data.table)
my_files_i <- unlist(snakemake@input)
my_data_i <- lapply(my_files_i, function(x) {dat=fread(x, header=F, quote= "",col.names=c("isomir",x)); return (dat)})
isomir_readcount <- Reduce(function(...) merge(..., all = T, by = c("isomir")), my_data_i)
isomir_readcount <- setnames(isomir_readcount,c("isomir", snakemake@params$samples))
isomir_readcount[is.na(isomir_readcount)] <- 0
write.table(isomir_readcount,snakemake@output[[1]],row.names=F,col.names=T, quote=F, sep = '\t')
