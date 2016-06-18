

lambda_out <- read.table("../sfa_outputs/Deng2014PPS/sqrt_deng_counts_lambda.out")
f_out <- t(read.table("../sfa_outputs/Deng2014PPS/sqrt_deng_counts_F.out"))

indices_mat <- SFA.ExtractTopFeatures(f_out, top_features = 100, options="min")

library(devtools)

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()

deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

top_gene_names <- t(apply(indices_mat, 1, function(x) return (deng.gene_names[x])));
