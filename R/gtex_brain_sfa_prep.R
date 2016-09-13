

###############  Create the SFA script for the GTEx Brain data ######################


library(GTExV6Brain)
library(ggplot2)
library(CountClust)
counts <- exprs(GTExV6Brain)
meta_data <- pData(GTExV6Brain)
gene_names <- rownames(counts)

voom_counts <- t(limma::voom(counts)$E);
write.table(voom_counts, 
            file="../../sfa_seq/sfa_inputs/brain_gtex_voom_counts.txt", 
            row.names=FALSE,col.names=FALSE, quote=FALSE, sep="\t")

write.table(t(voom_counts), 
            file="../../sfa_seq/sfa_inputs/brain_gtex_voom_counts_transpose.txt", 
            row.names=FALSE,col.names=FALSE, quote=FALSE, sep="\t")

./sfa_mac -gen ../../sfa_inputs/brain_gtex_voom_counts.txt -g 1259 -n 16069 -k 6 -iter 50 -r 800 -o ../../sfa_outputs/gtex_brain_sparse_loadings


###  Plot the GTEx Brain FactorGG plots

ll <- read.table("../sfa_outputs/GTEXBrain2013/gtex_brain_sparse_loadings_lambda.out")

omega <- ll

tissue_labels <- meta_data[ ,3]

tissue_levels <- unique(tissue_labels);

tissue_labels[grep("accumbens", tissue_labels)] = "Brain - N. accumbens"
tissue_labels[grep("Caudate", tissue_labels)] = "Brain - Caudate"
tissue_labels[grep("Putamen", tissue_labels)] = "Brain - Putamen"
tissue_labels[grep("cingulate", tissue_labels)] = "Brain - Anterior cortex (BA24)"

brain_labels <- tissue_labels[grep("Brain", tissue_labels)]

# assign tissue labels
rownames(omega) <- paste0("X", 1:length(brain_labels))
annotation <- data.frame(
  sample_id = paste0("X", 1:length(brain_labels)),
  label = factor(brain_labels,
                        levels = rev(c("Brain - Spinal cord (cervical c-1)",
                                       "Brain - Substantia nigra",
                                       "Brain - Amygdala",
                                       "Brain - Hypothalamus",
                                       "Brain - Hippocampus",
                                       "Brain - Putamen",
                                       "Brain - Caudate",
                                       "Brain - N. accumbens",
                                       "Brain - Frontal Cortex (BA9)",
                                       "Brain - Anterior cortex (BA24)",
                                       "Brain - Cortex",
                                       "Brain - Cerebellum",
                                       "Brain - Cerebellar Hemisphere") ) ) )

# define colors of the clusers
cols <- c("blue", "darkgoldenrod1", "cyan", "red")

cols <- c("blue", "darkgoldenrod1", "cyan", "red", "green", "yellow")

# pdf("plots/gtex-figures/brain-barplot.pdf",
#     height = 4, width = 4)
library(flashr)
FactorGGStack(loadings = omega,
              annotation = annotation,
              palette = cols,
              yaxis_label = "Tissue type",
              order_sample = TRUE,
              figure_title = "SFA Loadings Structure Plot (Sparse Loadings case)",
              legend_labels = NULL,
              scale=TRUE,
              axis_tick = list(axis_ticks_length = .1,
                               axis_ticks_lwd_y = .1,
                               axis_ticks_lwd_x = .1,
                               axis_label_size = 7,
                               axis_label_face = "bold"))

#dev.off()


ll <- read.table("../sfa_outputs/GTEXBrain2013_transpose/gtex_brain_sparse_factors_F.out")
omega <- t(ll)
cols <- c("blue", "darkgoldenrod1", "cyan", "red", "green", "yellow")

# pdf("plots/gtex-figures/brain-barplot.pdf",
#     height = 4, width = 4)
library(flashr)
FactorGGStack(loadings = omega,
              annotation = annotation,
              palette = cols,
              yaxis_label = "Tissue type",
              order_sample = TRUE,
              figure_title = "SFA Loadings Structure Plot (Sparse Factor case)",
              legend_labels = NULL,
              scale=TRUE,
              axis_tick = list(axis_ticks_length = .1,
                               axis_ticks_lwd_y = .1,
                               axis_ticks_lwd_x = .1,
                               axis_label_size = 7,
                               axis_label_face = "bold"))
