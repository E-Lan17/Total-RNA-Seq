##START :
library(sleuth)
library(data.table)
library(grid)
library(gridExtra)
library(dplyr)
 
# set the number of available cores to 4
options(mc.cores = 4L)

# get t2g
tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}
t2g <- tx2gene()
t2g

# get samples
sample_id <- dir(file.path("/home/XXX/Bureau/", "Quant"))
kal_dirs <- file.path("/home/XXX/Bureau", "Quant", sample_id)
s2c <- read.table(file.path("/home/XXX/Bureau", "/Meta", "samples.txt"), header = TRUE, stringsAsFactors=FALSE)

s2c <- dplyr::select(s2c, sample = sample, condition)
#s2c <- dplyr::select(s2c, sample = sample, condition, genotype)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c

# Clean h5 (remove .N version transcript)
library(rhdf5)
files <- list.files("/home/XXX/Bureau/Quant", pattern=".h5", recursive=TRUE, full.names=TRUE)
for (currentFile in files) {
  oldids <- h5read(currentFile, "/aux/ids")
  newids <- gsub("\\..*","", oldids)
  h5write(newids, currentFile, "/aux/ids")
}

# create sleuth object
#so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, extra_bootstrap_summary = TRUE, 
#read_bootstrap_tpm = TRUE, normalize = TRUE, 
#transformation_function = function(x) log2(x + 0.5), num_cores = max(4))
  #aggregation_column = 'ens_gene'))
#summarydata$sampletype <- relevel(summarydata$sampletype, ref = "MOV10_knockdown")
# pour changer le referentiel

so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, extra_bootstrap_summary = TRUE, 
read_bootstrap_tpm = TRUE, normalize = TRUE, aggregation_column = 'ens_gene')

# Fit models
so = sleuth_fit(so, ~condition + reps, 'full')
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~reps, "reduced")
so = sleuth_fit(so, ~1, "reduced")
so = sleuth_lrt(so,"reduced", "full")
#lrt good for sgnifiance WT good for the value of the fold change
#so = sleuth_wt(so,"conditionP1",'full')
#so = sleuth_wt(so,"conditionW8",'full')

# Store results
res_lrt= sleuth_results(so, test = "reduced:full", test_type = "lrt", show_all = TRUE)
res_lrt= sleuth_results(so, test = "reduced:full", test_type = "lrt", show_all = FALSE)
results_lrt_significant <- res_lrt$target_id[which(res_lrt$qval <= 0.05)]
head(res_lrt)

#res_wald <- sleuth_results(so, 'conditionP1','full', test_type = 'wt')
#results_wald_significant <- res_wald$target_id[which(res_wald$qval <= 0.05)]
#head(res_wald)
table = kallisto_table(so, use_filtered = FALSE, normalized = TRUE,
                       include_covariates = TRUE)
head(table)
table_1= subset(table, target_id == "ENSMUST00000109085")
write.csv(table_1, file = "/home/XXX/Bureau/ENSMUST00000109085.csv", row.names=FALSE)
models(so)

# Show plot
plot_bootstrap(so, "ENSMUST00000109087", units = "tpm", color_by = "condition")
plot_bootstrap(so, "ENSMUST00000109085", units = "est_counts", color_by = "condition")

sleuth_live(so)

#Lower divergence values represent samples that are more similar to each other
png("/home/XXX/Bureau/DTE_sleuth_sample_heatmap.png", width=7, height=7, units = "in", res = 300)
plot_sample_heatmap(so, use_filtered = FALSE, units = "tpm",color_high = "#214eff",
color_mid = "#FFC300", color_low = "#ff0000")
dev.off()

# plot_pca(so, color_by ="condition")
png("/home/XXX/Bureau/DTE_sleuth_PCA.png", width=7, height=7, units = "in", res = 300)
plot_pca(so, color_by = "condition", text_labels = TRUE, units = "tpm")
dev.off()

png("/home/XXX/Bureau/DTE_sleuth_bout.png", width=7, height=7, units = "in", res = 300)
plot_bootstrap(so, "ENSMUST00000117757", units = "tpm",color_by = "condition")
dev.off()



png("/home/XXX/Bureau/DGE_sleuth_transcript_heatmap.png", width=7, height=7, units = "in", res = 300)
plot_transcript_heatmap(so ,
                        transcripts = subset(res_lrt[order(res_lrt$qval),],
                        qval < 0.05)$target_id[1:40], 
                        units = "tpm", cluster_transcripts = TRUE,  
                        color_high = "#214eff", color_mid = "#FFC300", color_low = "#ff0000")
dev.off()

sig_transcripts <- res_lrt %>% filter(target_id == "ENSMUST00000109597" | target_id == "ENSMUST00000109087" )
transcripts = subset(res_lrt[order(res_lrt$target_id),],target_id == "ENSMUST00000109597 ENSMUST00000109087")
transcripts = subset(res_lrt[order(res_lrt$ext_gene),],ext_gene  == "Gnas" )
transcripts = subset(res_lrt[order(res_lrt$qval),],qval < 0.05)$target_id[1:40]

png("/home/XXX/Bureau/DTE_sleuth_transcript_heatmap.png", width=7, height=7, units = "in", res = 300)
plot_transcript_heatmap(so,
                        transcripts, 
                        units = "tpm", cluster_transcripts = TRUE, use_filtered = FALSE,
                        color_high = "#214eff", color_mid = "#FFC300", color_low = "#ff0000")
dev.off()


plot_transcript_heatmap(
  so,
  transcripts = head(table, n = 10)$target_id,
  labels_row = table$gene_name[0:10]
)

ulimit -s
Cstack_info()
