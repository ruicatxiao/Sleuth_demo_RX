# Set working directory one level above where all the kallisto h5 are located
# Folder outline: level 1 (WD)---level 2(Contains all individual runs)---level 3(Kallisto)---.h5 files
# need to create a run info table
##      sample condition
## 1 SRR493366  scramble
## 2 SRR493367  scramble
## 3 SRR493368  scramble
## 4 SRR493369   HOXA1KD
## 5 SRR493370   HOXA1KD
## 6 SRR493371   HOXA1KD


# sample prep
sample_id <- dir(file.path(".", "results"))
kal_dirs <- file.path(".", "results", sample_id, "kallisto")
s2c <- read.table(file.path(".", "results", "metadata", "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)


# running sleuth
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')

models(so)





# Check model

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

#Top 20 differentially expressed genes
head(sleuth_significant, 20)

# Take a single gene and show its expression level comparison
plot_bootstrap(so, "ENST00000263734", units = "est_counts", color_by = "condition")



# Total visulization
sleuth_live(so)