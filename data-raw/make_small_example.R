# Create a smaller example dataset for faster vignette builds
#
# Original: 105K rows (26,456 sites × 4 contrasts)
# Target:   ~6,000 sites × 2 contrasts
#
# Biased sampling: keeps ALL sites that overlap with PTMsigDB (needed for
# PTM-SEA vignette to produce results with min_size=10), plus random
# additional sites up to the target.
#
# Usage: Rscript data-raw/make_small_example.R

set.seed(42)

load("data/combined_test_diff_example.rda")

cat("Original dimensions:", nrow(combined_test_diff_example), "x",
    ncol(combined_test_diff_example), "\n")
cat("Original contrasts:", paste(unique(combined_test_diff_example$contrast),
    collapse = ", "), "\n")

# Keep only 2 contrasts
keep_contrasts <- c("KO_vs_WT", "KO_vs_WT_at_Early")
df <- combined_test_diff_example[combined_test_diff_example$contrast %in% keep_contrasts, ]

all_sites <- unique(df$site)
cat("Sites in 2 contrasts:", length(all_sites), "\n")

# Find sites that overlap with PTMsigDB (must keep all of these)
bundled_zip <- system.file("extdata", "ptmsigdb_kinase.rds.zip", package = "prophosqua")
temp_dir <- tempdir()
unzip(bundled_zip, exdir = temp_dir)
ptmsigdb_file <- file.path(temp_dir, "ptmsigdb_filtered_KINASE_15mer.rds")
pathways <- readRDS(ptmsigdb_file)
ptmsigdb_ids <- unique(gsub(";[ud]$", "", unlist(pathways)))

all_site_ids <- paste0(
  vapply(trimws(toupper(df$SequenceWindow)),
         function(x) prophosqua:::trim_flanking_seq(x, trim_to = 15),
         character(1)),
  "-p"
)
ptmsigdb_sites <- unique(df$site[all_site_ids %in% ptmsigdb_ids])
cat("PTMsigDB-matching sites (must keep):", length(ptmsigdb_sites), "\n")

# Sample additional random sites to reach target
n_target <- 6000
remaining_sites <- setdiff(all_sites, ptmsigdb_sites)
n_extra <- max(0, n_target - length(ptmsigdb_sites))
extra_sites <- sample(remaining_sites, min(n_extra, length(remaining_sites)))

sampled_sites <- union(ptmsigdb_sites, extra_sites)
cat("Total sampled sites:", length(sampled_sites), "\n")

combined_test_diff_example <- df[df$site %in% sampled_sites, ]

cat("Final dimensions:", nrow(combined_test_diff_example), "x",
    ncol(combined_test_diff_example), "\n")
cat("Final contrasts:", paste(unique(combined_test_diff_example$contrast),
    collapse = ", "), "\n")
cat("Sites per contrast:\n")
print(table(combined_test_diff_example$contrast))

save(combined_test_diff_example, file = "data/combined_test_diff_example.rda",
     compress = "xz")

cat("Saved to data/combined_test_diff_example.rda\n")
cat("File size:", round(file.size("data/combined_test_diff_example.rda") / 1024, 1), "KB\n")
