# Hand-made DEA result tables for the DPA/DPU table logic: the smallest site and
# protein results with the columns the pairing reads.

# One row per protein x contrast, with the columns test_diff() needs: the effect
# size, its standard error and its degrees of freedom.
protein_dea_table <- function(protein_ids = c("P1", "P2", "P3"), contrasts = c("a_vs_b", "c_vs_b")) {
  grid <- expand.grid(
    protein_Id = protein_ids,
    contrast = contrasts,
    stringsAsFactors = FALSE
  )
  n <- nrow(grid)
  data.frame(
    protein_Id = grid$protein_Id,
    contrast = grid$contrast,
    description = paste("protein", grid$protein_Id),
    protein_length = rep(300L, n),
    gene_name = paste0("GENE", sub("^P", "", grid$protein_Id)),
    diff = seq(-1, 1, length.out = n),
    std.error = rep(0.2, n),
    df = rep(6, n),
    std.error.unmoderated = rep(0.25, n),
    df.unmoderated = rep(4, n),
    statistic = seq(-3, 3, length.out = n),
    FDR = rep(0.01, n),
    estimate_type = rep("observed", n),
    stringsAsFactors = FALSE
  )
}

# The phospho side carries the same columns plus the site identifier and the
# site annotation. A PTM reader keys its row annotation on protein and site, so
# modAA, posInProtein and SequenceWindow reach every annotated sheet of the DEA
# output, the diff_exp_analysis one included. `sites` names one site per protein
# in `protein_ids`.
site_dea_table <- function(protein_ids = c("P1", "P2"), contrasts = c("a_vs_b", "c_vs_b")) {
  tab <- protein_dea_table(protein_ids, contrasts)
  tab$site <- paste0(tab$protein_Id, "~S10")
  tab$posInProtein <- 10L
  tab$modAA <- "S"
  tab$SequenceWindow <- "AAAAAAASAAAAAAA"
  tab$diff <- tab$diff + 0.5
  tab
}

# A pair in the shape the pair computation reads.
dea_result_pair <- function(site = site_dea_table(), protein = protein_dea_table()) {
  list(
    site = list(differential_results = site),
    protein = list(differential_results = protein)
  )
}
