# Synthetic prolfquapp DEA output directories.
#
# The compute functions read a DEA run the way prolfquapp writes it: a
# Results_WU_* subdirectory holding a DE_*.xlsx with a diff_exp_analysis sheet
# and, for a phospho run, a normalized_abundances sheet carrying the site
# annotation. These builders produce the smallest output with that shape.

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
    statistic = seq(-3, 3, length.out = n),
    FDR = rep(0.01, n),
    estimate_type = rep("observed", n),
    stringsAsFactors = FALSE
  )
}

# The phospho side carries the same columns plus the site identifier. `sites`
# names one site per protein in `protein_ids`.
site_dea_table <- function(protein_ids = c("P1", "P2"), contrasts = c("a_vs_b", "c_vs_b")) {
  tab <- protein_dea_table(protein_ids, contrasts)
  tab$site <- paste0(tab$protein_Id, "~S10")
  tab$diff <- tab$diff + 0.5
  tab
}

site_annotation_table <- function(protein_ids = c("P1", "P2")) {
  data.frame(
    site = paste0(protein_ids, "~S10"),
    posInProtein = rep(10L, length(protein_ids)),
    modAA = rep("S", length(protein_ids)),
    SequenceWindow = rep("AAAAAAASAAAAAAA", length(protein_ids)),
    protein_Id = protein_ids,
    gene_name = paste0("GENE", sub("^P", "", protein_ids)),
    protein_length = rep(300L, length(protein_ids)),
    stringsAsFactors = FALSE
  )
}

make_dea_output <- function(diff_exp, normalized = NULL, name = "DEA_fixture") {
  dea_dir <- tempfile(pattern = name)
  results <- file.path(dea_dir, "Results_WU_fixture")
  dir.create(results, recursive = TRUE)

  sheets <- list(diff_exp_analysis = diff_exp)
  if (!is.null(normalized)) {
    sheets$normalized_abundances <- normalized
  }
  writexl::write_xlsx(sheets, file.path(results, "DE_fixture.xlsx"))

  dea_dir
}

# A matched pair: two of the three proteins carry a tested site, so exactly one
# protein has no site and one site has a protein.
make_dea_pair <- function() {
  list(
    phospho = make_dea_output(
      site_dea_table(),
      site_annotation_table(),
      name = "DEA_phospho"
    ),
    protein = make_dea_output(protein_dea_table(), name = "DEA_protein")
  )
}
