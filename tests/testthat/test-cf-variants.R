cf_pair <- function() {
  paths <- anndata_pair_fixture()
  list(pair = .read_dea_pair(paths$dirs$phospho, paths$dirs$protein), paths = paths)
}

cf_from_pair <- function(pair) suppressWarnings(.compute_cf_dea_from_pair(pair))

test_that("CF returns both imputed variants on the site axes", {
  fixture <- cf_pair()
  cf <- cf_from_pair(fixture$pair)

  expect_named(cf$variants, c("correct_first_protein_imputed", "correct_first_site_protein_imputed"))
  columns <- c("site", "contrast", "diff.site", "FDR.site", "estimate_type", "n_protein_imputed")
  expect_true(all(columns %in% names(cf$variants$correct_first_protein_imputed$results)))
  expect_named(cf$variants$correct_first_site_protein_imputed, "abundances")
  for (variant in cf$variants) {
    expect_equal(dim(variant$abundances), c(nrow(fixture$pair$site$obs), nrow(fixture$pair$site$var)))
  }
  # The fixture's proteome is complete, so imputing it changes nothing.
  expect_equal(
    cf$variants$correct_first_protein_imputed$results$diff.site,
    cf$results$diff.site
  )
})

test_that("a sample whose protein is missing stays in the protein-imputed variant only", {
  fixture <- cf_pair()
  pair <- fixture$pair
  protein <- pair$protein
  key <- protein$sample_key
  target_protein <- protein$normalized_abundances$protein_Id[[1]]
  target_sample <- protein$normalized_abundances[[key]][[1]]
  knocked_out <- protein$normalized_abundances$protein_Id == target_protein &
    protein$normalized_abundances[[key]] == target_sample
  pair$protein$normalized_abundances$normalized_abundance[knocked_out] <- NA

  cf <- cf_from_pair(pair)
  site_rows <- pair$site$normalized_abundances$protein_Id == target_protein &
    pair$site$normalized_abundances[[pair$site$sample_key]] == target_sample &
    !is.na(pair$site$normalized_abundances$normalized_abundance)
  site_row <- pair$site$normalized_abundances[which(site_rows)[1], ]

  current <- cf$wide_data[cf$wide_data$site == site_row$site, target_sample, drop = TRUE]
  expect_true(length(current) == 0L || is.na(current))

  imputed_protein <- pair$protein$imputed_abundances$normalized_abundance[
    pair$protein$imputed_abundances$protein_Id == target_protein &
      pair$protein$imputed_abundances[[key]] == target_sample
  ]
  observed_median <- stats::median(
    pair$protein$normalized_abundances$normalized_abundance[
      pair$protein$normalized_abundances[[key]] == target_sample
    ],
    na.rm = TRUE
  )
  expect_equal(
    cf$variants$correct_first_protein_imputed$abundances[target_sample, site_row$site],
    site_row$normalized_abundance - imputed_protein + observed_median
  )

  protein_imputed <- cf$variants$correct_first_protein_imputed$results
  protein_imputed <- protein_imputed[protein_imputed$site == site_row$site, ]
  expect_true(all(protein_imputed$n_protein_imputed == 1L))
  expect_true(all(cf$results$n_protein_imputed == 0L))
})

test_that("the CorrectFirst axis leaves out sites whose protein the proteome did not quantify", {
  fixture <- cf_pair()
  pair <- fixture$pair
  dropped <- pair$protein$normalized_abundances$protein_Id[[1]]
  for (table in c("normalized_abundances", "imputed_abundances")) {
    long <- pair$protein[[table]]
    pair$protein[[table]] <- long[long$protein_Id != dropped, ]
  }
  site_long <- pair$site$normalized_abundances
  kept <- unique(site_long$site[site_long$protein_Id != dropped])
  expect_lt(length(kept), nrow(pair$site$var))

  cf <- cf_from_pair(pair)
  for (variant in cf$variants) {
    expect_setequal(colnames(variant$abundances), kept)
  }
})

test_that("the site-imputed variant is the lm-filled site minus the imputed protein, never fitted", {
  fixture <- cf_pair()
  pair <- fixture$pair
  cf <- cf_from_pair(pair)
  expect_named(cf$variants$correct_first_site_protein_imputed, "abundances")

  site <- pair$site
  key <- site$sample_key
  fitted_sites <- site$imputation$site[site$imputation$route == "fitted"]
  gaps <- site$normalized_abundances[
    is.na(site$normalized_abundances$normalized_abundance) &
      site$normalized_abundances$site %in% fitted_sites,
  ]
  expect_gt(nrow(gaps), 0)
  gap <- gaps[1, ]
  sample <- as.character(gap[[key]])
  at <- function(long, sample_key, rows) long$normalized_abundance[rows & long[[sample_key]] == sample]
  filled_site <- at(site$imputed_abundances, key, site$imputed_abundances$site == gap$site)
  protein <- pair$protein
  imputed_protein <- at(
    protein$imputed_abundances,
    protein$sample_key,
    protein$imputed_abundances$protein_Id == gap$protein_Id
  )
  observed_median <- stats::median(
    protein$normalized_abundances$normalized_abundance[protein$normalized_abundances[[protein$sample_key]] == sample],
    na.rm = TRUE
  )
  expect_equal(
    cf$variants$correct_first_site_protein_imputed$abundances[sample, gap$site],
    filled_site - imputed_protein + observed_median
  )
})

test_that("the site-imputed variant fills sites from their own fit only", {
  fixture <- cf_pair()
  pair <- fixture$pair
  gaps <- stats::aggregate(
    is.na(normalized_abundance) ~ site,
    data = pair$site$normalized_abundances,
    FUN = sum
  )
  fitted_sites <- pair$site$imputation$site[pair$site$imputation$route == "fitted"]
  gapped <- intersect(gaps$site[gaps[[2]] > 0], fitted_sites)
  skip_if(length(gapped) < 2L, "fixture has fewer than two fitted sites with gaps")
  pair$site$imputation$route[pair$site$imputation$site == gapped[[1]]] <- "lod_refit"

  filled <- .lm_filled_sites(pair$site)
  cells <- function(long, id) long$normalized_abundance[long$site == id]
  expect_equal(cells(filled, gapped[[1]]), cells(pair$site$normalized_abundances, gapped[[1]]))
  expect_equal(cells(filled, gapped[[2]]), cells(pair$site$imputed_abundances, gapped[[2]]))
  expect_false(anyNA(cells(filled, gapped[[2]])))
})

test_that("CF fits the model formula the DEAs recorded", {
  fixture <- cf_pair()
  pair <- fixture$pair
  samples <- pair$site$obs[[pair$site$sample_key]]
  batch <- stats::setNames(rep(c("x", "y"), length.out = length(samples)), samples)
  add_batch <- function(long, key) {
    long$batch <- unname(batch[as.character(long[[key]])])
    long
  }
  for (side in c("site", "protein")) {
    key <- pair[[side]]$sample_key
    pair[[side]]$normalized_abundances <- add_batch(pair[[side]]$normalized_abundances, key)
    pair[[side]]$imputed_abundances <- add_batch(pair[[side]]$imputed_abundances, key)
    pair[[side]]$formula <- "normalized_abundance ~ G_ + batch"
  }

  cf <- cf_from_pair(pair)
  fitted <- cf$variants$correct_first_protein_imputed$results
  fitted <- fitted[fitted$estimate_type == "observed", ]
  site <- fitted$site[[1]]
  contrast <- fitted$contrast[[1]]
  abundances <- cf$variants$correct_first_protein_imputed$abundances
  design <- pair$site$obs
  data <- data.frame(
    y = abundances[design[[pair$site$sample_key]], site],
    G_ = design$G_,
    batch = unname(batch[as.character(design[[pair$site$sample_key]])])
  )
  fit <- stats::lm(y ~ G_ + batch, data = data)
  # The contrast is written "G_<first> - G_<second>".
  levels_in_contrast <- sub("^G_", "", trimws(strsplit(cf$contrasts[[contrast]], "-", fixed = TRUE)[[1]]))
  at <- function(group) stats::predict(fit, newdata = data.frame(G_ = group, batch = "x"))
  expected <- unname(at(levels_in_contrast[[1]]) - at(levels_in_contrast[[2]]))
  estimated <- fitted$diff.site[fitted$site == site & fitted$contrast == contrast]
  expect_equal(estimated, expected, tolerance = 1e-8)
})

test_that("CF accepts the same model written differently", {
  fixture <- cf_pair()
  pair <- fixture$pair
  pair$site$formula <- "normalized_abundance ~ G_+batch"
  pair$protein$formula <- "normalized_abundance ~ batch + G_"
  expect_identical(.cf_model_string(pair), "~ G_ + batch")
})

test_that("the current CF is the observed site minus observed protein plus the sample's protein median", {
  fixture <- cf_pair()
  pair <- fixture$pair
  cf <- cf_from_pair(pair)
  expect_identical(.cf_model_string(pair), "~ G_")

  site <- pair$site$normalized_abundances
  protein <- pair$protein$normalized_abundances
  key <- pair$site$sample_key
  median_s <- tapply(
    protein$normalized_abundance[!grepl("^rev_", protein$protein_Id)],
    protein[[pair$protein$sample_key]][!grepl("^rev_", protein$protein_Id)],
    stats::median,
    na.rm = TRUE
  )
  joined <- merge(
    site[, c(key, "protein_Id", "site", "G_", "normalized_abundance")],
    stats::setNames(
      protein[, c(pair$protein$sample_key, "protein_Id", "normalized_abundance")],
      c(key, "protein_Id", "protein_abundance")
    ),
    by = c(key, "protein_Id")
  )
  joined$usage <- joined$normalized_abundance - joined$protein_abundance + unname(median_s[as.character(joined[[key]])])

  samples <- as.character(pair$site$obs[[key]])
  by_hand <- tidyr::pivot_wider(
    joined[, c("site", key, "usage")],
    names_from = tidyselect::all_of(key),
    values_from = "usage"
  )
  expect_equal(
    .cf_abundances(cf$wide_data, samples, pair$site$var$site),
    .cf_abundances(by_hand, samples, pair$site$var$site)
  )

  observed <- cf$results[cf$results$estimate_type == "observed", ]
  expect_gt(nrow(observed), 0)
  for (i in seq_len(nrow(observed))) {
    groups <- sub("^G_", "", trimws(strsplit(cf$contrasts[[observed$contrast[[i]]]], "-", fixed = TRUE)[[1]]))
    rows <- joined[joined$site == observed$site[[i]], ]
    expect_equal(
      observed$diff.site[[i]],
      mean(rows$usage[rows$G_ == groups[[1]]], na.rm = TRUE) -
        mean(rows$usage[rows$G_ == groups[[2]]], na.rm = TRUE),
      tolerance = 1e-10
    )
  }
})

test_that("CF refuses DEAs that used different models", {
  fixture <- cf_pair()
  pair <- fixture$pair
  pair$protein$formula <- "normalized_abundance ~ G_ + batch"
  expect_error(cf_from_pair(pair), "different models")
})

test_that("CF needs the imputedData layer of both DEAs", {
  fixture <- cf_pair()
  for (side in c("site", "protein")) {
    pair <- fixture$pair
    pair[[side]]$imputed_abundances <- NULL
    expect_error(cf_from_pair(pair), "no imputedData layer")
  }
})

test_that("CF asks for lm_impute when a DEA carries no imputation at all", {
  fixture <- cf_pair()
  for (side in c("site", "protein")) {
    pair <- fixture$pair
    pair[[side]]$imputed_abundances <- NULL
    pair[[side]]$imputation <- NULL
    expect_error(cf_from_pair(pair), "lm_impute model")
  }
})

test_that("the CF variants survive the MuData round trip", {
  fixture <- ptm_result_fixture()
  statistics <- read_ptm_h5mu(fixture$output, PTM_statistics)
  variants <- statistics$get_cf()$variants
  expect_named(variants, c("correct_first_protein_imputed", "correct_first_site_protein_imputed"))

  rewritten <- tempfile(fileext = ".h5mu")
  statistics$write_h5mu(rewritten)
  again <- read_ptm_h5mu(rewritten, PTM_statistics)$get_cf()$variants
  expect_equal(again, variants)

  container <- prolfquapp::read_h5mu(fixture$output)
  cf <- container$modalities$enriched_CF
  enriched_var <- as.data.frame(container$modalities$enriched$var)
  expect_equal(as.data.frame(cf$var), enriched_var[rownames(enriched_var) %in% rownames(cf$var), , drop = FALSE])
  expect_equal(dim(variants$correct_first_site_protein_imputed$abundances), c(nrow(cf$obs), nrow(cf$var)))
  expect_setequal(cf$layers_keys(), names(variants))
  expect_true(all(c("correct_first__a_vs_b", "correct_first_protein_imputed__a_vs_b") %in% cf$varm_keys()))
  expect_false(any(grepl("^correct_first_site_protein_imputed", cf$varm_keys())))
})

test_that("the paired-input stage carries no CF variant", {
  fixture <- ptm_result_fixture()
  inputs <- read_ptm_h5mu(fixture$output, PTM_statistics)$get_inputs()
  enriched <- inputs$get_enriched()
  expect_false(any(
    c("correct_first_protein_imputed", "correct_first_site_protein_imputed") %in% enriched$layers_keys()
  ))
  expect_false(any(grepl("^correct_first", enriched$varm_keys())))
})
