# The GSEAResult JSON export: string_gsea structure, canonical window ids.

fake_gsea_class <- methods::setClass(
  "FakeGseaResult",
  representation = methods::representation(
    geneList = "numeric",
    result = "data.frame",
    geneSets = "list"
  )
)

make_fake_gsea <- function() {
  gene_list <- c(
    "AAAAAAASAAAAAAA-p" = 2.5,
    "BBBBBBBSBBBBBBB-p" = 1.5,
    "CCCCCCCSCCCCCCC-p" = -2.0
  )
  fake_gsea_class(
    geneList = gene_list,
    result = data.frame(
      ID = "KINASE-PSP_CDK2",
      Description = "KINASE-PSP_CDK2",
      NES = -1.9,
      p.adjust = 0.003,
      core_enrichment = "CCCCCCCSCCCCCCC-p",
      stringsAsFactors = FALSE
    ),
    geneSets = list(
      "KINASE-PSP_CDK2" = c(
        "AAAAAAASAAAAAAA-p",
        "CCCCCCCSCCCCCCC-p",
        "ZZZZZZZSZZZZZZZ-p" # Not in the ranked data: full set, not mapped.
      )
    )
  )
}

test_that("gsea_result_data builds pool and terms with canonical window ids", {
  results <- list(A_vs_B = make_fake_gsea())

  out <- gsea_result_data(results, category = "PTM-SEA")

  contrast <- out$data$A_vs_B
  expect_equal(contrast$contrast, "A_vs_B")
  expect_equal(
    names(contrast$gene_pool),
    c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB", "CCCCCCCSCCCCCCC")
  )
  first <- contrast$gene_pool[["AAAAAAASAAAAAAA"]]
  expect_equal(first$input_label, "AAAAAAASAAAAAAA-p")
  expect_equal(first$input_value, 2.5)
  expect_equal(first$rank, 1)

  term <- contrast$categories[["PTM-SEA"]]$terms[[1]]
  expect_equal(term$term_id, "KINASE-PSP_CDK2")
  expect_equal(term$category, "PTM-SEA")
  expect_equal(term$direction, "bottom")
  expect_equal(term$fdr, 0.003)
  # Full mapped set: intersection of the gene set with the ranked data.
  expect_setequal(unlist(term$gene_ids), c("AAAAAAASAAAAAAA", "CCCCCCCSCCCCCCC"))
  expect_equal(term$genes_mapped, 2)
  expect_equal(term$genes_in_set, 3)
  expect_equal(unlist(term$leading_edge_ids), "CCCCCCCSCCCCCCC")

  expect_equal(out$rank_lists$A_vs_B$entries[["BBBBBBBSBBBBBBB-p"]], 1.5)
})

test_that("mea_gsea_result_data maps leading substrates onto the pool", {
  ranks <- list(
    A_vs_B = c("AAAAAAASAAAAAAA" = 2.5, "BBBBBBBSBBBBBBB" = 1.5, "CCCCCCCSCCCCCCC" = -2.0)
  )
  mea_clean <- data.frame(
    contrast = "A_vs_B",
    kinase = "CDK2",
    NES = 2.1,
    FDR = 0.01,
    Leading.substrates = "AAAAAAAsAAAAAAA;BBBBBBBsBBBBBBB",
    stringsAsFactors = FALSE
  )
  term2gene <- data.frame(
    term = rep("CDK2", 3),
    # Kinase-library writes the phospho residue lower-case.
    gene = c("AAAAAAAsAAAAAAA", "BBBBBBBsBBBBBBB", "ZZZZZZZsZZZZZZZ"),
    stringsAsFactors = FALSE
  )

  out <- mea_gsea_result_data(mea_clean, ranks, term2gene)

  term <- out$data$A_vs_B$categories$MEA$terms[[1]]
  expect_equal(term$term_id, "CDK2")
  expect_equal(term$direction, "top")
  expect_equal(term$method, "mea")
  expect_setequal(unlist(term$gene_ids), c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB"))
  expect_setequal(unlist(term$leading_edge_ids), c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB"))
  expect_equal(term$genes_in_set, 3)
  expect_equal(names(out$data$A_vs_B$gene_pool), names(ranks$A_vs_B))
})

test_that("read_mea_ranks reads the per-contrast rnk files", {
  dir <- file.path(tempdir(), "mea_rnk_test")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  # The header row is what prep_kinaselib.R writes; reading it as data would
  # add a bogus rank entry and coerce every statistic to a string.
  writeLines(
    c(
      "SequenceWindow\tstatistic.site",
      "AAAAAAASAAAAAAA\t2.5",
      "CCCCCCCSCCCCCCC\t-2.0"
    ),
    file.path(dir, "DPA_MEA_A_vs_B.rnk")
  )

  ranks <- read_mea_ranks(dir)

  expect_equal(names(ranks), "A_vs_B")
  expect_equal(ranks$A_vs_B, c("AAAAAAASAAAAAAA" = 2.5, "CCCCCCCSCCCCCCC" = -2.0))
  expect_type(ranks$A_vs_B, "double")
  empty_dir <- file.path(tempdir(), "mea_rnk_empty")
  dir.create(empty_dir, showWarnings = FALSE)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)
  expect_error(read_mea_ranks(empty_dir), "No .*rnk files")
})

test_that("write_gsea_result_json round-trips and keeps empty pools as objects", {
  results <- list(A_vs_B = make_fake_gsea())
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)

  write_gsea_result_json(gsea_result_data(results, category = "PTM-SEA"), path)

  parsed <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expect_named(parsed, c("data", "rank_lists"))
  contrast <- parsed$data$A_vs_B
  expect_equal(contrast$contrast, "A_vs_B")
  pool_ids <- names(contrast$gene_pool)
  term <- contrast$categories[["PTM-SEA"]]$terms[[1]]
  # The contract stringdbpy enforces: every term member resolves in the pool.
  expect_true(all(unlist(term$gene_ids) %in% pool_ids))
  expect_true(all(unlist(term$leading_edge_ids) %in% unlist(term$gene_ids)))
})
