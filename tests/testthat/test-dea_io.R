make_dea_dir <- function(files = character(), name = "DEA_test") {
  dea_dir <- file.path(tempfile(pattern = name))
  results <- file.path(dea_dir, "Results_WU_test")
  dir.create(results, recursive = TRUE)
  if (length(files) > 0) {
    file.create(file.path(results, files))
  }
  list(dea_dir = dea_dir, results = results)
}

test_that("get_dea_xlsx prefers the DE_ prefixed workbook", {
  d <- make_dea_dir(c("zzz_other.xlsx", "DE_run.xlsx"))
  expect_equal(basename(suppressMessages(get_dea_xlsx(d$dea_dir))), "DE_run.xlsx")
})

test_that("get_dea_xlsx falls back to a workbook without the prefix", {
  d <- make_dea_dir("summary.xlsx")
  expect_equal(basename(suppressMessages(get_dea_xlsx(d$dea_dir))), "summary.xlsx")
})

test_that("get_dea_xlsx searches the directory itself when Results_WU_ is absent", {
  dea_dir <- tempfile(pattern = "DEA_flat")
  dir.create(dea_dir)
  file.create(file.path(dea_dir, "DE_flat.xlsx"))
  expect_equal(basename(suppressMessages(get_dea_xlsx(dea_dir))), "DE_flat.xlsx")
})

test_that("get_dea_xlsx errors when there is no workbook", {
  d <- make_dea_dir()
  expect_error(get_dea_xlsx(d$dea_dir), "No Excel file found")
})

test_that("get_dea_parquet and get_dea_yaml resolve their files", {
  d <- make_dea_dir(c("lfqdata_normalized.parquet", "lfqdata.yaml"))
  expect_equal(basename(get_dea_parquet(d$dea_dir)), "lfqdata_normalized.parquet")
  expect_equal(basename(get_dea_yaml(d$dea_dir)), "lfqdata.yaml")
})

test_that("get_dea_file names the missing thing in its error", {
  d <- make_dea_dir()
  expect_error(get_dea_parquet(d$dea_dir), "No parquet file found")
  expect_error(get_dea_yaml(d$dea_dir), "No yaml file found")
})

test_that("get_sample_name_column reads the declared column", {
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines("sample_name: sampleName", yaml_file)
  expect_equal(get_sample_name_column(yaml_file), "sampleName")
})

test_that("get_sample_name_column rejects a missing or empty declaration", {
  empty <- tempfile(fileext = ".yaml")
  writeLines("other_key: 1", empty)
  expect_error(get_sample_name_column(empty), "No valid sample_name")

  blank <- tempfile(fileext = ".yaml")
  writeLines("sample_name: ''", blank)
  expect_error(get_sample_name_column(blank), "No valid sample_name")
})

test_that("get_dea_sample_name_column reads through the DEA directory", {
  d <- make_dea_dir()
  writeLines("sample_name: sampleName", file.path(d$results, "lfqdata.yaml"))
  expect_equal(get_dea_sample_name_column(d$dea_dir), "sampleName")
})

test_that("canonicalize_dea_sample_column renames the declared column", {
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines("sample_name: sampleName", yaml_file)
  data <- data.frame(sampleName = c("a", "b"), abundance = c(1, 2))

  out <- canonicalize_dea_sample_column(data, yaml_file)
  expect_named(out, c("Name", "abundance"))
  expect_equal(out$Name, c("a", "b"))
})

test_that("canonicalize_dea_sample_column is a no-op when already canonical", {
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines("sample_name: Name", yaml_file)
  data <- data.frame(Name = "a", abundance = 1)
  expect_identical(canonicalize_dea_sample_column(data, yaml_file), data)
})

test_that("canonicalize_dea_sample_column refuses to clobber an existing column", {
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines("sample_name: sampleName", yaml_file)
  data <- data.frame(sampleName = "a", Name = "b")
  expect_error(canonicalize_dea_sample_column(data, yaml_file), "both columns are present")
})

test_that("canonicalize_dea_sample_column errors when the declared column is absent", {
  yaml_file <- tempfile(fileext = ".yaml")
  writeLines("sample_name: sampleName", yaml_file)
  expect_error(
    canonicalize_dea_sample_column(data.frame(other = 1), yaml_file),
    "is absent from normalized DEA data"
  )
})

test_that("canonicalize_uniprot_ids takes the accession and leaves bare ids alone", {
  data <- data.frame(protein_Id = c("sp|P12345|PROT_HUMAN", "Q67890"))
  expect_equal(canonicalize_uniprot_ids(data)$protein_Id, c("P12345", "Q67890"))
})

test_that("canonicalize_uniprot_ids rejects a mapping that is not one-to-one", {
  data <- data.frame(protein_Id = c("sp|P12345|A_HUMAN", "tr|P12345|B_HUMAN"))
  expect_error(canonicalize_uniprot_ids(data), "not one-to-one")
})

test_that("canonicalize_uniprot_ids keeps repeated rows of the same identifier", {
  data <- data.frame(protein_Id = rep("sp|P12345|PROT_HUMAN", 3))
  expect_equal(canonicalize_uniprot_ids(data)$protein_Id, rep("P12345", 3))
})

test_that("canonicalize_uniprot_ids errors on a missing column", {
  expect_error(canonicalize_uniprot_ids(data.frame(x = 1)), "column is absent")
})

test_that("get_dea_ptm_site_info reads metadata from the abundance sheet", {
  skip_if_not_installed("writexl")
  d <- make_dea_dir()
  writexl::write_xlsx(
    list(
      normalized_abundances = data.frame(
        site = c("P12345_S10~PEP", "P12345_T20~PEP"),
        posInProtein = c(10L, 20L),
        modAA = c("S", "T"),
        SequenceWindow = c("AAAAAAASAAAAAAA", "AAAAAAATAAAAAAA"),
        protein_Id = "P12345",
        gene_name = "GENE",
        protein_length = 100L
      )
    ),
    file.path(d$results, "DE_sites.xlsx")
  )

  info <- suppressMessages(get_dea_ptm_site_info(d$dea_dir))
  expect_equal(nrow(info), 2)
  expect_true(all(
    c("site", "posInProtein", "modAA", "SequenceWindow", "protein_Id", "gene_name", "protein_length") %in% names(info)
  ))
  expect_equal(info$modAA, c("S", "T"))
})

test_that("get_dea_ptm_site_info rejects metadata that is not unique by site", {
  skip_if_not_installed("writexl")
  d <- make_dea_dir()
  writexl::write_xlsx(
    list(
      normalized_abundances = data.frame(
        site = c("P12345_S10~PEP", "P12345_S10~PEP"),
        posInProtein = c(10L, 11L),
        modAA = c("S", "S"),
        SequenceWindow = c("AAAAAAASAAAAAAA", "CCCCCCCSCCCCCCC"),
        protein_Id = "P12345",
        gene_name = "GENE",
        protein_length = 100L
      )
    ),
    file.path(d$results, "DE_sites.xlsx")
  )

  expect_error(
    suppressMessages(get_dea_ptm_site_info(d$dea_dir)),
    "not unique by site"
  )
})
