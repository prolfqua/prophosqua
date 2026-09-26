make_dea_dir <- function(files = character(), name = "DEA_test") {
  dea_dir <- file.path(tempfile(pattern = name))
  results <- file.path(dea_dir, "Results_WU_test")
  dir.create(results, recursive = TRUE)
  if (length(files) > 0) {
    file.create(file.path(results, files))
  }
  list(dea_dir = dea_dir, results = results)
}

test_that("get_dea_file finds the file under Results_WU_ and names a missing one", {
  d <- make_dea_dir("AnnData.h5ad")
  expect_identical(get_dea_file(d$dea_dir, "AnnData.h5ad"), file.path(d$results, "AnnData.h5ad"))
  expect_error(get_dea_file(d$dea_dir, "other.h5ad"), "No other.h5ad found")
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
