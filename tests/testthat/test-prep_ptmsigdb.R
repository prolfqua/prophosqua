test_that("write_gmt writes one tab separated line per set", {
  gmt <- tempfile(fileext = ".gmt")
  write_gmt(list(A = c("S1", "S2"), B = "S3"), gmt)

  expect_equal(
    readLines(gmt),
    c("A\tNA\tS1\tS2", "B\tNA\tS3")
  )
})

test_that("write_gmt round-trips through the reader PTM-SEA uses", {
  skip_if_not_installed("fgsea")

  pathways <- list(KINASE_A = c("SITE1", "SITE2"), KINASE_B = "SITE3")
  gmt <- tempfile(fileext = ".gmt")
  write_gmt(pathways, gmt)

  expect_equal(fgsea::gmtPathways(gmt), pathways)
})

test_that("write_gmt returns the path it wrote", {
  gmt <- tempfile(fileext = ".gmt")
  expect_equal(write_gmt(list(A = "S1"), gmt), gmt)
})

test_that("report_ptmsigdb_sources counts only the sub-sources that occur", {
  merged <- stats::setNames(
    list("S1", "S2", "S3"),
    c("KINASE-PSP_AKT1", "KINASE-PSP_MTOR", "PATH-NP_TGFB")
  )
  filtered <- merged["KINASE-PSP_AKT1"]

  expect_message(
    report_ptmsigdb_sources(merged, filtered),
    "KINASE-PSP: 2 total, 1 kept"
  )
  expect_message(
    report_ptmsigdb_sources(merged, filtered),
    "PATH-NP: 1 total, 0 kept"
  )
})
