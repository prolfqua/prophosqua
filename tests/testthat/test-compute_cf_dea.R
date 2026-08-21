test_that("derive_contrasts takes explicit contrasts as given", {
  annot <- data.frame(
    Contrast = c("G_a - G_b", "G_c - G_b"),
    ContrastName = c("a_vs_b", "c_vs_b")
  )
  expect_equal(
    derive_contrasts(annot),
    c(a_vs_b = "G_a - G_b", c_vs_b = "G_c - G_b")
  )
})

test_that("derive_contrasts drops empty explicit rows", {
  annot <- data.frame(
    Contrast = c("G_a - G_b", NA, ""),
    ContrastName = c("a_vs_b", "empty", "blank")
  )
  expect_equal(derive_contrasts(annot), c(a_vs_b = "G_a - G_b"))
})

test_that("derive_contrasts pairs every treatment with every control", {
  annot <- data.frame(
    Group = c("trt1", "trt2", "ctrl", "trt1"),
    Control = c("T", "T", "C", "T")
  )
  expect_equal(
    derive_contrasts(annot),
    c(trt1_vs_ctrl = "G_trt1 - G_ctrl", trt2_vs_ctrl = "G_trt2 - G_ctrl")
  )
})

test_that("derive_contrasts matches the column names case-insensitively", {
  annot <- data.frame(GROUP = c("trt", "ctrl"), CONTROL = c("t", "c"))
  expect_equal(derive_contrasts(annot), c(trt_vs_ctrl = "G_trt - G_ctrl"))
})

test_that("derive_contrasts refuses an annotation that defines nothing", {
  # The columns are there but no group is marked as control or treatment: the
  # shape of a template nobody filled in.
  annot <- data.frame(Group = c("g1", "g2"), Control = c("", ""))
  expect_error(
    derive_contrasts(annot, "template.tsv"),
    "No contrasts could be derived from template.tsv"
  )
})

test_that("derive_contrasts refuses an annotation with neither layout", {
  expect_error(
    derive_contrasts(data.frame(sample = "s1", condition = "c1")),
    "must have either"
  )
})
