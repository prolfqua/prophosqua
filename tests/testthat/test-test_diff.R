make_dpu_row <- function(
  protein_id,
  contrast,
  diff,
  moderated_se,
  moderated_df,
  raw_se,
  raw_df
) {
  data.frame(
    protein_Id = protein_id,
    contrast = contrast,
    diff = diff,
    std.error = moderated_se,
    df = moderated_df,
    std.error.unmoderated = raw_se,
    df.unmoderated = raw_df
  )
}

test_that("test_diff unmoderated reproduces a known Welch test", {
  site <- make_dpu_row("P1", "a_vs_b", 2, 0.2, 100, 0.3, 4)
  protein <- make_dpu_row("P1", "a_vs_b", 0.5, 0.25, 120, 0.4, 6)

  result <- test_diff(
    site,
    protein,
    join_column = c("protein_Id", "contrast"),
    variant = "unmoderated"
  )

  expect_equal(result$diff_diff, 1.5)
  expect_equal(result$SE_I, 0.5)
  expect_equal(result$df_I, 9.933774834437084, tolerance = 1e-12)
  expect_equal(result$tstatistic_I, 3)
  expect_equal(result$pValue_I, 0.01343834065037431, tolerance = 1e-12)
})

test_that("test_diff moderated uses the reported moderated pair", {
  site <- make_dpu_row("P1", "a_vs_b", 2, 0.2, 100, 0.3, 4)
  protein <- make_dpu_row("P1", "a_vs_b", 0.5, 0.25, 120, 0.4, 6)

  result <- test_diff(
    site,
    protein,
    join_column = c("protein_Id", "contrast"),
    variant = "moderated"
  )

  expect_equal(result$SE_I, sqrt(0.2^2 + 0.25^2))
  expect_equal(result$tstatistic_I, result$diff_diff / result$SE_I)
})

test_that("unmoderated invalid df rows are untestable and excluded from BH", {
  site <- do.call(
    rbind,
    list(
      make_dpu_row("P1", "a_vs_b", 2.0, 0.2, 100, 0.3, 4),
      make_dpu_row("P2", "a_vs_b", 1.5, 0.2, 100, 0.3, 0),
      make_dpu_row("P3", "a_vs_b", 1.0, 0.2, 100, 0.3, 4)
    )
  )
  protein <- do.call(
    rbind,
    list(
      make_dpu_row("P1", "a_vs_b", 0.5, 0.25, 120, 0.4, 6),
      make_dpu_row("P2", "a_vs_b", 0.5, 0.25, 120, 0.4, 6),
      make_dpu_row("P3", "a_vs_b", 0.5, 0.25, 120, 0.4, Inf)
    )
  )

  result <- test_diff(
    site,
    protein,
    join_column = c("protein_Id", "contrast"),
    variant = "unmoderated"
  )
  testable <- result$protein_Id %in% "P1"

  expect_true(all(is.na(result$df_I[!testable])))
  expect_true(all(is.na(result$tstatistic_I[!testable])))
  expect_true(all(is.na(result$pValue_I[!testable])))
  expect_true(all(is.na(result$FDR_I[!testable])))
  expect_equal(result$FDR_I[testable], result$pValue_I[testable])
})

test_that("test_diff rejects old DEA output for both variants", {
  site <- data.frame(
    protein_Id = "P1",
    contrast = "a_vs_b",
    diff = 1,
    std.error = 0.2,
    df = 10
  )

  for (variant in c("moderated", "unmoderated")) {
    expect_error(
      test_diff(
        site,
        site,
        join_column = c("protein_Id", "contrast"),
        variant = variant
      ),
      "std.error.unmoderated.*df.unmoderated"
    )
  }
})
