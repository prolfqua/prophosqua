test_that("enrichment vignette fixture contains nine portable JSON documents", {
  path <- test_path("../../inst/extdata/ptm_results_example.h5mu")
  if (!file.exists(path)) {
    path <- system.file(
      "extdata",
      "ptm_results_example.h5mu",
      package = "prophosqua"
    )
  }
  expect_true(file.exists(path))
  result <- read_ptm_h5mu(path, PTM_results)
  documents <- result$get_enrichment_documents()
  expected <- c(
    "PTMSEA__DPA",
    "PTMSEA__DPU",
    "PTMSEA__CF",
    "KinaseGSEA__DPA",
    "KinaseGSEA__DPU",
    "KinaseGSEA__CF",
    "MEA__DPA",
    "MEA__DPU",
    "MEA__CF"
  )

  expect_named(documents, expected)
  expect_true(all(vapply(
    documents,
    function(document) identical(document$version, "1.2.0"),
    logical(1)
  )))

  category_names <- c(PTMSEA = "PTM-SEA", KinaseGSEA = "KinaseLib", MEA = "MEA")
  for (method in names(category_names)) {
    category <- category_names[[method]]
    for (analysis in c("DPA", "DPU", "CF")) {
      document <- result$get_enrichment_document(method, analysis)
      payload <- jsonlite::fromJSON(document$json, simplifyVector = FALSE)
      decoded <- protsea::decode_gsea_json(document$json)

      expect_named(decoded, names(payload$data))
      expect_length(decoded, 2L)
      for (contrast in names(decoded)) {
        object <- decoded[[contrast]][[category]]
        category_payload <- payload$data[[contrast]]$categories[[category]]
        rank_payload <- payload$rank_lists[[contrast]]$entries

        expect_s4_class(object, "gseaResult")
        expect_equal(nrow(object@result), length(category_payload$terms))
        expect_identical(names(object@geneList), names(rank_payload))
        expect_equal(unname(object@geneList), unname(unlist(rank_payload)))
        expect_gt(length(object@geneSets), 0L)
        expect_true(all(
          c("exponent", "minGSSize", "maxGSSize") %in%
            names(object@params)
        ))
        native <- category_payload$gsea_result
        for (term_id in object@result$ID) {
          reproduced <- expected_gsea_trace(
            object@geneList,
            object@geneSets[[term_id]],
            exponent = object@params$exponent
          )
          expect_equal(
            unlist(native$running_scores[[term_id]], use.names = FALSE),
            reproduced$runningScore,
            tolerance = 1e-12
          )
          expect_equal(
            unlist(native$hit_indices[[term_id]], use.names = FALSE),
            which(reproduced$position == 1L)
          )
        }
      }
      plot <- enrichplot::gseaplot2(
        decoded[[1L]][[category]],
        geneSetID = 1L,
        pvalue_table = FALSE
      )
      expect_s3_class(plot, "gglist")
    }
  }
})

test_that("enrichment vignette has the required three-level tab structure", {
  path <- test_path("../../vignettes/ptm_enrichment.qmd")
  if (!file.exists(path)) {
    path <- system.file("doc", "ptm_enrichment.qmd", package = "prophosqua")
  }
  skip_if(
    !nzchar(path) || !file.exists(path),
    "package installed without vignette sources"
  )
  source <- readLines(path, warn = FALSE)
  fence <- grepl("^```", source)
  outside_code <- cumsum(fence) %% 2L == 0L & !fence
  top_starts <- which(outside_code & grepl("^# ", source))
  top_tabs <- sub("^# ", "", source[top_starts])

  expect_identical(
    top_tabs,
    c("Overview", "DPA", "DPU", "CorrectFirst DPU")
  )

  overview_end <- top_starts[[2L]] - 1L
  overview_source <- source[top_starts[[1L]]:overview_end]
  overview_fence <- grepl("^```", overview_source)
  overview_outside <- cumsum(overview_fence) %% 2L == 0L & !overview_fence
  overview_tabs <- sub(
    "^## ",
    "",
    grep("^## ", overview_source[overview_outside], value = TRUE)
  )
  expect_identical(
    overview_tabs,
    c("Summary", "Report provenance", "R session info")
  )

  for (analysis in c("DPA", "DPU", "CorrectFirst DPU")) {
    start <- top_starts[[match(analysis, top_tabs)]]
    following <- top_starts[top_starts > start]
    end <- if (length(following)) following[[1L]] - 1L else length(source)
    analysis_source <- source[start:end]
    analysis_fence <- grepl("^```", analysis_source)
    analysis_outside <- cumsum(analysis_fence) %% 2L == 0L & !analysis_fence
    method_starts <- which(
      analysis_outside & grepl("^## ", analysis_source)
    )
    methods <- sub("^## ", "", analysis_source[method_starts])
    expect_identical(methods, c("PTM-SEA", "Kinase GSEA", "MEA"))

    for (i in seq_along(method_starts)) {
      method_start <- method_starts[[i]]
      method_end <- if (i < length(method_starts)) {
        method_starts[[i + 1L]] - 1L
      } else {
        length(analysis_source)
      }
      views <- sub(
        "^### ",
        "",
        grep(
          "^### ",
          analysis_source[method_start:method_end],
          value = TRUE
        )
      )
      expect_identical(
        views,
        c(
          "Summary",
          "Dot plot",
          "Heatmap",
          "Volcano",
          "Running score",
          "Rank distributions",
          "Gene-set network",
          "Term similarity",
          "Results"
        )
      )
    }
  }

  expect_true(any(grepl("get_enrichment_document", source, fixed = TRUE)))
  expect_true(any(grepl("decode_gsea_json", source, fixed = TRUE)))
  expect_true(any(grepl("gseaplot2", source, fixed = TRUE)))
  expect_true(any(grepl("ridgeplot", source, fixed = TRUE)))
  expect_true(any(grepl("cnetplot", source, fixed = TRUE)))
  expect_true(any(grepl("pairwise_termsim", source, fixed = TRUE)))
  expect_true(any(grepl("emapplot", source, fixed = TRUE)))
  expect_true(any(grepl("treeplot", source, fixed = TRUE)))
  expect_false(any(grepl("read_excel|readRDS|read_xlsx", source)))
})
