make_stage_gsea <- function(category, empty = FALSE) {
  genes <- c(
    "AAAAAAASAAAAAAA" = 2.5,
    "BBBBBBBSBBBBBBB" = 1.5,
    "CCCCCCCSCCCCCCC" = -2
  )
  result <- data.frame(
    ID = "CDK2",
    Description = "CDK2",
    setSize = 2L,
    enrichmentScore = -0.7,
    NES = -1.9,
    pvalue = 0.001,
    qvalues = 0.003,
    rank = 3L,
    leading_edge = "tags=50%, list=33%, signal=75%",
    p.adjust = 0.003,
    core_enrichment = "CCCCCCCSCCCCCCC",
    # As clusterProfiler::GSEA names its rows.
    row.names = "CDK2",
    stringsAsFactors = FALSE
  )
  if (empty) {
    result <- result[FALSE, , drop = FALSE]
  }
  methods::new(
    methods::getClass("gseaResult", where = asNamespace("DOSE")),
    params = list(exponent = 1.5),
    organism = "unknown",
    setType = category,
    keytype = "sequence",
    readable = FALSE,
    geneList = genes,
    result = result,
    geneSets = list(CDK2 = c("AAAAAAASAAAAAAA", "CCCCCCCSCCCCCCC"))
  )
}

expected_gsea_trace <- function(ranks, members, exponent) {
  hits <- names(ranks) %in% members
  weights <- abs(ranks)^exponent
  increments <- ifelse(
    hits,
    weights / sum(weights[hits]),
    -1 / sum(!hits)
  )
  data.frame(
    runningScore = cumsum(increments),
    position = as.integer(hits)
  )
}

# Final PTM results as a pipeline run leaves them: the statistics, the
# enrichment files of each analysis and the final MuData beside them.
make_ptm_enrichment_fixture <- function(empty = FALSE) {
  paths <- anndata_pair_fixture()
  parameters <- list(
    run_kinase = TRUE,
    kinaselib = list(kin_type = "ST", threshold = 90, permutations = 100),
    analyses = list(
      dpa = list(subdir = "DPA"),
      dpu = list(subdir = "DPU"),
      cf = list(subdir = "CF")
    )
  )
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(paths$site),
    anndataR::read_h5ad(paths$protein),
    parameters = parameters
  )
  statistics <- suppressWarnings(PTM_statistics$new(inputs))
  root <- tempfile("ptm_results_")
  dir.create(root)
  statistics_path <- file.path(root, "PTM_statistics.h5mu")
  statistics$write_h5mu(statistics_path)
  files <- .ptm_enrichment_files(parameters, root)
  stages <- list()
  for (analysis in c("DPA", "DPU", "CF")) {
    for (stage in make_ptm_enrichment_stages(analysis, statistics, empty)) {
      key <- .ptm_key(class(stage)[1L], analysis)
      .write_ptm_stage_file(stage, files[[key]], .ptm_file_sha256(statistics_path))
      stages[[key]] <- stage
    }
  }
  path <- file.path(root, "PTM_results.h5mu")
  assemble_ptm_results(statistics_path, path)
  list(
    final = read_ptm_h5mu(path, PTM_results),
    path = path,
    root = root,
    statistics_path = statistics_path,
    files = files,
    stages = stages
  )
}

make_ptm_enrichment_stages <- function(analysis, statistics, empty) {
  ptmsea <- PTMSEA$new(statistics, analysis, .ptmsea_result(list(a_vs_b = make_stage_gsea("PTM-SEA", empty))))

  rank_table <- data.frame(
    SequenceWindow = c("AAAAAAASAAAAAAA", "BBBBBBBSBBBBBBB", "CCCCCCCSCCCCCCC"),
    statistic.site = c(2.5, 1.5, -2)
  )
  kinase_inputs <- KinaseInputs$new(
    statistics,
    analysis,
    list(seqwindows = rank_table["SequenceWindow"], ranks = list(a_vs_b = rank_table))
  )
  assignments <- KinaseAssignments$new(
    kinase_inputs,
    analysis,
    list(term2gene = data.frame(term = "CDK2", gene = rank_table$SequenceWindow))
  )
  kinase <- KinaseGSEA$new(
    assignments,
    analysis,
    .kinasegsea_result(list(a_vs_b = make_stage_gsea("KinaseLib", empty)))
  )
  mea <- MEA$new(assignments, analysis, .mea_result(list(a_vs_b = make_stage_gsea("MEA", empty))))
  list(ptmsea, kinase_inputs, assignments, kinase, mea)
}

ptm_enrichment_fixture <- local({
  values <- list()
  function(empty = FALSE) {
    key <- as.character(empty)
    if (is.null(values[[key]])) {
      values[[key]] <<- make_ptm_enrichment_fixture(empty)
    }
    values[[key]]
  }
})
