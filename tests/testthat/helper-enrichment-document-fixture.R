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

make_ptm_enrichment_fixture <- function(empty = FALSE) {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(paths$site),
    anndataR::read_h5ad(paths$protein),
    parameters = list(
      run_kinase = TRUE,
      kinaselib = list(kin_type = "ST", threshold = 90, permutations = 100),
      analyses = list(
        dpa = list(subdir = "DPA"),
        dpu = list(subdir = "DPU"),
        cf = list(subdir = "CF")
      )
    )
  )
  statistics <- suppressWarnings(PTM_statistics$new(inputs))
  branches <- unlist(
    lapply(c("DPA", "DPU", "CF"), make_ptm_enrichment_branches, statistics = statistics, empty = empty),
    recursive = FALSE
  )
  final <- PTM_results$new(statistics, branches)
  path <- tempfile(fileext = ".h5mu")
  final$write_h5mu(path)
  list(final = final, path = path)
}

make_ptm_enrichment_branches <- function(analysis, statistics, empty) {
  ptm_gsea <- make_stage_gsea("PTM-SEA", empty)
  ptmsea_result <- list(results = list(a_vs_b = ptm_gsea), all_clean = data.frame())
  ptmsea <- PTMSEA$new(statistics, analysis, ptmsea_result)

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

  kinase_gsea <- make_stage_gsea("KinaseLib", empty)
  kinase_result <- list(
    gsea_results = list(a_vs_b = kinase_gsea),
    all_results = data.frame(),
    gsea_info = data.frame(value = 1)
  )
  kinase <- KinaseGSEA$new(assignments, analysis, kinase_result)

  mea_json <- protsea::gsea_result_json_text(
    protsea::gsea_result_data(list(a_vs_b = kinase_gsea), category = "MEA", method = "gseapy")
  )
  mea_clean <- data.frame(
    contrast = "a_vs_b",
    kinase = "CDK2",
    NES = 2.1,
    pvalue = 0.001,
    FDR = 0.01,
    n_leading = 2,
    set_size = 3,
    Leading.substrates = "AAAAAAAsAAAAAAA;BBBBBBBsBBBBBBB"
  )
  if (empty) {
    mea_clean <- mea_clean[FALSE, , drop = FALSE]
  }
  motif <- MotifEnrichment$new(
    assignments,
    analysis,
    list(mea_results = mea_clean, gsea_json = mea_json)
  )
  mea <- MEA$new(
    motif,
    analysis,
    list(mea_clean = mea_clean, summary_df = data.frame(contrast = "a_vs_b", total_kinases = nrow(mea_clean)))
  )
  list(ptmsea, kinase, mea)
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
