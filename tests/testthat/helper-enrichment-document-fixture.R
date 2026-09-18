make_stage_gsea <- function(category, empty = FALSE) {
  genes <- c(
    "AAAAAAASAAAAAAA-p" = 2.5,
    "BBBBBBBSBBBBBBB-p" = 1.5,
    "CCCCCCCSCCCCCCC-p" = -2
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
    core_enrichment = "CCCCCCCSCCCCCCC-p",
    stringsAsFactors = FALSE
  )
  if (empty) {
    result <- result[FALSE, , drop = FALSE]
  }
  methods::new(
    "gseaResult",
    params = list(exponent = 1.5),
    organism = "unknown",
    setType = category,
    keytype = "sequence",
    readable = FALSE,
    geneList = genes,
    result = result,
    geneSets = list(CDK2 = c("AAAAAAASAAAAAAA-p", "CCCCCCCSCCCCCCC-p"))
  )
}

make_ptm_enrichment_fixture <- function(empty = FALSE) {
  paths <- anndata_pair_fixture()
  inputs <- DEA_enriched_total$new(
    anndataR::read_h5ad(paths$site),
    anndataR::read_h5ad(paths$protein),
    parameters = list(
      run_kinase = TRUE,
      analyses = list(
        dpa = list(subdir = "DPA"),
        dpu = list(subdir = "DPU"),
        cf = list(subdir = "CF")
      )
    )
  )
  statistics <- inputs$build(DPA_DPU)$build(PTM_statistics, cf = suppressWarnings(inputs$build(CF)))
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
  ptmsea_result <- list(
    results = list(a_vs_b = ptm_gsea),
    ranks = list(a_vs_b = ptm_gsea@geneList),
    all_clean = data.frame(),
    pathways = ptm_gsea@geneSets,
    data_info = data.frame(value = 1),
    ptmsigdb_summary = data.frame(value = 1),
    overlap_stats = data.frame(value = 1),
    n_overlap = 1L,
    n_our_sites = 3L,
    prep_info = data.frame(value = 1),
    results_info = data.frame(value = 1),
    has_results = !empty,
    analysis_inputs = list(source = "fixture")
  )
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
    ranks = list(a_vs_b = kinase_gsea@geneList),
    all_results = data.frame(),
    term2gene = assignments$get_results()$term2gene,
    term2gene_df = assignments$get_results()$term2gene,
    n_our_sequences = 3L,
    n_overlap_seqs = 2L,
    data_info = data.frame(value = 1),
    kl_info = data.frame(value = 1),
    assignment_stats = data.frame(value = 1),
    kinase_stats = data.frame(value = 1),
    ranks_info = data.frame(value = 1),
    gsea_info = data.frame(value = 1),
    has_results = !empty,
    analysis_inputs = list(source = "fixture")
  )
  kinase <- KinaseGSEA$new(assignments, analysis, kinase_result)

  motif <- MotifEnrichment$new(assignments, analysis, list(mea_results = data.frame(value = 1)))
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
  mea <- MEA$new(
    motif,
    analysis,
    list(
      mea_clean = mea_clean,
      summary_df = data.frame(contrast = "a_vs_b", total_kinases = nrow(mea_clean)),
      n_files = 1L,
      has_results = !empty,
      analysis_inputs = list(source = "fixture")
    )
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
