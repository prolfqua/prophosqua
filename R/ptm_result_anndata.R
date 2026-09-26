# A PTM result table is stored as one varm data frame per contrast, one row per
# site of the modality, NA where the site has no result. The frame keeps its
# site key but not the rest of the site annotation, which reading it back joins
# on that key; uns keeps only the column order of each table.

.ptm_analysis_payload <- function(data, site_var, prefix, contrasts) {
  .require_columns(data, "contrast", paste(prefix, "results"))
  block_columns <- setdiff(names(data), c(setdiff(names(site_var), "site"), "contrast"))
  keys <- .ptm_key(prefix, contrasts)
  frames <- lapply(contrasts, function(contrast) {
    .ptm_result_frame(data[data$contrast == contrast, block_columns, drop = FALSE], site_var)
  })
  list(
    values = stats::setNames(frames, keys),
    order = stats::setNames(rep(list(names(data)), length(keys)), keys)
  )
}

# The key a PTM result, stage or document is stored under: `<method>__<name>`,
# with the name URL-encoded so that any contrast name is a valid key.
.ptm_key <- function(prefix, name) {
  paste0(prefix, "__", utils::URLencode(name, reserved = TRUE, repeated = TRUE))
}

# The results of one contrast on the modality's sites, in var order.
.ptm_result_frame <- function(data, site_var) {
  if (anyDuplicated(data$site)) {
    stop("PTM results have duplicate rows for one site.", call. = FALSE)
  }
  frame <- dplyr::left_join(data.frame(site = site_var$site), data, by = "site")
  frame$site[!frame$site %in% data$site] <- NA
  rownames(frame) <- rownames(site_var)
  frame
}

.ptm_method_table <- function(adata, method, contrasts) {
  var <- as.data.frame(adata$var, stringsAsFactors = FALSE)
  order <- adata$uns$prophosqua$varm_order
  tables <- lapply(contrasts, function(contrast) {
    key <- .ptm_key(method, contrast)
    block <- adata$varm[[key]]
    if (!is.data.frame(block) || nrow(block) != nrow(var)) {
      stop("PTM result table is missing or not on the site axis: ", key, call. = FALSE)
    }
    block <- block[!is.na(block$site), , drop = FALSE]
    rownames(block) <- NULL
    site_annotation <- var[unique(c("site", setdiff(names(var), names(block))))]
    table <- dplyr::left_join(block, site_annotation, by = "site")
    table$contrast <- contrast
    dplyr::as_tibble(table[, as.character(unlist(order[[key]], use.names = FALSE)), drop = FALSE])
  })
  dplyr::bind_rows(tables)
}
