# Portable report values stored inside MuData, never serialized R objects.
# The explicit type and ordered names preserve empty tables, named ranks,
# factors, and matrix dimensions across R/Python HDF5 implementations.
.pack_ptm_value <- function(value) {
  if (is.null(value)) {
    return(list(type = "null"))
  }
  if (is.data.frame(value)) {
    return(list(type = "data.frame", columns = .pack_ptm_value(as.list(value)), row_names = rownames(value)))
  }
  if (is.factor(value)) {
    return(list(
      type = "factor",
      values = as.character(value),
      missing = is.na(value),
      levels = levels(value),
      ordered = is.ordered(value)
    ))
  }
  if (isS4(value)) {
    if (!inherits(value, "gseaResult")) {
      stop("Unsupported report object: ", class(value)[1L])
    }
    slots <- stats::setNames(
      lapply(methods::slotNames(value), function(key) methods::slot(value, key)),
      methods::slotNames(value)
    )
    return(list(type = "gseaResult", slots = .pack_ptm_value(slots)))
  }
  if (is.list(value)) {
    items <- lapply(value, .pack_ptm_value)
    names(items) <- sprintf("item_%06d", seq_along(items))
    return(list(type = "list", names = names(value), items = items))
  }
  if (!is.atomic(value)) {
    stop("Unsupported MuData value: ", class(value)[1L])
  }
  list(
    type = "atomic",
    storage = typeof(value),
    values = unname(as.vector(value)),
    missing = unname(as.vector(is.na(value))),
    names = names(value),
    dimensions = dim(value),
    dimnames = .pack_ptm_value(dimnames(value))
  )
}

.unpack_ptm_value <- function(value) {
  readers <- list(
    null = function(x) NULL,
    data.frame = function(x) {
      result <- .unpack_ptm_value(x$columns)
      structure(
        result,
        class = "data.frame",
        row.names = if (length(x$row_names)) as.character(unlist(x$row_names, use.names = FALSE)) else integer()
      )
    },
    factor = function(x) {
      values <- as.character(unlist(x$values, use.names = FALSE))
      values[as.logical(unlist(x$missing, use.names = FALSE))] <- NA_character_
      factor(values, levels = as.character(unlist(x$levels, use.names = FALSE)), ordered = x$ordered)
    },
    gseaResult = function(x) {
      do.call(methods::new, c(list(Class = "gseaResult"), .unpack_ptm_value(x$slots)))
    },
    list = function(x) {
      result <- lapply(x$items[sort(names(x$items))], .unpack_ptm_value)
      names(result) <- unlist(x$names, use.names = FALSE)
      result
    },
    atomic = function(x) {
      result <- as.vector(unlist(x$values, use.names = FALSE), mode = x$storage)
      result[as.logical(unlist(x$missing, use.names = FALSE))] <- NA
      dimensions <- unlist(x$dimensions, use.names = FALSE)
      if (length(dimensions)) {
        dim(result) <- as.integer(dimensions)
        dimnames(result) <- .unpack_ptm_value(x$dimnames)
      } else {
        value_names <- unlist(x$names, use.names = FALSE)
        if (length(value_names)) {
          names(result) <- as.character(value_names)
        }
      }
      result
    }
  )
  reader <- readers[[value$type]]
  if (is.null(reader)) {
    stop("Unknown portable PTM value type: ", value$type)
  }
  reader(value)
}

.require_ptm_fields <- function(value, fields, label) {
  missing <- fields[!fields %in% names(value)]
  if (length(missing)) {
    stop(label, " is incomplete: ", paste(missing, collapse = ", "))
  }
  if (any(vapply(value[fields], is.null, logical(1)))) {
    stop(label, " contains a missing required value.")
  }
  invisible(value)
}

.copy_ptm_value <- function(value) .unpack_ptm_value(.pack_ptm_value(value))
