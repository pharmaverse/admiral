deprecate_warn <- if (requireNamespace("lifecycle", quietly = TRUE)) {
  getExportedValue("lifecycle", "deprecate_warn")
} else {
  function(...) invisible()
}
