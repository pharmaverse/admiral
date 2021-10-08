as_hms <- if ("as_hms" %in% getNamespaceExports("hms")) {
  getExportedValue("hms", "as_hms")
} else {
  fun <- getExportedValue("hms", "as.hms")
  function(x, tz = attr(x, "tzone"), ...) {
    fun(x, tz, ...)
  }
}
