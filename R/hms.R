as_hms <- if ("as_hms" %in% getNamespaceExports("hms")) {
  getExportedValue("hms", "as_hms")
} else {
  getExportedValue("hms", "as.hms")
}
