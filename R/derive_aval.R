derive_aval <- function(dataset) {
  select_col <- function(pattern) grep(pattern, colnames(dataset), value = TRUE)

  stres <- select_col("^[A-Z]{2}STRES$")
  stresn <- select_col("^[A-Z]{2}STRESN$")
  stresu <- select_col("^[A-Z]{2}STRESU$")

  mutate(
    dataset,
    AVALC = !!sym(stres),
    AVAL = !!sym(stresn),
    AVALU = !!sym(stresu)
  )
}
