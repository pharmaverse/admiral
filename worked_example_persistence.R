# =============================================================================
# Worked example: what preserves / drops the `admiral_df` tag  (issue #3160 POC)
#
# Companion to worked_example.R. That script shows what summary.admiral_df()
# reports; this one shows WHEN the tag survives a pipeline at all, i.e. which
# operations keep or silently drop
#   * the `admiral_df` class, and
#   * the `admiral_keys` attribute
# and they are not always dropped together.
#
# CONTENTS
#   1. Mechanism            -- dplyr_reconstruct() copies attributes
#   2. Survives             -- the core verbs, joins, subsetting, bind_*
#   3. Drops both           -- summarise(), pivot_*, nest(), argument order
#   4. Partial drops        -- class gone but attribute kept (and vice versa)
#   5. Stale keys           -- rename()/select() leave keys pointing at nothing
#   6. Safe patterns        -- tag last, re-tag, check before relying on it
#                              (incl. check_tag() over each partial state)
#
# HOW TO RUN
#   From the package root:   Rscript worked_example_persistence.R
#
# NOTE: results below were verified with dplyr 1.2.1 / tidyr 1.3.2 / R 4.5.3.
# The behaviour is a property of dplyr's `dplyr_reconstruct()` method, not of
# admiral, so it can change with dplyr; re-run this script to re-check.
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(cli)
  library(rlang)
})

devtools::load_all()

rule <- function(title) cli::cli_h1(title)

# a small tagged BDS dataset, keys declared (not inferred) so we can watch the
# `admiral_keys` attribute travel through a pipeline
advs <- as_admiral_df(
  tribble(
    ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
    "PILOT01", "01-701-1015", "DIABP",     51, "BASELINE",
    "PILOT01", "01-701-1015", "SYSBP",    121, "BASELINE",
    "PILOT01", "01-701-1028", "DIABP",     79, "BASELINE",
    "PILOT01", "01-701-1028", "SYSBP",    130, "BASELINE"
  ),
  keys = c("USUBJID", "PARAMCD")
)

# an untagged record to bind on: a genuinely plain tibble, carrying neither the
# class nor the attribute (note that `as_tibble(advs)` would NOT do, it drops
# the class but keeps `admiral_keys` -- see section 4)
plain_row <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
  "PILOT01", "01-701-1015", "MAP",       74, "BASELINE"
)

# reports what is left of the tag after an operation
show <- function(label, data) {
  has_class <- inherits(data, "admiral_df")
  keys <- attr(data, "admiral_keys")
  fmt <- function(ok) if (ok) cli::col_green("kept") else cli::col_red("DROPPED")
  cli::cli_text(
    "{.strong {sprintf('%-34s', label)}} class: {fmt(has_class)}   ",
    "admiral_keys: {fmt(!is.null(keys))}"
  )
  invisible(data)
}

# =============================================================================
# 1. Mechanism
# =============================================================================
rule("1. Why anything survives at all")

cli::cli_bullets(c(
  "*" = "Most {.pkg dplyr} verbs rebuild their result with
         {.fun dplyr::dplyr_reconstruct}, whose {.cls data.frame} method copies
         {.emph all} attributes from the template -- including a custom class
         and {.var admiral_keys}.",
  "*" = "The tag is therefore lost only when an operation builds a
         {.emph new} data frame from scratch, or when the template it
         reconstructs from is a different (untagged) object.",
  "i" = "This is the opposite of the common assumption that attributes never
         survive {.pkg dplyr}."
))

# =============================================================================
# 2. Operations that KEEP both the class and the attribute
# =============================================================================
rule("2. Survives: the everyday derivation pipeline")

show("mutate()", mutate(advs, AVALC = as.character(AVAL)))
show("filter()", filter(advs, AVAL > 60))
show("arrange()", arrange(advs, desc(AVAL)))
show("select()", select(advs, USUBJID, PARAMCD, AVAL))
show("slice()", slice(advs, 1:2))
show("distinct()", distinct(advs, USUBJID, PARAMCD, AVAL))
show("relocate()", relocate(advs, AVAL))
show("count()", count(advs, USUBJID))
show("head()", head(advs, 2))
show("[ (rows)", advs[1:2, ])
show("[ (cols)", advs[, 1:3])
show("[[<-", {
  tmp <- advs
  tmp[["SRCDOM"]] <- "VS"
  tmp
})
show("bind_cols()", bind_cols(advs, tibble(FLAG = rep("Y", 4))))
show("bind_rows(tagged, plain)", bind_rows(advs, plain_row))
show("left_join(tagged, plain)", left_join(advs, tibble(
  USUBJID = "01-701-1015", TRT01P = "Placebo"
), by = "USUBJID"))
show("anti_join(tagged, plain)", anti_join(advs, tibble(
  USUBJID = "01-701-1028"
), by = "USUBJID"))
show("tidyr::fill()", tidyr::fill(advs, AVAL))
show("saveRDS() / readRDS()", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(advs, f)
  readRDS(f)
})

cli::cli_alert_success(
  "A derivation built from these verbs keeps its tag end to end."
)

# =============================================================================
# 3. Operations that SILENTLY drop both
# =============================================================================
rule("3. Drops both: new frame built from scratch")

show("summarise()", summarise(advs, N = n(), .by = USUBJID))
show("group_by() |> ungroup()", ungroup(group_by(advs, USUBJID)))
show("group_modify()", group_modify(group_by(advs, USUBJID), ~ .x))
show("tidyr::pivot_wider()", pivot_wider(
  advs,
  names_from = PARAMCD, values_from = AVAL
))
show("tidyr::pivot_longer()", pivot_longer(advs, AVAL))
show("nest() |> unnest()", unnest(nest(advs, data = AVAL), data))
show("tidyr::complete()", complete(advs, USUBJID, PARAMCD))
show("as_tibble(lapply(x, ...))", as_tibble(lapply(advs, identity)))

cli::cli_h2("Argument order and join side decide the outcome")

# `bind_rows()` reconstructs from its FIRST argument
show("bind_rows(tagged, plain)", bind_rows(advs, plain_row))
show("bind_rows(plain, tagged)", bind_rows(plain_row, advs))

# joins reconstruct from `x`, so a tagged `y` is not enough
show("left_join(tagged, plain)", left_join(
  advs, tibble(USUBJID = "01-701-1015", TRT01P = "Placebo"),
  by = "USUBJID"
))
show("left_join(plain, tagged)", left_join(
  tibble(USUBJID = c("01-701-1015", "01-701-1028")), advs,
  by = "USUBJID"
))

cli::cli_alert_warning(
  "Same data, different argument order -- keep the tagged dataset on the left."
)

# =============================================================================
# 4. Partial drops -- the two are NOT always lost together
# =============================================================================
rule("4. Partial drops: class and attribute part company")

cli::cli_h2("Class gone, attribute kept")
show("group_by()", group_by(advs, USUBJID))
show("rowwise()", rowwise(advs))
show("as_tibble()", as_tibble(advs))
show("as.data.frame()", as.data.frame(advs))
show("tidyr::drop_na()", drop_na(advs))

grouped <- group_by(advs, USUBJID)
cli::cli_bullets(c(
  "!" = "{.fun group_by} replaces the class vector outright:
         {.val {class(grouped)}}.",
  " " = "{.code inherits(x, 'admiral_df')} is already {.val FALSE}
         {.emph inside} the grouped block, so {.fun summary} dispatches to
         {.fun summary.data.frame} and the diagnostic is silently lost --
         {.fun ungroup} then removes {.var admiral_keys} as well.",
  " " = "{.fun as_tibble} / {.fun as.data.frame} leave the opposite: an
         orphaned {.var admiral_keys} with no class to dispatch on."
))

# the orphaned-attribute case, made concrete
plain <- as_tibble(advs)
cli::cli_text(
  "as_tibble(advs): class {.val {class(plain)}}, ",
  "but attr = {.val {attr(plain, 'admiral_keys')}}"
)

# =============================================================================
# 5. Stale keys -- preserved, and therefore able to go wrong
# =============================================================================
rule("5. Stale keys: rename() and select() keep the OLD names")

# the keys are stored as plain variable names, so nothing updates them
renamed <- rename(advs, SUBJID = USUBJID)
cli::cli_bullets(c(
  " " = "columns after rename(): {.val {colnames(renamed)}}",
  "!" = "'admiral_keys' still says: {.val {attr(renamed, 'admiral_keys')}}"
))

cli::cli_h2("summary() warns instead of quietly weakening the grain check")

# partially stale: only PARAMCD is left, so the check is weaker than declared
cli::cli_text("Partially stale (USUBJID renamed away):")
res_partial <- summary(renamed)
print(res_partial)

# fully stale: no declared key is left at all, no grain check is possible
cli::cli_text("Fully stale (both key variables gone):")
res_none <- summary(select(renamed, -PARAMCD))
print(res_none)

cli::cli_bullets(c(
  "i" = "Without the warning, the first case would silently check the grain on
         {.var PARAMCD} alone and report the verdict as a declared one -- here
         that invents duplicates, and with different data it would just as
         easily hide real ones.",
  "i" = "Inferred keys never trigger this warning: they are derived from the
         columns actually present."
))

# =============================================================================
# 6. Safe patterns
# =============================================================================
rule("6. Safe patterns")

cli::cli_bullets(c(
  "v" = "Tag {.emph last}: run {.fun as_admiral_df} / {.fun set_admiral_keys}
         after the reshaping is done, not before.",
  "v" = "Re-tag after {.fun summarise}, {.fun ungroup}, {.fun pivot_wider} and
         friends -- {.fun as_admiral_df} is idempotent and cheap.",
  "v" = "Keep the tagged dataset as the {.emph first} argument of
         {.fun bind_rows} and the {.arg x} of a join.",
  "v" = "Re-tag after {.fun rename} if a key variable was renamed; the
         attribute is preserved but not rewritten."
))

cli::cli_h2("Re-tagging restores everything")

reshaped <- advs %>%
  group_by(USUBJID, PARAMCD) %>%
  summarise(AVAL = mean(AVAL), .groups = "drop")
show("after group_by |> summarise", reshaped)
show("after re-tagging", as_admiral_df(
  reshaped,
  keys = c("USUBJID", "PARAMCD")
))

cli::cli_h2("Checking the tag before relying on it")

check_tag <- function(data) {
  has_class <- inherits(data, "admiral_df")
  keys <- attr(data, "admiral_keys")
  stale <- setdiff(keys, colnames(data))

  # `keys` has three distinct states and they behave differently in summary():
  #   NULL           -- no attribute at all, summary() infers the grain
  #   character(0)   -- an attribute IS set, so the inference fallback is
  #                     skipped and the grain is not checked at all, silently
  #   names          -- checked as declared, but only for the names still present
  # note `all(NULL %in% cols)` is TRUE (all(logical(0))), so a bare `all()` here
  # would report "keys present" for a dataset that declares no keys whatsoever
  declared <- if (is.null(keys)) {
    "none -- summary() infers the grain"
  } else if (length(keys) == 0) {
    cli::col_red("empty -- attribute set but no keys, so no grain check at all")
  } else {
    paste(keys, collapse = ", ")
  }

  status <- if (length(keys) == 0) {
    "n/a"
  } else if (length(stale) == 0) {
    cli::col_green("all present")
  } else if (length(stale) < length(keys)) {
    cli::col_red(paste0(
      "STALE: ", paste(stale, collapse = ", "),
      " missing; grain checked on ",
      paste(intersect(keys, colnames(data)), collapse = ", "), " only"
    ))
  } else {
    cli::col_red("STALE: no declared key left, grain not checked")
  }

  # cli_text()/cli_bullets() collapse runs of spaces, so padded labels come out
  # unaligned there; cli_verbatim() passes the line through untouched. It does
  # not interpret inline markup either, hence the explicit cli::col_*() above.
  field <- function(label, value) {
    cli::cli_verbatim(paste0("  ", sprintf("%-18s", label), value))
  }
  field("is an admiral_df:", if (has_class) {
    cli::col_green("TRUE")
  } else {
    cli::col_red("FALSE")
  })
  field("declared keys:", declared)
  field("keys present:", status)

  # the two halves are independent: an orphaned attribute is the dangerous case,
  # because re-tagging with `as_admiral_df(x)` keeps it (`keys = NULL` leaves the
  # attribute untouched) and promotes it back to "declared". This one is a
  # cli_bullets() message rather than a field, so that it wraps to the console.
  if (!has_class && !is.null(keys)) {
    cli::cli_bullets(c("!" = paste(
      "orphaned {.var admiral_keys}: {.fun summary} dispatches to",
      "{.fun summary.data.frame}, so the diagnostic is silently gone --",
      "re-tag with explicit {.arg keys} rather than {.code as_admiral_df(x)}."
    )))
  }
  invisible(data)
}

# the fully tagged reference case: both halves present, nothing stale
cli::cli_text("Fully tagged:")
check_tag(advs)

cli::cli_h2("Partial coverage: each half can be missing on its own")

# class kept, no attribute -- the benign case. summary() falls back to
# infer_admiral_keys(), so the grain is still checked, just marked "(inferred)".
# Beware: `all(NULL %in% colnames(data))` is TRUE, so a naive presence check
# reports "keys present" for a dataset that declares no keys at all.
cli::cli_text("Class only (no {.var admiral_keys}) -- summary() infers:")
check_tag(as_admiral_df(reshaped))

# class kept, attribute set but EMPTY. `is.null(keys)` is FALSE, so summary()
# skips the inference fallback; `length(keys) > 0` is FALSE, so it also skips
# the grain check and prints no Grain line -- no check, no warning, no trace.
# set_admiral_keys() produces exactly this when the spec defines no `key_seq`.
no_keys <- advs
attr(no_keys, "admiral_keys") <- character(0)
cli::cli_text("Class only (empty keys) -- summary() checks nothing:")
check_tag(no_keys)

# attribute kept, class gone -- the dangerous case. `plain` came from
# as_tibble(advs) in section 4; summary() now dispatches to summary.data.frame()
# and the admiral diagnostic disappears without an error.
cli::cli_text("Keys only (orphaned attribute):")
check_tag(plain)

cli::cli_h2("Partial coverage: the keys themselves can be half-present")

# `renamed` (section 5) renamed USUBJID away, so one declared key survives; the
# grain check silently narrows to PARAMCD unless something flags it
cli::cli_text("Partially stale keys:")
check_tag(renamed)

# and with PARAMCD dropped too, no declared key is left at all
cli::cli_text("Fully stale keys:")
check_tag(select(renamed, -PARAMCD))

cli::cli_bullets(c(
  "i" = "The class and {.var admiral_keys} travel independently, so
         {.code inherits(x, 'admiral_df')} alone is not evidence that the grain
         will be checked -- and {.var admiral_keys} alone is not evidence that
         {.fun summary} will dispatch to the admiral method.",
  "i" = "{.code length(keys) == 0} needs splitting into {.code NULL} (infer) and
         {.code character(0)} (check nothing): they print the same but behave
         differently in {.fun summary.admiral_df}.",
  "v" = "Re-tag with explicit {.arg keys} after any of these --
         {.code as_admiral_df(x)} on its own restores the class but keeps
         whatever stale attribute was riding along."
))

cli::cli_h1("Done")
