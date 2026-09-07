# Run `summary.admiral_df()` over every dataset in {pharmaverseadam}
#
# A breadth test for the #3160 summary feature: 30 real ADaM datasets across
# several therapeutic areas, none of them written with this feature in mind.
# It exercises the parts that are hardest to cover with hand-built fixtures --
# `get_admiral_df_type()` on unusual structures, `infer_admiral_keys()` when no
# keys are declared, and graceful degradation when expected variables are absent.
#
# Every non-ADSL dataset is summarized twice: once alone, and once with
# `adsl =` so the coverage facts and the orphan check run. The ADSL is chosen
# by `USUBJID` overlap rather than by name, because the therapeutic-area
# datasets do not all share a subject universe -- and where the overlap is low,
# that is itself the finding.
#
# Nothing here asserts: it prints, and the roll-up table at the end is the
# thing to read. Run with:
#   Rscript worked_example_pharmaverseadam.R

library(admiral)
library(pharmaverseadam)
library(dplyr)
library(tibble)

# plain-text cli output so the captured log is readable without ANSI codes
options(cli.unicode = FALSE, cli.num_colors = 1, width = 88)

ds_names <- data(package = "pharmaverseadam")$results[, "Item"]
adsl_names <- grep("^adsl", ds_names, value = TRUE)

load_ds <- function(nm) {
  e <- new.env()
  data(list = nm, package = "pharmaverseadam", envir = e)
  e[[nm]]
}

# the ADSL a dataset actually belongs to, by subject overlap; `overlap` is the
# share of the dataset's subjects the candidate accounts for, so 1 means every
# subject is covered and 0 means the two share no subjects at all
best_adsl <- function(ds) {
  if (!"USUBJID" %in% names(ds)) {
    return(NULL)
  }
  subj <- unique(ds$USUBJID)
  overlap <- vapply(
    adsl_names,
    function(a) mean(subj %in% load_ds(a)$USUBJID),
    numeric(1)
  )
  list(name = names(which.max(overlap)), overlap = unname(max(overlap)))
}

rule <- function(...) {
  cat("\n", strrep("=", 88), "\n", sprintf(...), "\n", strrep("=", 88), "\n\n",
    sep = ""
  )
}

# one row per (dataset, with/without adsl) for the roll-up at the end
log <- list()
record <- function(name, adsl_arg, s) {
  checks <- c(
    if (isTRUE(s$bds$multiple_baselines > 0)) "multiple_baselines",
    if (isTRUE(s$occds$occ_flag_dups > 0)) "occ_flag_dups",
    if (isTRUE(s$vs_adsl$n_orphans > 0)) "n_orphans"
  )
  log[[length(log) + 1]] <<- tibble(
    dataset = name,
    adsl = adsl_arg,
    type = s$type %||% NA_character_,
    obs = s$n_obs %||% NA_integer_,
    subj = s$n_subjects %||% NA_integer_,
    # no `keys` element at all means the structure check could not run
    key_src = s$key_source %||% "not checked",
    dup_keys = s$n_duplicate_keys %||% NA_integer_,
    orphans = s$vs_adsl$n_orphans %||% NA_integer_,
    adsl_cov = if (is.null(s$vs_adsl)) {
      NA_character_
    } else {
      sprintf("%d/%d", s$vs_adsl$n_common, s$vs_adsl$n_adsl)
    },
    failed = if (length(checks) > 0) paste(checks, collapse = ", ") else "",
    keys = paste(s$keys, collapse = "+")
  )
}

for (nm in ds_names) {
  ds <- load_ds(nm)
  rule("%s  --  %d obs, %d cols", nm, nrow(ds), ncol(ds))

  adm <- as_admiral_df(ds)

  cat("--- summary(", nm, ") ---\n\n", sep = "")
  s <- try(print(summary(adm)), silent = TRUE)
  if (!inherits(s, "try-error")) {
    record(nm, "none", s)
  } else {
    cat("ERROR: ", conditionMessage(attr(s, "condition")), "\n", sep = "")
  }

  # ADSL is the subject-level source of truth for everything else, so it is
  # only the non-ADSL datasets that get the second pass
  if (get_admiral_df_type(ds) == "ADSL") {
    next
  }

  pair <- best_adsl(ds)
  if (is.null(pair)) {
    cat("\n(no USUBJID: the `adsl =` pass is skipped)\n")
    next
  }

  cat(
    "\n--- summary(", nm, ", adsl = ", pair$name, ")  [",
    sprintf("%.1f%%", 100 * pair$overlap), " of subjects covered] ---\n\n",
    sep = ""
  )
  s2 <- try(print(summary(adm, adsl = load_ds(pair$name))), silent = TRUE)
  if (!inherits(s2, "try-error")) {
    record(nm, pair$name, s2)
  } else {
    cat("ERROR: ", conditionMessage(attr(s2, "condition")), "\n", sep = "")
  }
}

rule("Roll-up")
options(width = 200)
print(bind_rows(log), n = Inf, width = Inf)
