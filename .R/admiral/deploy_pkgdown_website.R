devtools::document()
devtools::install(upgrade = "never")
if (dir.exists(".git/worktrees")) {
  unlink(".git/worktrees", recursive = TRUE)
}
pkgdown::deploy_to_branch()
