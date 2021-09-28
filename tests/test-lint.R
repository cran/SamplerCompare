bad_linters <- c("commented_code_linter", "object_name_linter",
                 "cyclocomp_linter")
if (Sys.getenv("NOT_CRAN") == "true") {
  linter_mask <- !(names(lintr::default_linters) %in% bad_linters)
  lintr::expect_lint_free(linters = lintr::default_linters[linter_mask])
}
