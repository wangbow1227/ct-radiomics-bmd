required <- c(
  "dplyr", "tibble", "purrr",
  "caret",
  "rpart",
  "randomForest",
  "glmnet",
  "Boruta",
  "pROC",
  "logistf",
  "dcurves",
  "ggplot2", "gridExtra"
)

install_missing <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) == 0) {
    message("All packages already installed: ", paste(pkgs, collapse = ", "))
  } else {
    message("Installing: ", paste(missing, collapse = ", "))
    install.packages(missing)
  }
}

install_missing(required)
