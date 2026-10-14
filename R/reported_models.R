required_pkgs <- c("dplyr", "tibble", "caret", "pROC")
missing <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
invisible(lapply(required_pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

SEED <- 12345
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(".")
}
script_dir  <- get_script_dir(); project_dir <- dirname(script_dir)
data_dir    <- file.path(project_dir, "data")
out_dir     <- file.path(project_dir, "output", "reported_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_cohort <- function(name) {
  path <- file.path(data_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Expected analysis-ready data file not found: ", path,
         ". See data/README.md for the required files and their structure.")
  }
  readRDS(path)
}

META_COLS <- c("ID", "Type", "Type_new1", "Type_new2", "Age", "Sex2", "BMI")
train_raw <- read_cohort("train")
test_raw  <- read_cohort("test")
feat_cols <- setdiff(colnames(train_raw), META_COLS)
internal  <- rbind(train_raw, test_raw)
internal  <- internal[order(internal$ID, method = "radix"), ]
full_scaled <- bind_cols(internal %>% select(-all_of(feat_cols)),
                         as.data.frame(scale(as.matrix(internal[, feat_cols]))))
train <- full_scaled[full_scaled$ID %in% train_raw$ID, ]
test  <- full_scaled[full_scaled$ID %in% test_raw$ID, ]

fs_dir  <- file.path(data_dir, "feature_sets")

MODELS <- list(
  list(set = "lowBMD_Boruta", resp = "Type_new2",
       label = "Normal vs. low BMD, Boruta + logistic regression"),
  list(set = "lowBMD_LASSO",  resp = "Type_new2",
       label = "Normal vs. low BMD, LASSO + logistic regression"),
  list(set = "osteo_Boruta",  resp = "Type_new1",
       label = "Osteoporosis, Boruta + logistic regression"),
  list(set = "osteo_LASSO",   resp = "Type_new1",
       label = "Osteoporosis, LASSO + logistic regression")
)

rows <- list()
for (m in MODELS) {
  fpath <- file.path(fs_dir, paste0(m$set, ".txt"))
  if (!file.exists(fpath)) {
    stop("Expected feature-list file not found: ", fpath,
         ". The six feature-list files are supplied together with the data files; see data/README.md.")
  }
  ft <- readLines(fpath)
  ft <- ft[nzchar(ft)]
  miss <- setdiff(ft, colnames(train))
  if (length(miss)) stop("features absent from the data: ", paste(miss, collapse = ", "))
  fml <- as.formula(paste(m$resp, "~ Age + Sex2 +", paste(ft, collapse = " + ")))
  fit <- suppressWarnings(glm(fml, data = train, family = binomial))
  auc_of <- function(d) as.numeric(suppressMessages(
    auc(roc(d[[m$resp]], predict(fit, d, type = "response")))))
  rows[[length(rows) + 1]] <- data.frame(
    Model = m$label, N_features = length(ft),
    Train_AUC = round(auc_of(train), 3), Test_AUC = round(auc_of(test), 3),
    Aliased_terms = sum(is.na(coef(fit))),
    Max_abs_coefficient = round(max(abs(coef(fit)), na.rm = TRUE), 1))
  cat(sprintf("%-50s n=%2d  train AUC %.3f  test AUC %.3f\n",
              m$label, length(ft), auc_of(train), auc_of(test)))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "reported_models_performance.csv"), row.names = FALSE)

ov <- function(a, b) length(intersect(
  readLines(file.path(fs_dir, paste0(a, ".txt"))),
  readLines(file.path(fs_dir, paste0(b, ".txt")))))
cat(sprintf("\nFeatures shared by the two osteoporosis selections: %d\n",
            ov("osteo_Boruta", "osteo_LASSO")))
cat(sprintf("Features shared by the two low BMD selections:      %d\n",
            ov("lowBMD_Boruta", "lowBMD_LASSO")))
cat("\nWritten to:", out_dir, "\n")
