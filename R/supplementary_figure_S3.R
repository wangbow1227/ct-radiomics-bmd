required_pkgs <- c("dplyr", "tibble", "caret", "randomForest")
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
out_dir     <- file.path(project_dir, "output", "supplementary_figure_S3")
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

feature_set_dir <- file.path(project_dir, "data", "feature_sets")
read_feature_set <- function(name) {
  path <- file.path(feature_set_dir, paste0(name, ".txt"))
  if (!file.exists(path)) {
    stop("Expected feature-list file not found: ", path,
         ". The six feature-list files are supplied together with the data files; see data/README.md.")
  }
  v <- trimws(readLines(path))
  v[nzchar(v)]
}

PANELS <- list(
  list(panel = "A", resp = "Type_new2", task = "Normal vs. low BMD",                selection = "None",   feats = feat_cols),
  list(panel = "B", resp = "Type_new2", task = "Normal vs. low BMD",                selection = "Boruta", feats = read_feature_set("lowBMD_Boruta")),
  list(panel = "C", resp = "Type_new1", task = "Non-osteoporosis vs. osteoporosis", selection = "None",   feats = feat_cols),
  list(panel = "D", resp = "Type_new1", task = "Non-osteoporosis vs. osteoporosis", selection = "Boruta", feats = read_feature_set("osteo_Boruta"))
)

importance_rows <- list()
for (p in PANELS) {
  d <- train %>%
    mutate(.y = factor(.data[[p$resp]], levels = c(0, 1))) %>%
    select(.y, Age, Sex2, BMI, all_of(p$feats))
  set.seed(SEED)
  rf <- randomForest(.y ~ ., data = d, importance = TRUE)

  png(file.path(out_dir, paste0("Supplementary_Figure_S3", p$panel, ".png")),
      width = 7, height = 8, units = "in", res = 600)
  varImpPlot(rf, type = 2, main = "", n.var = 30)
  dev.off()

  gini <- sort(importance(rf, type = 2)[, 1], decreasing = TRUE)
  importance_rows[[p$panel]] <- tibble(Panel = p$panel, Task = p$task, Preselection = p$selection,
                                       Rank = seq_along(gini), Feature = names(gini),
                                       MeanDecreaseGini = unname(gini))
  cat(sprintf("S3%s top 5: %s\n", p$panel, paste(names(gini)[1:5], collapse = ", ")))
}

importance_all <- bind_rows(importance_rows)
write.csv(importance_all, file.path(out_dir, "Supplementary_Figure_S3_importance.csv"), row.names = FALSE)

top10 <- importance_all %>% filter(Rank <= 10) %>% split(.$Panel) %>% lapply(function(t) t$Feature)
cat(sprintf("\nTop-10 overlap, A vs B: %d of 10\n", length(intersect(top10$A, top10$B))))
cat(sprintf("Top-10 overlap, C vs D: %d of 10\n", length(intersect(top10$C, top10$D))))
cat(sprintf("\nOutputs written to %s\n", out_dir))
