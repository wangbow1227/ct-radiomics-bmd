required_pkgs <- c("dplyr", "tibble", "caret", "rpart", "randomForest", "glmnet", "pROC", "Boruta")
missing <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
invisible(lapply(required_pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

N_SPLITS <- 30
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(".")
}
script_dir  <- get_script_dir(); project_dir <- dirname(script_dir)
data_dir    <- file.path(project_dir, "data")
out_dir     <- file.path(project_dir, "output", "stability")
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
model_df <- internal %>%
  mutate(Type = as.factor(Type)) %>%
  select(Type, Type_new1, Type_new2, Age, Sex2, BMI, all_of(feat_cols))
cat(sprintf("Cohort n = %d\n\n", nrow(model_df)))

TASKS <- list(
  list(resp = "Type_new2", label = "Normal vs. Low BMD"),
  list(resp = "Type_new1", label = "Non-osteoporosis vs. Osteoporosis")
)

auc_rows  <- list()
feat_rows <- list()

for (s in seq_len(N_SPLITS)) {
  seed_s <- 1000 + s
  set.seed(seed_s)
  idx   <- createDataPartition(model_df$Type, p = 0.7, list = FALSE)
  tr_raw <- model_df[idx, ]
  te_raw <- model_df[-idx, ]

  mu <- sapply(model_df[, feat_cols], mean)
  sdv <- sapply(model_df[, feat_cols], sd); sdv[sdv == 0] <- 1
  scl <- function(df) {
    m <- sweep(as.matrix(df[, feat_cols]), 2, mu, "-")
    as.data.frame(sweep(m, 2, sdv, "/"))
  }
  tr <- bind_cols(tr_raw %>% select(-all_of(feat_cols)), scl(tr_raw))
  te <- bind_cols(te_raw %>% select(-all_of(feat_cols)), scl(te_raw))

  for (tk in TASKS) {
    y_tr <- factor(tr[[tk$resp]], levels = c(0, 1))
    y_te <- factor(te[[tk$resp]], levels = c(0, 1))
    X_tr <- tr[, feat_cols]

    set.seed(seed_s)
    bor <- Boruta(x = X_tr, y = y_tr, pValue = 0.01, maxRuns = 100, doTrace = 0)
    bor_feats <- getSelectedAttributes(TentativeRoughFix(bor), withTentative = FALSE)
    if (length(bor_feats) == 0) bor_feats <- feat_cols[1:5]

    set.seed(seed_s)
    cvfit <- cv.glmnet(as.matrix(X_tr), y_tr, family = "binomial", alpha = 1)
    co <- coef(cvfit, s = "lambda.min")
    las_feats <- setdiff(rownames(co)[as.numeric(co) != 0], "(Intercept)")
    if (length(las_feats) == 0) las_feats <- feat_cols[1:5]

    for (fs in list(list(n = "Boruta", f = bor_feats), list(n = "LASSO", f = las_feats))) {
      for (ft in fs$f) {
        feat_rows[[length(feat_rows) + 1]] <-
          tibble(Task = tk$label, Method = fs$n, Feature = ft, Split = s)
      }
    }

    fit_eval <- function(model_name, method, feats) {
      fml <- as.formula(paste("y ~ Age + Sex2 + BMI +", paste(feats, collapse = " + ")))
      d_tr <- cbind(tr[, c("Age", "Sex2", "BMI", feats)], y = y_tr)
      d_te <- cbind(te[, c("Age", "Sex2", "BMI", feats)], y = y_te)
      set.seed(seed_s)
      prob <- switch(method,
        decision_tree = predict(rpart(fml, data = d_tr, method = "class"), d_te, type = "prob")[, "1"],
        random_forest = predict(randomForest(fml, data = d_tr), d_te, type = "prob")[, "1"],
        logistic      = predict(glm(fml, data = d_tr, family = "binomial"), d_te, type = "response"))
      r <- pROC::roc(y_te, prob, levels = c("0", "1"), direction = "<", quiet = TRUE)
      tibble(Task = tk$label, Model = model_name, Split = s,
             AUC = as.numeric(pROC::auc(r)), n_features = length(feats))
    }

    auc_rows[[length(auc_rows) + 1]] <- fit_eval("Boruta + decision tree",       "decision_tree", bor_feats)
    auc_rows[[length(auc_rows) + 1]] <- fit_eval("Boruta + random forest",       "random_forest", bor_feats)
    auc_rows[[length(auc_rows) + 1]] <- fit_eval("Boruta + logistic regression", "logistic",      bor_feats)
    auc_rows[[length(auc_rows) + 1]] <- fit_eval("LASSO + logistic regression",  "logistic",      las_feats)
  }
  cat(sprintf("split %2d/%d done\n", s, N_SPLITS))
}

auc_df <- bind_rows(auc_rows)
write.csv(auc_df, file.path(out_dir, "stability_auc_by_split.csv"), row.names = FALSE)

summary_df <- auc_df %>%
  group_by(Task, Model) %>%
  summarise(n_splits = n(),
            Median_AUC = round(median(AUC), 3),
            Min_AUC    = round(min(AUC), 3),
            Max_AUC    = round(max(AUC), 3),
            IQR        = paste0(round(quantile(AUC, .25), 3), "-", round(quantile(AUC, .75), 3)),
            Median_n_features = round(median(n_features)),
            .groups = "drop") %>%
  arrange(Task, Model)
write.csv(summary_df, file.path(out_dir, "stability_auc_summary.csv"), row.names = FALSE)
cat("\n=== Supplementary Table S12 (stability of testing AUC) ===\n")
print(as.data.frame(summary_df), row.names = FALSE)

freq_df <- bind_rows(feat_rows) %>%
  count(Task, Method, Feature, name = "Times_selected") %>%
  mutate(Pct_of_resamples = round(100 * Times_selected / N_SPLITS)) %>%
  arrange(Task, Method, desc(Times_selected))
write.csv(freq_df, file.path(out_dir, "stability_feature_frequency.csv"), row.names = FALSE)

cat("\n=== Supplementary Table S13 (features selected in >=80% of resamples) ===\n")
print(as.data.frame(freq_df %>% filter(Pct_of_resamples >= 80)), row.names = FALSE)
cat("\nDone. See:", out_dir, "\n")
