required_pkgs <- c("dplyr", "tibble", "purrr", "caret", "rpart", "randomForest", "pROC", "logistf")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Install them first, e.g. install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
}
invisible(lapply(required_pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

SEED <- 12345
BORUTA_PVALUE  <- 0.01
BORUTA_MAXRUNS <- 100

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

script_dir   <- get_script_dir()
project_dir  <- dirname(script_dir)
data_dir     <- file.path(project_dir, "data")
out_dir      <- file.path(project_dir, "output", "main_analysis")
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
cat(sprintf("Internal cohort n = %d (train %d / test %d), radiomic features = %d\n",
            nrow(internal), nrow(train_raw), nrow(test_raw), length(feat_cols)))

feat_mean <- sapply(internal[, feat_cols], mean)
feat_sd   <- sapply(internal[, feat_cols], sd)
feat_sd[feat_sd == 0] <- 1

scale_features <- function(df) {
  scaled <- sweep(as.matrix(df[, feat_cols]), 2, feat_mean, "-")
  scaled <- sweep(scaled, 2, feat_sd, "/")
  as.data.frame(scaled)
}

train <- bind_cols(train_raw %>% select(-all_of(feat_cols)), scale_features(train_raw))
test  <- bind_cols(test_raw  %>% select(-all_of(feat_cols)), scale_features(test_raw))

table1_data <- internal %>%
  mutate(
    LowBMD       = factor(Type_new2, levels = c(0, 1), labels = c("Normal", "Low BMD")),
    Osteoporosis = factor(Type_new1, levels = c(0, 1), labels = c("No", "Yes")),
    HU           = original_firstorder_Mean
  )

summarise_by_group <- function(df, group_var) {
  df %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      N = n(),
      `Age, mean` = round(mean(Age), 2), `Age, SD` = round(sd(Age), 2),
      `BMI, mean` = round(mean(BMI), 2), `BMI, SD` = round(sd(BMI), 2),
      `Female, %` = round(mean(Sex2) * 100),
      `HU, mean` = round(mean(HU), 2), `HU, SD` = round(sd(HU), 2),
      .groups = "drop"
    )
}
write.csv(summarise_by_group(table1_data, "LowBMD"), file.path(out_dir, "Table1_NormalVsLowBMD.csv"), row.names = FALSE)
write.csv(summarise_by_group(table1_data, "Osteoporosis"), file.path(out_dir, "Table1_NonosteoVsOsteo.csv"), row.names = FALSE)

select_features <- function(response_var) {
  y <- factor(train[[response_var]], levels = c(0, 1))
  X <- train[, feat_cols]
  set.seed(SEED)
  b <- Boruta::Boruta(x = X, y = y, pValue = BORUTA_PVALUE,
                      maxRuns = BORUTA_MAXRUNS, doTrace = 0)
  bf <- Boruta::getSelectedAttributes(Boruta::TentativeRoughFix(b), withTentative = FALSE)
  set.seed(SEED)
  cv <- glmnet::cv.glmnet(as.matrix(X), y, family = "binomial", alpha = 1)
  co <- coef(cv, s = "lambda.min")
  lf <- setdiff(rownames(co)[as.numeric(co) != 0], "(Intercept)")
  list(boruta = bf, lasso = lf, lambda = cv$lambda.min)
}

FEATURE_SOURCE <- Sys.getenv("BMD_FEATURES", unset = "reported")
if (!FEATURE_SOURCE %in% c("reported", "select")) {
  stop("BMD_FEATURES must be 'reported' or 'select', got: ", FEATURE_SOURCE)
}
if (FEATURE_SOURCE == "reported") {
  fs_dir <- file.path(data_dir, "feature_sets")
  rdset <- function(n) {
    path <- file.path(fs_dir, paste0(n, ".txt"))
    if (!file.exists(path)) {
      stop("Expected feature-list file not found: ", path,
           ". The six feature-list files are supplied together with the data files; see data/README.md.")
    }
    v <- readLines(path)
    v <- v[nzchar(v)]
    miss <- setdiff(v, feat_cols)
    if (length(miss)) stop("features absent from the data: ", paste(miss, collapse = ", "))
    v
  }
  BMD_BORUTA_FEATS   <- rdset("lowBMD_Boruta")
  BMD_LASSO_FEATS    <- rdset("lowBMD_LASSO")
  OSTEO_BORUTA_FEATS <- rdset("osteo_Boruta")
  OSTEO_LASSO_FEATS  <- rdset("osteo_LASSO")
  LAMBDA_BMD <- NA_real_; LAMBDA_OSTEO <- NA_real_
  cat("Feature sets: data/feature_sets (Supplementary Table S3)\n")
} else {
  sel_bmd <- select_features("Type_new2")
  sel_op  <- select_features("Type_new1")
  BMD_BORUTA_FEATS   <- sel_bmd$boruta
  BMD_LASSO_FEATS    <- sel_bmd$lasso
  OSTEO_BORUTA_FEATS <- sel_op$boruta
  OSTEO_LASSO_FEATS  <- sel_op$lasso
  LAMBDA_BMD   <- sel_bmd$lambda
  LAMBDA_OSTEO <- sel_op$lambda
  cat(sprintf("Feature selection (Boruta pValue %.3g, maxRuns %d, TentativeRoughFix; LASSO lambda.min)\n",
              BORUTA_PVALUE, BORUTA_MAXRUNS))
}
cat(sprintf("  low BMD      Boruta %2d   LASSO %2d\n",
            length(BMD_BORUTA_FEATS), length(BMD_LASSO_FEATS)))
cat(sprintf("  osteoporosis Boruta %2d   LASSO %2d\n",
            length(OSTEO_BORUTA_FEATS), length(OSTEO_LASSO_FEATS)))

make_formula <- function(feats) {
  as.formula(paste(".y ~ Age + Sex2 + BMI +", paste(feats, collapse = " + ")))
}

hanley_mcneil <- function(auc, n_pos, n_neg) {
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt(max(0, (auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg)))
  c(max(0, auc - 1.96 * se), min(1, auc + 1.96 * se))
}

compute_metrics <- function(truth, prob, pred_class) {
  truth <- factor(truth, levels = c("0", "1"))
  pred_class <- factor(pred_class, levels = c("0", "1"))
  cm <- caret::confusionMatrix(pred_class, truth, positive = "1")
  roc_obj <- pROC::roc(truth, prob, levels = c("0", "1"), direction = "<", quiet = TRUE)
  tibble(
    Sensitivity = unname(cm$byClass["Sensitivity"]) * 100,
    Specificity = unname(cm$byClass["Specificity"]) * 100,
    Accuracy    = unname(cm$overall["Accuracy"]) * 100,
    AUC         = as.numeric(pROC::auc(roc_obj))
  )
}

fit_predict <- function(method, formula, train_df, test_df) {
  if (method == "decision_tree") {
    set.seed(SEED)
    fit <- rpart::rpart(formula, data = train_df, method = "class")
    prob_train  <- predict(fit, train_df, type = "prob")[, "1"]
    prob_test   <- predict(fit, test_df,  type = "prob")[, "1"]
    class_train <- predict(fit, train_df, type = "class")
    class_test  <- predict(fit, test_df,  type = "class")
  } else if (method == "random_forest") {
    set.seed(SEED)
    fit <- randomForest::randomForest(formula, data = train_df)
    prob_train  <- predict(fit, train_df, type = "prob")[, "1"]
    prob_test   <- predict(fit, test_df,  type = "prob")[, "1"]
    class_train <- predict(fit, train_df, type = "class")
    class_test  <- predict(fit, test_df,  type = "class")
  } else if (method == "logistic_regression") {
    fit <- glm(formula, data = train_df, family = "binomial")
    prob_train  <- predict(fit, train_df, type = "response")
    prob_test   <- predict(fit, test_df,  type = "response")
    class_train <- factor(ifelse(prob_train > 0.5, "1", "0"), levels = c("0", "1"))
    class_test  <- factor(ifelse(prob_test  > 0.5, "1", "0"), levels = c("0", "1"))
  } else stop("unknown method: ", method)
  list(fit = fit, prob_train = prob_train, prob_test = prob_test,
       class_train = class_train, class_test = class_test)
}

run_task <- function(response_var, task_label, boruta_feats, lasso_feats, boruta_note, lasso_note) {
  cat(sprintf("\n---- %s ----\n", task_label))
  cat(sprintf("Boruta feature set: %d features (%s)\n", length(boruta_feats), boruta_note))
  cat(sprintf("LASSO  feature set: %d features (%s)\n", length(lasso_feats), lasso_note))

  train_df <- train %>% mutate(.y = factor(.data[[response_var]], levels = c(0, 1))) %>%
    select(.y, Age, Sex2, BMI, all_of(feat_cols))
  test_df <- test %>% mutate(.y = factor(.data[[response_var]], levels = c(0, 1))) %>%
    select(.y, Age, Sex2, BMI, all_of(feat_cols))

  strategies <- list(
    list(name = "No preselection + decision tree", method = "decision_tree", feats = feat_cols),
    list(name = "Boruta + decision tree",          method = "decision_tree", feats = boruta_feats),
    list(name = "No preselection + random forest", method = "random_forest", feats = feat_cols),
    list(name = "Boruta + random forest",          method = "random_forest", feats = boruta_feats),
    list(name = "Boruta + logistic regression",    method = "logistic_regression", feats = boruta_feats),
    list(name = "LASSO + logistic regression",     method = "logistic_regression", feats = lasso_feats)
  )

  results <- map_dfr(strategies, function(s) {
    fm <- make_formula(s$feats)
    res <- fit_predict(s$method, fm, train_df, test_df)
    m_train <- compute_metrics(train_df$.y, res$prob_train, res$class_train)
    m_test  <- compute_metrics(test_df$.y,  res$prob_test,  res$class_test)
    bind_rows(
      m_train %>% mutate(Strategy = s$name, Cohort = "Training", n_features = length(s$feats)),
      m_test  %>% mutate(Strategy = s$name, Cohort = "Testing",  n_features = length(s$feats))
    ) %>% mutate(Task = task_label, .before = 1) %>%
      mutate(.prob_train = list(res$prob_train), .prob_test = list(res$prob_test),
             .fit = list(res$fit), .method = s$method, .feats = list(s$feats))
  })

  attr(results, "train_y") <- train_df$.y
  attr(results, "test_y")  <- test_df$.y
  results
}

write_feature_sets <- function() {
  sets <- list(
    list(Task = "Normal vs. Low BMD",                Method = "Boruta", Features = BMD_BORUTA_FEATS),
    list(Task = "Normal vs. Low BMD",                Method = "LASSO",  Features = BMD_LASSO_FEATS),
    list(Task = "Non-osteoporosis vs. Osteoporosis", Method = "Boruta", Features = OSTEO_BORUTA_FEATS),
    list(Task = "Non-osteoporosis vs. Osteoporosis", Method = "LASSO",  Features = OSTEO_LASSO_FEATS))
  out <- map_dfr(sets, function(s)
    tibble(Task = s$Task, Method = s$Method, n_features = length(s$Features), Feature = s$Features))
  write.csv(out, file.path(out_dir, "selected_feature_sets.csv"), row.names = FALSE)
  cat("\n=== Selected feature sets ===\n")
  print(as.data.frame(out %>% count(Task, Method, n_features)), row.names = FALSE)
}
write_feature_sets()

prune_collinear <- function(feats, thresh = 0.999) {
  M <- as.matrix(train[, feats])
  cc <- abs(cor(M))
  keep <- character(0); dropped <- character(0)
  for (f in feats) {
    if (length(keep) == 0 || max(cc[f, keep]) < thresh) keep <- c(keep, f)
    else dropped <- c(dropped, f)
  }
  list(keep = keep, dropped = dropped)
}

write_aor <- function(response_var, task_label, feats, method_label) {
  stopifnot(requireNamespace("logistf", quietly = TRUE))
  pr <- prune_collinear(feats)
  cat(sprintf("  %s / %s: %d features, %d collinear dropped, %d modelled\n",
              task_label, method_label, length(feats), length(pr$dropped), length(pr$keep)))
  train_df <- train %>% mutate(.y = factor(.data[[response_var]], levels = c(0, 1)))
  fml <- as.formula(paste(".y ~ Age + Sex2 + BMI +", paste(pr$keep, collapse = " + ")))
  fit <- logistf::logistf(fml, data = train_df,
                          control = logistf::logistf.control(maxit = 500))
  est <- tibble(
    Task = task_label, Method = method_label, Variable = names(fit$coefficients),
    aOR      = round(exp(fit$coefficients), 2),
    CI_lower = round(exp(fit$ci.lower), 2),
    CI_upper = round(exp(fit$ci.upper), 2),
    P_value  = signif(fit$prob, 3)
  ) %>% filter(Variable != "(Intercept)")
  if (length(pr$dropped) > 0) {
    est <- bind_rows(est, tibble(
      Task = task_label, Method = method_label, Variable = pr$dropped,
      aOR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, P_value = NA_real_))
  }
  est
}

boruta_note <- if (FEATURE_SOURCE == "reported") "data/feature_sets" else sprintf("pValue %.3g, maxRuns %d", BORUTA_PVALUE, BORUTA_MAXRUNS)
results_bmd <- run_task("Type_new2", "Normal vs. Low BMD",
                         BMD_BORUTA_FEATS, BMD_LASSO_FEATS, boruta_note,
                         if (is.na(LAMBDA_BMD)) "data/feature_sets" else sprintf("lambda.min %.8f", LAMBDA_BMD))
results_op  <- run_task("Type_new1", "Non-osteoporosis vs. Osteoporosis",
                         OSTEO_BORUTA_FEATS, OSTEO_LASSO_FEATS, boruta_note,
                         if (is.na(LAMBDA_OSTEO)) "data/feature_sets" else sprintf("lambda.min %.8f", LAMBDA_OSTEO))

aor_table <- bind_rows(
  write_aor("Type_new2", "Normal vs. Low BMD",                BMD_BORUTA_FEATS,   "Boruta"),
  write_aor("Type_new2", "Normal vs. Low BMD",                BMD_LASSO_FEATS,    "LASSO"),
  write_aor("Type_new1", "Non-osteoporosis vs. Osteoporosis", OSTEO_BORUTA_FEATS, "Boruta"),
  write_aor("Type_new1", "Non-osteoporosis vs. Osteoporosis", OSTEO_LASSO_FEATS,  "LASSO"))
write.csv(aor_table, file.path(out_dir, "adjusted_odds_ratios.csv"), row.names = FALSE)
cat("\n=== Adjusted odds ratios ===\n")
print(as.data.frame(aor_table %>% count(Task, Method, name = "n_terms")), row.names = FALSE)
cat("aOR ranges by task/method:\n")
print(as.data.frame(aor_table %>% group_by(Task, Method) %>%
  summarise(aOR_min = min(aOR, na.rm = TRUE), aOR_max = max(aOR, na.rm = TRUE),
            n_significant = sum(P_value < 0.05, na.rm = TRUE), .groups = "drop")), row.names = FALSE)

all_results <- bind_rows(results_bmd, results_op)

results_table <- all_results %>%
  select(Task, Strategy, Cohort, n_features, Sensitivity, Specificity, Accuracy, AUC) %>%
  mutate(across(c(Sensitivity, Specificity, Accuracy), ~ round(., 2)), AUC = round(AUC, 3))

cat("\n=== Table 3 / Table 4 performance ===\n")
print(as.data.frame(results_table), row.names = FALSE)
write.csv(results_table, file.path(out_dir, "Table3_4.csv"), row.names = FALSE)

plot_roc_panel <- function(task_results, y_truth, cohort_name, panel_label) {
  cols <- c("black", "blue", "red", "green4", "hotpink", "gold4")
  order_names <- c("No preselection + decision tree", "No preselection + random forest",
                    "Boruta + decision tree", "Boruta + random forest",
                    "Boruta + logistic regression", "LASSO + logistic regression")
  legend_labels <- c("Decision tree", "Random forest", "Boruta + decision tree",
                      "Boruta + random forest", "Boruta + logistic regression",
                      "LASSO + logistic regression")
  probs <- if (cohort_name == "Training") task_results$.prob_train else task_results$.prob_test
  names(probs) <- task_results$Strategy
  plot(NULL, xlim = c(1, 0), ylim = c(0, 1), xlab = "Specificity", ylab = "Sensitivity", main = panel_label)
  abline(0, 1, col = "grey")
  aucs <- c()
  for (i in seq_along(order_names)) {
    prob <- probs[[order_names[i]]]
    r <- pROC::roc(y_truth, prob, levels = c("0", "1"), direction = "<", quiet = TRUE)
    lines(r$specificities, r$sensitivities, col = cols[i], lwd = 2)
    aucs[i] <- round(as.numeric(pROC::auc(r)), 3)
  }
  legend("bottomright", legend = paste0(legend_labels, " (AUC ", sprintf("%.3f", aucs), ")"),
         col = cols, lwd = 2, cex = 0.6, box.lty = 0)
}

png(file.path(out_dir, "Figure2.png"), width = 9, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2))
plot_roc_panel(results_bmd, attr(results_bmd, "train_y"), "Training", "A) Normal vs Low BMD - Training")
plot_roc_panel(results_bmd, attr(results_bmd, "test_y"),  "Testing",  "B) Normal vs Low BMD - Testing")
plot_roc_panel(results_op,  attr(results_op, "train_y"),  "Training", "C) Osteoporosis - Training")
plot_roc_panel(results_op,  attr(results_op, "test_y"),   "Testing",  "D) Osteoporosis - Testing")
dev.off()
cat(sprintf("\nFigure 2 written to %s\n", file.path(out_dir, "Figure2.png")))

external_path <- file.path(data_dir, "external.rds")

if (!file.exists(external_path)) {
  cat("\nEXTERNAL VALIDATION COHORT: NOT FOUND at", external_path, "\n")
} else {

  cat("\n=== External validation cohort ===\n")
  L1L4_e <- readRDS(external_path)
  cat(sprintf("External cohort n = %d\n", nrow(L1L4_e)))

  scale_external <- function(df) as.data.frame(scale(as.matrix(df[, feat_cols])))
  external <- bind_cols(
    L1L4_e %>% select(ID, Age, Sex2, BMI, Type_new1, Type_new2),
    scale_external(L1L4_e)
  )

  score_external <- function(task_results, response_var, task_label) {
    fitted_strategies <- task_results %>% filter(Cohort == "Training")
    truth_all <- factor(external[[response_var]], levels = c(0, 1))
    map_dfr(seq_len(nrow(fitted_strategies)), function(i) {
      row <- fitted_strategies[i, ]
      fit <- row$.fit[[1]]
      if (row$.method == "logistic_regression") {
        prob <- predict(fit, newdata = external, type = "response")
        pred_class <- factor(ifelse(prob > 0.5, "1", "0"), levels = c("0", "1"))
      } else {
        prob <- predict(fit, newdata = external, type = "prob")[, "1"]
        pred_class <- predict(fit, newdata = external, type = "class")
      }
      compute_metrics(truth_all, prob, pred_class) %>%
        mutate(Task = task_label, Strategy = row$Strategy, Cohort = "External", .before = 1)
    })
  }

  ext_bmd <- score_external(results_bmd, "Type_new2", "Normal vs. Low BMD")
  ext_op  <- score_external(results_op,  "Type_new1", "Non-osteoporosis vs. Osteoporosis")

  ext_results <- bind_rows(ext_bmd, ext_op) %>%
    mutate(across(c(Sensitivity, Specificity, Accuracy), ~ round(., 2)), AUC = round(AUC, 3))

  cat("=== External-validation performance ===\n")
  print(as.data.frame(ext_results), row.names = FALSE)
  write.csv(ext_results, file.path(out_dir, "Table3_4_external.csv"), row.names = FALSE)

  auc_ci_external <- function(task_results, response_var, task_label) {
    fitted_strategies <- task_results %>% filter(Cohort == "Training")
    truth_all <- factor(external[[response_var]], levels = c(0, 1))
    map_dfr(seq_len(nrow(fitted_strategies)), function(i) {
      row <- fitted_strategies[i, ]
      fit <- row$.fit[[1]]
      prob <- if (row$.method == "logistic_regression") {
        predict(fit, newdata = external, type = "response")
      } else {
        predict(fit, newdata = external, type = "prob")[, "1"]
      }
      roc_obj <- pROC::roc(truth_all, prob, levels = c("0", "1"), direction = "<", quiet = TRUE)
      auc_v <- as.numeric(pROC::auc(roc_obj))
      ci <- hanley_mcneil(auc_v, sum(truth_all == "1"), sum(truth_all == "0"))
      tibble(Task = task_label, Strategy = row$Strategy,
             AUC      = round(auc_v, 3),
             CI_lower = round(ci[1], 3),
             CI_upper = round(ci[2], 3),
             CI_width = round(ci[2] - ci[1], 3))
    })
  }
  auc_ci_table <- bind_rows(
    auc_ci_external(results_bmd, "Type_new2", "Normal vs. Low BMD"),
    auc_ci_external(results_op,  "Type_new1", "Non-osteoporosis vs. Osteoporosis"))
  cat("\n=== Hanley-McNeil 95% CIs for the 12 external-validation AUCs ===\n")
  print(as.data.frame(auc_ci_table), row.names = FALSE)
  write.csv(auc_ci_table, file.path(out_dir, "external_auc_ci.csv"), row.names = FALSE)

  summarise_range <- function(df, task, cohort) {
    s <- df %>% filter(Task == task)
    sprintf("%s | %s | AUC %.3f-%.3f | Sens %.2f-%.2f | Spec %.2f-%.2f | Acc %.2f-%.2f",
            task, cohort, min(s$AUC), max(s$AUC), min(s$Sensitivity), max(s$Sensitivity),
            min(s$Specificity), max(s$Specificity), min(s$Accuracy), max(s$Accuracy))
  }
  cat("\n=== External-cohort summary ranges ===\n")
  cat(summarise_range(ext_results, "Normal vs. Low BMD", "External"), "\n")
  cat(summarise_range(ext_results, "Non-osteoporosis vs. Osteoporosis", "External"), "\n")

  cat(sprintf("\nExternal AUC range: %.3f-%.3f (Normal vs. Low BMD), %.3f-%.3f (Osteoporosis)\n",
              min(ext_bmd$AUC), max(ext_bmd$AUC), min(ext_op$AUC), max(ext_op$AUC)))
}

cat("\nDone. See:", out_dir, "\n")

get_cm <- function(task_results, y_truth, cohort_name, strategy_name) {
  row <- task_results[task_results$Strategy == strategy_name, ]
  fit  <- row$.fit[[1]]
  method <- row$.method[1]
  pred_class <- if (method == "decision_tree") {
    if (cohort_name == "Training") predict(fit, train %>% mutate(.y = factor(y_truth, levels = c(0, 1))), type = "class")
    else predict(fit, test %>% mutate(.y = factor(y_truth, levels = c(0, 1))), type = "class")
  } else {
    if (cohort_name == "Training") predict(fit, train, type = "class")
    else predict(fit, test, type = "class")
  }
  truth <- factor(y_truth, levels = c("0", "1"))
  pred_class <- factor(pred_class, levels = c("0", "1"))
  cm <- caret::confusionMatrix(pred_class, truth, positive = "1")
  tb <- cm$table
  tibble(
    Cohort = cohort_name, Strategy = strategy_name,
    TN = tb["0", "0"], FP = tb["1", "0"], FN = tb["0", "1"], TP = tb["1", "1"]
  )
}

cm_rows <- list()
for (strat in c("Boruta + decision tree", "Boruta + random forest")) {
  cm_rows[[length(cm_rows) + 1]] <- get_cm(results_bmd, attr(results_bmd, "train_y"), "Training", strat) %>% mutate(Task = "Normal vs. Low BMD")
  cm_rows[[length(cm_rows) + 1]] <- get_cm(results_bmd, attr(results_bmd, "test_y"),  "Testing",  strat) %>% mutate(Task = "Normal vs. Low BMD")
  cm_rows[[length(cm_rows) + 1]] <- get_cm(results_op,  attr(results_op, "train_y"),  "Training", strat) %>% mutate(Task = "Non-osteoporosis vs. Osteoporosis")
  cm_rows[[length(cm_rows) + 1]] <- get_cm(results_op,  attr(results_op, "test_y"),   "Testing",  strat) %>% mutate(Task = "Non-osteoporosis vs. Osteoporosis")
}

cm_table <- bind_rows(cm_rows) %>% select(Task, Strategy, Cohort, TN, FP, FN, TP)
print(as.data.frame(cm_table), row.names = FALSE)
write.csv(cm_table, file.path(out_dir, "confusion_matrices_train_test.csv"), row.names = FALSE)

if (file.exists(external_path)) {
  get_cm_ext <- function(task_results, response_var, strategy_name, task_label) {
    row <- task_results[task_results$Strategy == strategy_name & task_results$Cohort == "Training", ]
    fit <- row$.fit[[1]]
    method <- row$.method[1]
    truth <- factor(external[[response_var]], levels = c(0, 1))
    if (method == "logistic_regression") {
      prob <- predict(fit, newdata = external, type = "response")
      pred_class <- factor(ifelse(prob > 0.5, "1", "0"), levels = c("0", "1"))
    } else {
      pred_class <- predict(fit, newdata = external, type = "class")
    }
    pred_class <- factor(pred_class, levels = c("0", "1"))
    cm <- caret::confusionMatrix(pred_class, truth, positive = "1")
    tb <- cm$table
    tibble(Task = task_label, Strategy = strategy_name, Cohort = "External",
           TN = tb["0", "0"], FP = tb["1", "0"], FN = tb["0", "1"], TP = tb["1", "1"])
  }
  ext_cm_rows <- list()
  for (strat in c("Boruta + decision tree", "Boruta + random forest")) {
    ext_cm_rows[[length(ext_cm_rows) + 1]] <- get_cm_ext(results_bmd, "Type_new2", strat, "Normal vs. Low BMD")
    ext_cm_rows[[length(ext_cm_rows) + 1]] <- get_cm_ext(results_op,  "Type_new1", strat, "Non-osteoporosis vs. Osteoporosis")
  }
  ext_cm_table <- bind_rows(ext_cm_rows) %>% select(Task, Strategy, Cohort, TN, FP, FN, TP)
  print(as.data.frame(ext_cm_table), row.names = FALSE)
  write.csv(ext_cm_table, file.path(out_dir, "confusion_matrices_external.csv"), row.names = FALSE)
}

if (file.exists(external_path)) {
  row_rf_op <- results_op[results_op$Strategy == "Boruta + random forest" & results_op$Cohort == "Training", ]
  fit_rf_op <- row_rf_op$.fit[[1]]
  pred_class_ext <- predict(fit_rf_op, newdata = external, type = "class")
  prob_ext <- predict(fit_rf_op, newdata = external, type = "prob")[, "1"]
  truth_ext <- factor(external$Type_new1, levels = c(0, 1))
  misclass_idx <- which(as.character(pred_class_ext) != as.character(truth_ext))
  misclass_df <- tibble(
    ID = external$ID[misclass_idx],
    True_label = ifelse(truth_ext[misclass_idx] == "1", "Osteoporosis", "Non-osteoporosis"),
    Predicted_label = ifelse(pred_class_ext[misclass_idx] == "1", "Osteoporosis", "Non-osteoporosis"),
    Predicted_prob_osteoporosis = round(prob_ext[misclass_idx], 3),
    Age = external$Age[misclass_idx], BMI = round(external$BMI[misclass_idx], 1)
  )
  print(as.data.frame(misclass_df), row.names = FALSE)
  write.csv(misclass_df, file.path(out_dir, "misclassified_external_cases.csv"), row.names = FALSE)
}

fit_and_eval_baseline <- function(formula, response_var, train_df, test_df, external_df) {
  fit <- glm(formula, data = train_df, family = "binomial")
  eval_one <- function(df, cohort) {
    prob <- predict(fit, newdata = df, type = "response")
    pred_class <- factor(ifelse(prob > 0.5, "1", "0"), levels = c("0", "1"))
    truth <- factor(df[[response_var]], levels = c(0, 1))
    cm <- caret::confusionMatrix(pred_class, truth, positive = "1")
    roc_obj <- pROC::roc(truth, prob, levels = c("0", "1"), direction = "<", quiet = TRUE)
    auc_v <- as.numeric(pROC::auc(roc_obj))
    ci <- hanley_mcneil(auc_v, sum(truth == "1"), sum(truth == "0"))
    tibble(
      Sensitivity = unname(cm$byClass["Sensitivity"]) * 100,
      Specificity = unname(cm$byClass["Specificity"]) * 100,
      Accuracy = unname(cm$overall["Accuracy"]) * 100,
      AUC = auc_v, AUC_lo = ci[1], AUC_hi = ci[2],
      Cohort = cohort
    )
  }
  bind_rows(eval_one(train_df, "Training"), eval_one(test_df, "Testing"), eval_one(external_df, "External"))
}

run_baseline_task <- function(response_var, task_label) {
  train_df <- train_raw %>% mutate(.y = .data[[response_var]], HU = original_firstorder_Mean)
  test_df  <- test_raw  %>% mutate(.y = .data[[response_var]], HU = original_firstorder_Mean)
  ext_df   <- L1L4_e    %>% mutate(.y = .data[[response_var]], HU = original_firstorder_Mean)

  baseline_models <- list(
    "Clinical only (Age+Sex+BMI)" = as.formula(".y ~ Age + Sex2 + BMI"),
    "HU only" = as.formula(".y ~ HU"),
    "Clinical + HU" = as.formula(".y ~ Age + Sex2 + BMI + HU")
  )

  map_dfr(names(baseline_models), function(model_name) {
    res <- fit_and_eval_baseline(baseline_models[[model_name]], ".y", train_df, test_df, ext_df)
    res %>% mutate(Task = task_label, Model = model_name, .before = 1)
  })
}

baseline_results <- bind_rows(
  run_baseline_task("Type_new2", "Normal vs. Low BMD"),
  run_baseline_task("Type_new1", "Non-osteoporosis vs. Osteoporosis")
) %>%
  mutate(across(c(Sensitivity, Specificity, Accuracy), ~ round(., 3)),
         across(c(AUC, AUC_lo, AUC_hi), ~ round(., 3))) %>%
  select(Task, Model, Sensitivity, Specificity, Accuracy, AUC, AUC_lo, AUC_hi, Cohort)

print(as.data.frame(baseline_results), row.names = FALSE)
write.csv(baseline_results, file.path(out_dir, "baseline_model_comparison.csv"), row.names = FALSE)

dca_pkgs <- c("dcurves", "ggplot2", "gridExtra")
missing_dca_pkgs <- dca_pkgs[!dca_pkgs %in% rownames(installed.packages())]
if (length(missing_dca_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_dca_pkgs, collapse = ", "),
       ". Install them first, e.g. install.packages(c(",
       paste(sprintf('"%s"', missing_dca_pkgs), collapse = ", "), "))")
}
invisible(lapply(dca_pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

strategy_names <- c("Boruta + decision tree", "Boruta + random forest",
                     "Boruta + logistic regression", "LASSO + logistic regression")

get_test_probs <- function(task_results) {
  out <- map(strategy_names, function(nm) {
    row <- task_results[task_results$Strategy == nm & task_results$Cohort == "Testing", ]
    row$.prob_test[[1]]
  })
  names(out) <- strategy_names
  out
}
probs_bmd_test <- get_test_probs(results_bmd)
probs_op_test  <- get_test_probs(results_op)

get_ext_probs <- function(task_results) {
  fitted <- task_results %>% filter(Cohort == "Training")
  out <- map(strategy_names, function(nm) {
    row <- fitted[fitted$Strategy == nm, ]
    fit <- row$.fit[[1]]; method <- row$.method[1]
    if (method == "logistic_regression") predict(fit, newdata = external, type = "response")
    else predict(fit, newdata = external, type = "prob")[, "1"]
  })
  names(out) <- strategy_names
  out
}
probs_bmd_ext <- get_ext_probs(results_bmd)
probs_op_ext  <- get_ext_probs(results_op)

truth_bmd_test <- attr(results_bmd, "test_y")
truth_op_test  <- attr(results_op,  "test_y")
truth_bmd_ext  <- factor(external$Type_new2, levels = c(0, 1))
truth_op_ext   <- factor(external$Type_new1, levels = c(0, 1))

calib_one <- function(truth, prob, strategy, task, cohort) {
  y <- as.numeric(as.character(truth))
  brier <- mean((prob - y)^2)
  eps <- 1e-6
  p_clamped <- pmin(pmax(prob, eps), 1 - eps)
  logit_p <- log(p_clamped / (1 - p_clamped))
  m_slope <- tryCatch(glm(y ~ logit_p, family = binomial()), error = function(e) NULL)
  m_int   <- tryCatch(glm(y ~ offset(logit_p), family = binomial()), error = function(e) NULL)
  tibble(Task = task, Strategy = strategy, Cohort = cohort,
         Brier = round(brier, 3),
         Calib_intercept = if (!is.null(m_int) && abs(coef(m_int)[1]) < 50)
                             round(unname(coef(m_int)[1]), 3) else NA_real_,
         Calib_slope = if (!is.null(m_slope)) round(unname(coef(m_slope)[2]), 3) else NA_real_)
}

calib_rows <- list()
for (nm in strategy_names) {
  calib_rows[[length(calib_rows) + 1]] <- calib_one(truth_bmd_test, probs_bmd_test[[nm]], nm, "Normal vs. Low BMD", "Testing")
  calib_rows[[length(calib_rows) + 1]] <- calib_one(truth_bmd_ext,  probs_bmd_ext[[nm]],  nm, "Normal vs. Low BMD", "External")
  calib_rows[[length(calib_rows) + 1]] <- calib_one(truth_op_test,  probs_op_test[[nm]],  nm, "Non-osteoporosis vs. Osteoporosis", "Testing")
  calib_rows[[length(calib_rows) + 1]] <- calib_one(truth_op_ext,   probs_op_ext[[nm]],   nm, "Non-osteoporosis vs. Osteoporosis", "External")
}
calib_table <- bind_rows(calib_rows)
write.csv(calib_table, file.path(out_dir, "calibration_metrics.csv"), row.names = FALSE)
cat("\n=== Calibration (Brier / intercept / slope) ===\n")
print(as.data.frame(calib_table), row.names = FALSE)

plot_calib <- function(truth, prob_list, main_title, filename) {
  png(filename, width = 6, height = 6, units = "in", res = 600)
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "Predicted probability", ylab = "Observed proportion", main = main_title)
  abline(0, 1, col = "grey", lty = 2)
  cols <- c("red", "green4", "hotpink", "gold4")
  y <- as.numeric(as.character(truth))
  for (i in seq_along(prob_list)) {
    prob <- prob_list[[i]]
    bins <- cut(prob, breaks = unique(quantile(prob, probs = seq(0, 1, 0.2), na.rm = TRUE)), include.lowest = TRUE)
    obs <- tapply(y, bins, mean); pred <- tapply(prob, bins, mean)
    lines(pred, obs, col = cols[i], type = "b", lwd = 2, pch = 16)
  }
  legend("bottomright", legend = names(prob_list), col = cols, lwd = 2, pch = 16, cex = 0.7, box.lty = 0)
  dev.off()
}
plot_calib(truth_bmd_test, probs_bmd_test, "A) Normal vs Low BMD - Testing", file.path(out_dir, "calibration_plot_lowbmd_testing.png"))
plot_calib(truth_bmd_ext,  probs_bmd_ext,  "B) Normal vs Low BMD - External", file.path(out_dir, "calibration_plot_lowbmd_external.png"))
plot_calib(truth_op_test,  probs_op_test,  "C) Osteoporosis - Testing", file.path(out_dir, "calibration_plot_osteoporosis_testing.png"))
plot_calib(truth_op_ext,   probs_op_ext,   "D) Osteoporosis - External", file.path(out_dir, "calibration_plot_osteoporosis_external.png"))

dca_one <- function(truth, prob_list, thresholds = seq(0, 0.5, 0.01)) {
  y <- as.numeric(as.character(truth))
  df <- data.frame(y = y)
  for (nm in names(prob_list)) df[[make.names(nm)]] <- prob_list[[nm]]
  preds <- make.names(names(prob_list))
  fm <- as.formula(paste("y ~", paste(preds, collapse = " + ")))
  label_list <- setNames(as.list(names(prob_list)), preds)
  dca(fm, data = df, thresholds = thresholds, label = label_list)
}

dca_bmd_test <- dca_one(truth_bmd_test, probs_bmd_test)
dca_bmd_ext  <- dca_one(truth_bmd_ext,  probs_bmd_ext)
dca_op_test  <- dca_one(truth_op_test,  probs_op_test)
dca_op_ext   <- dca_one(truth_op_ext,   probs_op_ext)

save_dca_csv <- function(dca_obj, filename) {
  write.csv(as_tibble(dca_obj$dca), filename, row.names = FALSE)
}
save_dca_csv(dca_bmd_test, file.path(out_dir, "dca_net_benefit_lowbmd_testing.csv"))
save_dca_csv(dca_bmd_ext,  file.path(out_dir, "dca_net_benefit_lowbmd_external.csv"))
save_dca_csv(dca_op_test,  file.path(out_dir, "dca_net_benefit_osteoporosis_testing.csv"))
save_dca_csv(dca_op_ext,   file.path(out_dir, "dca_net_benefit_osteoporosis_external.csv"))

p1 <- plot(dca_bmd_test) + ggtitle("A) Normal vs Low BMD - Testing")
p2 <- plot(dca_bmd_ext)  + ggtitle("B) Normal vs Low BMD - External")
p3 <- plot(dca_op_test)  + ggtitle("C) Osteoporosis - Testing")
p4 <- plot(dca_op_ext)   + ggtitle("D) Osteoporosis - External")

png(file.path(out_dir, "dca_plots.png"), width = 11, height = 10, units = "in", res = 600)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

cat("\nDCA and calibration outputs written to", out_dir, "\n")
