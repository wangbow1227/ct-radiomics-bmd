required_pkgs <- c("dplyr", "tibble", "caret", "glmnet", "Boruta")
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
out_dir    <- file.path(project_dir, "output", "selection_figures")
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

TASKS <- list(
  list(resp = "Type_new2", key = "lowBMD",
       label = "Normal vs. low BMD",        lambda = 0.02008833,
       bor_p = 0.01, bor_runs = 100),
  list(resp = "Type_new1", key = "osteoporosis",
       label = "Osteoporosis",              lambda = 0.01622183,
       bor_p = 0.01, bor_runs = 100)
)

for (tk in TASKS) {
  y  <- factor(train[[tk$resp]], levels = c(0, 1))
  Xm <- as.matrix(train[, feat_cols])

  fit <- glmnet(Xm, y, family = "binomial", alpha = 1, nlambda = 300)
  n_at_lambda <- sum(as.numeric(coef(fit, s = tk$lambda))[-1] != 0)
  set.seed(SEED)
  cvfit <- cv.glmnet(Xm, y, family = "binomial", alpha = 1)

  png(file.path(out_dir, sprintf("LASSO_coefficient_path_%s.png", tk$key)),
      width = 7, height = 5.5, units = "in", res = 600)
  plot(fit, xvar = "lambda", label = FALSE)
  abline(v = log(tk$lambda), col = "red", lty = 2, lwd = 2)
  title(sprintf("LASSO coefficient paths - %s", tk$label), line = 2.5, cex.main = 1.0)
  mtext(sprintf("red dashed line: lambda = %.8f (%d features retained)", tk$lambda, n_at_lambda),
        side = 3, line = 0.4, cex = 0.8)
  dev.off()

  png(file.path(out_dir, sprintf("LASSO_cv_curve_%s.png", tk$key)),
      width = 7, height = 5.5, units = "in", res = 600)
  plot(cvfit)
  abline(v = log(tk$lambda), col = "red", lty = 2, lwd = 2)
  title(sprintf("LASSO cross-validation - %s", tk$label), line = 2.5, cex.main = 1.0)
  dev.off()

  lam <- exp(seq(log(max(fit$lambda)), log(min(fit$lambda)), length.out = 2000))
  nz  <- apply(as.matrix(coef(fit, s = lam)), 2, function(v) sum(v[-1] != 0))
  png(file.path(out_dir, sprintf("LASSO_features_vs_lambda_%s.png", tk$key)),
      width = 7, height = 5.5, units = "in", res = 600)
  plot(log(lam), nz, type = "l", lwd = 2,
       xlab = expression(log(lambda)), ylab = "Number of retained features",
       main = sprintf("Retained features vs lambda - %s", tk$label), cex.main = 1.0)
  abline(v = log(tk$lambda), col = "red", lty = 2, lwd = 2)
  abline(h = n_at_lambda, col = "grey50", lty = 3)
  legend("topright", bty = "n", cex = 0.8,
         legend = sprintf("lambda = %.8f -> %d features", tk$lambda, n_at_lambda))
  dev.off()

  set.seed(SEED)
  b <- Boruta(x = train[, feat_cols], y = y, pValue = tk$bor_p,
              maxRuns = tk$bor_runs, doTrace = 0)
  n_conf <- length(getSelectedAttributes(TentativeRoughFix(b), withTentative = FALSE))
  png(file.path(out_dir, sprintf("Boruta_importance_%s.png", tk$key)),
      width = 11, height = 6, units = "in", res = 600)
  par(mar = c(9, 4, 4, 1))
  plot(b, xlab = "", ylab = "Importance (Z-score)", las = 2, cex.axis = 0.45,
       main = sprintf("Boruta importance - %s", tk$label), cex.main = 1.0)
  mtext(sprintf("pValue = %g, maxRuns = %d; %d features confirmed",
                tk$bor_p, tk$bor_runs, n_conf), side = 3, line = 0.3, cex = 0.8)
  dev.off()

  cat(sprintf("%-34s LASSO lambda %.8f -> %d features; Boruta -> %d confirmed\n",
              tk$label, tk$lambda, n_at_lambda, n_conf))
}

cat("\nFigures written to:", out_dir, "\n")
