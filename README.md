# CT Radiomics-Based Classification of Bone Mineral Density

Analysis code for:

> Wang B, Lai Y-H, Lai C-S, Wang C-K, Tai T-W, Tsai Y-S.
> **Abdominal Computed Tomography Radiomics-Based Classification of Bone Mineral Density.**
> *BMC Medical Imaging.*

The study evaluates whether radiomic features extracted from routine non-contrast
abdominal CT can classify DXA-derived bone mineral density (BMD) categories, as an
opportunistic screening adjunct — not a replacement — for DXA.

---

## Repository contents

| File | Purpose |
|---|---|
| `R/main_analysis.R` | Complete analysis: Table 1, Tables 3 and 4 including external validation, Figure 2, confusion matrices, clinical and HU comparison models, calibration and decision-curve analysis |
| `R/reported_models.R` | Fits the four logistic regression models from the feature sets in `data/feature_sets/` and reports their AUCs |
| `R/stability_analysis.R` | Repeats splitting, feature selection and model fitting over 30 stratified 70:30 splits (Supplementary Tables S12 and S13) |
| `R/supplementary_figure_S3.R` | Random forest variable importance (Supplementary Figure S3) |
| `R/selection_figures.R` | LASSO coefficient paths, cross-validation curves, retained features against the LASSO penalty, and Boruta importance |
| `install_packages.R` | Installs the R packages used by the scripts |
| `data/README.md` | Structure of the externally supplied data and feature-list files |

---

## Methods implemented

**Cohort and features.** The 474 radiomic features are averaged across the four vertebral
levels (L1–L4); a patient is included only if a usable record exists at all four levels. The
scripts read three analysis-ready data files — one row per patient, with the averaged,
unscaled features — described under Data availability below.

**Standardisation.** Radiomic features are standardised to z-scores using the mean and
standard deviation of the full internal cohort. The external cohort is standardised to its
own mean and standard deviation.

**Train/test split.** A single 70:30 split, stratified on the three-level BMD category
(normal / osteopenia / osteoporosis), drawn once with `caret::createDataPartition` after
`set.seed(12345)`. The split is at patient level — 375 patients for training and 159 for
testing — and is supplied as the ready-made `train.rds` and `test.rds` files, so the
scripts never re-split the data.

**Feature selection.** Boruta with the package defaults (`pValue = 0.01`,
`maxRuns = 100`), with tentative features resolved by `TentativeRoughFix()`, and LASSO
logistic regression (`glmnet`, `alpha = 1`) with the penalty chosen by cross-validation
(λ = 0.02008833 for normal vs. low BMD, λ = 0.01622183 for osteoporosis). Both methods
involve random resampling, so the features in each model are fixed in six feature-list
files (their contents equal Supplementary Table S3), which the scripts read from
`data/feature_sets/`. Like the data files, they are supplied separately and are not part
of this repository.

| Task | Boruta | LASSO |
|---|---|---|
| Normal vs. low BMD | 28 | 20 |
| Non-osteoporosis vs. osteoporosis | 27 | 13 |

Six features are common to the two osteoporosis feature sets.

**Models.** Logistic regression, decision tree (`rpart`) and random forest
(`randomForest`), with default parameters for the decision tree and random forest. Every
model includes age, sex and BMI together with the radiomic features. A
predicted-probability threshold of 0.5 is applied in the training, testing and external
cohorts. AUCs are computed with `pROC`, with 95% confidence intervals by the method of
Hanley and McNeil.

**Stability analysis.** Splitting, standardisation, Boruta and LASSO selection and model
fitting are repeated over 30 independent stratified 70:30 splits per task, recording the
testing AUC and how often each feature is selected.

---

## Data availability

**No patient data is included in this repository, and none should be committed to it.**

The datasets analysed contain individual-level clinical, DXA and radiomic data from
patients at National Cheng Kung University Hospital. They were collected under IRB
approval `A-ER-111-221` with a waiver of consent for retrospective analysis, and are
not publicly redistributable. They are available from the corresponding author on
reasonable request, subject to institutional data-sharing approval.

To run the analysis, place the three analysis-ready data files in `data/`. Each holds one
row per patient — a study ID, the BMD category and the two binary outcomes, age, sex, BMI,
and the 474 L1–L4-averaged (unscaled) radiomic features. The exact structure is described
in [`data/README.md`](data/README.md). These files are supplied by the corresponding
author together with this repository; they are not part of it.

```
data/
├── train.rds           # training cohort, n = 375
├── test.rds            # testing cohort, n = 159
├── external.rds        # external validation cohort, n = 43
└── feature_sets/       # six feature-list text files (Supplementary Table S3)
```

This repository contains code only: no patient data, no feature-list files, and no fitted
model objects are committed, and `.gitignore` excludes everything under `data/` except its
README to keep it that way.

---

## Requirements

Developed and tested on **R 4.3.1**.

```r
source("install_packages.R")
```

| Package | Version used | Purpose |
|---|---|---|
| dplyr, tibble, purrr | 1.2.1 / 3.3.1 / 1.2.2 | data manipulation |
| caret | 7.0-1 | stratified splitting, confusion matrices |
| rpart | 4.1.27 | decision tree |
| randomForest | 4.7-1.2 | random forest |
| glmnet | 4.1-9 | LASSO feature selection |
| Boruta | 8.0.0 | Boruta feature selection |
| logistf | 1.26.1 | adjusted odds ratios (Firth-penalised logistic regression) |
| pROC | 1.19.1 | ROC curves and AUC |
| dcurves | 0.5.0 | decision-curve analysis |
| ggplot2, gridExtra | 4.0.3 / 2.3.1 | figure panels |

---

## Running the analysis

```bash
Rscript R/main_analysis.R
```

By default the models are built from the feature sets in `data/feature_sets/`. To run
Boruta and LASSO selection on the training cohort instead:

```bash
BMD_FEATURES=select Rscript R/main_analysis.R
```

The other scripts are run the same way:

```bash
Rscript R/reported_models.R
Rscript R/stability_analysis.R
Rscript R/supplementary_figure_S3.R
Rscript R/selection_figures.R
```

`stability_analysis.R` takes roughly 20–40 minutes, because Boruta is run on all 474
features for every split and task.

---

## Outputs

Each script writes to its own folder under `output/`.

| Output | Corresponds to |
|---|---|
| `main_analysis/Table1_NormalVsLowBMD.csv`, `main_analysis/Table1_NonosteoVsOsteo.csv` | Table 1 |
| `main_analysis/Table3_4.csv` | Tables 3 and 4, training and testing cohorts |
| `main_analysis/Table3_4_external.csv`, `main_analysis/external_auc_ci.csv` | Tables 3 and 4, external validation cohort, with AUC 95% confidence intervals |
| `main_analysis/Figure2.png` | Figure 2 |
| `main_analysis/confusion_matrices_train_test.csv`, `main_analysis/confusion_matrices_external.csv` | Supplementary Table S9 |
| `main_analysis/baseline_model_comparison.csv` | Supplementary Table S10 (clinical-only, HU-only, clinical + HU) |
| `main_analysis/calibration_metrics.csv` | Supplementary Table S11 |
| `main_analysis/calibration_plot_*.png` | Supplementary Figure S5 |
| `main_analysis/dca_net_benefit_*.csv`, `main_analysis/dca_plots.png` | Supplementary Figure S6 |
| `main_analysis/misclassified_external_cases.csv` | misclassified external cases described in the Results |
| `main_analysis/selected_feature_sets.csv` | features used in each model |
| `main_analysis/adjusted_odds_ratios.csv` | adjusted odds ratios for each feature set |
| `stability/` | Supplementary Tables S12 and S13 |
| `supplementary_figure_S3/` | Supplementary Figure S3 |
| `reported_models/reported_models_performance.csv` | AUCs of the four logistic regression models |
| `selection_figures/` | LASSO and Boruta figures |

---

## License

Released under the MIT License — see [`LICENSE`](LICENSE).

The license covers the analysis code only. It does not cover the patient data, which is
not distributed here, nor the manuscript text and figures.

---

## Funding

Supported by the National Science and Technology Council, Taiwan — *AI-Powered Cancer
Databases Integrating Clinical Information, Digital Pathology, Medical Imaging and
Genomics for Precision Medicine*, Grant No. NSTC 112-2321-B-006-022.

## Contact

Correspondence: Dr. Yi-Shan Tsai, MD — Department of Diagnostic Radiology,
National Cheng Kung University Hospital, College of Medicine,
National Cheng Kung University, Tainan, Taiwan.
