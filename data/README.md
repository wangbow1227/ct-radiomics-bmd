# Data

**No patient data is stored in this repository, and none should be committed to it.**
`.gitignore` excludes every data file in this folder for that reason.

The analysis reads three analysis-ready R data files and six feature-list text files,
all supplied by the corresponding author together with this repository:

```
data/
├── train.rds           # training cohort, n = 375
├── test.rds            # testing cohort, n = 159
├── external.rds        # external validation cohort, n = 43
└── feature_sets/
    ├── lowBMD_Boruta.txt    # 28 features
    ├── lowBMD_LASSO.txt     # 20 features
    ├── osteo_Boruta.txt     # 27 features
    ├── osteo_LASSO.txt      # 13 features
    ├── Type3_Boruta.txt     # exploratory three-level model
    └── Type3_LASSO.txt      # exploratory three-level model
```

## Structure

Each file is a data frame with one row per patient and the same columns:

| Column | Contents |
|---|---|
| `ID` | study-specific pseudonymous patient identifier |
| `Type` | BMD category: 1 = normal, 2 = osteopenia, 3 = osteoporosis |
| `Type_new1` | non-osteoporosis (0) vs. osteoporosis (1) |
| `Type_new2` | normal (0) vs. low BMD (1) |
| `Age` | age in years |
| `Sex2` | 1 = female, 0 = male |
| `BMI` | body mass index |
| 474 feature columns | radiomic features (PyRadiomics names), averaged across the four vertebral levels (L1–L4), unscaled |

The training/testing split was drawn once, stratified on the three-level BMD category,
with `caret::createDataPartition` after `set.seed(12345)`; it is fixed in these files, and
the scripts never re-split the data. Feature standardisation is performed by the scripts:
the internal cohort (train + test) is standardised to z-scores using the mean and standard
deviation of the full internal cohort, and the external cohort is standardised to its own
mean and standard deviation.

## Availability

These files contain individual-level patient data collected at National Cheng Kung
University Hospital under IRB approval `A-ER-111-221`, with consent waived for
retrospective analysis. They are not publicly redistributable and must not be uploaded
here. They are available from the corresponding author on reasonable request, subject to
institutional data-sharing approval.

The `feature_sets/` files contain no patient data: each lists the names of the radiomic
features in one reported model, identical to Supplementary Table S3 of the manuscript.
They are kept out of the repository together with the data files; the repository itself
contains code only.
