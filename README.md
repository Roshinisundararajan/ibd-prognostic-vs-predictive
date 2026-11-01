# Prognostic or treatment-predictive? Baseline mucosal transcriptomics in IBD

A placebo-controlled evaluation of baseline mucosal transcriptomic signatures for
prognosis and biologic selection in inflammatory bowel disease.

**AI in Medicine course project, Department of Biotechnology, IIT Madras, 2026.**
Roshini S (BT26S013) and Deepikka SK (BT25D400).

---

## The question

Published transcriptomic signatures have been proposed for choosing between biologics in
IBD, but nearly all were derived in single treated arms, where a marker of *who recovers*
cannot be separated from a marker of *who benefits from the drug*. This project tests those
signatures against a randomised placebo arm, which is the design that can tell the two apart.

Four objectives:

1. **Prediction.** Does baseline expression improve on seven routine clinical variables?
2. **Prognostic versus predictive.** Are published signatures prognostic, or do they modify
   the treatment effect?
3. **Shared response biology.** Do infliximab and vedolizumab share response-associated
   biology?
4. **External transfer and pooling.** Does a locked model transfer to independent cohorts?

## Headline results

| Result | Value |
|---|---|
| Seven clinical variables, out of fold | AUROC 0.704 (0.634 to 0.771) |
| Clinical + expression | AUROC 0.731 (0.675 to 0.786), difference +0.027 (−0.034 to +0.091) |
| Published signatures, prognostic | 9 of 13 |
| Published signatures, treatment-predictive | 0 of 13, bounded below 2.39-fold per SD |
| Hallmark pathways | 1 of 50 prognostic (OXPHOS, OR 2.63), 0 predictive |
| Average treatment effect (AIPW) | +7.0 percentage points |
| Individual benefit ranking (AUTOC) | negative in all 8 learner × feature-set combinations |
| Cross-drug agreement after scan-day adjustment | odds-ratio ratio 1.03 (0.21 to 5.09) |
| Locked model, external transfer | AUROC 0.666 (GSE207022), 0.663 (GSE73661) |
| **Compact prognostic model, held out** | **0.723 → 0.751**, Brier 0.0953 → 0.0933, both intervals exclude zero |
| Feature-selection leakage, measured deliberately | inflates AUROC by +0.073 |

The compact prognostic model is the constructive result and is implemented as a
[bedside calculator](calculator/). The leakage figure is larger than the apparent
contribution of the transcriptome, which is why every number above comes from a pipeline
with all learned steps confined to training folds.

## Cohorts

| Accession | Cohort | n | Events | Role |
|---|---|---|---|---|
| [GSE206285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206285) | UC, ustekinumab 364 / placebo 186 (UNIFI induction) | 550 | 63 clinical remissions at week 8 | Discovery |
| [GSE73661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73661) | UC, infliximab 23 / vedolizumab 40 | 66 | 18 endoscopic healing | Cross-drug, second transfer target |
| [GSE207022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207022) | Crohn's, ustekinumab 63 / placebo 38 | 101 | 12 mucosal healing at week 8 | External validation, pooling |

Cohorts were kept separate throughout. No expression matrices were merged and no
cross-study batch correction was applied.

## Repository layout

```
analysis/        R scripts, run from this directory; they read data/ and write out/
  data/          gene sets and signature lists (small, committed)
results/
  discovery/     GSE206285 outputs (objectives i and ii)
  external/      GSE207022 and GSE73661 outputs (objectives iii and iv)
  causal_ml/     meta-learner, AIPW and class-imbalance outputs
calculator/      self-contained HTML implementation of the compact model
docs/            full results write-up and the final report
```

`results/` holds the committed copies of what the scripts produced. `analysis/out/` is
where a fresh run writes, and is not tracked.

Two result files appear in both `results/discovery/` and `results/external/` with different
contents (`05_metrics.csv`, `06_interaction_results.csv`). They come from different runs and
are deliberately kept apart rather than overwritten.

## Reproducing the analysis

Requires R 4.3 or later.

```r
install.packages(c("data.table", "matrixStats", "glmnet", "ranger", "gbm",
                   "logistf", "msigdbr", "survival"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "singscore"))
```

Then, from inside `analysis/`:

```bash
cd analysis
mkdir -p out
Rscript 00_setup.R                  # check packages
Rscript 00a_fetch_hallmark_genesets.R
Rscript 01_download_annotate.R      # downloads from GEO; slow, several GB
Rscript 02_build_cohort.R
Rscript 04_models_cv.R              # objective (i), model ladder
Rscript 05_evaluate.R
Rscript 06_signatures_interaction.R # objective (ii), 13 signatures
Rscript 10_treatment_benefit.R      # objective (ii), causal ML
Rscript 12_prog_cv.R                # compact prognostic model
Rscript 08_crossdrug_final.R        # objective (iii)
Rscript 07_replicate.R              # objective (iv)
```

Scripts use paths relative to the working directory, so run them from `analysis/`, not from
`analysis/scripts`. `analysis/README_RUNBOOK.md` has the detailed order and expected runtimes.

### A note on the raw data

The GEO series matrices, the platform annotation and the derived `.rds` cohort objects are
**not** in this repository. They total roughly 320 MB, which is past what GitHub is meant to
hold, and every one of them is regenerated by `01_download_annotate.R` and
`02_build_cohort.R` directly from GEO. See `analysis/data/README.md`.

## Methodological notes

- All learned preprocessing (imputation, scaling, univariate filtering, penalty tuning) is
  fitted inside training folds only. Two leaky pipelines were run **deliberately**, to
  measure how much each shortcut inflates performance, and are excluded from every
  substantive comparison.
- A permuted-outcome control reruns the whole pipeline on shuffled labels and returns
  0.479, 0.491 and 0.473 against 0.730 with real labels.
- Interaction models use Firth-penalised logistic regression, with Benjamini-Hochberg
  correction applied across the thirteen signatures and the fifty pathways separately.
- Every null interaction is reported with its minimum detectable effect, so a negative
  finding is a stated bound rather than an absence.
- In GSE73661, array scan day is confounded with both drug and response. 948 genes pass
  FDR 0.10 in the infliximab arm before adjustment and none afterwards, so those genes must
  not be reported as a response signature.

## Citation of the data

Data were obtained from NCBI GEO and are the property of their original depositors. Please
cite the original trials and series: Sands et al., *NEJM* 2019 (GSE206285); Arijs et al.,
*Gut* 2018 (GSE73661); Feagan et al., *NEJM* 2016 and the UNITI programme (GSE207022).
Full reference list in `docs/`.

## Licence

Code in this repository is released under the MIT licence (see `LICENSE`). The report and
figures in `docs/` are coursework by the named authors. Replace or remove `LICENSE` if you
would rather not grant reuse rights.
