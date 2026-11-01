# UC transcriptomic signatures — 2-day runbook (revised 17 Oct 2026)

Yes, two days works — but not the two days the first version of this runbook
described. This revision exists because I reviewed those scripts properly and
found five bugs, two of which would have produced plausible-looking numbers
that were wrong. Those are the expensive kind: a script that crashes costs an
hour, a script that silently scrambles your clinical data costs the project.

Every script here has been executed end to end against a synthetic dataset
with the same shape as GSE206285 (550 patients, 50 healthy controls, ~13%
remission, missing clinical values, two signature formats). They run.

---

## What changed, and why it mattered

**1. Clinical covariates were silently misaligned.** `model.matrix()` drops
rows containing `NA`. CRP and calprotectin have missing values, so the
clinical block came out shorter than the expression matrix, and `clin_mm[tr, ]`
then indexed *different patients* from `y[tr]`. No error, no warning, wrong
answer. Clinical dummy coding is now done by hand with NAs preserved
(`build_clin_mm`), imputed inside each training fold, and there is an
assertion that stops the run if the row counts ever disagree.

**2. The clinical baseline was the wrong variables.** The old script looked
for clinical data in the series matrix, where GSE206285 keeps only age and
sex. The seven routine variables — CRP, calprotectin, Mayo score, disease
duration, prior anti-TNF, and the two neutrophil scores — are in a *separate
supplementary file*. Script 01 now downloads it. Without it, "does the
transcriptome beat clinical variables?" reduces to "does it beat age and
sex?", which nobody is asking.

**3. Confidence intervals were roughly 40% too narrow.** With three CV
repeats each patient appears three times in the predictions file. The
bootstrap resampled *rows*, treating one patient as three independent people,
and the DeLong test comparing models had the same flaw. Script 05 now
averages predictions to one per patient first and resamples **patients**.
Model-vs-model comparison is a paired bootstrap on the same draws.

**4. Healthy controls were being labelled "treated".** `ifelse(is_placebo,
"placebo", "treated")` assigns everything that is not placebo to the drug arm
— including the 50 non-randomised control biopsies. They survived only by
accident, because their outcome was missing. Arm assignment is now explicit
and anything matching neither arm is dropped and counted.

**5. Outcome mapping would have failed on the first run.** The old defaults
looked for `"yes"`/`"no"`; GSE206285 writes `Y`/`N`. The mapping is now
pre-filled from the actual parsed field names, and if a field name has
changed the script stops and prints the fields that do exist.

Smaller fixes: the univariate filter divided by the overall SD instead of the
pooled within-group SD, which penalises exactly the genes that separate the
groups; calibration deciles crashed on tied predictions; the filter was
vectorised, which took script 04 from tens of minutes to a few.

---

## Two things were added

**A deliberate leak, for the report.** `enet_leaky` picks its genes using the
whole dataset, test folds included. It exists to be wrong. Script 05 prints
the AUROC gap between it and the honest model — that number is your evidence
that in-fold refitting is not a formality. It makes a good slide and it is the
kind of thing that distinguishes a careful project from a competent one.

**An endpoint audit, which runs before any interaction test.** A
treatment-by-biomarker interaction can only exist on an endpoint where the
drug separates from placebo at all. Script 06 checks each endpoint first and
tells you whether the predictive question is even well posed on it. As far as
I can tell no published analysis of this dataset reports this check, and if
the remission and healing endpoints behave differently, that contrast is the
most interesting thing you will find in two days.

---

## Second revision, same day — three more fixes

**Signature scoring was running on the filtered matrix.** Script `02` drops the
bottom quarter of genes by variance before saving, and script `06` then scored
signatures on what was left. Rank-based methods rank a patient's genes against
each other, so removing a quarter of the genome shifts every rank. `02` now
saves `expr_all` beside `expr`, and `06` scores on `expr_all`. If you add xCell
or ssGSEA later, they go on `expr_all` too — xCell on a filtered matrix returns
numbers that look fine and mean nothing.

**The prognostic/predictive rule forced exclusivity.** A marker can be both:
prognostic in the placebo arm *and* modified by treatment. The rule is now a
2×2 on two BH-adjusted q values — prognostic only, predictive only, both, or
null. The primary decision is BH q at 0.10 on each axis, fixed in advance; the
confidence interval is reported beside it but does not make the call.

**The placebo-arm coefficient is now named as such.** With placebo as the
reference level of `arm`, the signature coefficient is its effect *in the
placebo arm*, not a pooled effect. The output column is `OR_placebo_arm` so
nobody misreads the table.

One thing not fixed in code, because it is a decision rather than a bug: if you
build the five-representation grid, do not pick the winner and report its
number as the headline. Choosing the best of five on the same sixty-three
events inflates it. Either select inside the outer fold, or pre-specify pathway
scores as primary before looking and report the rest as secondary.

## The one risk that is still real

**The signature gene lists.** Nothing in script 06 runs without them, and the
old runbook's estimate of "20 minutes" is wrong — supplementary tables arrive
as inconsistent XLSX and PDF, and gene symbols from 2017–2019 papers need
checking against current HGNC names. Budget **three hours**, start it on Day 1
morning in parallel with everything else, and do not let it slip to Day 2.

Format, in `data/signatures/`:

- `<name>.txt` — one gene symbol per line, if the signature has no direction
- `<name>_up.txt` **and** `<name>_dn.txt` — if it does

Use the pair whenever the paper gives direction. Scoring an up/down signature
as if every gene pointed the same way produces a different signature, and
quietly weakens exactly the thing you are testing.

---

## Revised day plan

### Day 1, 09:00 — Roshini
- `00_setup.R`. Start it and leave it. Shorter list than before: pROC, PRROC
  and sva are gone.
- `01_download_annotate.R`. Check two things in the output: that the
  supplementary clinical table was found, and that the printed field names
  include remission and treatment. If the clinical table did not download,
  fetch it from the GEO page by hand now — do not defer this.

### Day 1, 10:30 — Deepikka, in parallel, all morning
- The signature gene lists. This is the critical path for objective 2.
- Script 06 depends only on script 02, so it can run before 04 finishes.

### Day 1, 11:00 — Roshini
- `02_build_cohort.R`. Read the printed tables: arm mapping, outcome by arm,
  exclusions, clinical missingness. If any of them look wrong, stop here —
  everything downstream is wasted work if the label is wrong.
- `03_qc_batch.R`. Fifteen minutes. The randomisation check (PC1 ~ arm)
  should be null; that null is a result worth one line in the report.

### Day 1, afternoon
- `04_models_cv.R`. Now a few minutes rather than an hour, so you can afford
  `N_REPEATS <- 3`. Say in the report which value you used.
- `05_evaluate.R`. Seconds.
- Read the metrics table and the model comparison before doing anything else.

### Day 2, morning — both
- `06_signatures_interaction.R` once the gene lists land.
- Write the results paragraph straight from `out/05_metrics.csv`,
  `out/05_model_comparison.csv` and `out/06_interaction_results.csv`.

### Day 2, afternoon — both
- Report: methods from the proposal, results from `out/`, limitations.
- Slides: one results figure (`05_roc.png` or `05_calibration.png`), one
  table, and the leakage number.

---

## Output files

| File | From | Use |
|---|---|---|
| `01_pheno_<GSE>.csv` | 01 | parsed metadata, one row per sample |
| `01_clinical_<GSE>.csv` | 01 | the supplementary clinical table |
| `01_expr_<GSE>.rds` | 01 | gene-level expression matrix |
| `02_cohort_<GSE>.rds` | 02 | expression + outcome + arm + clinical, aligned |
| `02_exclusions_<GSE>.csv` | 02 | CONSORT-style record of what was dropped |
| `03_pca.png`, `03_qc_boxplot.png`, `03_qc_notes.txt` | 03 | QC for the report |
| `04_cv_predictions.csv` | 04 | out-of-fold predictions, every model |
| `04_stability.csv` | 04 | how often each gene was selected |
| `05_metrics.csv` | 05 | AUROC, PR-AUC, Brier, calibration, ECE + CIs |
| `05_model_comparison.csv` | 05 | paired ΔAUROC against clinical alone |
| `05_roc.png`, `05_calibration.png`, `05_decision_curve.png` | 05 | main figures |
| `06_endpoint_audit.csv` | 06 | does the drug separate, endpoint by endpoint |
| `06_signature_scores.csv` | 06 | per-patient signature scores |
| `06_interaction_results.csv` | 06 | interaction OR, CI, BH-adjusted p |
| `07_replication.csv` | 07 | locked model applied to the second cohort |

---

## Batch correction

Do not ComBat the discovery cohort. One platform, one trial — correcting
within it removes signal and buys nothing.

Use ComBat only to draw a picture of two cohorts on one PCA. Never fit a model
on ComBat output that mixes cohorts, and never pass the outcome into `mod`:
that writes the thing you are trying to predict into the data.

For transfer across platforms, do not merge matrices. Convert to within-sample
ranks and score per patient. Ranks are platform-independent, so the model
transfers with no correction at all. That is what script 07 does.

---

## Cut list — state these as scoped-out future work

- Cross-drug comparison on GSE73661 (objective iii)
- Two-stage inverse-variance meta-analysis (objective iv)
- Causal forest, uplift modelling, per-patient treatment benefit
- Isotonic recalibration — deliberately replaced by two-parameter Platt
  scaling, because ~60 events cannot support a stable isotonic fit. Say that;
  it is a defensible methodological choice, not an omission.
- The web calculator
- Replication (script 07) if Day 2 is tight

Scoping these out is normal and costs nothing. Implying they were done is the
only thing that costs marks.

---

## Honest note on the numbers

At roughly one patient in eight reaching remission, expect a modest AUROC and
expect the transcriptome to add little over clinical variables. That is a
result, and your proposal is already framed to report it that way.

Report the calibration slope and the smallest detectable interaction odds
ratio next to every null. "We found nothing" and "we ruled out effects larger
than this" are the same experiment written two ways, and only one of them is
a finding.

---

## Before you submit

Delete `_make_test_data.R` and anything in `out/` or `data/signatures/` with
`fake` or `test` in the name. Those files exist only so the pipeline could be
proven to run before you had real data in it.
