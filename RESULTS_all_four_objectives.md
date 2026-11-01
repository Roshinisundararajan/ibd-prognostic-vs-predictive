# UC biologic selection — results for all four objectives

Run 20 Oct 2026. Discovery results (objectives i and ii) are Roshini's, from the
handoff pack. Objectives (iii) and (iv) were run here on the real data.

**What ran in a reduced form.** Three packages could not be installed in the
environment that produced the new results — limma, logistf and singscore
(Bioconductor is unreachable from it). So: the cross-drug statistics are
ordinary t rather than moderated t; the GSE207022 interaction models are plain
logistic regression rather than Firth-penalised; and Hallmark scoring used the
built-in rank score. The last of these is harmless — for an up-only gene set
both engines are linear in the mean within-sample rank, and every score is
standardised, so the per-SD odds ratios are identical. The first two are not
harmless and are flagged where they matter. Re-run both on the laptop before
the figures are final.

---

## Objective (i) — prediction model (discovery, GSE206285)

550 UC patients, 63 remissions (11.5%), 186 placebo and 364 treated.

| Model | AUROC (95% CI) | PR-AUC | Brier | Calibration slope |
|---|---|---|---|---|
| Arm only | 0.505 (0.427–0.579) | 0.109 | 0.101 | 0.65 |
| Seven clinical variables | 0.704 (0.634–0.771) | 0.220 | 0.097 | 0.71 |
| Arm + signature scores | 0.656 (0.588–0.721) | 0.180 | 0.101 | 0.61 |
| Elastic net, expression | 0.709 (0.653–0.762) | 0.198 | 0.103 | 0.58 |
| **Elastic net, clinical + expression** | **0.731 (0.675–0.786)** | 0.219 | 0.100 | 0.64 |
| Random forest | 0.730 (0.675–0.785) | 0.196 | 0.112 | 0.78 |
| *Leaky model (deliberate)* | *0.782 (0.731–0.832)* | *0.256* | *0.095* | *0.93* |

The transcriptome adds about 0.027 AUROC over seven routine clinical variables,
with heavily overlapping intervals. The honest statement is that it does not
measurably beat the clinical baseline at this sample size.

**The leakage demonstration is the most quotable number here.** Selecting genes
on the whole dataset instead of inside each training fold moves AUROC from
0.731 to 0.782 — a gain of 0.051, roughly twice the entire apparent benefit of
adding the transcriptome. A paper reporting 0.78 from a leaky pipeline would
look like a better result than the honest 0.73, and this is how that happens.

---

## Objective (ii) — prognostic versus predictive (discovery)

### The endpoint audit, run first

| Endpoint | n | events | OR (treated vs placebo) | 95% CI | p | arms separate |
|---|---|---|---|---|---|---|
| Clinical remission, week 8 | 550 | 63 | 1.91 | 1.03–3.56 | 0.041 | **yes** |
| Mucosal healing, week 8 | 542 | 83 | 1.08 | 0.66–1.77 | 0.77 | no |

The drug separates from placebo on remission and not on healing. That contrast
is worth stating plainly: the predictive question is well posed on remission
and not on healing, so a null interaction on healing means something different
from a null on remission.

### Interaction tests

Eleven Linggi signatures: **8 prognostic only, 3 null, 0 predictive.** Fifty
Hallmark pathways: **1 prognostic only (oxidative phosphorylation), 49 null, 0
predictive.**

Smallest interaction odds ratio the design could have detected, per standard
deviation: **2.40** for the signatures, **2.46** for the pathways. So the result
is not "we found nothing" — it is "we ruled out treatment-by-biomarker
interactions larger than roughly 2.4-fold per SD."

This is the project's central finding. Every signature that looks informative
in a treated arm is tracking who recovers, not who benefits from this drug.

---

## Objective (iii) — is response biology shared across biologics?

GSE73661 baseline, one sample per patient: infliximab 23 (8 responders),
vedolizumab 40 (9 responders). 10,000 permutations.

| Level | Model | rho | permutation null 95% | p | IFX q<0.10 | VDZ q<0.10 |
|---|---|---|---|---|---|---|
| Gene (17,480) | unadjusted | +0.482 | −0.487 to +0.485 | 0.053 | 1064 | 0 |
| Gene | adjusted for scan day | +0.319 | −0.500 to +0.493 | 0.242 | 0 | 0 |
| Hallmark (50) | unadjusted | +0.675 | −0.722 to +0.719 | 0.084 | 16 | 0 |
| Hallmark | adjusted for scan day | +0.542 | −0.734 to +0.735 | 0.240 | 0 | 0 |

**Read the null band before the p-value.** It is roughly ±0.49 at gene level and
±0.72 at pathway level. On simulated data with independent features the same
test gives a band of ±0.025. The difference is correlation between features:
real genes move in blocks, so the effective number of independent features is
in the dozens, and the correlation between two t-vectors swings enormously
under the null. Any parametric test of this correlation would have called
rho = 0.48 overwhelming evidence of shared response biology. It is not
distinguishable from noise.

The earlier instability is gone: at 2,000 permutations the unadjusted p moved
between 0.046 and 0.059 across runs. At 10,000 it is 0.053 and stable.

**Conclusion.** Response programmes point the same way in both drugs but not
beyond what the null allows, and nothing survives correction for scan day. At
23 and 40 patients this is an underpowered result, not evidence of
drug-specific biology.

### The infliximab signature is scan day

| drug | n | responders | scan days | response vs scan day (Fisher) | R² of response on scan day | residual df |
|---|---|---|---|---|---|---|
| IFX | 23 | 8 | 8 | 0.098 | 0.54 | 14 |
| VDZ | 40 | 9 | 15 | 0.298 | 0.41 | 24 |

Unadjusted, infliximab has 1,064 genes and 16 Hallmark pathways at BH q < 0.10,
and vedolizumab has none of either. After adjusting for scan day: **zero and
zero**, and not one of the 1,064 genes survives.

Scan day explains 54% of the variance in IFX response status and the batches
nearly separate the outcome — the 2012-11-08 arrays are 3 responders and no
non-responders, the July and early-October arrays are all non-responders. Two
readings are defensible and the report should give both: those features are
batch artefact, or the IFX arm cannot separate response from scan day at all.
The design does not distinguish them. What is certain is that 1,064 unadjusted
hits must not be reported as a response signature.

This agrees with the project's own earlier check from the other direction —
`results/03_batch_clustering/batch_confounding_tests.csv` has scan day versus
drug at p = 5e-05 in this cohort, and GSE207022 clean at 0.88.

The adjustment spends 7 of 22 degrees of freedom in the IFX arm, so it removes
real signal along with the artefact. That is a limitation of the cohort, not of
the method.

---

## Objective (iv) — replication

### Model transfer — this one works

The discovery model was locked on all 550 patients, converted to within-sample
ranks, and applied unchanged. No ComBat, no merged matrices.

| Target | shared genes | model size | n | events | AUROC (95% CI) | Brier | Calibration slope |
|---|---|---|---|---|---|---|---|
| GSE207022 (Crohn's, ustekinumab) | 12,571 | 40 genes | 101 | 12 | **0.666 (0.534–0.793)** | 0.111 | 1.65 |
| GSE73661 (UC, infliximab + vedolizumab) | 10,212 | 26 genes | 66 | 18 | **0.663 (0.513–0.798)** | 0.231 | 1.38 |

Both intervals exclude 0.50, though the second only just. The signal transfers
across array platform, across disease (UC to Crohn's), across endpoint, and —
in GSE73661 — across drug class, to patients given infliximab or vedolizumab
rather than ustekinumab.

That last point is the interesting one. A model trained on ustekinumab-treated
and placebo UC patients ranks responders to two entirely different mechanisms
at roughly the same accuracy it achieves in its own cohort. That is what a
prognostic signal looks like, and it corroborates objective (ii) from a
completely different direction.

Both calibration slopes are above 1, meaning the predicted probabilities are
too compressed for the new cohorts: the ordering carries more information than
the probabilities do. Report the slope rather than hiding it. (Note: script 07
prints a sentence about slopes "well below 1" regardless of the value — ignore
that line, it is a bug in the narration, not in the number.)

### Pooled interaction estimate — ran, and should be reported as uninformative

GSE207022 does have a randomised placebo arm (38 placebo, 63 ustekinumab among
patients with an outcome), so the interaction test and the pooling both ran.

All 11 signatures and all 50 pathways: **null** in GSE207022, and **null**
pooled across the two trials. No heterogeneity worth interpreting (I² = 0 for
9 of 11 signatures).

The reason to distrust the pooled number is not the result but the precision:

- GSE207022's endpoint audit fails — OR 1.94 (0.49–7.68), Fisher p = 0.53. The
  drug does not demonstrably separate from placebo there, so an interaction
  tests modification of an effect that was never established in that cohort.
- With 12 events, the smallest interaction odds ratio detectable in GSE207022
  is **20.8 per SD**, against 2.4 on discovery. Nothing biological lives there.
- Inverse-variance weighting gives GSE207022 **7.7%** of the pooled estimate.
  The pooled result is the discovery result with a rounding error attached.
- The two estimates are not estimating the same quantity anyway: UC versus
  Crohn's, clinical remission versus mucosal healing.

Write it as: the pooled estimate was computed, it is null, and it is dominated
by the discovery cohort — the replication trial contributes under a tenth of
the weight and cannot meaningfully confirm or refute an interaction.

**One extra caveat on these specific numbers.** The GSE207022 interaction
models ran without Firth penalisation, because logistf was unavailable. With 12
events and an interaction term, unpenalised logistic regression drifts toward
separation, which is why some intervals run past 40. Re-run with logistf
installed; the calls will not change but the intervals will tighten.

---

## What the four objectives say together

A pre-treatment biopsy carries a real signal about who will reach remission.
That signal is prognostic: it is no better than routine clinical variables at
this sample size, it shows no interaction with treatment on either endpoint
where it was tested, and it transfers to patients given two other drug classes
in a different disease. Nothing in these data supports using a baseline mucosal
transcriptome to choose between biologics, and the design bounds that negative
at interactions larger than about 2.4-fold per standard deviation.

The methodological results stand on their own: 0.051 AUROC of leakage
inflation, and 1,064 differential genes that evaporate when scan day is
accounted for.

---

## Still outstanding

- Re-run objective (iii) with limma, and the GSE207022 interactions with
  logistf. Direction will hold; digits will move.
- Only Linggi signatures were scored. West 2017, Martin 2019 and Czarnewski
  2019 are cited in the proposal as benchmarked signatures and their gene
  lists are not in `data/signatures/`.
- `10_treatment_benefit.R` (per-patient treatment benefit) was not run.
- The web calculator.

## Files

In `results/full_run/`: `05_metrics.csv`, `06_endpoint_audit.csv`,
`06_interaction_results.csv`, `06b_hallmark_interaction.csv` (discovery);
`06_endpoint_audit_GSE207022.csv`, `06_interaction_results_GSE207022.csv`,
`06b_hallmark_interaction_GSE207022.csv`; `07_replication_GSE207022.csv`,
`07_replication_GSE73661.csv`; `08_crossdrug_summary.csv`,
`08_gene_stats_unadjusted.csv`, `08_gene_stats_adjusted.csv`,
`08_hallmark_stats.csv`, `08_confounding_within_drug.csv`; `09_meta.csv`.
Scripts in `scripts/`: `02_build_repl_cohorts.R`, `08_crossdrug_final.R`.
