# Week-8 remission estimator

A single self-contained HTML page, no dependencies, no network calls. Open `index.html`
in any browser, or view the hosted copy:
Open `index.html` in any browser, or launch the interactive calculator here: [UC Week-8 Remission Estimator](https://roshinisundararajan.github.io/ibd-prognostic-vs-predictive/calculator/).

## What it implements

The four-parameter clinical prognostic model from the report (section 7.7), fitted on all
550 discovery patients with uniform shrinkage of 0.917 from 500 bootstrap resamples and
the intercept re-estimated afterwards. The number of predictors was fixed in advance by
Riley's sample-size criteria, which at 63 events support roughly four to six parameters.

```
logit(p) = 2.677 + 0.738·ustekinumab − 0.381·Mayo − 0.581·antiTNF − 0.177·log2(calprotectin + 1)
```

| Predictor | Odds ratio |
|---|---|
| Ustekinumab versus placebo | 2.09 |
| Mayo score, per point | 0.68 |
| Prior anti-TNF exposure | 0.56 |
| Faecal calprotectin, per doubling | 0.84 |

Optimism-corrected AUROC 0.72, calibration slope 0.92. Under ten-times-repeated five-fold
cross-validation with every step refitted inside the training fold, AUROC 0.72 and
calibration slope 0.89.

## Two things it deliberately does not do

It carries **no biopsy term**. The principal-component feature that raises discrimination
to 0.75 needs a full expression array projected onto the discovery cohort's PCA basis,
which is not a point-of-care input.

It applies the **same treatment odds ratio to every patient**. The project found no
treatment-predictive signal, so the difference the page displays is the trial's average
benefit conditioned on that patient's prognosis, not a personalised benefit estimate. The
page says so on its face.

## Status

Research use only. Derived from one induction trial, not externally validated, and not a
clinical decision tool.
