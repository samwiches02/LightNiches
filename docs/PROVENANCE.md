# Provenance: what came from master_script.R and what is new

Every function and block in this repository carries one of three tags:

- **[YOURS]**: code from Sam's `master_script.R`, logic unchanged.
- **[YOURS, EDITED]**: code from `master_script.R`, with each change listed under the tag.
- **[NEW]**: written by Claude (Anthropic) in October 2026.

Comments without a `[Claude]` marker are Sam's original comments, copied from
`master_script.R`. Comments marked `[Claude]` were added by Claude.

## Summary by file

| File | Mostly | Notes |
|---|---|---|
| `config.R` | yours | Your values and comments; new sections for paths, optics, reference pool, statistics |
| `R/functions_inputs.R` | yours, edited | Tides, capture tide, Kd lookups, hourly table, geom; silent fallbacks now stop |
| `R/functions_raytrace.R` | yours, edited | `prepare_dem()`, ambient, extraction, centroid are yours; azimuth fix, refraction, `lambert = FALSE`, mid-hour sun and the orientation test are new |
| `R/functions_models.R` | new | `reference_kd()`, `spectral_atten()`, `metrics()` and the daytime-integration are yours; ladder, CV, within-site metrics, bootstrap and logger-number test are new |
| `R/utils.R` | new | `make_key()` is yours |
| `scripts/01_inputs.R` | yours, edited | Your steps, reorganised as a stage |
| `scripts/02_raytrace.R` | yours, edited | Your `run_pipeline()` loops, minus Kd and prediction |
| `scripts/03_models.R` | new | Replaces the toggle-and-rerun workflow |
| `scripts/04_figures.R` | mixed | Fig 3, Fig 7, S1, S6 adapted from yours; Figs 4–6 new |
| `run_all.R` | new | |

## Removed from master_script.R

Altitude- and depth-dependent ambient weighting (`w_ambient_combined()` and its
constants), all in-sample fitting grids (`fit_alt_weights()`,
`fit_combined_weights()`, `fit_ambient_weights()`, per-site best W),
`KDPAR_OVERRIDE`, `kd_implied`, `val_strat`, `find_reference_loggers()` and the
hand-pasted `REFERENCE_LOGGER` list, the GT19 azimuth diagnostic, the DEM
verification plots, the run-summary printout, and `metrics_log.csv` appending.

The original `master_script.R` should be kept unchanged in `archive/` as the
record of how the thesis results were produced.
