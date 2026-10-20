# LightNiches

Predicts underwater PAR at every pixel of a photogrammetric reef DEM by combining
surface PAR, refracted sun position, tide-corrected depth, a diffuse attenuation
coefficient (Kd), and ray-traced direct and ambient shading. Validated against
HOBO light loggers at Jiigurru (Lizard Island), Great Barrier Reef.

## Quick start

```r
# once, in the project folder
install.packages("renv")
renv::restore()        # installs the exact package versions in renv.lock
```

```sh
Rscript run_all.R      # runs every stage in order, each in a fresh R session
```

The first run of stage 03 writes `output/reference_pool.csv` and stops on
purpose, before any model is scored. Check the pool, then resume:

```sh
Rscript run_all.R 03
```

## Stages

| Stage | Script | Reads | Writes | Time |
|---|---|---|---|---|
| 01 | `scripts/01_inputs.R` | `data/` | `output/hourly_inputs.rds`, `capture_tide.csv`, `kd_satellite.csv` | seconds |
| 02 | `scripts/02_raytrace.R` | stage 01 output, `data/dems/` | `output/cache.rds`, `output/raytrace/`, `output/dems/` | hours |
| 03 | `scripts/03_models.R` | `output/cache.rds` | reference pool, eval rows, predictions, all metrics | minutes |
| 04 | `scripts/04_figures.R` | stage 02–03 output | `output/figures/` | minutes |

Stages talk to each other only through files in `output/`. Stage 02 saves each
site as it finishes and skips sites already done, so an interrupted run resumes
where it stopped (`FORCE_RAYTRACE <- TRUE` in `config.R` redoes them).

Each stage writes `output/manifest_<stage>.txt` with the run time, git commit,
every setting in `config.R`, and the R session.

## Where things are

```
config.R                 every setting (the only place any is assigned)
run_all.R                runs the stages
R/utils.R                keys, paths, manifests
R/functions_inputs.R     tides, capture-tide correction, Kd lookups, hourly table
R/functions_raytrace.R   DEM preparation, sun angles, shading, extraction
R/functions_models.R     reference Kd, model ladder, cross-validation, metrics
scripts/01..04           the four stages
data/                    raw inputs (not in git; see below)
output/                  everything generated (not in git)
```

## Data (`data/`)

| File | Source |
|---|---|
| `env_long_depth_kd_longlat.csv` | logger lux, surface PAR, tide, logger elevation (dpHB); built by `read_environment.R` from data collected by Brambilla et al. |
| `geom_new.csv` | logger positions annotated on the DEMs |
| `kdpar_lookup.csv`, `kd_raw_per_wavelength.csv` | satellite Kd(λ), NESP TWQ 5.3 MODIS product, via `kdpar_full.R` |
| `tides/tides_Lizard{2017,2018,2019}.csv` | Queensland Government tide records |
| `dems/<SITE><yy>.tif` | photogrammetry DEMs, e.g. `CB19.tif` |

Data availability: [repository, DOI — to add].

## Models

| Model | Components |
|---|---|
| M0 | surface PAR only |
| M1 | + constant Kd (0.2 m⁻¹, a priori) |
| M2s / M2r | satellite Kd / reference-logger Kd instead of constant |
| M3s / M3r | + direct shading with the refracted sun |
| M4s / M4r | + ambient shading, one weight chosen by leave-one-event-out CV |

Diagnostics: `M3r_unrefracted` (in-air sun), `M4r_siteKd` (site-constant
reference Kd), `M4r_buffer_2cm` / `_10cm` (buffer sensitivity), `S_spectralKd`
(band-by-band Kd, supplement).

All models are scored on the same rows (`output/eval_keys.csv`), which exclude
every logger in the reference pool.

## Citation

[paper citation and Zenodo DOI — to add]
