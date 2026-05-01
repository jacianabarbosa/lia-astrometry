[Leia em Português do Brasil](README.pt-BR.md)

# Lia

**Lia** is a Python pipeline for pre-screening moving-object candidates in IASC FITS image sequences, with Gaia DR3 astrometric refinement and support for manual validation in Astrometrica.

Lia does not replace Astrometrica, MPC validation, or the official IASC workflow. It is designed as a pre-screening and prioritization layer to reduce manual inspection effort and generate auditable outputs for subsequent human validation.

> The name **Lia** is used as a project name and personal dedication by the author. It is not an acronym.

**Author:** Jaciana Barbosa
**Repository:** `lia-astrometry`
**Version:** `v1.5.0`
**License:** MIT

## Overview

IASC/Pan-STARRS practice and campaign sets usually contain four FITS frames of the same field. The scientific task is to find sources that move coherently while most stars and galaxies remain fixed. Lia automates the first screening step: it detects point-like sources, builds candidate tracks, rejects obvious artifacts, ranks the remaining candidates, and writes outputs that can be reviewed manually in Astrometrica.

The output is intentionally conservative. A high-ranked candidate is not a confirmed asteroid discovery. It is a candidate worth measuring or rejecting in the normal Astrometrica workflow.

## Example Output

![Example candidates](images/example_candidates.png)

Visualization produced by Lia for an IASC sequence. Each row is one candidate; columns are the four frames in chronological order. Colored circles mark the measured centroid and arrows show the motion direction.

## Positioning

Lia is a pre-screening pipeline. It prioritizes candidates; it does not confirm discoveries, submit observations, or replace the final measurement process.

Final astrometric measurements must be performed or reviewed in Astrometrica. MPC-style reports produced by Lia are draft or auxiliary outputs only. Official campaign submission should continue to follow the standard Astrometrica/IASC workflow.

The triage score is heuristic. It is useful for ranking and auditing, but it is not a calibrated likelihood estimate.

## Features

- FITS loading with `float64` image handling for numerical stability.
- Mid-exposure JD/MJD extraction from `DATE-OBS` or `MJD-OBS` plus exposure time.
- Local sky estimation with `photutils.Background2D`.
- Saturation masking for Pan-STARRS pixels flagged near `65535`, typical of saturated, defective, or off-CCD regions, preventing spurious detections and inflated local background RMS estimates.
- Source detection with `DAOStarFinder`.
- Sub-pixel centroid refinement using a Moffat PSF model with Gaussian fallback.
- FWHM and morphology checks to reject hot pixels, cosmic-ray-like detections, blends, and extended sources.
- Optional Gaia DR3 WCS refinement with two-stage cross-matching and residual tracking.
- Header-only WCS mode for internal Gaia-vs-header comparisons.
- Four-frame track construction in a local tangent-plane representation.
- Linear kinematic validation with residuals and `R^2`.
- SkyBot/IMCCE neighborhood checks for known Solar System objects.
- JSON, text report, draft MPC text, PNG cutouts, and CSV validation exports.

## Workflow

```text
IASC FITS set
-> Lia pre-screening
-> prioritized candidate list
-> manual inspection and measurement in Astrometrica
-> Astrometrica MPC report
-> IASC submission when appropriate
-> IASC operational feedback, when available
```

## Installation

```bash
git clone https://github.com/jacianabarbosa/lia-astrometry.git
cd lia-astrometry

python3 -m venv venv
source venv/bin/activate

venv/bin/python -m pip install --upgrade pip
venv/bin/python -m pip install -r requirements.txt
venv/bin/python -m pytest tests/test_basic.py -v
```

The current test suite has 92 tests covering public identity, source detection, WCS handling, Gaia refinement, kinematic scoring, saturation masking, header-mode behavior, JSON traceability, and CSV metric export.

## Usage

Basic run:

```bash
venv/bin/python src/detector.py --images images/XY14_p10 --output results/XY14_p10
```

Adjust the detection threshold:

```bash
venv/bin/python src/detector.py --sigma 5.0
venv/bin/python src/detector.py --sigma 6.5
```

Skip the SkyBot/MPC neighborhood query:

```bash
venv/bin/python src/detector.py --no-mpc
```

Observer metadata can be passed directly or through environment variables:

```bash
venv/bin/python src/detector.py --observer "Jaciana Barbosa" --email "name@example.com"
export LIA_OBS="Jaciana Barbosa"
export LIA_EMAIL="name@example.com"
```

## Command-Line Options

```text
usage: detector.py [-h] [--images IMAGES] [--output OUTPUT] [--sigma SIGMA]
                   [--no-mpc] [--observer OBSERVER] [--email EMAIL]
                   [--wcs-mode {gaia,header}]

Lia - pre-screen moving-object candidates in IASC FITS sequences

options:
  -h, --help            show this help message and exit
  --images IMAGES       input folder containing exactly four FITS frames
  --output OUTPUT       output folder for reports, JSON, logs, and plots
  --sigma SIGMA         local-SNR detection threshold, default 5.5
  --no-mpc              skip SkyBot/IMCCE neighborhood checks
  --observer OBSERVER   observer name for draft auxiliary reports
  --email EMAIL         observer email for draft auxiliary reports
  --wcs-mode {gaia,header}
                        WCS mode: Gaia DR3 refinement or reconstructed header WCS
```

## Running With Gaia And Header-Only WCS Modes

Gaia DR3 refinement is the default mode:

```bash
venv/bin/python src/detector.py --images fits/set01 --output results/set01_gaia --wcs-mode gaia
```

Header-only mode disables only Gaia refinement. Background estimation, PSF fitting, tracking, scoring, filters, thresholds, and output generation remain the same:

```bash
venv/bin/python src/detector.py --images fits/set01 --output results/set01_header --wcs-mode header
```

Technical note: Lia does not call `astropy.wcs.WCS(header)` directly for Pan-STARRS frames. These FITS headers can include proprietary distortion keywords such as `PCA1X0Y2` and `PCA2X0Y2`, which the standard parser may interpret as PC matrix entries and fail with singular matrices. Lia reconstructs the header WCS from the primary astrometric keywords (`CTYPE`, `CRVAL`, `CRPIX`, `CD`/`CDELT`, and `CROTA`). Both modes start from this reconstruction; only the Gaia DR3 refinement step differs.

This allows an internal comparison using the same code path:

```text
candidate recovery
ranking changes
RA/Dec differences
angular separation from Astrometrica measurements
Gaia RMS and number of Gaia matches
cases where Gaia improves, does not change, or worsens preliminary coordinates
```

Only the final Astrometrica report should be submitted to IASC. The Gaia/header comparison is methodological and internal; do not submit duplicate reports.

## Outputs

For each processed set, Lia writes:

- `*_candidates.json`: structured run metadata, candidates, tracks, scores, morphology, WCS status, Gaia status, SkyBot status, and manual-validation placeholders.
- `*_report.txt`: readable operational report for candidate review.
- `*_MPC_report.txt`: draft MPC-style auxiliary output for review only.
- `*_candidates.png`: visual cutouts for top candidates.
- `*_pipeline.log`: scientific audit trail for the run.

The JSON contains `run_metadata` with `run_id`, `pipeline_name`, `pipeline_version`, `repository`, `input_set`, `timestamp_execution_utc`, `wcs_mode`, and `sigma`.

Each candidate includes `candidate_id`, `triage_class`, `heuristic_score`, `score_components`, `decision_reasons`, `track`, `motion`, `photometry`, `morphology`, and manual-review fields for later Astrometrica/IASC evaluation. Portuguese legacy fields may remain during the transition for backward compatibility, but the English schema is the public target.

## Exporting Metrics

Lia includes a CSV exporter for validation studies:

```bash
venv/bin/python tools/export_metrics.py results/ --out validation_metrics --recursive
```

This produces:

- `validation_metrics_sets.csv`: one row per processed set.
- `validation_metrics_candidates.csv`: one row per candidate.

The CSV files include empty columns such as `measured_in_astrometrica`, `included_in_astrometrica_mpc`, `iasc_feedback`, `manual_classification`, and `notes` where campaign results can be added after human review.

## Score And Flags

The triage score ranges from 0 to 12 and combines:

| Component | Max | Rationale |
|---|---:|---|
| Linearity | 3 | four-frame track residual around a linear trajectory |
| Velocity consistency | 2 | uniformity of frame-to-frame steps |
| Photometric stability | 2 | relative flux stability across frames |
| Point-source morphology | 2 | rejection of diffuse, blended, or artifact-like sources |
| Morphology consistency | 1 | stability of source shape across frames |
| Elongation | 1 | penalty for trailed or extended detections |
| Typical motion range | 1 | total displacement consistent with short IASC main-belt candidate sequences |

Total: 12. The score is a heuristic for ranking, not a calibrated likelihood.

Classes:

| Score | Class | Meaning |
|---:|---|---|
| 8-12 | `STRONG` | inspect first |
| 6-7 | `MODERATE` | plausible candidate |
| 4-5 | `WEAK` | low-priority candidate |
| 0-3 | `DISCARDED` | likely artifact or rejected source |

The score is a reproducible triage heuristic, not a calibrated statistical confidence value.

## Validation Plan

### Main Operational Evaluation

Run Lia in Gaia mode on IASC sets:

```text
IASC FITS set
-> Lia with Gaia DR3
-> prioritized candidate list
-> manual inspection in Astrometrica
-> Astrometrica MPC report
-> IASC submission
-> IASC operational feedback, when available
```

Primary metrics:

- reduction of candidate search space;
- top-1/top-3/top-5 usefulness;
- candidates measurable in Astrometrica;
- candidates included in an Astrometrica MPC report;
- IASC feedback, when available;
- false positives and failure modes.

### Internal Gaia Comparison

Run the same sets in both modes:

```bash
venv/bin/python src/detector.py --wcs-mode gaia
venv/bin/python src/detector.py --wcs-mode header
```

Compare recovery, ranking, Gaia RMS, match counts, and angular differences relative to Astrometrica measurements. This comparison is internal; only the final Astrometrica-generated report should enter the official IASC workflow.

### Future Statistical Calibration

The heuristic score is a ranking tool. Future work can calibrate it against validated IASC outcomes using logistic regression, ROC analysis, or beta calibration to estimate calibrated real-candidate confidence. This requires a sufficient number of positive and negative outcomes from the Astrometrica/IASC workflow.

## Limitations

- Lia does not replace Astrometrica.
- The triage score is heuristic, not a calibrated statistical confidence value.
- Gaia DR3 refinement may fail, may be unnecessary, or may not improve every field.
- Header WCS may already be sufficient in some IASC frames.
- PSF fitting can fail for low-SNR, blended, saturated, or edge sources.
- `Background2D` parameters require empirical validation on campaign data.
- False-positive rejection can reject real candidates under degraded seeing or severe artifacts.
- IASC feedback is operational validation, not universal ground truth.
- Lia does not perform orbit determination, digest2 scoring, shift-and-stack, CNN classification, or automatic MPC/IASC submission.

## Project Structure

```text
lia-astrometry/
├── README.md
├── docs/
│   └── methodology.md
├── src/
│   └── detector.py
├── tests/
│   └── test_basic.py
├── tools/
│   └── export_metrics.py
├── requirements.txt
└── pytest.ini
```

## Methodology

The full technical methodology is maintained in [docs/methodology.md](docs/methodology.md). It covers FITS metadata handling, mid-exposure timing, Background2D, PSF fitting, Gaia DR3 cross-matching, WCS residuals, kinematic validation, false-positive rejection, output schema, and planned validation.

## References

- Gaia Collaboration et al. (2023), *Gaia Data Release 3*, Astronomy & Astrophysics, 674, A1.
- Stetson, P. B. (1987), *DAOPHOT: A Computer Program for Crowded-Field Stellar Photometry*, Publications of the Astronomical Society of the Pacific, 99, 191.
- Bradley et al. (2024), *astropy/photutils: source detection and photometry tools*.
- Astropy Collaboration et al. (2022), *The Astropy Project: sustaining and growing a community-oriented open-source project*.
- IMCCE SkyBot cone-search documentation, used through `astroquery.imcce.Skybot`.
- Smullen et al. (2025), *TRIPP: TRansient Image Processing Pipeline*, arXiv:2501.18142.
- International Astronomical Search Collaboration (IASC) and AIsteroid legacy detection workflow materials.
- Minor Planet Center astrometry and observation-format documentation.
- IASC public campaign materials and Astrometrica workflow guidance.

## License

MIT License. See [LICENSE](LICENSE).

## Changelog

### v1.5.0

- Renamed the project from previous internal/public names to **Lia**.
- Updated repository identity to `lia-astrometry`.
- Updated documentation for international publication readiness.
- Added or revised scientific traceability fields for validation with IASC/Astrometrica workflows.
- Clarified that Lia is a pre-screening pipeline and does not replace Astrometrica or official MPC/IASC validation.
- Prepared documentation for internal comparison between Gaia DR3 refinement and header-only WCS mode.

### Earlier History

Earlier internal builds developed the core detection, Background2D, PSF fitting, Gaia DR3 WCS refinement, SkyBot checks, audit logs, and draft MPC output. Lia is the current public identity for that work.
