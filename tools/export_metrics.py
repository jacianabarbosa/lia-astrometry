#!/usr/bin/env python3
"""
export_metrics.py — Lia
=======================
Read Lia JSON outputs and export validation metrics to CSV.

  validation_metrics_sets.csv       — one row per processed FITS set
  validation_metrics_candidates.csv — one row per candidate

Usage:
    python tools/export_metrics.py results/ --out validation_metrics

The CSV files include empty manual-review columns for later Astrometrica and
IASC campaign validation.

Fields exported per input set:
    run_id, input_set, pipeline_version, wcs_mode, sigma,
    n_sources_F1..F4, n_tracks_attempted, n_unique_candidates,
    n_rejected_false_positive, n_strong, n_moderate, n_weak, n_discarded,
    top1_candidate_id, top1_score, top1_class,
    top3_candidate_ids, top5_candidate_ids,
    gaia_status_F1..F4, gaia_rms_pre_F1..F4, gaia_rms_pos_F1..F4,
    gaia_n_matches_F1..F4, execution_time_s, data_obs, timestamp_execution_utc

Fields exported per candidate:
    candidate_id, run_id, input_set, wcs_mode, pipeline_version,
    rank, heuristic_score, score_max, triage_class, score_percent,
    s_linearity, s_velocity, s_photometry, s_morphology,
    s_morph_consistency, s_elongation, s_velocity_range,
    flags, ra_deg, dec_deg, vel_arcsec_min, move_total_px,
    linearity_residual_px, step_consistency_px, flux_cv, pointedness,
    mean_snr, mean_fwhm_px, median_position_sigma_arcsec,
    mpc_status, gaia_static_status, n_rejections, n_penalties,
    rejection_reasons, penalty_reasons, data_obs,
    --- manual columns (empty) ---
    measured_in_astrometrica, included_in_astrometrica_mpc, iasc_feedback,
    manual_classification, notes
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def _safe(val, default=""):
    if val is None:
        return default
    return val


def _join(lst, sep="; "):
    if not lst:
        return ""
    return sep.join(str(x) for x in lst if x is not None)


def process_json(path: Path) -> tuple[dict, list[dict]]:
    with open(path, encoding="utf-8") as f:
        dados = json.load(f)

    rm   = dados.get("run_metadata", {})
    mg   = dados.get("global_metrics") or dados.get("global_metrics", {})
    gaia = mg.get("gaia_refinement", [])

    def gaia_campo(idx, field, default=""):
        if idx < len(gaia):
            return _safe(gaia[idx].get(field), default)
        return default

    n_sources = mg.get("n_sources_by_frame") or mg.get("n_fontes_by_frame", [])

    set_row = {
        "run_id"             : _safe(rm.get("run_id") or mg.get("run_id")),
        "repository"         : _safe(rm.get("repository") or dados.get("repository")),
        "input_set"          : _safe(dados.get("input_set")),
        "pipeline_version"   : _safe(rm.get("pipeline_version", dados.get("pipeline", ""))),
        "wcs_mode"           : _safe(rm.get("wcs_mode") or mg.get("wcs_mode")),
        "sigma"     : _safe(mg.get("sigma")),
        "n_sources_F1"       : n_sources[0] if len(n_sources) > 0 else "",
        "n_sources_F2"       : n_sources[1] if len(n_sources) > 1 else "",
        "n_sources_F3"       : n_sources[2] if len(n_sources) > 2 else "",
        "n_sources_F4"       : n_sources[3] if len(n_sources) > 3 else "",
        "n_tracks_attempted" : _safe(mg.get("n_tracks_attempted")),
        "n_unique_candidates": _safe(mg.get("n_unique_candidates")),
        "n_rejected_false_positive"    : _safe(mg.get("n_rejected_false_positive")),
        "n_strong"           : _safe(mg.get("class_distribution", {}).get("STRONG")),
        "n_moderate"         : _safe(mg.get("class_distribution", {}).get("MODERATE")),
        "n_weak"             : _safe(mg.get("class_distribution", {}).get("WEAK")),
        "n_discarded"        : _safe(mg.get("class_distribution", {}).get("DISCARDED")),
        "top1_candidate_id"  : _safe(mg.get("top1_candidate_id")),
        "top3_candidate_ids" : _join(mg.get("top3_candidate_ids", [])),
        "top5_candidate_ids" : _join(mg.get("top5_candidate_ids", [])),
        "gaia_status_F1"     : gaia_campo(0, "status"),
        "gaia_status_F2"     : gaia_campo(1, "status"),
        "gaia_status_F3"     : gaia_campo(2, "status"),
        "gaia_status_F4"     : gaia_campo(3, "status"),
        "gaia_rms_pre_F1"    : gaia_campo(0, "rms_pre_arcsec"),
        "gaia_rms_pre_F2"    : gaia_campo(1, "rms_pre_arcsec"),
        "gaia_rms_pre_F3"    : gaia_campo(2, "rms_pre_arcsec"),
        "gaia_rms_pre_F4"    : gaia_campo(3, "rms_pre_arcsec"),
        "gaia_rms_pos_F1"    : gaia_campo(0, "rms_pos_arcsec"),
        "gaia_rms_pos_F2"    : gaia_campo(1, "rms_pos_arcsec"),
        "gaia_rms_pos_F3"    : gaia_campo(2, "rms_pos_arcsec"),
        "gaia_rms_pos_F4"    : gaia_campo(3, "rms_pos_arcsec"),
        "gaia_n_matches_F1"  : gaia_campo(0, "n_matches"),
        "gaia_n_matches_F2"  : gaia_campo(1, "n_matches"),
        "gaia_n_matches_F3"  : gaia_campo(2, "n_matches"),
        "gaia_n_matches_F4"  : gaia_campo(3, "n_matches"),
        "execution_time_s"   : _safe(mg.get("execution_time_s")),
        "data_obs"           : _safe(dados.get("data_obs")),
        "timestamp_execution_utc" : _safe(rm.get("timestamp_execution_utc")),
    }

    # Infer top candidate fields from the candidate table.
    cands = dados.get("candidates") or dados.get("candidates", [])
    top1 = cands[0] if cands else {}
    set_row["top1_score"]  = _safe(top1.get("score_total"))
    set_row["top1_class"] = _safe(top1.get("triage_class"))

    linhas_cand = []
    input_set    = dados.get("input_set") or dados.get("input_set", "")
    data_obs    = dados.get("data_obs", "")
    wcs_mode    = _safe(rm.get("wcs_mode") or mg.get("wcs_mode"))
    pv          = _safe(rm.get("pipeline_version", dados.get("pipeline", "")))
    run_id      = _safe(rm.get("run_id") or mg.get("run_id"))

    for c in cands:
        sc   = c.get("score_components", {})
        phot  = c.get("photometry", {})
        motion  = c.get("motion", {})
        morph= c.get("morphology", {})
        inc  = c.get("astrometric_uncertainty", {})
        pos1 = c.get("frame1_position") or c.get("posicao_frame1", {})
        mpc  = c.get("mpc", {})
        gs   = c.get("gaia_static", {})
        rd   = c.get("decision_reasons", {})
        manual = c.get("manual_validation", {})

        snrs = phot.get("snr_by_frame") or []
        snr_medio = ""
        if snrs:
            validos = [s for s in snrs if s is not None]
            if validos:
                snr_medio = round(sum(validos) / len(validos), 2)

        linhas_cand.append({
            "candidate_id"           : _safe(c.get("candidate_id")),
            "run_id"                 : run_id,
            "input_set"              : input_set,
            "wcs_mode"               : wcs_mode,
            "pipeline_version"       : pv,
            "data_obs"               : data_obs,
            "rank"                   : _safe(c.get("rank")),
            "heuristic_score"        : _safe(c.get("heuristic_score", c.get("score_total"))),
            "score"                  : _safe(c.get("score_total")),
            "score_max"              : _safe(sc.get("score_max")),
            "triage_class"           : _safe(c.get("triage_class")),
            "score_percent"          : _safe(c.get("score_percent")),
            "s_linearity"          : _safe(sc.get("linearity")),
            "s_velocity"           : _safe(sc.get("velocity")),
            "s_photometry"           : _safe(sc.get("photometry")),
            "s_morphology"           : _safe(sc.get("morphology")),
            "s_morph_consistency"    : _safe(sc.get("morph_consistency")),
            "s_elongation"           : _safe(sc.get("elongation")),
            "s_velocity_range"       : _safe(sc.get("velocity_range")),
            "flags"                  : _join(c.get("flags", [])),
            "ra_deg"                 : _safe(pos1.get("ra_deg")),
            "dec_deg"                : _safe(pos1.get("dec_deg")),
            "vel_arcsec_min"         : _safe(motion.get("rate_arcsec_min") or motion.get("vel_arcsec_min")),
            "move_total_px"          : _safe(motion.get("total_px")),
            "linearity_residual_px"  : _safe(motion.get("linearity_residual_px") or motion.get("linearidade_px")),
            "step_consistency_px"    : _safe(motion.get("step_consistency_px") or motion.get("vel_consistencia")),
            "flux_cv"              : _safe(phot.get("flux_cv")),
            "pointedness"            : _safe(phot.get("pointedness") or phot.get("pointness")),
            "mean_snr"               : snr_medio,
            "mean_fwhm_px"          : _safe(morph.get("mean_fwhm_px")),
            "median_position_sigma_arcsec": _safe(
                inc.get("median_position_sigma_arcsec")
                or inc.get("sigma_pos_mediana_arcsec")
            ),
            "mpc_status"             : _safe(mpc.get("status")),
            "gaia_static_status"     : _safe(gs.get("status")),
            "n_rejections"            : len(rd.get("rejections", [])),
            "n_penalties"         : len(rd.get("penalties", [])),
            "rejection_reasons"        : _join(rd.get("rejections", []), " | "),
            "penalty_reasons"     : _join(rd.get("penalties", []), " | "),
            # Colunas manuais — preencher após inspeção
            "measured_in_astrometrica": _safe(
                manual.get("measured_in_astrometrica")
                or manual.get("medido_astrometrica")
            ),
            "included_in_astrometrica_mpc": _safe(
                manual.get("included_in_astrometrica_mpc")
                or manual.get("entrou_mpc")
            ),
            "iasc_feedback"          : _safe(manual.get("iasc_feedback")),
            "manual_classification"  : _safe(
                manual.get("manual_classification")
                or manual.get("classificacao_manual")
            ),
            "notes"                  : _safe(manual.get("notes") or manual.get("observacoes")),
        })

    return set_row, linhas_cand


def main():
    parser = argparse.ArgumentParser(
        description="Lia — export JSON validation metrics to CSV"
    )
    parser.add_argument(
        "folder",
        type=str,
        help="Directory containing JSON files generated by Lia, for example results/",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="validation_metrics",
        help="Output CSV prefix (default: validation_metrics)",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search for JSON files recursively in subdirectories",
    )
    args = parser.parse_args()

    folder = Path(args.folder)
    if not folder.exists():
        print(f"[ERROR] Folder not found: {folder}", file=sys.stderr)
        sys.exit(1)

    padrao = "**/*_candidates.json" if args.recursive else "*_candidates.json"
    jsons  = sorted(folder.glob(padrao))
    if not jsons:
        legacy_pattern = "**/*_candidates.json" if args.recursive else "*_candidates.json"
        jsons = sorted(folder.glob(legacy_pattern))

    if not jsons:
        print(f"[WARN] No *_candidates.json files found in {folder}", file=sys.stderr)
        sys.exit(0)

    linhas_conj = []
    linhas_cand = []

    for j in jsons:
        try:
            lc, lcs = process_json(j)
            linhas_conj.append(lc)
            linhas_cand.extend(lcs)
        except Exception as e:
            print(f"[WARN] Error processing {j.name}: {e}", file=sys.stderr)

    if not linhas_conj:
        print("[WARN] No JSON file was processed successfully.", file=sys.stderr)
        sys.exit(0)

    # CSV por input_set
    sets_csv = Path(f"{args.out}_sets.csv")
    with open(sets_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(linhas_conj[0].keys()))
        writer.writeheader()
        writer.writerows(linhas_conj)
    print(f"[OK] {sets_csv}  ({len(linhas_conj)} sets)")

    # CSV por candidate
    if linhas_cand:
        candidates_csv = Path(f"{args.out}_candidates.csv")
        with open(candidates_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(linhas_cand[0].keys()))
            writer.writeheader()
            writer.writerows(linhas_cand)
        print(f"[OK] {candidates_csv}  ({len(linhas_cand)} candidates)")
    else:
        print("[WARN] No candidates found in the JSON files.")


if __name__ == "__main__":
    main()
