"""
=============================================================================
Lia — detector.py
=============================================================================
Pre-screening pipeline for moving-object candidates in IASC/Pan-STARRS FITS
image sequences.

Author  : Jaciana Barbosa
Version : 1.5.0
License : MIT

Description:
    Lia detects coherent moving-object candidates across four FITS frames,
    estimates preliminary astrometric coordinates, checks nearby known Solar
    System objects through SkyBot/IMCCE, assigns a heuristic triage score, and
    writes auditable outputs for subsequent human validation in Astrometrica.

    Lia is a pre-screening and prioritization layer. It does not replace
    Astrometrica, MPC validation, the official IASC workflow, or human
    judgement.

Usage:
    python detector.py --images ../images/ --output ../results/
    python detector.py --wcs-mode gaia    # Gaia DR3 refinement (default)
    python detector.py --wcs-mode header  # header-only WCS mode

Dependencies:
    pip install -r requirements.txt

Changelog v1.5.0:
    - Public project identity changed to Lia; repository identity changed to
      lia-astrometry.
    - Added --wcs-mode {gaia,header} for internal Gaia-vs-header comparisons.
    - Added run_id, candidate_id, run_metadata, WCS mode, mid-exposure MJD/JD,
      astrometric uncertainty estimates, and manual validation placeholders.
    - Clarified the heuristic score as a triage score, not a calibrated
      confidence estimate.
    - Added CSV metric export support for validation studies.
=============================================================================
"""

import os
import sys
import json
import time
import uuid
import argparse
import logging
import warnings
import multiprocessing as mp
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord, Angle
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from astropy.modeling import fitting
from astropy.modeling.models import Gaussian2D, Moffat2D
import astropy.units as u

from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder

warnings.filterwarnings("ignore")

PROJECT_NAME = "Lia"
REPOSITORY_NAME = "lia-astrometry"
PIPELINE_VERSION = "1.5.0"
VERSION = PIPELINE_VERSION

# ─────────────────────────────────────────────
# THRESHOLDS CENTRALIZADOS
# Edite aqui para tuning futuro sem caçar valores pelo código.
# ─────────────────────────────────────────────
class T:
    # Detecção
    DETECTION_SIGMA       = 5.5       # sigma local de detecção (sobrescrito por --sigma)
    FWHM_PX              = 3.0       # FWHM assumido para o finder
    BACKGROUND_BOX       = (50, 50)  # malha local do céu para Background2D
    BACKGROUND_FILTER    = (3, 3)    # suaviza a malha sem apagar gradientes reais
    PSF_FIT_RADIUS_PX      = 8         # janela local para ajuste PSF sub-pixel
    PSF_MODEL           = "moffat"  # Moffat modela melhor asas de seeing atmosférico
    PSF_MAX_ITER         = 100
    FWHM_MIN_PX          = 1.2       # abaixo disso tende a hot pixel/raio cósmico
    FWHM_MAX_PX          = 8.0       # acima disso tende a galáxia, blend ou trail
    SHARPNESS_MIN        = 0.2
    SHARPNESS_MAX        = 1.0
    ROUNDNESS_MAX        = 0.8       # |roundness| máximo no finder

    # Associação entre frames
    MATCH_RADIUS_PX        = 20.0      # janela de casamento mútuo [px]

    # Filtro de motion (pós-deriva)
    MOVE_MIN_PX          = 2.0       # deslocamento mínimo total [px]
    MOVE_MAX_PX          = 120.0     # deslocamento máximo total [px]
    MIN_RESIDUAL_PX       = 1.8       # resíduo mínimo em relação à deriva

    # Refinamento WCS via Gaia DR3
    GAIA_QUERY_RADIUS_ARCMIN  = 6.0    # raio do cone search em torno do CRVAL
    GAIA_MAG_LIMIT            = 19.0   # G_mag máximo para evitar sources faint demais
    GAIA_MAG_MIN            = 12.0   # evitar saturadas (Pan-STARRS satura cedo)
    GAIA_MIN_MATCHES        = 8      # estrelas casadas mínimas para refinar (passada fina)
    GAIA_MATCH_RADIUS_ARCSEC  = 2.0    # tolerância de cross-match por estrela (passada fina)
    GAIA_DIFF_MIN_ARCSEC    = 0.3    # só tenta refinar se RMS_pre > este valor
    GAIA_DIFF_MAX_ARCSEC    = 0.9    # aceita result Gaia se RMS_pos < este (bruto ou afim)
    GAIA_TIMEOUT_S          = 30     # timeout da query
    GAIA_COARSE_RADIUS_ARCSEC  = 15.0   # raio da passada bruta (WCS Pan-STARRS tem ~5-7 arcsec erro)
    GAIA_MIN_COARSE_MATCHES  = 5      # mínimo de matches para estimar offset grosseiro
    GAIA_OFFSET_BIN_ARCSEC  = 1.5     # 2D voting bin for coarse offset
    GAIA_REFINEMENT_RADIUS_ARCSEC = 2.5     # cluster radius around the coarse-vote peak

    # False-positive rejection
    EDGE_MARGIN_PX      = 30.0      # minimum distance to frame edge
    HOT_PIXEL_POINTNESS_MAX   = 0.70      # pointness > threshold → hot pixel suspect
    HOT_PIXEL_SNR_MIN    = 3.0       # minimum SNR to avoid hot-pixel rejection
    CHAOTIC_FLUX_CV       = 1.5       # CV above this value means chaotic flux
    MIN_SNR           = 2.0       # minimum SNR in at least 3 of 4 frames

    # Score: linearity (0–3)
    LIN_EXCELENTE        = 0.8       # residual < → 3 pts
    LIN_GOOD              = 1.5       # residual < → 2 pts
    LIN_MARGINAL         = 2.5       # residual < → 1 pt; senão 0

    # Score: velocity/consistência (0–2)
    VEL_UNIFORM         = 3.0       # σ passos < → 2 pts
    VEL_MODERATE         = 5.0       # σ passos < → 1 pt; senão 0

    # Score: photometry / CV (0–2)
    PHOT_STABLE          = 0.25      # CV < → 2 pts
    PHOT_MODERATE         = 0.50      # CV < → 1 pt; senão 0  (era 0.45)

    # Score: morphology / pointness (0–2)
    MOR_POINTLIKE          = 0.12      # pont > → 2 pts  (era 0.10)
    MOR_MARGINAL         = 0.06      # pont > → 1 pt; senão 0  (era 0.05)

    # Score: elongation (0–1, novo componente)
    ELON_GOOD             = 1.6       # elong < → 1 pt; senão 0

    # Score: consistência morfológica entre frames (0–1, novo)
    MORPH_CONSIST_MAX_STD  = 0.08      # std(pontual_by_frame) < → 1 pt

    # Score: faixa de velocity típica de MBA (0–1)
    VEL_RANGE_MIN_PX     = 4.0
    VEL_RANGE_MAX_PX     = 50.0

    # Classification final
    SCORE_STRONG          = 8         # score >= → STRONG
    SCORE_MODERATE       = 6         # score >= → MODERATE
    SCORE_WEAK          = 4         # score >= → WEAK; senão DISCARDED

    # Mascaramento de pixels inválidos (Pan-STARRS marca saturados com 65535)
    SAT_PIXEL_VALUE      = 65535      # valor que indica saturação
    SAT_PIXEL_TOLERANCE = 35         # margem (mascara também 65500..65535)
    MIN_VALID_PIXEL     = 1          # abaixo disso é "buraco" / pixel ruim
    FRAME_SAT_WARNING_PCT  = 5.0        # >% saturados → avisa qualidade ruim
    FRAME_SAT_REJECT_PCT = 25.0      # >% saturados → aborta o pipeline

    # Deduplicação
    DEDUP_RADIUS_PX        = 25.0

    # MPC: minimum score for inclusion in the 80-column report
    MPC_SCORE_MIN        = 6

    # Drift: threshold for a stable source
    STABLE_DRIFT_PX    = 2.0
    MIN_STABLE_DRIFT_SOURCES  = 5

    # Astrometry: uncertainty and Gaia static-source rejection
    GAIA_STATIC_MATCH_ARCSEC = 1.5
    GAIA_STATIC_MIN_FRAMES   = 3


# ─────────────────────────────────────────────
# Cores para terminal + logging
# ─────────────────────────────────────────────
class Color:
    GREEN    = "\033[92m"
    YELLOW  = "\033[93m"
    RED = "\033[91m"
    BLUE     = "\033[94m"
    BOLD  = "\033[1m"
    RESET    = "\033[0m"

_logger = logging.getLogger(PROJECT_NAME)

def log(msg, level="INFO"):
    cores = {"INFO": Color.BLUE, "OK": Color.GREEN,
             "WARN": Color.YELLOW, "ERROR": Color.RED}
    cor = cores.get(level, Color.RESET)
    print(f"{cor}[{level}]{Color.RESET} {msg}")
    nivel_py = {"INFO": logging.INFO, "OK": logging.INFO,
                "WARN": logging.WARNING, "ERROR": logging.ERROR}.get(level, logging.INFO)
    _logger.log(nivel_py, msg)


def configure_file_logger(output_dir: Path, input_set_name: str) -> Path:
    log_file = output_dir / f"{input_set_name}_pipeline.log"
    handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    _logger.setLevel(logging.INFO)
    for h in list(_logger.handlers):
        _logger.removeHandler(h)
    _logger.addHandler(handler)
    return log_file


# ─────────────────────────────────────────────
# 1. CARREGAMENTO DAS IMAGENS
# ─────────────────────────────────────────────
def _build_panstarrs_wcs(header: fits.Header) -> WCS:
    """
    Constrói WCS a partir das chaves primárias do header Pan-STARRS.

    Por que não usar `WCS(header)` direto:
        Os FITS do Pan-STARRS incluem chaves proprietárias de distorção
        polinomial (PCA1X0Y2, PCA2X0Y2, …) que o parser do astropy
        interpreta como matriz PCi_j e encontra valores próximos de zero,
        resultando em SingularMatrixError. Solução: construir o WCS
        manualmente apenas com CTYPE/CRVAL/CRPIX/CDELT.
    """
    if not all(k in header for k in ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2",
                                     "CRPIX1", "CRPIX2")):
        raise ValueError("Header FITS does not contain the minimum WCS keywords "
                         "(CTYPE/CRVAL/CRPIX).")

    w = WCS(naxis=2)
    w.wcs.ctype = [header["CTYPE1"], header["CTYPE2"]]
    w.wcs.crval = [float(header["CRVAL1"]), float(header["CRVAL2"])]
    w.wcs.crpix = [float(header["CRPIX1"]), float(header["CRPIX2"])]

    if all(k in header for k in ("CD1_1", "CD1_2", "CD2_1", "CD2_2")):
        w.wcs.cd = [[float(header["CD1_1"]), float(header["CD1_2"])],
                    [float(header["CD2_1"]), float(header["CD2_2"])]]
    elif "CDELT1" in header and "CDELT2" in header:
        w.wcs.cdelt = [float(header["CDELT1"]), float(header["CDELT2"])]
        if "CROTA1" in header or "CROTA2" in header:
            crota = float(header.get("CROTA2", header.get("CROTA1", 0.0)))
            w.wcs.crota = [crota, crota]
    else:
        raise ValueError("Header FITS sem CD matrix nem CDELT — "
                         "escala do pixel desconhecida.")

    if "EQUINOX" in header:
        w.wcs.equinox = float(header["EQUINOX"])
    if "RADESYS" in header:
        w.wcs.radesys = str(header["RADESYS"]).strip()
    return w


# ─────────────────────────────────────────────
# 1b. REFINAMENTO ASTROMÉTRICO VIA GAIA DR3
# ─────────────────────────────────────────────

def _query_gaia_dr3(ra_centro: float, dec_centro: float,
                        raio_arcmin: float = None,
                        mag_lim: float = None,
                        mag_min: float = None):
    """
    Consulta o catálogo Gaia DR3 via TAP do ESA.

    Física/astrometria:
        Gaia DR3 (epoch 2016.0) é o catálogo astrométrico mais preciso
        disponível, com erros típicos < 1 mas em RA/Dec. Limites de
        magnitude evitam estrelas saturadas no Pan-STARRS (G < 12) e
        sources faint cuja detecção falha no detector (G > 19).
        RUWE < 1.4 filtra sources com astrometria degradada (binárias
        não resolvidas, extensões reais).

    Retorna astropy.Table com colunas ra, dec, pmra, pmdec,
    phot_g_mean_mag, ref_epoch. Retorna None se falhar ou < GAIA_MIN_MATCHES.
    """
    try:
        from astroquery.gaia import Gaia
        Gaia.ROW_LIMIT = 5000
        if raio_arcmin is None:
            raio_arcmin = T.GAIA_QUERY_RADIUS_ARCMIN
        if mag_lim is None:
            mag_lim = T.GAIA_MAG_LIMIT
        if mag_min is None:
            mag_min = T.GAIA_MAG_MIN

        query = f"""
        SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, ref_epoch
        FROM gaiadr3.gaia_source
        WHERE 1=CONTAINS(POINT('ICRS', ra, dec),
                         CIRCLE('ICRS', {ra_centro}, {dec_centro},
                                {raio_arcmin / 60.0}))
          AND phot_g_mean_mag BETWEEN {mag_min} AND {mag_lim}
          AND ruwe < 1.4
        """
        job = Gaia.launch_job_async(query, dump_to_file=False,
                                    verbose=False)
        tabela = job.get_results()
        if len(tabela) < T.GAIA_MIN_MATCHES:
            return None
        return tabela
    except Exception:
        return None


def _apply_gaia_proper_motion(tabela_gaia, data_obs_jd: float):
    """
    Propaga posições Gaia DR3 (epoch 2016.0) para a época da observação.

    Física:
        Estrelas de alto motion próprio podem se deslocar dezenas de
        mas em poucos anos. Sem corrigir, o cross-match falha exatamente
        para as estrelas mais brilhantes e próximas — as mais úteis para
        calibrate the WCS because they have lower astrometric uncertainty.
    """
    from astropy.time import Time
    t_gaia = Time(2016.0, format="jyear")
    t_obs  = Time(data_obs_jd, format="jd")
    dt_yr  = float((t_obs - t_gaia).to(u.year).value)

    pmra  = np.array(tabela_gaia["pmra"].filled(0.0),  dtype=np.float64)
    pmdec = np.array(tabela_gaia["pmdec"].filled(0.0), dtype=np.float64)
    ra    = np.array(tabela_gaia["ra"],  dtype=np.float64)
    dec   = np.array(tabela_gaia["dec"], dtype=np.float64)

    ra_corr  = ra  + pmra  * dt_yr / 3.6e6 / np.cos(np.radians(dec))
    dec_corr = dec + pmdec * dt_yr / 3.6e6

    tabela_gaia = tabela_gaia.copy()
    tabela_gaia["ra"]  = ra_corr
    tabela_gaia["dec"] = dec_corr
    return tabela_gaia


def _delta_ra_degrees(ra_a, ra_b):
    """Menor diferença angular RA_a - RA_b em graus, com wrap 0/360."""
    return (np.asarray(ra_a, dtype=np.float64) - np.asarray(ra_b, dtype=np.float64) + 180.0) % 360.0 - 180.0


def _cross_match_gaia(deteccoes_radec: np.ndarray,
                      gaia_radec: np.ndarray,
                      raio_arcsec: float = None):
    """
    Casa detecções (RA/Dec via WCS) com estrelas Gaia usando KDTree.

    Retorna (idx_det, idx_gaia) — índices alinhados das detecções e
    estrelas casadas dentro de raio_arcsec.

    Física:
        A aproximação esférica local (projeção plana sobre (RA*cosDec, Dec))
        é válida para fields de ~10 arcmin onde o erro de curvatura é
        < 0.01 arcsec — muito abaixo da precisão alvo.
    """
    from scipy.spatial import cKDTree
    if raio_arcsec is None:
        raio_arcsec = T.GAIA_MATCH_RADIUS_ARCSEC

    ra_med  = float(np.median(gaia_radec[:, 0]))
    dec_med = float(np.median(gaia_radec[:, 1]))
    cosdec  = np.cos(np.radians(dec_med))

    det_xy = np.column_stack([
        _delta_ra_degrees(deteccoes_radec[:, 0], ra_med) * cosdec * 3600.0,
        (deteccoes_radec[:, 1] - dec_med) * 3600.0,
    ])
    gaia_xy = np.column_stack([
        _delta_ra_degrees(gaia_radec[:, 0], ra_med) * cosdec * 3600.0,
        (gaia_radec[:, 1] - dec_med) * 3600.0,
    ])

    tree = cKDTree(gaia_xy)
    dist, idx = tree.query(det_xy, k=1, distance_upper_bound=raio_arcsec)
    valid = np.isfinite(dist) & (dist < raio_arcsec)
    return np.where(valid)[0], idx[valid]


def _estimate_coarse_offset(
    deteccoes_radec: np.ndarray,
    gaia_radec: np.ndarray,
    raio_arcsec: float = None,
    retornar_indices: bool = False,
) -> tuple:
    """
    Passada bruta do cross-match Gaia: raio grande para tolerar o erro
    inicial de 5-7 arcsec do WCS Pan-STARRS.

    Retorna (offset_ra_arcsec, offset_dec_arcsec, n_matches), onde os
    offsets são det - gaia (positivo = detecção deslocada para leste/norte).
    Com retornar_indices=True, retorna também os índices das detecções e
    estrelas Gaia que participaram do pico de translação.

    Física:
        O WCS Pan-STARRS pode ter erro sistemático de vários arcsec. Em field
        denso, vizinho mais próximo direto é enviesado: a estrela Gaia errada
        pode estar mais perto que a correspondente verdadeira. Por isso esta
        etapa usa votação 2D de todos os pares dentro do raio bruto e escolhe
        o pico de translação coerente antes de calcular a mediana robusta.
    """
    from scipy.spatial import cKDTree
    if raio_arcsec is None:
        raio_arcsec = T.GAIA_COARSE_RADIUS_ARCSEC

    ra_med  = float(np.median(gaia_radec[:, 0]))
    dec_med = float(np.median(gaia_radec[:, 1]))
    cosdec  = np.cos(np.radians(dec_med))

    det_xy = np.column_stack([
        _delta_ra_degrees(deteccoes_radec[:, 0], ra_med) * cosdec * 3600.0,
        (deteccoes_radec[:, 1] - dec_med) * 3600.0,
    ])
    gaia_xy = np.column_stack([
        _delta_ra_degrees(gaia_radec[:, 0], ra_med) * cosdec * 3600.0,
        (gaia_radec[:, 1] - dec_med) * 3600.0,
    ])

    tree = cKDTree(gaia_xy)
    pares_det = []
    pares_gaia = []
    deltas = []
    for i_det, p in enumerate(det_xy):
        for i_gaia in tree.query_ball_point(p, r=raio_arcsec):
            delta = p - gaia_xy[i_gaia]
            if abs(delta[0]) <= raio_arcsec and abs(delta[1]) <= raio_arcsec:
                pares_det.append(i_det)
                pares_gaia.append(i_gaia)
                deltas.append(delta)

    if not deltas:
        vazio = (np.array([], dtype=int), np.array([], dtype=int))
        return (0.0, 0.0, 0, *vazio) if retornar_indices else (0.0, 0.0, 0)

    deltas = np.asarray(deltas, dtype=np.float64)
    pares_det = np.asarray(pares_det, dtype=int)
    pares_gaia = np.asarray(pares_gaia, dtype=int)

    bin_arcsec = float(T.GAIA_OFFSET_BIN_ARCSEC)
    bins = np.arange(-raio_arcsec, raio_arcsec + bin_arcsec, bin_arcsec)
    hist, x_edges, y_edges = np.histogram2d(deltas[:, 0], deltas[:, 1],
                                            bins=[bins, bins])
    peak = np.unravel_index(np.argmax(hist), hist.shape)
    centro = np.array([
        0.5 * (x_edges[peak[0]] + x_edges[peak[0] + 1]),
        0.5 * (y_edges[peak[1]] + y_edges[peak[1] + 1]),
    ])

    raio_refino = float(T.GAIA_REFINEMENT_RADIUS_ARCSEC)
    dist_peak = np.hypot(deltas[:, 0] - centro[0], deltas[:, 1] - centro[1])
    cluster = dist_peak <= raio_refino

    # Garante pareamento 1:1 dentro do pico: para cada detecção, mantém o par
    # mais próximo do centro do pico. Isso reduz blends e duplicatas Gaia.
    candidates = np.where(cluster)[0]
    if candidates.size:
        ordem = candidates[np.argsort(dist_peak[candidates])]
        usados_det = set()
        usados_gaia = set()
        escolhidos = []
        for k in ordem:
            d = int(pares_det[k])
            g = int(pares_gaia[k])
            if d in usados_det or g in usados_gaia:
                continue
            usados_det.add(d)
            usados_gaia.add(g)
            escolhidos.append(k)
        escolhidos = np.asarray(escolhidos, dtype=int)
    else:
        escolhidos = np.array([], dtype=int)

    n_val = int(len(escolhidos))
    if n_val < T.GAIA_MIN_COARSE_MATCHES:
        idx_det_vazio = np.array([], dtype=int)
        idx_gaia_vazio = np.array([], dtype=int)
        return ((0.0, 0.0, n_val, idx_det_vazio, idx_gaia_vazio)
                if retornar_indices else (0.0, 0.0, n_val))

    res_ra = deltas[escolhidos, 0]
    res_dec = deltas[escolhidos, 1]
    idx_det = pares_det[escolhidos]
    idx_gaia = pares_gaia[escolhidos]
    off_ra = float(np.median(res_ra))
    off_dec = float(np.median(res_dec))

    if retornar_indices:
        return off_ra, off_dec, n_val, idx_det, idx_gaia
    return off_ra, off_dec, n_val


def _refine_wcs_gaia(frame: dict, sigma: float = None) -> dict:
    """
    Refina o WCS de um frame contra Gaia DR3 em duas passadas.

    Passada bruta (raio 15 arcsec):
        Estima offset puro de translação — suficiente para cobrir o erro
        sistemático inicial de 5-7 arcsec do WCS Pan-STARRS. Aplica o
        offset ao CRVAL antes da passada fina.

    Passada fina (raio 2 arcsec com WCS pré-corrigido):
        Cross-match preciso + ajuste afim (CD matrix) minimizando resíduos.

    Se a passada fina falhar, o offset bruto é preservado no WCS com status
    "gaia_coarse_offset_only" — melhor do que descartar a correção inteira.

    Status possíveis:
        "gaia_refinado"              — ajuste afim completo aceito
        "gaia_coarse_offset_only"   — só translação; passada fina falhou
        "gaia_skipped_low_pre_rms"     — WCS Pan-STARRS já era bom
        "gaia_network_failed"           — Gaia não respondeu
        "gaia_match_failed"          — estrelas insuficientes no field
        "gaia_high_post_rms"       — ajuste não convergiu para < limiar
        "invalid_wcs"               — frame não tem WCS válido
    """
    if not frame.get("wcs_ok"):
        frame["wcs_status"]          = "invalid_wcs"
        frame["wcs_rms_pre_arcsec"]  = None
        frame["wcs_rms_pos_arcsec"]  = None
        frame["gaia_n_matches"]      = 0
        frame["gaia_coarse_matches"] = 0
        frame["gaia_refined_matches"]  = 0
        frame["gaia_coarse_offset_arcsec"] = [0.0, 0.0]
        return frame

    wcs_original = frame["wcs"]
    if sigma is None:
        sigma = T.DETECTION_SIGMA

    # ── Detectar sources para cross-match ──────────────────────────────────
    data = np.asarray(frame["data"], dtype=np.float64)
    mascara = _invalid_pixel_mask(data)
    try:
        from photutils.background import Background2D, MedianBackground
        from photutils.detection import DAOStarFinder
        from astropy.stats import sigma_clipped_stats
        bkg = Background2D(data, box_size=T.BACKGROUND_BOX,
                           filter_size=T.BACKGROUND_FILTER,
                           bkg_estimator=MedianBackground(),
                           exclude_percentile=25.0, mask=mascara)
        bg  = np.asarray(bkg.background, dtype=np.float64)
        rms = np.asarray(bkg.background_rms, dtype=np.float64)
    except Exception:
        _, med, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
        bg  = np.full(data.shape, float(med))
        rms = np.full(data.shape, max(float(std), 1e-6))

    rms = np.where(np.isfinite(rms) & (rms > 0), rms, np.nanmedian(rms))
    rms = np.where(np.isfinite(rms) & (rms > 0), rms, 1.0)
    det_img = (data - bg) / rms
    det_img[mascara] = 0.0

    finder = DAOStarFinder(fwhm=T.FWHM_PX, threshold=4.0,
                           sharpness_range=(T.SHARPNESS_MIN, T.SHARPNESS_MAX),
                           roundness_range=(-T.ROUNDNESS_MAX, T.ROUNDNESS_MAX),
                           exclude_border=True)
    tabela_det = finder(det_img, mask=mascara)
    if tabela_det is None or len(tabela_det) < T.GAIA_MIN_COARSE_MATCHES:
        frame["wcs_status"]               = "gaia_match_failed"
        frame["wcs_rms_pre_arcsec"]       = None
        frame["wcs_rms_pos_arcsec"]       = None
        frame["gaia_n_matches"]           = 0
        frame["gaia_coarse_matches"]     = 0
        frame["gaia_refined_matches"]      = 0
        frame["gaia_coarse_offset_arcsec"] = [0.0, 0.0]
        return frame

    x_col = "x_centroid" if "x_centroid" in tabela_det.colnames else "xcentroid"
    y_col = "y_centroid" if "y_centroid" in tabela_det.colnames else "ycentroid"
    px_det = np.column_stack([tabela_det[x_col], tabela_det[y_col]])

    try:
        sky = wcs_original.all_pix2world(px_det, 0)
        det_radec = sky.astype(np.float64)
    except Exception:
        frame["wcs_status"]               = "gaia_match_failed"
        frame["wcs_rms_pre_arcsec"]       = None
        frame["wcs_rms_pos_arcsec"]       = None
        frame["gaia_n_matches"]           = 0
        frame["gaia_coarse_matches"]     = 0
        frame["gaia_refined_matches"]      = 0
        frame["gaia_coarse_offset_arcsec"] = [0.0, 0.0]
        return frame

    # ── Consultar Gaia ────────────────────────────────────────────────────
    ra_c  = float(wcs_original.wcs.crval[0])
    dec_c = float(wcs_original.wcs.crval[1])
    tabela_gaia = _query_gaia_dr3(ra_c, dec_c)
    if tabela_gaia is None:
        frame["wcs_status"]               = "gaia_network_failed"
        frame["wcs_rms_pre_arcsec"]       = None
        frame["wcs_rms_pos_arcsec"]       = None
        frame["gaia_n_matches"]           = 0
        frame["gaia_coarse_matches"]     = 0
        frame["gaia_refined_matches"]      = 0
        frame["gaia_coarse_offset_arcsec"] = [0.0, 0.0]
        return frame

    tabela_gaia = _apply_gaia_proper_motion(tabela_gaia, frame["jd"])
    frame["gaia_catalogo_refinado"] = tabela_gaia
    gaia_radec  = np.column_stack([
        np.array(tabela_gaia["ra"],  dtype=np.float64),
        np.array(tabela_gaia["dec"], dtype=np.float64),
    ])

    # ── PASSADA BRUTA: offset grosseiro de translação ─────────────────────
    off_ra, off_dec, n_bruto, idx_det_bruto, idx_gaia_bruto = _estimate_coarse_offset(
        det_radec, gaia_radec, retornar_indices=True
    )
    frame["gaia_coarse_matches"]     = n_bruto
    frame["gaia_coarse_offset_arcsec"] = [round(off_ra, 3), round(off_dec, 3)]

    log(f"  Gaia bruto [{frame['file']}]: N={n_bruto}, "
        f"offset=({off_ra:+.2f}, {off_dec:+.2f}) arcsec", "INFO")

    if n_bruto < T.GAIA_MIN_COARSE_MATCHES:
        frame["wcs_status"]         = "gaia_match_failed"
        frame["wcs_rms_pre_arcsec"] = None
        frame["wcs_rms_pos_arcsec"] = None
        frame["gaia_n_matches"]     = 0
        frame["gaia_refined_matches"] = 0
        return frame

    # Aplicar offset bruto ao CRVAL em duas iterações de refinamento.
    # Física: a mediana da passada bruta tem bias residual de ~0.1–0.2 arcsec
    # porque as detecções no primeiro passo ainda incluem o offset sistemático
    # do WCS Pan-STARRS. Uma segunda iteração com raio de busca médio (5 arcsec)
    # remove esse bias residual e converge para o erro puro de centroide.
    cosdec_bruto = np.cos(np.radians(dec_c))
    wcs_bruto = wcs_original.deepcopy()
    wcs_bruto.wcs.crval = [
        ra_c  - off_ra  / cosdec_bruto / 3600.0,
        dec_c - off_dec / 3600.0,
    ]
    wcs_bruto.wcs.set()

    # Segunda iteração: refinamento do offset com raio intermediário
    try:
        det_rd_iter = wcs_bruto.all_pix2world(px_det, 0).astype(np.float64)
        off_ra2, off_dec2, n_iter2 = _estimate_coarse_offset(
            det_rd_iter, gaia_radec, raio_arcsec=T.GAIA_MATCH_RADIUS_ARCSEC * 3
        )
        if n_iter2 >= T.GAIA_MIN_COARSE_MATCHES and (abs(off_ra2) > 0.05 or abs(off_dec2) > 0.05):
            cosdec2 = np.cos(np.radians(float(wcs_bruto.wcs.crval[1])))
            wcs_bruto.wcs.crval = [
                float(wcs_bruto.wcs.crval[0]) - off_ra2  / cosdec2 / 3600.0,
                float(wcs_bruto.wcs.crval[1]) - off_dec2 / 3600.0,
            ]
            wcs_bruto.wcs.set()
    except Exception:
        pass

    # Re-converter pixels para RA/Dec com WCS bruto corrigido
    try:
        sky_bruto = wcs_bruto.all_pix2world(px_det, 0)
        det_radec_corr = sky_bruto.astype(np.float64)
    except Exception:
        # Bruto aplicado mesmo que a passada fina falhe
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_coarse_offset_only"
        frame["wcs_rms_pre_arcsec"] = None
        frame["wcs_rms_pos_arcsec"] = None
        frame["gaia_n_matches"]   = n_bruto
        frame["gaia_refined_matches"] = 0
        return frame

    # ── PASSADA FINA: cross-match estreito + ajuste afim ─────────────────
    idx_det, idx_gaia = _cross_match_gaia(det_radec_corr, gaia_radec)
    n_fino = len(idx_det)
    frame["gaia_refined_matches"] = n_fino

    log(f"  Gaia fino  [{frame['file']}]: N={n_fino}", "INFO")

    if n_fino >= T.GAIA_MIN_MATCHES:
        idx_det_ajuste = idx_det
        idx_gaia_ajuste = idx_gaia
        match_origem = "fino"
    elif n_bruto >= T.GAIA_MIN_MATCHES:
        # Se a passada fina ainda fica curta, os pares do pico bruto são mais
        # informativos que descartar tudo: eles já representam uma translação
        # coerente do field e permitem estimar/validar a solução afim.
        idx_det_ajuste = idx_det_bruto
        idx_gaia_ajuste = idx_gaia_bruto
        match_origem = "bruto"
        log(f"  Gaia ajuste[{frame['file']}]: usando {len(idx_det_ajuste)} "
            f"pares do pico bruto", "WARN")
    else:
        idx_det_ajuste = np.array([], dtype=int)
        idx_gaia_ajuste = np.array([], dtype=int)
        match_origem = "nenhum"

    # Calcular RMS_pre usando o WCS bruto
    det_match  = det_radec_corr[idx_det_ajuste]
    gaia_match = gaia_radec[idx_gaia_ajuste]
    n_ajuste = len(idx_det_ajuste)
    frame["gaia_match_origem"] = match_origem
    frame["gaia_n_matches_ajuste"] = n_ajuste

    if n_ajuste >= T.GAIA_MIN_MATCHES:
        cosdec_med = np.cos(np.radians(float(np.median(gaia_match[:, 1]))))
        d_ra  = _delta_ra_degrees(det_match[:, 0], gaia_match[:, 0]) * cosdec_med * 3600.0
        d_dec = (det_match[:, 1] - gaia_match[:, 1]) * 3600.0
        rms_pre = float(np.sqrt(np.mean(d_ra**2 + d_dec**2)))
        frame["wcs_rms_pre_arcsec"] = round(rms_pre, 4)
        frame["gaia_n_matches"]     = n_ajuste

        if rms_pre <= T.GAIA_DIFF_MIN_ARCSEC:
            frame["wcs"]              = wcs_bruto
            frame["wcs_status"]       = "gaia_skipped_low_pre_rms"
            frame["wcs_rms_pos_arcsec"] = round(rms_pre, 4)
            frame["gaia_n_matches_usados"] = n_ajuste
            return frame

        # Se o RMS_pre já está dentro do critério de aceite, o WCS bruto
        # (translação pura) é suficiente — o ajuste afim não tem signal de
        # rotação/escala para aprender com erros de centroide aleatórios.
        if rms_pre < T.GAIA_DIFF_MAX_ARCSEC:
            frame["wcs"]              = wcs_bruto
            frame["wcs_status"]       = "gaia_refinado"
            frame["wcs_rms_pos_arcsec"] = round(rms_pre, 4)
            frame["gaia_n_matches_usados"] = n_ajuste
            return frame
    else:
        # Passada fina insuficiente — aceita offset bruto como melhor result
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_coarse_offset_only"
        frame["wcs_rms_pre_arcsec"] = None
        frame["wcs_rms_pos_arcsec"] = None
        frame["gaia_n_matches"]   = n_bruto
        return frame

    # ── Ajuste afim via lstsq com verificação de condicionamento ─────────
    # Física: se os resíduos pós-passada-bruta são uniformes (sem gradiente
    # espacial), o lstsq não tem signal suficiente para estimar rotação/escala
    # e produz coeficientes instáveis que pioram o WCS. Verificamos o número
    # de condição da matriz A; se > 1e6 (mal condicionada), não há informação
    # suficiente para o ajuste afim — mantemos o WCS da passada bruta.
    px_match    = px_det[idx_det_ajuste]
    ra_med_fit  = float(np.mean(gaia_match[:, 0]))
    dec_med_fit = float(np.mean(gaia_match[:, 1]))
    cosdec_fit  = np.cos(np.radians(dec_med_fit))
    target_x    = _delta_ra_degrees(gaia_match[:, 0], ra_med_fit) * cosdec_fit * 3600.0
    target_y    = (gaia_match[:, 1] - dec_med_fit)              * 3600.0

    # Centralizar pixels para melhorar condicionamento numérico
    px_cx = float(np.mean(px_match[:, 0]))
    px_cy = float(np.mean(px_match[:, 1]))
    px_scale = float(np.std(px_match)) if float(np.std(px_match)) > 0 else 1.0
    px_norm = np.column_stack([
        (px_match[:, 0] - px_cx) / px_scale,
        (px_match[:, 1] - px_cy) / px_scale,
        np.ones(len(px_match)),
    ])

    cond = np.linalg.cond(px_norm)
    if cond > 1e6:
        # Matriz mal condicionada — ajuste afim não converge; mantém bruto
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_coarse_offset_only"
        frame["wcs_rms_pos_arcsec"] = None
        frame["gaia_n_matches_usados"] = n_fino
        return frame

    mascara_fit = np.ones(len(px_match), dtype=bool)
    coef_x = coef_y = None
    try:
        for _ in range(4):
            if int(mascara_fit.sum()) < T.GAIA_MIN_MATCHES:
                break
            A = px_norm[mascara_fit]
            coef_x, _, _, _ = np.linalg.lstsq(A, target_x[mascara_fit], rcond=None)
            coef_y, _, _, _ = np.linalg.lstsq(A, target_y[mascara_fit], rcond=None)

            pred_x = px_norm[:, 0]*coef_x[0] + px_norm[:, 1]*coef_x[1] + coef_x[2]
            pred_y = px_norm[:, 0]*coef_y[0] + px_norm[:, 1]*coef_y[1] + coef_y[2]
            resid = np.hypot(pred_x - target_x, pred_y - target_y)
            resid_ok = resid[mascara_fit]
            med = float(np.median(resid_ok))
            mad = float(np.median(np.abs(resid_ok - med)))
            limite = max(T.GAIA_DIFF_MAX_ARCSEC * 2.0, med + 3.0 * 1.4826 * mad)
            nova = resid <= limite
            if int(nova.sum()) < T.GAIA_MIN_MATCHES or np.array_equal(nova, mascara_fit):
                break
            mascara_fit = nova
    except np.linalg.LinAlgError:
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_coarse_offset_only"
        frame["wcs_rms_pos_arcsec"] = None
        return frame

    if coef_x is None or coef_y is None or int(mascara_fit.sum()) < T.GAIA_MIN_MATCHES:
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_coarse_offset_only"
        frame["wcs_rms_pos_arcsec"] = None
        return frame

    frame["gaia_n_matches_usados"] = int(mascara_fit.sum())
    px_match_eval   = px_match[mascara_fit]
    gaia_match_eval = gaia_match[mascara_fit]

    # Reconstruir CD matrix a partir dos coeficientes normalizados.
    # Física: a projeção TAN é não linear para separações grandes —
    # reconstruir o CRVAL analiticamente a partir do ajuste afim introduz
    # erro de curvatura quando CRPIX está longe do centroide dos matches.
    # Solução: manter o CRVAL já refinado pela passada bruta (que é uma
    # translação pura, sem esse problema) e atualizar apenas a CD matrix.
    cd1_1 = float(coef_x[0]) / px_scale / cosdec_fit / 3600.0
    cd1_2 = float(coef_x[1]) / px_scale / cosdec_fit / 3600.0
    cd2_1 = float(coef_y[0]) / px_scale               / 3600.0
    cd2_2 = float(coef_y[1]) / px_scale               / 3600.0

    crpix_x = float(wcs_original.wcs.crpix[0])
    crpix_y = float(wcs_original.wcs.crpix[1])

    # CRVAL: usa o da passada bruta (já corrigido em translação) e apenas
    # ajusta o deslocamento residual medido pelo term constante do lstsq.
    # O term coef_x[2]/coef_y[2] representa o offset residual (em arcsec
    # no sistema local) do centroide dos pixels matched para o centroide Gaia.
    crval_bruto = list(wcs_bruto.wcs.crval)
    new_ra  = crval_bruto[0]
    new_dec = crval_bruto[1]

    wcs_novo = wcs_bruto.deepcopy()
    wcs_novo.wcs.crval = [new_ra, new_dec]
    wcs_novo.wcs.crpix = [crpix_x, crpix_y]
    wcs_novo.wcs.cd    = [[cd1_1, cd1_2], [cd2_1, cd2_2]]
    wcs_novo.wcs.set()

    # ── Avaliar RMS_pos ───────────────────────────────────────────────────
    try:
        sky_novo = wcs_novo.all_pix2world(px_match_eval, 0)
        d_ra_n   = _delta_ra_degrees(sky_novo[:, 0], gaia_match_eval[:, 0]) * cosdec_fit * 3600.0
        d_dec_n  = (sky_novo[:, 1] - gaia_match_eval[:, 1]) * 3600.0
        rms_pos  = float(np.sqrt(np.mean(d_ra_n**2 + d_dec_n**2)))
    except Exception:
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_coarse_offset_only"
        frame["wcs_rms_pos_arcsec"] = None
        return frame

    frame["wcs_rms_pos_arcsec"] = round(rms_pos, 4)
    log(f"  Gaia RMS   [{frame['file']}]: pre={frame['wcs_rms_pre_arcsec']} "
        f"pos={rms_pos:.4f} arcsec", "INFO")

    if rms_pos < T.GAIA_DIFF_MAX_ARCSEC:
        frame["wcs"]        = wcs_novo
        frame["wcs_status"] = "gaia_refinado"
    else:
        frame["wcs"]        = wcs_bruto
        frame["wcs_status"] = "gaia_coarse_offset_only"

    return frame


def _find_science_hdu(hdul: fits.HDUList) -> int:
    for i, h in enumerate(hdul):
        if h.data is not None and h.data.ndim == 2:
            return i
    raise ValueError("No HDU with a 2D image was found in the FITS file.")


def _extract_observation_times(header: fits.Header, file: str) -> dict:
    """Extract start and mid-exposure times, accepting DATE-OBS or MJD-OBS."""
    exptime = float(header.get("EXPTIME", header.get("EXPOSURE", 45.0)))
    date_obs = str(header.get("DATE-OBS", "")).strip()
    mjd_obs = header.get("MJD-OBS")

    try:
        if date_obs:
            t_inicio = Time(date_obs, format="isot", scale="utc")
            date_obs_inicio = date_obs
        elif mjd_obs is not None:
            t_inicio = Time(float(mjd_obs), format="mjd", scale="utc")
            date_obs_inicio = t_inicio.isot
        else:
            raise ValueError("DATE-OBS/MJD-OBS missing")
    except Exception:
        log(
            f"Invalid or missing timestamp in {file}: "
            f"DATE-OBS={date_obs!r}, MJD-OBS={mjd_obs!r}",
            "ERROR",
        )
        sys.exit(1)

    t_mid = t_inicio + (exptime / 2.0) * u.second
    return {
        "date_obs_inicio": date_obs_inicio,
        "date_obs_mid": t_mid.isot,
        "jd_inicio": float(t_inicio.jd),
        "jd_mid": float(t_mid.jd),
        "mjd_inicio": float(t_inicio.mjd),
        "mjd_mid": float(t_mid.mjd),
        "exptime": exptime,
    }


def load_fits(folder: Path, wcs_mode: str = "gaia") -> list[dict]:
    """
    Carrega os 4 files FITS da folder, ordenados por timestamp.

    wcs_mode:
        "gaia"   — refina o WCS contra Gaia DR3 (default).
        "header" — uses only o WCS/header Pan-STARRS construído; desativa
                   completamente a consulta Gaia mantendo todos os demais
                   parâmetros do pipeline idênticos (Background2D, PSF, score,
                   filters e thresholds).
    """
    files = sorted(folder.glob("*.fits")) + sorted(folder.glob("*.fit"))
    if len(files) != 4:
        log(f"Expected exactly 4 FITS files, found {len(files)}.", "ERROR")
        sys.exit(1)

    frames = []
    for file_path in files:
        with fits.open(file_path) as hdul:
            idx = _find_science_hdu(hdul)
            header = hdul[0].header.copy()
            if idx != 0:
                for k, v in hdul[idx].header.items():
                    if k and k not in ("SIMPLE", "BITPIX", "NAXIS",
                                       "EXTEND", "COMMENT", "HISTORY",
                                       "XTENSION", "PCOUNT", "GCOUNT"):
                        header[k] = v
            data = hdul[idx].data.astype(np.float64)

            try:
                wcs = _build_panstarrs_wcs(header)
                wcs_ok = True
            except Exception as e:
                log(f"Invalid WCS in {file_path.name}: {e}", "ERROR")
                wcs = None
                wcs_ok = False

            tempos = _extract_observation_times(header, file_path.name)

            frame_dict = {
                "file" : file_path.name,
                "data"    : data,
                "header"  : header,
                "wcs"     : wcs,
                "wcs_ok"  : wcs_ok,
                "date_obs": tempos["date_obs_mid"],
                "date_obs_inicio": tempos["date_obs_inicio"],
                "date_obs_mid": tempos["date_obs_mid"],
                "jd"      : tempos["jd_mid"],
                "jd_inicio": tempos["jd_inicio"],
                "jd_mid"  : tempos["jd_mid"],
                "mjd_inicio": tempos["mjd_inicio"],
                "mjd_mid" : tempos["mjd_mid"],
                "exptime" : tempos["exptime"],
            }

            if wcs_ok:
                # Measure invalid-pixel fraction before heavily saturated
                # frames contaminate source detection.
                mascara_frame = _invalid_pixel_mask(data)
                sat_pct = float(mascara_frame.mean() * 100.0)
                frame_dict["sat_pct"] = sat_pct
                if sat_pct > T.FRAME_SAT_REJECT_PCT:
                    log(
                        f"Frame {file_path.name} tem {sat_pct:.1f}% de pixels "
                        f"saturados/inválidos (limite: {T.FRAME_SAT_REJECT_PCT}%). "
                        f"Baixe um frame substituto e reexecute o pipeline.",
                        "ERROR",
                    )
                    sys.exit(1)
                elif sat_pct > T.FRAME_SAT_WARNING_PCT:
                    log(
                        f"Frame {file_path.name}: {sat_pct:.1f}% de pixels saturados "
                        f"— qualidade reduzida, mascaramento ativo.",
                        "WARN",
                    )
            else:
                frame_dict["sat_pct"] = 0.0

            frames.append(frame_dict)
        log(f"Carregado: {file_path.name}  ({frame_dict['date_obs']})  "
            f"[{data.shape[1]}x{data.shape[0]}]  WCS={'OK' if wcs_ok else 'FALHOU'}")

    if not any(f["wcs_ok"] for f in frames):
        log("No frame has a valid WCS. Cannot calculate RA/Dec.", "ERROR")
        sys.exit(1)

    frames.sort(key=lambda x: x["jd"])

    if wcs_mode == "gaia":
        # Refinamento astrométrico via Gaia DR3 — etapa separada, após WCS Pan-STARRS.
        # Cada frame é processado independentemente; fallback transparente se Gaia falhar.
        log("Refinando WCS contra Gaia DR3 (wcs_mode=gaia)...")
        for i, f in enumerate(frames):
            if f["wcs_ok"]:
                frames[i] = _refine_wcs_gaia(f)
                status  = frames[i].get("wcs_status", "?")
                rms_pre = frames[i].get("wcs_rms_pre_arcsec")
                rms_pos = frames[i].get("wcs_rms_pos_arcsec")
                n_match = frames[i].get("gaia_n_matches", 0)
                pre_str = f"{rms_pre:.3f}" if rms_pre is not None else "N/A"
                pos_str = f"{rms_pos:.3f}" if rms_pos is not None else "N/A"
                log(f"  WCS {frames[i]['file']}: {status} "
                    f"(RMS pre={pre_str} arcsec, pos={pos_str} arcsec, N={n_match})")
            else:
                frames[i]["wcs_status"]         = "invalid_wcs"
                frames[i]["wcs_rms_pre_arcsec"] = None
                frames[i]["wcs_rms_pos_arcsec"] = None
                frames[i]["gaia_n_matches"]     = 0
    else:
        # Modo header: apenas WCS Pan-STARRS construído, sem Gaia.
        # Registra status "gaia_disabled" em todos os frames para rastreabilidade.
        log("WCS em mode header — Gaia DR3 desativado (wcs_mode=header).", "WARN")
        for i, f in enumerate(frames):
            frames[i]["wcs_status"]               = "gaia_disabled"
            frames[i]["wcs_rms_pre_arcsec"]       = None
            frames[i]["wcs_rms_pos_arcsec"]       = None
            frames[i]["gaia_n_matches"]           = 0
            frames[i]["gaia_coarse_matches"]     = 0
            frames[i]["gaia_refined_matches"]      = 0
            frames[i]["gaia_coarse_offset_arcsec"] = [0.0, 0.0]
            frames[i]["gaia_catalogo_refinado"]   = None

    return frames


# ─────────────────────────────────────────────
# 2. DETECÇÃO DE FONTES EM CADA FRAME
# ─────────────────────────────────────────────
def _invalid_pixel_mask(data: np.ndarray) -> np.ndarray:
    """
    Retorna máscara booleana True para pixels que NÃO devem ser usados
    em detecção, photometry ou estimativa de background.

    Pan-STARRS marca pixels saturados, defeituosos ou off-CCD com 65535
    (limite do uint16) e ocasionalmente com 0 ou valores muito baixos.
    Esses pixels não carregam signal físico válido; passá-los ao
    Background2D infla a RMS local, e ao DAOStarFinder gera detecções
    espúrias que poluem o casamento entre frames.
    """
    sat_min = float(T.SAT_PIXEL_VALUE - T.SAT_PIXEL_TOLERANCE)
    return (
        (data >= sat_min)
        | (data < T.MIN_VALID_PIXEL)
        | ~np.isfinite(data)
    )


def _estimate_local_background(img: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Estima background e RMS locais com `Background2D`.

    Física/astrometria:
        Imagens reais raramente têm céu plano: vinheta, Lua, halos de estrelas
        brilhantes e gradientes do detector mudam o ruído em escalas de dezenas
        de pixels. Um threshold global superestima sources em regiões limpas e
        subestima ruído em bordas. A malha 50x50 mede o céu em escala maior que
        a PSF estelar, enquanto o filtro 3x3 suaviza células ruidosas sem apagar
        gradientes lentos. Assim, o critério n-sigma passa a ser local.
    """
    data = np.asarray(img, dtype=np.float64)
    # Pan-STARRS marca saturação com 65535. Sem essa máscara, o Background2D
    # infla a RMS em torno de detritos saturados (ex.: rastro de satélite),
    # elevando o threshold local e produzindo clusters de "sources" espúrias.
    mascara = _invalid_pixel_mask(data)
    try:
        bkg = Background2D(
            data,
            box_size=T.BACKGROUND_BOX,
            filter_size=T.BACKGROUND_FILTER,
            bkg_estimator=MedianBackground(),
            exclude_percentile=25.0,
            mask=mascara,
        )
        background = np.asarray(bkg.background, dtype=np.float64)
        rms = np.asarray(bkg.background_rms, dtype=np.float64)
    except Exception:
        data_valida = np.where(mascara, np.nan, data)
        _, median, std = sigma_clipped_stats(
            data_valida, sigma=3.0, maxiters=5, mask=mascara
        )
        background = np.full(data.shape, float(median), dtype=np.float64)
        rms = np.full(data.shape, max(float(std), 1e-6), dtype=np.float64)

    rms = np.where(np.isfinite(rms) & (rms > 0), rms, np.nanmedian(rms))
    rms = np.where(np.isfinite(rms) & (rms > 0), rms, 1.0).astype(np.float64)
    return background.astype(np.float64), rms


def _float_cutout(img: np.ndarray, y: float, x: float,
                   raio: int) -> tuple[np.ndarray, int, int]:
    """
    Extrai recorte usando índices inteiros apenas para limites de array.
    A coordenada científica `x,y` permanece float64 e é usada pelo modelo PSF.
    """
    ny, nx = img.shape
    yc = int(np.rint(float(y)))
    xc = int(np.rint(float(x)))
    y0 = max(0, yc - raio); y1 = min(ny, yc + raio + 1)
    x0 = max(0, xc - raio); x1 = min(nx, xc + raio + 1)
    return img[y0:y1, x0:x1].astype(np.float64), y0, x0


def _fwhm_moffat(gamma: float, alpha: float) -> float:
    """FWHM de uma Moffat circular: 2*gamma*sqrt(2^(1/alpha)-1)."""
    if gamma <= 0 or alpha <= 0:
        return float("nan")
    return float(2.0 * gamma * np.sqrt(2.0 ** (1.0 / alpha) - 1.0))


def _fit_subpixel_psf(img_sub: np.ndarray, x0: float, y0: float,
                          fwhm_px: float = 3.0,
                          mascara_invalida: np.ndarray | None = None) -> dict | None:
    """
    Refina uma detecção por ajuste PSF e retorna centro sub-pixel.

    Física/astrometria:
        O centroide simples é uma média ponderada dos pixels e se desloca com
        ruído, background residual, pixels saturados e blends. O ajuste PSF
        compara a source a um modelo contínuo, obtendo o centro do perfil óptico
        em coordenadas sub-pixel. Moffat é preferido em seeing ruim porque suas
        asas de potência representam espalhamento atmosférico melhor que uma
        Gaussiana pura; em fields bem comportados a Gaussiana fica como fallback
        estável.
    """
    cut, y_origin, x_origin = _float_cutout(img_sub, y0, x0, T.PSF_FIT_RADIUS_PX)
    if cut.size < 16 or min(cut.shape) < 5:
        return None

    # If the fitting window contains >30% saturated/invalid pixels, the PSF
    # model would be dominated by artifacts, so skip the fit.
    if mascara_invalida is not None:
        cut_mask, _, _ = _float_cutout(
            mascara_invalida.astype(np.float64), y0, x0, T.PSF_FIT_RADIUS_PX
        )
        if cut_mask.mean() > 0.30:
            return None

    yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
    cx0 = np.float64(x0 - x_origin)
    cy0 = np.float64(y0 - y_origin)
    peak = np.float64(np.nanmax(cut))
    if not np.isfinite(peak) or peak <= 0:
        return None

    fitter = fitting.LevMarLSQFitter()
    fitted = None
    model_name = T.PSF_MODEL.lower()

    try:
        if model_name == "moffat":
            gamma0 = max(float(fwhm_px) / 2.0, 0.8)
            model = Moffat2D(amplitude=peak, x_0=cx0, y_0=cy0,
                             gamma=gamma0, alpha=2.5)
            model.amplitude.bounds = (0.0, None)
            model.x_0.bounds = (cx0 - 3.0, cx0 + 3.0)
            model.y_0.bounds = (cy0 - 3.0, cy0 + 3.0)
            model.gamma.bounds = (0.4, 8.0)
            model.alpha.bounds = (1.1, 8.0)
            fitted = fitter(model, xx, yy, cut, maxiter=T.PSF_MAX_ITER)
            fwhm_fit = _fwhm_moffat(float(fitted.gamma.value),
                                    float(fitted.alpha.value))
        else:
            raise ValueError("usar fallback gaussiano")
    except Exception:
        try:
            sigma0 = max(float(fwhm_px) / 2.355, 0.5)
            model = Gaussian2D(amplitude=peak, x_mean=cx0, y_mean=cy0,
                               x_stddev=sigma0, y_stddev=sigma0)
            model.amplitude.bounds = (0.0, None)
            model.x_mean.bounds = (cx0 - 3.0, cx0 + 3.0)
            model.y_mean.bounds = (cy0 - 3.0, cy0 + 3.0)
            model.x_stddev.bounds = (0.4, 6.0)
            model.y_stddev.bounds = (0.4, 6.0)
            fitted = fitter(model, xx, yy, cut, maxiter=T.PSF_MAX_ITER)
            fwhm_fit = 2.355 * float(np.sqrt(
                fitted.x_stddev.value * fitted.y_stddev.value
            ))
            model_name = "gaussian"
        except Exception:
            return None

    if fitted is None:
        return None

    if isinstance(fitted, Moffat2D):
        x_fit = float(fitted.x_0.value + x_origin)
        y_fit = float(fitted.y_0.value + y_origin)
        amplitude = float(fitted.amplitude.value)
    else:
        x_fit = float(fitted.x_mean.value + x_origin)
        y_fit = float(fitted.y_mean.value + y_origin)
        amplitude = float(fitted.amplitude.value)

    if not (np.isfinite(x_fit) and np.isfinite(y_fit)
            and np.isfinite(fwhm_fit) and np.isfinite(amplitude)):
        return None

    modelo_2d = fitted(xx, yy)
    resid = cut - modelo_2d
    resid_rms = float(np.sqrt(np.mean(resid ** 2)))
    fluxo = float(np.sum(np.clip(modelo_2d, 0.0, None)))

    return {
        "x": np.float64(x_fit),
        "y": np.float64(y_fit),
        "flux": np.float64(fluxo),
        "fwhm": np.float64(fwhm_fit),
        "amplitude": np.float64(amplitude),
        "resid_rms": np.float64(resid_rms),
        "modelo": model_name,
    }


def detect_sources(img: np.ndarray, sigma: float = 5.5,
                    fwhm_px: float = 3.0) -> list[tuple]:
    """
    Detecta sources pontuais com background local + refinamento PSF.

    Retorna lista de (y, x, flux, tamanho), preservando o contrato público.
    As coordenadas vêm do ajuste PSF em float64, sem arredondamento científico.
    """
    data = np.asarray(img, dtype=np.float64)
    # Pan-STARRS marca saturação com 65535. Sem essa máscara, clusters saturados
    # aparecem como picos acima do threshold e entram no casamento entre frames
    # gerando tracks com cinemática absurda que elimina também tracks reais.
    mascara = _invalid_pixel_mask(data)
    background, rms = _estimate_local_background(data)
    img_sub = (data - background).astype(np.float64)
    detection_img = (img_sub / rms).astype(np.float64)
    # Zerar pixels mascarados evita que (65535 - bg) / rms produza picos
    # artificiais nos buracos antes de o DAOStarFinder receber a image.
    detection_img[mascara] = 0.0

    # DAOStarFinder fica como gerador de sementes em uma image já corrigida
    # por background local. A medição final não é o centroide do DAO, e sim o
    # centro do modelo PSF ajustado no recorte da source.
    finder = DAOStarFinder(
        fwhm=fwhm_px,
        threshold=float(sigma),
        sharpness_range=(T.SHARPNESS_MIN, T.SHARPNESS_MAX),
        roundness_range=(-T.ROUNDNESS_MAX, T.ROUNDNESS_MAX),
        exclude_border=True,
    )
    tabela = finder(detection_img, mask=mascara)
    if tabela is None or len(tabela) == 0:
        return []

    x_col = "x_centroid" if "x_centroid" in tabela.colnames else "xcentroid"
    y_col = "y_centroid" if "y_centroid" in tabela.colnames else "ycentroid"

    sources = []
    for row in tabela:
        y0 = np.float64(row[y_col])
        x0 = np.float64(row[x_col])
        psf = _fit_subpixel_psf(img_sub, x0=x0, y0=y0, fwhm_px=fwhm_px,
                                    mascara_invalida=mascara)
        if psf is None:
            continue

        # FWHM físico: raios cósmicos/hot pixels tendem a FWHM sub-pixel,
        # enquanto blends, galáxias e tracks têm FWHM grande demais para a PSF
        # estelar do input_set. Esse filtro protege a astrometria antes do WCS.
        fwhm_fit = float(psf["fwhm"])
        if not (T.FWHM_MIN_PX <= fwhm_fit <= T.FWHM_MAX_PX):
            continue

        y = float(np.float64(psf["y"]))
        x = float(np.float64(psf["x"]))
        flux = float(np.float64(psf["flux"]))
        tam = float(np.pi * (fwhm_fit / 2.0) ** 2)
        sources.append((y, x, flux, tam))
    return sources


def _detect_sources_worker(args: tuple[int, np.ndarray, float, float]) -> tuple[int, list[tuple]]:
    """Worker top-level para multiprocessing com `spawn` no macOS/Apple Silicon."""
    idx, data, sigma, fwhm_px = args
    return idx, detect_sources(data, sigma=sigma, fwhm_px=fwhm_px)


def _detect_sources_frames_parallel(frames: list[dict], sigma: float) -> list[list[tuple]]:
    """
    Processa os 4 frames em paralelo.

    Em Apple Silicon M2, quatro images FITS independentes são uma carga
    naturalmente paralelizável: background local e PSF fitting não compartilham
    estado entre frames. Limitamos a 4 processos para casar com a cadência do
    IASC e evitar overhead de cópia maior que o ganho computacional.
    """
    tarefas = [
        (i, np.asarray(f["data"], dtype=np.float64), float(sigma), float(T.FWHM_PX))
        for i, f in enumerate(frames)
    ]
    n_proc = min(4, len(tarefas), os.cpu_count() or 1)
    if n_proc <= 1:
        return [detect_sources(f["data"], sigma=sigma) for f in frames]

    try:
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=n_proc) as pool:
            results = pool.map(_detect_sources_worker, tarefas)
        results.sort(key=lambda item: item[0])
        return [sources for _, sources in results]
    except Exception as e:
        log(f"Multiprocessing indisponível; usando processamento serial: {e}", "WARN")
        return [detect_sources(f["data"], sigma=sigma) for f in frames]


# ─────────────────────────────────────────────
# 3. IDENTIFICAÇÃO DE OBJETOS EM MOVIMENTO
# ─────────────────────────────────────────────
def _mutual_match(fontes_a: list[tuple], fontes_b: list[tuple],
                 radius_px: float = None) -> dict[int, int]:
    """
    Casamento mútuo de vizinho mais próximo entre duas listas de sources.
    Retorna {idx_a: idx_b} onde o casamento é recíproco e dentro de radius_px.
    """
    if radius_px is None:
        radius_px = T.MATCH_RADIUS_PX
    if not fontes_a or not fontes_b:
        return {}

    pa = np.array([(f[0], f[1]) for f in fontes_a])
    pb = np.array([(f[0], f[1]) for f in fontes_b])

    d2 = ((pa[:, None, :] - pb[None, :, :]) ** 2).sum(axis=2)

    melhor_a_para_b = np.argmin(d2, axis=1)
    melhor_b_para_a = np.argmin(d2, axis=0)

    pares = {}
    for i, j in enumerate(melhor_a_para_b):
        if melhor_b_para_a[j] == i and d2[i, j] < radius_px ** 2:
            pares[i] = int(j)
    return pares


def _pixel_scale_arcsec(wcs: WCS | None) -> float:
    """Escala média do pixel em arcsec/pixel a partir do WCS."""
    if wcs is None:
        return 0.26
    try:
        escalas = proj_plane_pixel_scales(wcs) * 3600.0
        escalas = np.asarray(escalas, dtype=np.float64)
        escalas = escalas[np.isfinite(escalas) & (escalas > 0)]
        if escalas.size:
            return float(np.median(escalas))
    except Exception:
        pass
    return 0.26


def _sources_to_sky_plane(sources: list[tuple], frame: dict,
                               ra_ref: float, dec_ref: float) -> np.ndarray:
    """
    Projeta sources detectadas para um plano tangente local em arcsec.
    Usa o WCS refinado do frame; se o WCS falhar, retorna coordenadas em pixel.
    """
    if not sources:
        return np.empty((0, 2), dtype=np.float64)

    pix = np.array([(f[1], f[0]) for f in sources], dtype=np.float64)
    if not frame.get("wcs_ok") or frame.get("wcs") is None:
        return np.array([(f[0], f[1]) for f in sources], dtype=np.float64)

    try:
        radec = frame["wcs"].all_pix2world(pix, 0).astype(np.float64)
        cosdec = np.cos(np.radians(float(dec_ref)))
        return np.column_stack([
            _delta_ra_degrees(radec[:, 0], ra_ref) * cosdec * 3600.0,
            (radec[:, 1] - dec_ref) * 3600.0,
        ])
    except Exception:
        return np.array([(f[0], f[1]) for f in sources], dtype=np.float64)


def _circular_ra_mean_degrees(ras: list[float]) -> float:
    """Média circular de RA em graus, correta em fields cruzando 0/360."""
    ang = np.radians(np.asarray(ras, dtype=np.float64))
    ra = np.degrees(np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang))))
    return float(ra % 360.0)


def _mutual_sky_match(fontes_a: list[tuple], fontes_b: list[tuple],
                         frame_a: dict, frame_b: dict,
                         radius_px: float = None) -> dict[int, int]:
    """
    Casamento mútuo em coordenadas celestes projetadas pelo WCS refinado.
    O raio operacional continua definido em pixels, convertido para arcsec.
    """
    if radius_px is None:
        radius_px = T.MATCH_RADIUS_PX
    if not fontes_a or not fontes_b:
        return {}

    try:
        ra_ref = _circular_ra_mean_degrees([
            float(frame_a["wcs"].wcs.crval[0]),
            float(frame_b["wcs"].wcs.crval[0]),
        ])
        dec_ref = float(np.mean([
            float(frame_a["wcs"].wcs.crval[1]),
            float(frame_b["wcs"].wcs.crval[1]),
        ]))
        escala = float(np.mean([
            _pixel_scale_arcsec(frame_a.get("wcs")),
            _pixel_scale_arcsec(frame_b.get("wcs")),
        ]))
        pa = _sources_to_sky_plane(fontes_a, frame_a, ra_ref, dec_ref)
        pb = _sources_to_sky_plane(fontes_b, frame_b, ra_ref, dec_ref)
        raio = float(radius_px) * escala
    except Exception:
        return _mutual_match(fontes_a, fontes_b, radius_px=radius_px)

    if len(pa) == 0 or len(pb) == 0:
        return {}

    d2 = ((pa[:, None, :] - pb[None, :, :]) ** 2).sum(axis=2)
    melhor_a_para_b = np.argmin(d2, axis=1)
    melhor_b_para_a = np.argmin(d2, axis=0)

    pares = {}
    for i, j in enumerate(melhor_a_para_b):
        if melhor_b_para_a[j] == i and d2[i, j] < raio ** 2:
            pares[i] = int(j)
    return pares



def find_movers(frames: list[dict], sigma: float = 5.5) -> tuple[list[dict], dict]:
    """
    Compara sources entre os 4 frames e identifica objects com motion
    residual após remover a deriva do field.

    Com WCS refinado por Gaia, o casamento mútuo opera sobre coordenadas
    celestes projetadas em plano tangente local, não apenas em pixels brutos.

    Retorna (movers, global_metrics).
    """
    todas_fontes = _detect_sources_frames_parallel(frames, sigma=sigma)
    n_fontes_by_frame = []
    for i, f in enumerate(frames):
        sources = todas_fontes[i]
        n_fontes_by_frame.append(len(sources))
        log(f"  Frame {i+1} ({f['date_obs'][:19]}): {len(sources)} sources")

    pares_12 = _mutual_sky_match(todas_fontes[0], todas_fontes[1],
                                    frames[0], frames[1])
    pares_23 = _mutual_sky_match(todas_fontes[1], todas_fontes[2],
                                    frames[1], frames[2])
    pares_34 = _mutual_sky_match(todas_fontes[2], todas_fontes[3],
                                    frames[2], frames[3])

    raw_candidates = []
    n_tracks_attempted = 0

    for i0, i1 in pares_12.items():
        if i1 not in pares_23:
            continue
        i2 = pares_23[i1]
        if i2 not in pares_34:
            continue
        i3 = pares_34[i2]
        n_tracks_attempted += 1

        f0 = todas_fontes[0][i0]
        f1 = todas_fontes[1][i1]
        f2 = todas_fontes[2][i2]
        f3 = todas_fontes[3][i3]
        track = [(f0[0], f0[1], f0[2]),
                  (f1[0], f1[1], f1[2]),
                  (f2[0], f2[1], f2[2]),
                  (f3[0], f3[1], f3[2])]

        pos  = np.array([(t[0], t[1]) for t in track], dtype=np.float64)
        move = float(np.hypot(pos[-1, 0] - pos[0, 0], pos[-1, 1] - pos[0, 1]))
        raw_candidates.append({
            "track"    : track,
            "move_total": move,
            "fluxes"   : [t[2] for t in track],
        })

    deriva_dx, deriva_dy = 0.0, 0.0
    n_estaveis = 0
    if raw_candidates:
        vetores = np.array([
            (c["track"][-1][1] - c["track"][0][1],
             c["track"][-1][0] - c["track"][0][0])
            for c in raw_candidates
        ])
        deriva_dx = float(np.median(vetores[:, 0]))
        deriva_dy = float(np.median(vetores[:, 1]))
        res_inicial = np.hypot(vetores[:, 0] - deriva_dx,
                               vetores[:, 1] - deriva_dy)
        mascara_estaveis = res_inicial < T.STABLE_DRIFT_PX
        if mascara_estaveis.sum() >= T.MIN_STABLE_DRIFT_SOURCES:
            deriva_dx = float(np.median(vetores[mascara_estaveis, 0]))
            deriva_dy = float(np.median(vetores[mascara_estaveis, 1]))
        n_estaveis = int(mascara_estaveis.sum())
        log(f"Deriva instrumental estimada (dx, dy): "
            f"({deriva_dx:+.2f}, {deriva_dy:+.2f}) px  "
            f"[{n_estaveis} estrelas usadas]")

        movers = []
        for c, v in zip(raw_candidates, vetores):
            res = float(np.hypot(v[0] - deriva_dx, v[1] - deriva_dy))
            if T.MOVE_MIN_PX < c["move_total"] < T.MOVE_MAX_PX and res > T.MIN_RESIDUAL_PX:
                c["residual"] = res
                movers.append(c)
    else:
        movers = []

    metrics = {
        "n_fontes_by_frame"        : n_fontes_by_frame,
        "n_tracks_attempted"        : n_tracks_attempted,
        "drift_dx_px"              : round(deriva_dx, 3),
        "drift_dy_px"              : round(deriva_dy, 3),
        "n_estaveis_deriva"         : n_estaveis,
        "sigma"            : sigma,
        "background_model"          : "photutils.Background2D",
        "background_box_size"       : list(T.BACKGROUND_BOX),
        "background_filter_size"    : list(T.BACKGROUND_FILTER),
        "psf_model"                 : T.PSF_MODEL,
        "psf_fwhm_range_px"         : [T.FWHM_MIN_PX, T.FWHM_MAX_PX],
        "multiprocessing_processes" : min(4, len(frames), os.cpu_count() or 1),
        "saturacao_pct_by_frame"   : [round(f.get("sat_pct", 0.0), 2) for f in frames],
        "gaia_refinement"          : [
            {
                "frame"               : f["file"],
                "status"              : f.get("wcs_status", "invalid_wcs"),
                "rms_pre_arcsec"      : f.get("wcs_rms_pre_arcsec"),
                "rms_pos_arcsec"      : f.get("wcs_rms_pos_arcsec"),
                "n_matches"           : f.get("gaia_n_matches", 0),
                "coarse_matches"     : f.get("gaia_coarse_matches", 0),
                "coarse_offset_arcsec" : f.get("gaia_coarse_offset_arcsec", [0.0, 0.0]),
                "n_matches_fino"      : f.get("gaia_refined_matches", 0),
                "n_matches_ajuste"    : f.get("gaia_n_matches_ajuste", 0),
                "n_matches_usados"    : f.get("gaia_n_matches_usados", 0),
                "match_origem"        : f.get("gaia_match_origem"),
            }
            for f in frames
        ],
    }

    return movers, metrics


# ─────────────────────────────────────────────
# 4. ANÁLISE E SCORE DE CADA CANDIDATO
# ─────────────────────────────────────────────
def _is_near_edge(tx: float, ty: float, shape: tuple,
                   margem_px: float = None) -> bool:
    """Retorna True se o centroide está a menos de margem_px da borda."""
    if margem_px is None:
        margem_px = T.EDGE_MARGIN_PX
    ny, nx = shape
    return tx < margem_px or tx > (nx - margem_px) \
        or ty < margem_px or ty > (ny - margem_px)


def _frame_morphology(img: np.ndarray, ty: float, tx: float,
                          r: int = 20) -> dict:
    """
    Calcula métricas morfológicas simples a partir do recorte local.

    Retorna:
        pointlike     : razão pico/total (0–1; alta = ponto, baixa = estendida)
        elongation  : razão eixo maior/menor via PCA dos pixels acima do limiar
                      (1.0 = circular; >1.6 = estendida/artefato de trail)
        fwhm_est_px : estimativa de FWHM via limiar de meia potência [px]
        valid      : False se recorte degenerado (borda, fluxo nulo)
    """
    ny, nx = img.shape
    # Inteiros são usados somente para delimitar o recorte em memória. O centro
    # astrométrico `tx,ty` continua float64 em todos os cálculos de distância.
    yc = int(np.rint(float(ty)))
    xc = int(np.rint(float(tx)))
    y0 = max(0, yc - r); y1 = min(ny, yc + r + 1)
    x0 = max(0, xc - r); x1 = min(nx, xc + r + 1)
    cut = img[y0:y1, x0:x1].astype(np.float64)

    if cut.size < 4:
        return {"pointlike": 0.0, "elongation": 1.0, "fwhm_est_px": None, "valid": False}

    bg = float(np.percentile(cut, 30))
    cut_sub = np.clip(cut - bg, 0.0, None)
    total = float(cut_sub.sum())

    if total <= 0.0:
        return {"pointlike": 0.0, "elongation": 1.0, "fwhm_est_px": None, "valid": False}

    pointlike = float(cut_sub.max() / total)

    # Elongation via PCA sobre pixels acima de 10% do pico
    peak = float(cut_sub.max())
    limiar = 0.10 * peak
    yy, xx = np.where(cut_sub > limiar)
    elongation = 1.0
    if len(yy) >= 5:
        weights = cut_sub[yy, xx]
        cy = float(np.average(yy, weights=weights))
        cx = float(np.average(xx, weights=weights))
        dy = yy - cy; dx = xx - cx
        # Tensor de inércia ponderado
        Ixx = float(np.sum(weights * dx * dx)) / weights.sum()
        Iyy = float(np.sum(weights * dy * dy)) / weights.sum()
        Ixy = float(np.sum(weights * dx * dy)) / weights.sum()
        disc = np.sqrt(max(0.0, ((Ixx - Iyy) / 2) ** 2 + Ixy ** 2))
        lam1 = (Ixx + Iyy) / 2 + disc
        lam2 = (Ixx + Iyy) / 2 - disc
        if lam2 > 0.01:
            elongation = float(np.sqrt(lam1 / lam2))
        else:
            elongation = float(np.sqrt(lam1 / 0.01))

    # FWHM estimada: raio onde a intensidade cai à metade do pico
    half = peak / 2.0
    cy_cut = float(ty) - y0; cx_cut = float(tx) - x0
    yg, xg = np.ogrid[:cut_sub.shape[0], :cut_sub.shape[1]]
    dist_map = np.sqrt((yg - cy_cut) ** 2 + (xg - cx_cut) ** 2)
    mask_acima = cut_sub >= half
    fwhm_est = None
    if mask_acima.any():
        raios_acima = dist_map[mask_acima]
        fwhm_est = float(2.0 * raios_acima.max())

    psf_fit = _fit_subpixel_psf(cut_sub, x0=cx_cut, y0=cy_cut,
                                    fwhm_px=T.FWHM_PX)
    fwhm_psf = float(psf_fit["fwhm"]) if psf_fit is not None else None
    psf_model = psf_fit["modelo"] if psf_fit is not None else None

    return {
        "pointlike"    : pointlike,
        "elongation" : elongation,
        "fwhm_est_px": fwhm_psf if fwhm_psf is not None else fwhm_est,
        "psf_model" : psf_model,
        "valid"     : True,
    }


def _check_false_positive(c: dict, frames: list[dict],
                               snrs: list[float],
                               morphology_metrics: list[dict],
                               fluxes: list[float]) -> list[str]:
    """
    Avalia heurísticas de falso positivo e retorna lista de razões de rejeição.
    Lista vazia significa candidate sem evidência de ser artefato.

    Heurísticas:
    - Borda do frame (qualquer posição na track)
    - SNR insuficiente em mais de 1 frame
    - Hot-pixel profile (extreme pointness + weak SNR)
    - Brilho completamente caótico (CV >> limiar)
    - Morfologia estendida grave e inconsistente (elongation muito alta)
    """
    reasons = []
    track   = c["track"]
    image_shape = frames[0]["data"].shape

    # Borda: rejeição (não apenas flag)
    frames_na_borda = [
        fi for fi in range(4)
        if _is_near_edge(track[fi][1], track[fi][0], image_shape)
    ]
    if frames_na_borda:
        reasons.append(
            f"posição na borda em frame(s) {[f+1 for f in frames_na_borda]} "
            f"(margem {T.EDGE_MARGIN_PX:.0f} px)"
        )

    # SNR insuficiente
    snrs_validos = [s for s in snrs if s is not None]
    n_snr_baixo = sum(1 for s in snrs_validos if s < T.MIN_SNR)
    if n_snr_baixo >= 2:
        reasons.append(
            f"SNR insuficiente em {n_snr_baixo}/4 frames "
            f"(mínimo {T.MIN_SNR:.1f}; valores={[round(s,1) for s in snrs_validos]})"
        )

    # Hot pixel: pointness muito alta + SNR borderline
    pointlike_by_frame = [m["pointlike"] for m in morphology_metrics if m["valid"]]
    pont_media = float(np.mean(pointlike_by_frame)) if pointlike_by_frame else 0.0
    snr_medio  = float(np.mean(snrs_validos)) if snrs_validos else 0.0
    if pont_media > T.HOT_PIXEL_POINTNESS_MAX and snr_medio < T.HOT_PIXEL_SNR_MIN:
        reasons.append(
            f"perfil suspeito de hot pixel: pointness={pont_media:.3f} "
            f"(>{T.HOT_PIXEL_POINTNESS_MAX}) com SNR médio={snr_medio:.1f} "
            f"(<{T.HOT_PIXEL_SNR_MIN})"
        )

    # Brilho caótico
    flux_array = np.array(fluxes, dtype=np.float64)
    if flux_array.mean() > 0:
        cv = float(flux_array.std() / flux_array.mean())
        if cv > T.CHAOTIC_FLUX_CV:
            reasons.append(
                f"flux caótico: CV={cv:.2f} (>{T.CHAOTIC_FLUX_CV}) — "
                "inconsistente com objeto real"
            )

    # Elongation grave (artefato de trail CCD, raio cósmico, etc.)
    elongs = [m["elongation"] for m in morphology_metrics if m["valid"]]
    if elongs:
        elong_max = float(max(elongs))
        elong_med = float(np.mean(elongs))
        # Só rejeita se TODOS os frames mostram elongation grave
        if elong_med > 3.0 and elong_max > 4.0:
            reasons.append(
                f"morphology sistematicamente estendida: "
                f"elongation média={elong_med:.2f}, máx={elong_max:.2f} "
                "(possível artefato, trail CCD ou galáxia de background)"
            )

    # FWHM incompatível com a PSF estelar: a posição MPC precisa vir do centro
    # de uma source pointlike. FWHM sub-pixel é assinatura clássica de raio cósmico
    # ou hot pixel; FWHM muito largo desloca o centro por blend/galáxia/trail.
    fwhms = [float(m["fwhm_est_px"]) for m in morphology_metrics
             if m["valid"] and m.get("fwhm_est_px") is not None
             and np.isfinite(m["fwhm_est_px"])]
    if fwhms:
        fwhm_med = float(np.median(fwhms))
        if fwhm_med < T.FWHM_MIN_PX or fwhm_med > T.FWHM_MAX_PX:
            reasons.append(
                f"FWHM incompatível com PSF estelar: mediana={fwhm_med:.2f} px "
                f"(faixa aceita {T.FWHM_MIN_PX:.1f}-{T.FWHM_MAX_PX:.1f} px)"
            )

    return reasons


def _calculate_local_snr(img: np.ndarray, ty: float, tx: float,
                        signal_radius: int = 4, raio_bg: int = 12) -> float:
    """
    SNR local simples: signal = soma de pixels no signal_radius após subtrair
    background estimado em annulus entre signal_radius e raio_bg.
    """
    ny, nx = img.shape
    yc = int(np.rint(float(ty)))
    xc = int(np.rint(float(tx)))
    y0s = max(0, yc - signal_radius); y1s = min(ny, yc + signal_radius + 1)
    x0s = max(0, xc - signal_radius); x1s = min(nx, xc + signal_radius + 1)
    corte_sinal = img[y0s:y1s, x0s:x1s]

    y0b = max(0, yc - raio_bg); y1b = min(ny, yc + raio_bg + 1)
    x0b = max(0, xc - raio_bg); x1b = min(nx, xc + raio_bg + 1)
    corte_bg = img[y0b:y1b, x0b:x1b]

    bg_med  = float(np.median(corte_bg))
    bg_std  = float(np.std(corte_bg))
    signal   = float(np.sum(corte_sinal - bg_med))
    n_pix   = corte_sinal.size
    noise   = bg_std * np.sqrt(n_pix) if n_pix > 0 else 1.0
    return float(signal / noise) if noise > 0 else 0.0


def _calculate_astrometric_uncertainties(
    frames: list[dict],
    morphology_metrics: list[dict],
    snrs: list[float],
) -> dict:
    """Estimativa por frame: erro centroidal FWHM/SNR combinado ao RMS Gaia."""
    by_frame = []
    sigma_total = []
    for fi, frame in enumerate(frames):
        snr = float(snrs[fi]) if fi < len(snrs) and snrs[fi] is not None else 0.0
        morph = morphology_metrics[fi] if fi < len(morphology_metrics) else {}
        fwhm_px = morph.get("fwhm_est_px") or T.FWHM_PX
        escala = _pixel_scale_arcsec(frame.get("wcs"))
        if escala is None:
            escala = 1.0

        sigma_centroid = None
        if snr > 0 and np.isfinite(snr):
            sigma_centroid = float(escala * float(fwhm_px) / (2.355 * snr))

        sigma_wcs = frame.get("wcs_rms_pos_arcsec")
        sigma_wcs = float(sigma_wcs) if sigma_wcs is not None else None

        componentes = [v for v in (sigma_centroid, sigma_wcs)
                       if v is not None and np.isfinite(v)]
        sigma_pos = float(np.sqrt(np.sum(np.square(componentes)))) if componentes else None
        if sigma_pos is not None:
            sigma_total.append(sigma_pos)

        by_frame.append({
            "frame_index": fi,
            "sigma_centroid_arcsec": round(sigma_centroid, 4)
                                      if sigma_centroid is not None else None,
            "sigma_wcs_arcsec": round(sigma_wcs, 4) if sigma_wcs is not None else None,
            "sigma_pos_arcsec": round(sigma_pos, 4) if sigma_pos is not None else None,
            "metodo": "FWHM/(2.355*SNR) combinado em quadratura com RMS Gaia",
        })

    return {
        "by_frame": by_frame,
        "sigma_pos_mediana_arcsec": round(float(np.median(sigma_total)), 4)
                                    if sigma_total else None,
    }


def _evaluate_gaia_static_source(c: dict, frames: list[dict]) -> dict:
    """Marca tracks compatíveis com a mesma source Gaia em vários frames."""
    matches = []
    source_ids = []
    for fi, frame in enumerate(frames):
        tabela = frame.get("gaia_catalogo_refinado")
        if tabela is None or not frame.get("wcs_ok") or len(tabela) == 0:
            continue

        try:
            ra, dec = pixel_to_radec(c["track"][fi][1], c["track"][fi][0], frame["wcs"])
            gaia_radec = np.column_stack([
                np.array(tabela["ra"], dtype=np.float64),
                np.array(tabela["dec"], dtype=np.float64),
            ])
            idx_det, idx_gaia = _cross_match_gaia(
                np.array([[ra, dec]], dtype=np.float64),
                gaia_radec,
                raio_arcsec=T.GAIA_STATIC_MATCH_ARCSEC,
            )
            if len(idx_det) == 0:
                continue

            gi = int(idx_gaia[0])
            sid = str(tabela["source_id"][gi]) if "source_id" in tabela.colnames else str(gi)
            source_ids.append(sid)
            coord_c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
            coord_g = SkyCoord(ra=float(tabela["ra"][gi]) * u.deg,
                               dec=float(tabela["dec"][gi]) * u.deg)
            matches.append({
                "frame_index": fi,
                "source_id": sid,
                "dist_arcsec": round(float(coord_c.separation(coord_g).to(u.arcsec).value), 4),
                "g_mag": round(float(tabela["phot_g_mean_mag"][gi]), 3)
                         if "phot_g_mean_mag" in tabela.colnames else None,
            })
        except Exception:
            continue

    mesma_fonte = False
    if source_ids:
        _, contagens = np.unique(source_ids, return_counts=True)
        mesma_fonte = int(np.max(contagens)) >= T.GAIA_STATIC_MIN_FRAMES

    return {
        "status": "fonte_estatica_gaia" if mesma_fonte else "no_match_estatico",
        "n_matches": len(matches),
        "matches": matches,
    }


def analyze_candidate(c: dict, frames: list[dict]) -> dict:
    """
    Calcula métricas, score decomposto, morphology robusta, flags diagnósticas
    e razões de decisão para um candidate.

    Score (0–12):
        linearity        : 0–3  (resíduo da regressão linear na track)
        velocity         : 0–2  (consistência entre passos consecutivos)
        photometry         : 0–2  (coeficiente de variação do flux)
        morphology         : 0–2  (pointness média do perfil)
        consistencia_morfo : 0–1  (estabilidade morfológica entre frames)
        elongation         : 0–1  (compacidade; penaliza sources estendidas)
        velocity_range          : 0–1  (velocity total na faixa típica de MBA)

    Candidatos com razões de rejeição recebem triage_class=DISCARDED independente do score.
    """
    track    = c["track"]
    pos       = np.array([(t[1], t[0]) for t in track], dtype=np.float64)  # (x, y)
    t_idx     = np.arange(4, dtype=np.float64)
    image_shape = frames[0]["data"].shape

    penalty_reasons: list[str] = []   # penaliza score mas não rejeita
    rejection_reasons:    list[str] = []   # rejeição direta → DISCARDED

    # ── Linearidade ──────────────────────────────────────────────────────────
    px = np.polyfit(t_idx, pos[:, 0], 1)
    py = np.polyfit(t_idx, pos[:, 1], 1)
    res_x = pos[:, 0] - np.polyval(px, t_idx)
    res_y = pos[:, 1] - np.polyval(py, t_idx)
    lin = float(np.mean(np.sqrt(res_x ** 2 + res_y ** 2)))
    ss_res = float(np.sum(res_x ** 2 + res_y ** 2))
    pos_centered = pos - np.mean(pos, axis=0)
    ss_tot = float(np.sum(pos_centered[:, 0] ** 2 + pos_centered[:, 1] ** 2))
    linear_r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else None

    # ── Velocidade / consistência de passos ──────────────────────────────────
    passos_x = np.diff(pos[:, 0])
    passos_y = np.diff(pos[:, 1])
    vel_cons = float(np.std(passos_x) + np.std(passos_y))

    # Velocidade angular em arcsec/min
    vel_arcsec_min = None
    if frames[0]["wcs_ok"] and frames[-1]["wcs_ok"]:
        try:
            ra0, dec0 = frames[0]["wcs"].all_pix2world(
                track[0][1], track[0][0], 0)
            ra3, dec3 = frames[-1]["wcs"].all_pix2world(
                track[-1][1], track[-1][0], 0)
            c0 = SkyCoord(ra=float(ra0) * u.deg, dec=float(dec0) * u.deg)
            c3 = SkyCoord(ra=float(ra3) * u.deg, dec=float(dec3) * u.deg)
            sep_arcsec = float(c0.separation(c3).to(u.arcsec).value)
            dt_min = (frames[-1]["jd"] - frames[0]["jd"]) * 1440.0
            if dt_min > 0:
                vel_arcsec_min = round(sep_arcsec / dt_min, 4)
        except Exception:
            pass

    # ── Fotometria e morphology robusta por frame ─────────────────────────────
    fluxes       : list[float] = []
    snrs          : list[float] = []
    morphology_metrics: list[dict]  = []

    for fi, frame in enumerate(frames):
        ty, tx = track[fi][0], track[fi][1]
        morph = _frame_morphology(frame["data"], ty, tx, r=20)
        morphology_metrics.append(morph)
        fluxes.append(morph["pointlike"] * 1.0)   # placeholder; recalcular abaixo

        # Brilho: soma do recorte subtraído de background
        r = 20
        yc = int(np.rint(float(ty)))
        xc = int(np.rint(float(tx)))
        y0_ = max(0, yc - r); y1_ = min(frame["data"].shape[0], yc + r + 1)
        x0_ = max(0, xc - r); x1_ = min(frame["data"].shape[1], xc + r + 1)
        cut = frame["data"][y0_:y1_, x0_:x1_].astype(np.float64)
        bg  = np.percentile(cut, 30)
        fluxes[-1] = float(np.sum(np.clip(cut - bg, 0, None)))

        snrs.append(_calculate_local_snr(frame["data"], ty, tx))

    flux_array = np.array(fluxes, dtype=np.float64)
    flux_cv   = float(np.std(flux_array) / (np.mean(flux_array) + 1e-10)) \
                  if np.all(flux_array > 0) else float("inf")

    pointlike_by_frame = [m["pointlike"] for m in morphology_metrics if m["valid"]]
    pont           = float(np.mean(pointlike_by_frame)) if pointlike_by_frame else 0.0
    pont_std       = float(np.std(pointlike_by_frame))  if len(pointlike_by_frame) > 1 else 0.0

    elongations_by_frame   = [m["elongation"] for m in morphology_metrics if m["valid"]]
    mean_elongation    = float(np.mean(elongations_by_frame)) if elongations_by_frame else 1.0

    fwhms_by_frame    = [m["fwhm_est_px"] for m in morphology_metrics
                      if m["valid"] and m["fwhm_est_px"] is not None]
    fwhm_medio     = float(np.mean(fwhms_by_frame)) if fwhms_by_frame else None

    # ── Verificação de falso positivo ────────────────────────────────────────
    rejection_reasons = _check_false_positive(
        c, frames, snrs, morphology_metrics, fluxes
    )
    astrometric_uncertainty = _calculate_astrometric_uncertainties(
        frames, morphology_metrics, snrs
    )
    gaia_static = _evaluate_gaia_static_source(c, frames)
    if gaia_static["status"] == "fonte_estatica_gaia":
        rejection_reasons.append(
            f"compatível com source estática Gaia em {gaia_static['n_matches']}/4 frames"
        )

    # ── Score por componente ─────────────────────────────────────────────────
    move = c["move_total"]

    if   lin < T.LIN_EXCELENTE:  s_lin = 3
    elif lin < T.LIN_GOOD:        s_lin = 2
    elif lin < T.LIN_MARGINAL:   s_lin = 1
    else:
        s_lin = 0
        penalty_reasons.append(f"linearity ruim (resíduo={lin:.2f} px)")

    if   vel_cons < T.VEL_UNIFORM:  s_velocity = 2
    elif vel_cons < T.VEL_MODERATE:  s_velocity = 1
    else:
        s_velocity = 0
        penalty_reasons.append(f"velocity irregular (σ={vel_cons:.2f} px)")

    if   flux_cv < T.PHOT_STABLE:  s_phot = 2
    elif flux_cv < T.PHOT_MODERATE: s_phot = 1
    else:
        s_phot = 0
        penalty_reasons.append(f"flux instável (CV={flux_cv:.2f})")

    if   pont > T.MOR_POINTLIKE:   s_morph = 2
    elif pont > T.MOR_MARGINAL:  s_morph = 1
    else:
        s_morph = 0
        penalty_reasons.append(f"morphology estendida (pointness={pont:.4f})")

    # Consistência morfológica entre frames (novo)
    s_morph_consistency = 1 if pont_std < T.MORPH_CONSIST_MAX_STD else 0
    if s_morph_consistency == 0:
        penalty_reasons.append(
            f"morphology inconsistente entre frames (std pointness={pont_std:.3f})"
        )

    # Elongation (novo)
    s_elongation = 1 if mean_elongation < T.ELON_GOOD else 0
    if s_elongation == 0:
        penalty_reasons.append(
            f"elongation elevada (média={mean_elongation:.2f}, limiar={T.ELON_GOOD})"
        )

    s_range = 1 if T.VEL_RANGE_MIN_PX < move < T.VEL_RANGE_MAX_PX else 0

    raw_score = s_lin + s_velocity + s_phot + s_morph + s_morph_consistency + s_elongation + s_range
    score_max   = 3 + 2 + 2 + 2 + 1 + 1 + 1   # = 12

    # Candidatos com razão de rejeição → DISCARDED independente de score
    if rejection_reasons:
        score   = 0
        triage_class  = "DISCARDED"
        cli_color = Color.RED
    else:
        score = raw_score
        if   score >= T.SCORE_STRONG:    triage_class, cli_color = "STRONG",    Color.GREEN
        elif score >= T.SCORE_MODERATE: triage_class, cli_color = "MODERATE", Color.YELLOW
        elif score >= T.SCORE_WEAK:    triage_class, cli_color = "WEAK",    Color.YELLOW
        else:                           triage_class, cli_color = "DISCARDED", Color.RED

    # Normalized heuristic score for display only; it is not calibrated confidence.
    score_percent = min(95, max(5, int(score / score_max * 100)))

    # ── Flags diagnósticas ───────────────────────────────────────────────────
    flags: list[str] = []

    if s_lin >= 2:   flags.append("LINEARIDADE_BOA")
    elif s_lin == 0: flags.append("LINEARIDADE_RUIM")

    if s_velocity >= 2:   flags.append("VELOCIDADE_CONSISTENTE")
    elif s_velocity == 0: flags.append("VELOCIDADE_IRREGULAR")

    if s_phot >= 2:   flags.append("BRILHO_ESTAVEL")
    elif s_phot == 0: flags.append("BRILHO_INSTAVEL")

    if s_morph >= 2:   flags.append("MORFOLOGIA_PONTUAL")
    elif s_morph == 0: flags.append("MORFOLOGIA_ESTENDIDA")

    if s_elongation == 0:         flags.append("ELONGACAO_ALTA")
    if s_morph_consistency == 0: flags.append("MORFO_INCONSISTENTE")

    borda = any(
        _is_near_edge(track[fi][1], track[fi][0], image_shape)
        for fi in range(4)
    )
    if borda:
        flags.append("EDGE_FRAME")

    if rejection_reasons:
        flags.append("REJECTED_FP")
    if gaia_static["status"] == "fonte_estatica_gaia":
        flags.append("GAIA_FONTE_ESTATICA")

    # ── Log de decisão ───────────────────────────────────────────────────────
    if rejection_reasons:
        for r in rejection_reasons:
            log(f"  REJEIÇÃO FP: {r}", "WARN")
    elif penalty_reasons:
        for r in penalty_reasons:
            log(f"  Penalização: {r}", "WARN")
    else:
        log(f"  Promovido: linearity={lin:.2f}px, vel_cons={vel_cons:.2f}, "
            f"CV={flux_cv:.2f}, pont={pont:.4f}", "OK")

    return {
        **c,
        "linearity"           : lin,
        "linear_r2"             : linear_r2,
        "vel_consistencia"      : vel_cons,
        "vel_arcsec_min"        : vel_arcsec_min,
        "fluxes"               : flux_array.tolist(),
        "flux_cv"             : flux_cv,
        "pointlike"               : pont,
        "pont_std"              : pont_std,
        "mean_elongation"           : mean_elongation,
        "mean_fwhm_px"         : fwhm_medio,
        "morphology_metrics"        : morphology_metrics,
        "snrs"                  : snrs,
        "astrometric_uncertainty": astrometric_uncertainty,
        "gaia_static"           : gaia_static,
        "score"                 : score,
        "raw_score"           : raw_score,
        "score_max"             : score_max,
        "score_linearity"     : s_lin,
        "score_velocity"      : s_velocity,
        "score_photometry"      : s_phot,
        "score_morphology"      : s_morph,
        "score_morph_consistency"   : s_morph_consistency,
        "score_elongation"      : s_elongation,
        "score_velocity_range"       : s_range,
        "score_percent"         : score_percent,
        "score_normalized"      : round(score / score_max, 4),
        "heuristic_score"       : score,
        "priority_score"        : score,
        "triage_class"                : triage_class,
        "triage_class"          : {
            "STRONG": "STRONG",
            "MODERATE": "MODERATE",
            "WEAK": "WEAK",
            "DISCARDED": "DISCARDED",
        }.get(triage_class, triage_class),
        "cli_color"               : cli_color,
        "flags"                 : flags,
        "penalty_reasons"    : penalty_reasons,
        "rejection_reasons"       : rejection_reasons,
    }


# ─────────────────────────────────────────────
# 5. CONVERSÃO PIXEL → RA/DEC
# ─────────────────────────────────────────────
def pixel_to_radec(x: float, y: float, wcs: WCS) -> tuple[float, float]:
    if wcs is None:
        raise ValueError("WCS ausente — impossível converter pixel em RA/Dec.")
    ra, dec = wcs.all_pix2world(x, y, 0)
    return float(ra), float(dec)


# ─────────────────────────────────────────────
# 6. CONSULTA AO MINOR PLANET CENTER
# ─────────────────────────────────────────────
def query_position_catalog(ra: float, dec: float, data_obs: str,
                                  raio_arcmin: float = 2.0) -> dict:
    """
    Consulta o SkyBot (IMCCE) para procurar correspondência posicional próxima.

    Retorna dict com schema normalizado:
        status: not_queried | query_failed | no_match |
                ambiguous_match  | likely_match
    """
    if not (0 <= ra < 360) or not (-90 <= dec <= 90):
        return {
            "status": "query_failed",
            "reason": f"coordenadas fora do domínio (RA={ra:.3f}, Dec={dec:.3f})",
            "match" : None,
        }

    try:
        from astroquery.imcce import Skybot

        coord    = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        epoch    = Time(data_obs, format="isot", scale="utc")
        result = Skybot.cone_search(coord, rad=raio_arcmin * u.arcmin,
                                       epoch=epoch)

        if result is not None and len(result) > 0:
            # Mais de um objeto no cone → ambiguous
            status = "ambiguous_match" if len(result) > 1 else "likely_match"
            obj    = result[0]

            num  = str(obj["Number"])  if "Number" in result.colnames else ""
            name = str(obj["Name"])    if "Name"   in result.colnames \
                   else str(obj.get("Designation", "desconhecido"))
            mag  = float(obj["V"])     if "V"      in result.colnames else None
            dist = float(obj["centerdist"].to(u.arcmin).value) \
                   if "centerdist" in result.colnames else None
            type = str(obj["Type"])    if "Type"   in result.colnames else None

            return {
                "status"  : status,
                "reason"  : None,
                "n_objects_in_cone": len(result),
                "match"   : {
                    "name_designation"     : f"{num} {name}".strip(),
                    "type"                : type,
                    "magnitude_v"         : round(mag, 2) if mag is not None else None,
                    "distancia_arcmin"    : round(dist, 4) if dist is not None else None,
                },
            }

        return {
            "status"  : "no_match",
            "reason"  : None,
            "n_objects_in_cone": 0,
            "match"   : None,
        }

    except ImportError:
        return {
            "status": "query_failed",
            "reason": "astroquery is not installed",
            "match" : None,
        }
    except Exception as e:
        msg = str(e)
        # SkyBot retorna "No table found" quando o cone não contém objects
        if "no table found" in msg.lower() or "no result" in msg.lower():
            return {
                "status"           : "no_match",
                "reason"           : None,
                "n_objects_in_cone": 0,
                "match"            : None,
            }
        return {
            "status": "query_failed",
            "reason": msg[:200],
            "match" : None,
        }


def _add_mpc_flag(c: dict) -> None:
    """Adiciona flag MPC ao candidate com base no result da consulta."""
    mpc = c.get("mpc", {})
    status = mpc.get("status", "not_queried")
    flags = c["flags"]

    if status == "no_match":
        flags.append("MPC_SEM_MATCH")
    elif status == "likely_match":
        flags.append("MPC_MATCH_PROVAVEL")
    elif status == "ambiguous_match":
        flags.append("MPC_MATCH_AMBIGUO")
    elif status == "query_failed":
        flags.append("MPC_CONSULTA_FALHOU")


# ─────────────────────────────────────────────
# 7. GERAÇÃO DO RELATÓRIO MPC
# ─────────────────────────────────────────────
def format_mpc_ra(ra_deg: float) -> str:
    ra_wrap = float(ra_deg) % 360.0
    return Angle(ra_wrap * u.deg).to_string(
        unit=u.hourangle, sep=" ", precision=2, pad=True,
    )


def format_mpc_dec(dec_deg: float) -> str:
    dec_clip = max(-90.0, min(90.0, float(dec_deg)))
    return Angle(dec_clip * u.deg).to_string(
        unit=u.deg, sep=" ", precision=1, pad=True, alwayssign=True,
    )


def format_mpc_date(jd: float) -> str:
    t  = Time(jd, format="jd", scale="utc")
    dt = t.to_datetime()
    frac = (dt.hour * 3600 + dt.minute * 60
            + dt.second + dt.microsecond / 1e6) / 86400.0
    return f"{dt.year:04d} {dt.month:02d} {dt.day + frac:08.5f}"


def generate_mpc_report(candidates: list[dict], frames: list[dict],
                        input_set_name: str, output_dir: Path,
                        observer: dict) -> Path:
    """
    Gera file de relatório no formato MPC 80-colunas.
    Apenas candidates com score >= 6 (MODERATE+) são incluídos.
    Magnitude deixada em branco — requer calibração fotométrica no Astrometrica.
    """
    lines  = []
    included = [c for c in candidates if c["score"] >= T.MPC_SCORE_MIN]

    for i, c in enumerate(included):
        track     = c["track"]
        designation = f"TMP{i+1:04d}"

        for fi, frame in enumerate(frames):
            if not frame["wcs_ok"]:
                continue
            tx, ty   = track[fi][1], track[fi][0]
            ra, dec  = pixel_to_radec(tx, ty, frame["wcs"])
            mpc_date = format_mpc_date(frame["jd"])
            ra_mpc   = format_mpc_ra(ra)
            dec_mpc  = format_mpc_dec(dec)

            row = (
                "     "
                + f"{designation:<7s}"
                + " "
                + " "
                + "C"
                + mpc_date
                + " "
                + ra_mpc
                + " "
                + dec_mpc
                + "         "
                + "     "
                + " "
                + "     "
                + "F51"
            )
            lines.append(row[:80].ljust(80))

    path   = output_dir / f"{input_set_name}_MPC_report.txt"
    observer_name  = observer.get("name",  "IASC Observer")
    observer_email = observer.get("email", "observer@example.com")
    with open(path, "w", encoding="ascii", errors="replace") as f:
        f.write("COD F51\n")
        f.write("OBS IASC Citizen Scientist\n")
        f.write(f"MEA {observer_name}\n")
        f.write("TEL 1.8-m f/4.4 Ritchey-Chretien + CCD\n")
        f.write(f"ACK {PROJECT_NAME} pipeline v{PIPELINE_VERSION}\n")
        f.write(f"AC2 {observer_email}\n")
        f.write("----- ------------------- ------------------ "
                "------------------ ----- ---\n")
        for row in lines:
            f.write(row + "\n")
        f.write("-----\n")

    return path


# ─────────────────────────────────────────────
# 8. VISUALIZAÇÃO
# ─────────────────────────────────────────────
def generate_visualization(candidates: list[dict], frames: list[dict],
                       output_dir: Path, input_set_name: str):
    """
    Gera PNG com os top-5 candidates (não-DISCARDED).
    Cada row: 4 recortes do candidate nos 4 frames.
    Inclui: rank, score, flags principais, vetor de trajetória.
    """
    top = [c for c in candidates if c["triage_class"] != "DISCARDED"][:5]
    if not top:
        log("No valid candidate for visualization.", "WARN")
        return

    triage_colors = {
        "STRONG"   : "#00e07a",
        "MODERATE": "#ffaa00",
        "WEAK"   : "#ff7733",
        "DISCARDED": "#ff2222",
    }
    times = [f["date_obs"][11:16] + " UT" for f in frames]
    n = len(top)
    fig, axes = plt.subplots(n, 4, figsize=(20, 4.5 * n), squeeze=False)
    fig.patch.set_facecolor("#0d0d1a")

    for ri, cand in enumerate(top):
        cor    = triage_colors[cand["triage_class"]]
        track = cand["track"]
        flags  = cand.get("flags", [])

        # Flags mais relevantes para exibir na image (máximo 3)
        flags_exibir = [f for f in flags
                        if f not in ("MPC_SEM_MATCH",)][:3]
        flags_str = "  ".join(flags_exibir) if flags_exibir else ""

        for fi, frame in enumerate(frames):
            ax = axes[ri][fi]
            ty, tx = track[fi][0], track[fi][1]
            r = 55
            y0 = max(0, int(ty) - r); y1 = min(frame["data"].shape[0], int(ty) + r)
            x0 = max(0, int(tx) - r); x1 = min(frame["data"].shape[1], int(tx) + r)
            cut = frame["data"][y0:y1, x0:x1]
            p1, p99 = np.percentile(cut, (1, 99.5))
            cut_n = np.clip((cut - p1) / (p99 - p1 + 1e-10), 0, 1)
            ax.imshow(cut_n, cmap="gray", origin="lower")

            # Círculo na posição atual
            cx_local = tx - x0
            cy_local = ty - y0
            ax.add_patch(plt.Circle((cx_local, cy_local), 9,
                                    color=cor, fill=False, lw=1.8))

            # Posições anteriores como cruzes e vetor de trajetória
            prev_coords = []
            for pfi in range(fi):
                ppy, ppx = track[pfi][0], track[pfi][1]
                px_l = ppx - x0
                py_l = ppy - y0
                ax.plot(px_l, py_l, "+", color=cor,
                        alpha=0.45, ms=7, mew=1.3)
                prev_coords.append((px_l, py_l))

            # Vetor de trajetória: seta do primeiro ao atual frame (visível no frame 3+)
            if fi >= 2 and len(prev_coords) >= 1:
                px_prev, py_prev = prev_coords[0]
                ax.annotate("",
                    xy=(cx_local, cy_local),
                    xytext=(px_prev, py_prev),
                    arrowprops=dict(
                        arrowstyle="->",
                        color=cor,
                        lw=1.2,
                        alpha=0.7,
                    ),
                )

            ax.set_facecolor("#0d0d1a")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_edgecolor(cor); sp.set_linewidth(1.5)

            # Label lateral esquerdo (apenas no frame 0)
            if fi == 0:
                mpc_status = cand.get("mpc", {}).get("status", "")
                mpc_label  = {
                    "no_match"      : "SkyBot: sem match",
                    "likely_match" : "SkyBot: match",
                    "ambiguous_match"  : "SkyBot: ambiguous",
                    "query_failed": "SkyBot: falhou",
                    "not_queried" : "",
                }.get(mpc_status, "")

                ra_str = ""
                if frame["wcs_ok"]:
                    ra_v, dec_v = pixel_to_radec(tx, ty, frame["wcs"])
                    ra_str = f"RA {ra_v:.4f}° Dec {dec_v:.4f}°"

                ylabel_lines = [
                    f"#{ri+1} {cand['triage_class']}",
                    f"Score {cand['score']}/{cand.get('score_max', 12)}  "
                    f"({cand['score_percent']}%)",
                    f"Lin={cand['linearity']:.2f}px  Vel={cand['vel_consistencia']:.2f}",
                ]
                if ra_str:
                    ylabel_lines.append(ra_str)
                if mpc_label:
                    ylabel_lines.append(mpc_label)

                ax.set_ylabel(
                    "\n".join(ylabel_lines),
                    color=cor, fontsize=6.5, labelpad=4,
                )

            # Título do frame (apenas na row 0)
            if ri == 0:
                ax.set_title(f"Frame {fi+1}\n{times[fi]}",
                             color="white", fontsize=9)

            # Score e rank no canto superior direito do frame 0 de cada candidate
            if fi == 0:
                ax.text(
                    0.97, 0.97,
                    f"#{ri+1}  {cand['score']}/{cand.get('score_max', 12)}",
                    transform=ax.transAxes,
                    color=cor, fontsize=8, fontweight="bold",
                    ha="right", va="top",
                    bbox=dict(boxstyle="round,pad=0.2",
                              facecolor="#0d0d1a", alpha=0.7,
                              edgecolor=cor, linewidth=0.8),
                )

            # Flags resumidas no canto inferior (frame 3)
            if fi == 3 and flags_str:
                ax.text(
                    0.02, 0.03,
                    flags_str,
                    transform=ax.transAxes,
                    color=cor, fontsize=5.5,
                    ha="left", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.2",
                              facecolor="#0d0d1a", alpha=0.65,
                              edgecolor=cor, linewidth=0.6),
                )

    plt.suptitle(
        f"{PROJECT_NAME} v{PIPELINE_VERSION} | {input_set_name} | "
        f"{frames[0]['date_obs'][:10]}\n"
        "● posição atual   + posição anterior   → vetor trajetória",
        color="white", fontsize=11, y=1.01,
    )
    plt.tight_layout(h_pad=0.8, w_pad=0.4)
    path = output_dir / f"{input_set_name}_candidates.png"
    plt.savefig(path, dpi=130, bbox_inches="tight", facecolor="#0d0d1a")
    plt.close()
    log(f"Visualização salva: {path}", "OK")


# ─────────────────────────────────────────────
# 9. RELATÓRIO TEXTO FINAL
# ─────────────────────────────────────────────
def _diagnostic_summary(c: dict) -> str:
    """
    Gera uma frase técnica curta explicando o ranking do candidate.
    Inclui razões de penalização e rejeição quando presentes.
    """
    parts = []

    # Rejection reasons têm precedência
    if c.get("rejection_reasons"):
        parts.append("REJEITADO: " + c["rejection_reasons"][0])
        if len(c["rejection_reasons"]) > 1:
            parts.append(f"+{len(c['rejection_reasons'])-1} razão(ões) adicionais")
        return "; ".join(parts) + "."

    lin = c["linearity"]
    if lin < T.LIN_EXCELENTE:
        parts.append(f"track linear (res={lin:.2f} px)")
    elif lin < T.LIN_MARGINAL:
        parts.append(f"linearity marginal (res={lin:.2f} px)")
    else:
        parts.append(f"linearity ruim (res={lin:.2f} px)")

    vel = c["vel_consistencia"]
    if vel < T.VEL_UNIFORM:
        parts.append("velocity uniforme")
    else:
        parts.append(f"velocity irregular (σ={vel:.2f})")

    cv = c["flux_cv"]
    if cv < T.PHOT_STABLE:
        parts.append("flux estável")
    elif cv < T.PHOT_MODERATE:
        parts.append(f"flux moderadamente variável (CV={cv:.2f})")
    else:
        parts.append(f"flux instável (CV={cv:.2f})")

    elong = c.get("mean_elongation", 1.0)
    if elong >= T.ELON_GOOD:
        parts.append(f"elongação elevada ({elong:.2f})")

    if c.get("penalty_reasons"):
        n = len(c["penalty_reasons"])
        parts.append(f"{n} penalização(ões) de score")

    mpc_status = c.get("mpc", {}).get("status", "not_queried")
    if mpc_status == "likely_match":
        name = c["mpc"].get("match", {}).get("name_designation", "?")
        dist = c["mpc"].get("match", {}).get("distancia_arcmin", "?")
        parts.append(f"match SkyBot: {name} a {dist}'")
    elif mpc_status == "no_match":
        parts.append("no positional match in SkyBot")
    elif mpc_status == "ambiguous_match":
        parts.append("múltiplos objects no cone SkyBot — verificar")

    if "EDGE_FRAME" in c.get("flags", []):
        parts.append("posição na borda do frame")

    return "; ".join(parts) + "."


def generate_text_report(candidates: list[dict], frames: list[dict],
                        input_set_name: str, output_dir: Path,
                        global_metrics: dict) -> Path:
    """Gera relatório operacional em text com priorização e diagnóstico."""
    path = output_dir / f"{input_set_name}_report.txt"
    lines  = []
    sep     = "=" * 72
    sub     = "-" * 72
    sub2    = "·" * 72

    wcs_mode = global_metrics.get("wcs_mode", "gaia")
    run_id   = global_metrics.get("run_id", "n/a")

    # ── Cabeçalho ──────────────────────────────────────────────────────────
    lines += [
        sep,
        "  Lia — Analysis Report",
        f"  Versão do pipeline   : {PIPELINE_VERSION}",
        f"  Data de processamento: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"  Conjunto de images  : {input_set_name}",
        f"  Modo WCS             : {wcs_mode.upper()}",
        f"  run_id               : {run_id}",
        f"  Observatório         : Pan-STARRS / IASC (código F51)",
        f"  Data de observação   : {frames[0]['date_obs'][:10]}",
        f"  Intervalo temporal   : {frames[0]['date_obs'][11:19]} → "
        f"{frames[-1]['date_obs'][11:19]} UTC",
        sep, "",
    ]

    # ── Métricas globais do input_set ───────────────────────────────────────
    mg = global_metrics
    lines += [
        "  PROCESSED SET METRICS",
        sub,
        f"  Modo WCS             : {wcs_mode.upper()}",
        f"  run_id               : {run_id}",
        f"  Frames válidos       : {mg.get('n_valid_frames', 4)} de 4",
        f"  Sigma de detecção    : {mg.get('sigma', '?')}",
        f"  Fontes detectadas    : " +
        "  |  ".join(
            f"F{i+1}={n}" for i, n in enumerate(mg.get("n_fontes_by_frame", []))
        ),
        f"  Trilhas tentadas     : {mg.get('n_tracks_attempted', '?')}",
        f"  Unique candidates    : {mg.get('n_unique_candidates', '?')}",
        f"  Rejeitados FP        : {mg.get('n_rejected_false_positive', '?')}",
        f"  Deriva do field      : dx={mg.get('drift_dx_px', 0):+.2f} px  "
        f"dy={mg.get('drift_dy_px', 0):+.2f} px  "
        f"[{mg.get('n_estaveis_deriva', '?')} estrelas]",
        f"  Tempo de execução    : {mg.get('execution_time_s', '?')} s",
        "", sub,
        f"  Total candidates    : {len(candidates)}",
        f"  STRONG    (8-12) : {sum(1 for c in candidates if c['triage_class']=='STRONG')}",
        f"  MODERATE (6-7)  : {sum(1 for c in candidates if c['triage_class']=='MODERATE')}",
        f"  WEAK    (4-5)  : {sum(1 for c in candidates if c['triage_class']=='WEAK')}",
        f"  DISCARDED (0-3)  : {sum(1 for c in candidates if c['triage_class']=='DISCARDED')}",
    ]

    # Score resumido
    scores_validos = [c["score"] for c in candidates if c["triage_class"] != "DISCARDED"]
    if scores_validos:
        lines += [
            f"  Score (not DISCARDED) : "
            f"min={min(scores_validos)}  max={max(scores_validos)}  "
            f"med={round(float(np.median(scores_validos)), 1)}",
        ]
    lines += ["", sep, ""]

    # ── Seção de priorização operacional ───────────────────────────────────
    lines += [
        "  OPERATIONAL PRIORITIZATION",
        "  (use as an Astrometrica inspection guide, not as a final decision)",
        sub,
    ]

    grupos = {
        "INSPECIONAR PRIMEIRO"      : [],
        "INSPECIONAR SE HOUVER TEMPO": [],
        "LIKELY ARTIFACT"         : [],
        "IGNORAR"                   : [],
    }
    for c in candidates:
        cl = c["triage_class"]
        mpc_status = c.get("mpc", {}).get("status", "not_queried")
        borda = "EDGE_FRAME" in c.get("flags", [])
        if cl == "STRONG" and mpc_status in ("no_match", "not_queried"):
            grupos["INSPECIONAR PRIMEIRO"].append(c)
        elif cl in ("STRONG", "MODERATE"):
            grupos["INSPECIONAR SE HOUVER TEMPO"].append(c)
        elif cl == "WEAK":
            grupos["LIKELY ARTIFACT"].append(c)
        else:
            grupos["IGNORAR"].append(c)

    rank_global = {c["rank"]: c for c in candidates if "rank" in c}

    for group_name, group_candidates in grupos.items():
        if not group_candidates:
            continue
        lines += [f"  ┌─ {group_name} ({len(group_candidates)}) "]
        for c in group_candidates:
            rank_label = f"#{c.get('rank', '?')}"
            flags_str  = "  ".join(c.get("flags", []))
            resumo     = _diagnostic_summary(c)
            lines += [
                f"  │  {rank_label}  score {c['score']}/{c.get('score_max', 12)}  {c['triage_class']}",
                f"  │     Lin={c['score_linearity']}/3  "
                f"Vel={c['score_velocity']}/2  "
                f"Fot={c['score_photometry']}/2  "
                f"Mor={c['score_morphology']}/2  "
                f"Faixa={c['score_velocity_range']}/1",
                f"  │     Flags: {flags_str if flags_str else '—'}",
                f"  │     {resumo}",
                f"  │",
            ]
        lines[-1] = lines[-1].replace("  │", "  └─")
        lines += [""]

    lines += [sep, ""]

    # ── Detalhamento por candidate ─────────────────────────────────────────
    lines += ["  DETALHAMENTO POR CANDIDATO", sub, ""]

    for c in candidates:
        if c["triage_class"] == "DISCARDED":
            continue
        track = c["track"]
        tx0, ty0 = track[0][1], track[0][0]
        if frames[0]["wcs_ok"]:
            ra, dec = pixel_to_radec(tx0, ty0, frames[0]["wcs"])
        else:
            ra, dec = float("nan"), float("nan")
        mpc = c.get("mpc", {"status": "not_queried"})
        mpc_status = mpc.get("status", "not_queried")

        # Status MPC legível
        if mpc_status == "likely_match":
            match = mpc.get("match", {})
            status_mpc_txt = (
                f"MATCH PROVÁVEL — {match.get('name_designation', '?')}  "
                f"mag={match.get('magnitude_v', '?')}  "
                f"dist={match.get('distancia_arcmin', '?')}'  "
                f"type={match.get('type', '?')}"
            )
        elif mpc_status == "ambiguous_match":
            status_mpc_txt = (
                f"MATCH AMBÍGUO — {mpc.get('n_objects_in_cone', '?')} "
                "objects no cone de busca"
            )
        elif mpc_status == "no_match":
            status_mpc_txt = "SEM CORRESPONDÊNCIA no cone de 2'"
        elif mpc_status == "query_failed":
            status_mpc_txt = f"CONSULTA FALHOU — {mpc.get('reason', '')}"
        else:
            status_mpc_txt = "Not queried (offline mode)"

        flags_str = "  ".join(c.get("flags", [])) or "—"

        # Razões de decisão
        razoes_rej = c.get("rejection_reasons", [])
        razoes_pen = c.get("penalty_reasons", [])
        score_max  = c.get("score_max", 12)

        lines += [
            f"  CANDIDATO #{c.get('rank', '?'):02}  [{c.get('candidate_id', '?')}]",
            sub,
            f"  Classification    : {c['triage_class']}  "
            f"(score {c['score']}/{score_max}  |  normalized {c['score_percent']}%)",
            f"  Flags            : {flags_str}",
            f"  Diagnostic      : {_diagnostic_summary(c)}",
        ]

        if razoes_rej:
            lines.append(f"  ── Rejection reasons ───────────────────────────────────")
            for r in razoes_rej:
                lines.append(f"  ✗  {r}")

        if razoes_pen:
            lines.append(f"  ── Penalty reasons ────────────────────────────────")
            for r in razoes_pen:
                lines.append(f"  ↓  {r}")

        lines += [
            "",
            f"  ── Decomposição do score ────────────────────────────────────",
            f"  Linearidade        : {c['score_linearity']}/3  "
            f"(resíduo {c['linearity']:.3f} px)",
            f"  Velocidade         : {c['score_velocity']}/2  "
            f"(σ passos {c['vel_consistencia']:.3f})",
            f"  Fotometria         : {c['score_photometry']}/2  "
            f"(CV flux {c['flux_cv']:.3f})",
            f"  Morfologia         : {c['score_morphology']}/2  "
            f"(pointness média {c['pointlike']:.4f})",
            f"  Consist. morfológ. : {c.get('score_morph_consistency', '?')}/1  "
            f"(std pointness {c.get('pont_std', 0.0):.4f})",
            f"  Elongation         : {c.get('score_elongation', '?')}/1  "
            f"(elongação média {c.get('mean_elongation', 1.0):.2f})",
            f"  Faixa vel. típica  : {c['score_velocity_range']}/1  "
            f"({c['move_total']:.1f} px total)",
            f"  Score total        : {c['score']}/{score_max}",
            "",
            f"  ── Position (Frame 1) ────────────────────────────────────────",
            f"  Pixel            : x={tx0:.1f}  y={ty0:.1f}",
            f"  Coordenadas      : RA = {ra:.6f}°  Dec = {dec:.6f}°",
            f"  RA  (HH MM SS)   : {format_mpc_ra(ra)}",
            f"  Dec (±DD MM SS)  : {format_mpc_dec(dec)}",
            "",
            f"  ── Movimento ────────────────────────────────────────────────",
            f"  Deslocamento total : {c['move_total']:.2f} px",
            f"  Direction (dx, dy)   : ({track[-1][1]-track[0][1]:.1f}, "
            f"{track[-1][0]-track[0][0]:.1f}) px",
        ]

        if c.get("vel_arcsec_min") is not None:
            lines.append(f"  Velocidade angular : {c['vel_arcsec_min']:.2f} arcsec/min")

        fwhm_str = (f"{c['mean_fwhm_px']:.1f} px"
                    if c.get("mean_fwhm_px") is not None else "n/d")
        lines += [
            "",
            f"  ── Fotometria e Morfologia ──────────────────────────────────",
            f"  Brilho por frame  : {[round(b, 1) for b in c['fluxes']]}",
            f"  SNR estimado      : {[round(s, 1) for s in c.get('snrs', [])]}",
            f"  Variation (CV)     : {c['flux_cv']:.3f}  "
            + ("✓ estável" if c["flux_cv"] < T.PHOT_STABLE else "! variável"),
            f"  Perfil pointlike    : {c['pointlike']:.4f}  "
            + ("✓ ponto" if c["pointlike"] > T.MOR_MARGINAL else "! estendido"),
            f"  Mean elongation   : {c.get('mean_elongation', 1.0):.2f}  "
            + ("✓" if c.get("mean_elongation", 1.0) < T.ELON_GOOD else "! estendida"),
            f"  FWHM estimada     : {fwhm_str}",
            "",
            f"  ── Status MPC/SkyBot ────────────────────────────────────────",
            f"  {status_mpc_txt}",
            "",
            f"  ── Trilha frame a frame ─────────────────────────────────────",
        ]
        for fi, t in enumerate(track):
            if frames[fi]["wcs_ok"]:
                ra_f, dec_f = pixel_to_radec(t[1], t[0], frames[fi]["wcs"])
                coord_str = f"RA={ra_f:.5f}° Dec={dec_f:.5f}°"
            else:
                coord_str = "RA=?  Dec=?  (WCS inválido)"
            snr_str = f"  SNR≈{c['snrs'][fi]:.1f}" if c.get("snrs") else ""
            lines.append(
                f"  Frame {fi+1} ({frames[fi]['date_obs'][11:19]} UTC)  "
                f"x={t[1]:.1f} y={t[0]:.1f}  flux={t[2]:.0f}  {coord_str}{snr_str}"
            )

        lines += ["", sub, ""]

    # ── Metodologia ────────────────────────────────────────────────────────
    lines += [
        sep,
        "  METODOLOGIA",
        sub,
        "  1. Background local via photutils.Background2D (box 50x50, filtro 3x3)",
        "  2. Point-source seeds in local-SNR image + Moffat/Gaussian PSF fit",
        "     para centro sub-pixel e FWHM compatível com PSF estelar",
        "  3. Rastreamento em plano celeste local via WCS por frame",
        "     + casamento mútuo recíproco entre frames consecutivos",
        "  4. Refinamento WCS por Gaia DR3 com fallback para Pan-STARRS",
        "  5. Rejeição de falsos positivos: borda, SNR baixo, hot pixel heurístico,",
        "     flux caótico, elongação/FWHM incompatível com source pointlike",
        f"  6. Score 0-{3+2+2+2+1+1+1}: linearity(3) + velocity(2) + photometry(2)",
        "     + morphology(2) + consist.morph(1) + elongação(1) + velocity_range(1)",
        "  7. Morfologia: pointness, elongação via PCA, FWHM PSF, consist. frames",
        "  8. Conversão pixel→RA/Dec via WCS Pan-STARRS refinado por Gaia DR3",
        "  9. Correspondência posicional: SkyBot/IMCCE, cone 2 arcmin",
        " 10. MPC 80-colunas com magnitude omitida (requer calibração Astrometrica)",
        "  NOTE: the score is a heuristic triage score, not calibrated confidence.",
        "        Thresholds are centralized in T.* for reproducible tuning.",
        "",
        "  LIMITAÇÕES",
        sub,
        "  - No absolute photometric calibration (blank magnitude in MPC draft)",
        "  - WCS usa Pan-STARRS como solução inicial e Gaia DR3 quando disponível",
        "  - Refined heuristic score: does not replace digest2 or CNN methods",
        "  - SkyBot indica vizinhança posicional — ausência de match ≠ novidade",
        "  - Validação humana no Astrometrica é obrigatória antes de enviar ao IASC",
        "  - Heurísticas de falso positivo podem rejeitar candidates legítimos em",
        "    fields com PSF degradada — revisar DISCARDED com flag REJECTED_FP",
        "",
        f"  Pipeline : {PROJECT_NAME} v{PIPELINE_VERSION}",
        "  Autora   : Jaciana Barbosa",
        sep,
    ]

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return path


# ─────────────────────────────────────────────
# 10. EXPORTAÇÃO JSON
# ─────────────────────────────────────────────
def _serialize_track(c: dict, frames: list[dict]) -> list[dict]:
    """Serializa a track completa frame a frame para o JSON."""
    track_json = []
    for fi, t in enumerate(c["track"]):
        ty, tx, flux = t[0], t[1], t[2]
        frame = frames[fi]

        ra_deg, dec_deg = None, None
        ra_fmt, dec_fmt = None, None
        if frame["wcs_ok"]:
            try:
                ra_deg, dec_deg = pixel_to_radec(tx, ty, frame["wcs"])
                ra_fmt  = format_mpc_ra(ra_deg)
                dec_fmt = format_mpc_dec(dec_deg)
                ra_deg  = round(ra_deg, 6)
                dec_deg = round(dec_deg, 6)
            except Exception:
                pass

        snr = c.get("snrs", [None, None, None, None])[fi]
        incertezas_frame = c.get("astrometric_uncertainty", {}).get("by_frame") or []
        incerteza_frame = incertezas_frame[fi] if fi < len(incertezas_frame) else {}

        track_json.append({
            "frame_index"   : fi,
            "timestamp_utc" : frame["date_obs"],
            "timestamp_inicio_utc": frame.get("date_obs_inicio", frame["date_obs"]),
            "jd"            : round(frame["jd"], 6),
            "jd_inicio"     : round(frame.get("jd_inicio", frame["jd"]), 6),
            "jd_mid"        : round(frame.get("jd_mid", frame["jd"]), 6),
            "mjd_inicio"    : round(frame.get("mjd_inicio", Time(frame["jd"], format="jd").mjd), 8),
            "mjd_mid"       : round(frame.get("mjd_mid", Time(frame["jd"], format="jd").mjd), 8),
            "x"             : round(tx, 2),
            "y"             : round(ty, 2),
            "ra_deg"        : ra_deg,
            "dec_deg"       : dec_deg,
            "ra_fmt"        : ra_fmt,
            "dec_fmt"       : dec_fmt,
            "flux"          : round(float(flux), 1),
            "snr"           : round(float(snr), 2) if snr is not None else None,
            "astrometric_uncertainty": incerteza_frame,
            "wcs_valid"    : frame["wcs_ok"],
        })
    return track_json


def export_json(candidates: list[dict], frames: list[dict],
                  input_set_name: str, output_dir: Path,
                  global_metrics: dict,
                  run_metadata: dict | None = None) -> Path:
    """
    Exporta JSON estruturado com:
    - run_metadata: run_id, pipeline_version, wcs_mode, timestamp, input_set
    - métricas globais do input_set
    - por candidate: candidate_id rastreável, score decomposto, track completa,
      flags, status MPC normalizado, colunas manuais vazias para avaliação
    """
    if run_metadata is None:
        run_metadata = {
            "run_id"           : "n/a",
            "pipeline_name"    : PROJECT_NAME,
            "pipeline_version" : PIPELINE_VERSION,
            "repository"       : REPOSITORY_NAME,
            "wcs_mode"         : global_metrics.get("wcs_mode", "gaia"),
            "timestamp_execution_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "timestamp_execution_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "input_set"        : str(output_dir),
            "input_set_name"    : input_set_name,
            "input_set_name"   : input_set_name,
            "sigma"            : global_metrics.get("sigma"),
            "observer"       : {},
        }

    # Cabeçalho global
    output = {
        "pipeline"       : f"{PROJECT_NAME} v{PIPELINE_VERSION}",
        "pipeline_name"  : PROJECT_NAME,
        "repository"     : REPOSITORY_NAME,
        "input_set"       : input_set_name,
        "input_set"      : input_set_name,
        "processado_em"  : datetime.now().isoformat(timespec="seconds"),
        "observatory"   : "Pan-STARRS / IASC (F51)",
        "data_obs"       : frames[0]["date_obs"][:10],
        "run_metadata"   : run_metadata,
        "global_metrics": global_metrics,
        "global_metrics" : global_metrics,
        "candidates"     : [],
        "candidates"     : [],
    }

    for c in candidates:
        tx0, ty0 = c["track"][0][1], c["track"][0][0]
        ra_r, dec_r, ra_fmt, dec_fmt = None, None, None, None
        if frames[0]["wcs_ok"]:
            try:
                ra_r, dec_r = pixel_to_radec(tx0, ty0, frames[0]["wcs"])
                ra_fmt  = format_mpc_ra(ra_r)
                dec_fmt = format_mpc_dec(dec_r)
                ra_r    = round(ra_r, 6)
                dec_r   = round(dec_r, 6)
            except Exception:
                pass

        entry = {
            # Identificação rastreável
            "candidate_id"      : c.get("candidate_id"),
            "run_id"            : run_metadata.get("run_id"),
            "wcs_mode"          : run_metadata.get("wcs_mode"),
            "rank"              : c.get("rank"),
            "triage_class"      : c.get("triage_class"),
            "score_total"       : c["score"],
            "heuristic_score"   : c.get("heuristic_score", c["score"]),
            "priority_score"    : c.get("priority_score", c["score"]),
            "score_normalized"  : c.get("score_normalized"),
            "score_percent"     : c["score_percent"],

            # Score decomposto
            "score_components" : {
                "linearity"       : c["score_linearity"],
                "velocity"        : c["score_velocity"],
                "photometry"        : c["score_photometry"],
                "morphology"        : c["score_morphology"],
                "morph_consistency": c.get("score_morph_consistency", 0),
                "elongation"        : c.get("score_elongation", 0),
                "velocity_range"         : c["score_velocity_range"],
                "maximos"           : {
                    "linearity": 3, "velocity": 2, "photometry": 2,
                    "morphology": 2, "morph_consistency": 1,
                    "elongation": 1, "velocity_range": 1,
                },
                "score_max"         : c.get("score_max", 12),
            },
            # Flags diagnósticas
            "flags"             : c.get("flags", []),

            # Razões de decisão auditáveis
            "decision_reasons"    : {
                "penalties" : c.get("penalty_reasons", []),
                "rejections"    : c.get("rejection_reasons", []),
            },
            "decision_reasons"  : {
                "penalties": c.get("penalty_reasons", []),
                "rejections": c.get("rejection_reasons", []),
            },

            # Position resumida (frame 1)
            "posicao_frame1"    : {
                "pixel_x" : round(tx0, 2),
                "pixel_y" : round(ty0, 2),
                "ra_deg"  : ra_r,
                "dec_deg" : dec_r,
                "ra_fmt"  : ra_fmt,
                "dec_fmt" : dec_fmt,
            },

            # Métricas de motion
            "motion"         : {
                "total_px"        : round(c["move_total"], 3),
                "residuo_deriva"  : round(c.get("residual", 0.0), 3),
                "linearidade_px"  : round(c["linearity"], 4),
                "linear_r2"       : round(c["linear_r2"], 6)
                                    if c.get("linear_r2") is not None else None,
                "vel_consistencia": round(c["vel_consistencia"], 4),
                "vel_arcsec_min"  : c.get("vel_arcsec_min"),
            },
            "motion"            : {
                "total_px": round(c["move_total"], 3),
                "field_drift_residual_px": round(c.get("residual", 0.0), 3),
                "linearity_residual_px": round(c["linearity"], 4),
                "linear_r2": round(c["linear_r2"], 6)
                             if c.get("linear_r2") is not None else None,
                "step_consistency_px": round(c["vel_consistencia"], 4),
                "rate_arcsec_min": c.get("vel_arcsec_min"),
            },

            # Fotometria e morphology
            "photometry"        : {
                "flux_by_frame": [round(b, 1) for b in c["fluxes"]],
                "flux_cv"       : round(c["flux_cv"], 4)
                                    if c["flux_cv"] != float("inf") else None,
                "pointness"    : round(c["pointlike"], 5),
                "pont_std"        : round(c.get("pont_std", 0.0), 5),
                "snr_by_frame"   : [round(s, 2) if s is not None else None
                                     for s in c.get("snrs", [])],
            },
            "photometry"        : {
                "flux_by_frame": [round(b, 1) for b in c["fluxes"]],
                "flux_cv": round(c["flux_cv"], 4)
                           if c["flux_cv"] != float("inf") else None,
                "pointedness": round(c["pointlike"], 5),
                "pointedness_std": round(c.get("pont_std", 0.0), 5),
                "snr_by_frame": [round(s, 2) if s is not None else None
                                 for s in c.get("snrs", [])],
            },

            "astrometric_uncertainty": c.get("astrometric_uncertainty", {}),
            "gaia_static"           : c.get("gaia_static", {
                "status": "nao_avaliado",
                "n_matches": 0,
                "matches": [],
            }),

            "morphology"        : {
                "mean_elongation": round(c.get("mean_elongation", 1.0), 4),
                "mean_fwhm_px"   : c.get("mean_fwhm_px"),
                "by_frame"       : [
                    {
                        "elongation" : round(m["elongation"], 4),
                        "pointlike"    : round(m["pointlike"], 5),
                        "fwhm_est_px": round(float(m["fwhm_est_px"]), 4)
                                       if m["fwhm_est_px"] is not None else None,
                        "psf_model" : m.get("psf_model"),
                        "valid"     : m["valid"],
                    }
                    for m in c.get("morphology_metrics", [])
                ],
            },
            # Trilha completa frame a frame
            "track"             : _serialize_track(c, frames),

            # Status MPC normalizado
            "mpc"               : c.get("mpc", {"status": "not_queried"}),

            # Colunas de avaliação manual — preencher após inspeção no Astrometrica/IASC
            "manual_validation"  : {
                "medido_astrometrica"  : None,
                "entrou_mpc"          : None,
                "feedback_iasc"       : None,
                "classificacao_manual": None,
                "observacoes"         : None,
            },
            "manual_validation" : {
                "measured_in_astrometrica"      : None,
                "included_in_astrometrica_mpc" : None,
                "iasc_feedback"                : None,
                "manual_classification"        : None,
                "notes"                        : None,
            },
        }
        output["candidates"].append(entry)
        output["candidates"].append(entry)

    json_file = output_dir / f"{input_set_name}_candidates.json"
    with open(json_file, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    return json_file


# ─────────────────────────────────────────────
# 11. PIPELINE PRINCIPAL
# ─────────────────────────────────────────────
def _resolve_observer(args) -> dict:
    name = (
        args.observer
        or os.environ.get("LIA_OBS")
        or "IASC Observer"
    )
    email = (
        args.email
        or os.environ.get("LIA_EMAIL")
        or "observer@example.com"
    )
    return {"name": name, "email": email}


def _deduplicate(candidates: list[dict], radius_px: float = None) -> list[dict]:
    if radius_px is None:
        radius_px = T.DEDUP_RADIUS_PX
    unique_candidates = []
    for c in candidates:
        y0, x0 = c["track"][0][0], c["track"][0][1]
        dup = any(
            np.hypot(x0 - u["track"][0][1], y0 - u["track"][0][0]) < radius_px
            for u in unique_candidates
        )
        if not dup:
            unique_candidates.append(c)
    return unique_candidates


def main():
    parser = argparse.ArgumentParser(
        description="Lia — pre-screen moving-object candidates in IASC FITS sequences"
    )
    parser.add_argument("--images",     type=str, default="images")
    parser.add_argument("--output",     type=str, default="results")
    parser.add_argument("--sigma",      type=float, default=5.5)
    parser.add_argument("--no-mpc", dest="no_mpc", action="store_true")
    parser.add_argument("--observer", type=str, default=None)
    parser.add_argument("--email",      type=str, default=None)
    parser.add_argument(
        "--wcs-mode",
        choices=["gaia", "header"],
        default="gaia",
        dest="wcs_mode",
        help=(
            "WCS refinement mode. "
            "'gaia' (default): refines against Gaia DR3. "
            "'header': uses only o WCS/header Pan-STARRS construído, "
            "without querying Gaia. All other parameters (Background2D, "
            "PSF, score, filters, thresholds) remain identical."
        ),
    )
    args = parser.parse_args()

    t_inicio = time.monotonic()
    run_id   = str(uuid.uuid4())

    images_dir = Path(args.images)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    initial_files = sorted(images_dir.glob("*.fits")) + sorted(images_dir.glob("*.fit"))
    input_set_name = (initial_files[0].name.split("_")[0].split(".")[0]
                     if initial_files else "pipeline")

    observer = _resolve_observer(args)

    print(f"\n{Color.BOLD}{'='*60}")
    print(f"  {PROJECT_NAME} v{PIPELINE_VERSION}")
    print(f"  Observer   : {observer['name']}  <{observer['email']}>")
    print(f"  WCS mode   : {args.wcs_mode.upper()}")
    print(f"  run_id     : {run_id}")
    print(f"{'='*60}{Color.RESET}\n")

    log_file = configure_file_logger(output_dir, input_set_name)
    log(f"Pipeline log: {log_file}")
    log(f"run_id: {run_id}  |  wcs_mode: {args.wcs_mode}")
    log("Carregando images FITS...")
    frames = load_fits(images_dir, wcs_mode=args.wcs_mode)
    input_set_name = frames[0]["file"].split("_")[0].split(".")[0]

    log("\nDetectando objects em motion...")
    movers, global_metrics = find_movers(frames, sigma=args.sigma)
    log(f"Raw candidates found: {len(movers)}", "OK")

    if not movers:
        log("No moving object detected.", "WARN")
        sys.exit(0)

    movers = _deduplicate(movers)
    log(f"After deduplication: {len(movers)} unique candidates", "OK")

    log("\nAnalisando candidates...")
    candidates = [analyze_candidate(m, frames) for m in movers]
    candidates.sort(key=lambda x: -x["score"])

    # Atribuir rank e candidate_id rastreável após ordenação
    mode_label = args.wcs_mode.upper()
    for i, c in enumerate(candidates):
        c["rank"]         = i + 1
        c["candidate_id"] = f"{input_set_name}_{mode_label}_C{i+1:03d}"

    if not args.no_mpc:
        log("\nConsultando SkyBot/IMCCE...")
        for c in candidates:
            if c["triage_class"] == "DISCARDED":
                c["mpc"] = {"status": "not_queried", "reason": None, "match": None}
                continue
            if not frames[0]["wcs_ok"]:
                c["mpc"] = {"status": "query_failed",
                            "reason": "WCS inválido", "match": None}
                continue
            tx, ty  = c["track"][0][1], c["track"][0][0]
            ra, dec = pixel_to_radec(tx, ty, frames[0]["wcs"])
            log(f"  #{c['rank']}  RA={ra:.4f} Dec={dec:.4f}...")
            c["mpc"] = query_position_catalog(ra, dec, frames[0]["date_obs"])
            status   = c["mpc"]["status"]
            if status == "likely_match":
                log(f"  → match: {c['mpc']['match']['name_designation']}", "WARN")
            elif status == "ambiguous_match":
                log(f"  → ambiguous ({c['mpc']['n_objects_in_cone']} objects)", "WARN")
            elif status == "no_match":
                log("  → no match", "OK")
            else:
                log(f"  → {status}: {c['mpc'].get('reason', '')}", "WARN")
    else:
        for c in candidates:
            c["mpc"] = {"status": "not_queried", "reason": None, "match": None}

    # Adicionar flags MPC após consulta
    for c in candidates:
        _add_mpc_flag(c)

    # Completar métricas globais
    t_fim = time.monotonic()
    n_rejected_false_positive = sum(
        1 for c in candidates if "REJECTED_FP" in c.get("flags", [])
    )
    global_metrics.update({
        "n_valid_frames"         : sum(1 for f in frames if f["wcs_ok"]),
        "n_unique_candidates"      : len(movers),
        "n_final_candidates"      : len(candidates),
        "n_rejected_false_positive"          : n_rejected_false_positive,
        "class_distribution"             : {
            "STRONG"   : sum(1 for c in candidates if c["triage_class"] == "STRONG"),
            "MODERATE": sum(1 for c in candidates if c["triage_class"] == "MODERATE"),
            "WEAK"   : sum(1 for c in candidates if c["triage_class"] == "WEAK"),
            "DISCARDED": sum(1 for c in candidates if c["triage_class"] == "DISCARDED"),
        },
        "top1_candidate_id"        : candidates[0].get("candidate_id") if candidates else None,
        "top3_candidate_ids"       : [c.get("candidate_id") for c in candidates[:3]],
        "top5_candidate_ids"       : [c.get("candidate_id") for c in candidates[:5]],
        "wcs_mode"                 : args.wcs_mode,
        "run_id"                   : run_id,
        "execution_time_s"         : round(t_fim - t_inicio, 1),
    })

    run_metadata = {
        "run_id"              : run_id,
        "pipeline_name"       : PROJECT_NAME,
        "pipeline_version"    : PIPELINE_VERSION,
        "repository"          : REPOSITORY_NAME,
        "wcs_mode"            : args.wcs_mode,
        "sigma"               : args.sigma,
        "timestamp_execution_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "timestamp_execution_utc"  : datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input_set"           : str(images_dir.resolve()),
        "input_set_name"      : input_set_name,
        "input_set_name"       : input_set_name,
        "observer"          : observer,
    }

    # Resumo terminal
    score_max_str = str(candidates[0].get("score_max", 12)) if candidates else "12"
    print(f"\n{Color.BOLD}{'─'*62}")
    print("  RESUMO DOS CANDIDATOS")
    print(f"{'─'*62}{Color.RESET}")
    print(f"  {'#':>3}  {'Score':>7}  {'Norm':>5}  {'Class':<10}  "
          f"{'MPC':<28}  Flags")
    print(f"  {'─'*3}  {'─'*7}  {'─'*5}  {'─'*10}  {'─'*28}  {'─'*24}")

    _flag_abrev = {
        "LINEARIDADE_BOA"       : "LIN✓",
        "LINEARIDADE_RUIM"      : "LIN✗",
        "VELOCIDADE_CONSISTENTE": "VEL✓",
        "VELOCIDADE_IRREGULAR"  : "VEL✗",
        "BRILHO_ESTAVEL"        : "PHO✓",
        "BRILHO_INSTAVEL"       : "PHO✗",
        "MORFOLOGIA_PONTUAL"    : "MOR✓",
        "MORFOLOGIA_ESTENDIDA"  : "MOR✗",
        "ELONGACAO_ALTA"        : "ELG✗",
        "MORFO_INCONSISTENTE"   : "MCI✗",
        "MPC_SEM_MATCH"         : "noMPC",
        "MPC_MATCH_PROVAVEL"    : "MPC!",
        "MPC_MATCH_AMBIGUO"     : "MPC?",
        "MPC_CONSULTA_FALHOU"   : "MPCx",
        "EDGE_FRAME"            : "EDGE",
        "REJECTED_FP"          : "FP✗",
    }

    for c in candidates:
        if c["triage_class"] == "DISCARDED":
            continue
        mpc     = c.get("mpc", {})
        mpc_str = (mpc.get("match", {}) or {}).get("name_designation") \
                  or mpc.get("status", "?")
        flags_short = " ".join(
            _flag_abrev.get(f, f) for f in c.get("flags", [])
        )
        sm = c.get("score_max", 12)
        print(
            f"  {Color.BOLD}{c['rank']:>3}{Color.RESET}  "
            f"{c['cli_color']}{c['score']:>4}/{sm}  "
            f"{c['score_percent']:>4}%  "
            f"{c['triage_class']:<10}{Color.RESET}  "
            f"{mpc_str[:28]:<28}  "
            f"{flags_short}"
        )

    # Gerar outputs
    log("\nGenerating text report...")
    text_file = generate_text_report(candidates, frames, input_set_name,
                                  output_dir, global_metrics)
    log(f"Report saved: {text_file}", "OK")

    log("Generating MPC report...")
    mpc_file = generate_mpc_report(candidates, frames, input_set_name,
                                  output_dir, observer)
    log(f"MPC report saved: {mpc_file}", "OK")

    log("Generating visualization...")
    generate_visualization(candidates, frames, output_dir, input_set_name)

    log("Exportando JSON...")
    json_file = export_json(candidates, frames, input_set_name,
                             output_dir, global_metrics, run_metadata)
    log(f"JSON saved: {json_file}", "OK")

    print(f"\n{Color.GREEN}{Color.BOLD}Analysis complete. "
          f"Resultados em: {output_dir}/{Color.RESET}")
    print(f"   → {text_file.name}")
    print(f"   → {mpc_file.name}")
    print(f"   → {input_set_name}_candidates.png")
    print(f"   → {json_file.name}")
    print(f"   → {log_file.name}")
    print()


if __name__ == "__main__":
    main()
