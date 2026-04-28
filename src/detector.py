"""
=============================================================================
TRITON — detector.py
=============================================================================
Pipeline de pré-triagem de asteroides em imagens FITS do IASC/Pan-STARRS.

Autora  : Jaciana Barbosa
Versão  : 1.4.2
Licença : MIT

Descrição:
    Detecta objetos em movimento em 4 frames FITS, calcula coordenadas
    astrométricas (RA/Dec), consulta catálogos de objetos conhecidos via
    SkyBot/IMCCE para verificar correspondência posicional próxima, gera
    score de probabilidade e produz relatório no formato MPC para envio
    ao IASC após validação no Astrometrica.

    Esta ferramenta é um auxiliar de pré-triagem. Não substitui a
    validação astrométrica no software Astrometrica nem o julgamento humano.

Uso:
    python detector.py --imagens ../imagens/ --output ../resultados/

Dependências:
    pip install -r requirements.txt

Changelog v1.4.2:
    - WCS Pan-STARRS preserva CROTA1/2 quando usa CDELT
    - WCS bruto pós-offset também preserva CROTA antes da passada fina Gaia
    - Diferenças de RA no cross-match/RMS usam wrap 0/360
    - Offset bruto usa votação 2D de translação e pareamento 1:1 no pico
    - Ajuste afim Gaia usa sigma-clipping para remover pares incoerentes
    - Log de arquivo é configurado antes do carregamento/refinamento Gaia

Changelog v1.4.1:
    - Cross-match Gaia em duas passadas para contornar erro inicial de 5-7 arcsec do WCS
      Pan-STARRS: passada bruta (15 arcsec) estima offset de translação e corrige CRVAL;
      passada fina (2 arcsec) com WCS pré-corrigido realiza ajuste afim completo
    - Mesmo que a passada fina falhe, o offset bruto é preservado no WCS (status
      "gaia_offset_bruto_apenas"), melhorando o alinhamento entre frames
    - Métricas por frame incluem n_matches_bruto, offset_bruto_arcsec, n_matches_fino
    - Nova função _estimar_offset_grosseiro isolada para testabilidade

Changelog v1.4.0:
    - Refinamento astrométrico via Gaia DR3 por frame
    - WCS Pan-STARRS construído primeiro, refinado contra Gaia depois
    - Movimento próprio Gaia (epoch 2016.0) propagado para data da observação
    - Ajuste afim (CD matrix) minimizando resíduos cross-match
    - Fallback transparente: mantém WCS Pan-STARRS se Gaia falhar
    - Removida heurística cinemática em pixels relativos
      (substituída por astrometria absoluta Gaia)
    - Status do refinamento por frame em métricas globais

Changelog v1.3.2:
    - Deriva instrumental estimada por par de frames (1→2, 2→3, 3→4)
      ANTES da validação cinemática
    - Trilha corrigida pela deriva acumulada usada como input da
      validação cinemática; trilha original preservada para
      fotometria e relatório MPC
    - Aviso explícito quando frame específico tem deslocamento
      sistemático (problema de pointing, tracking ou WCS)
    - Ordem das operações reorganizada: deriva → correção → validação

Changelog v1.3.1:
    - Mascaramento de pixels saturados/inválidos antes do
      Background2D, DAOStarFinder e ajuste PSF
    - Aviso de qualidade por frame; aborta se saturação >25%
    - Métricas globais incluem fração de saturação por frame

Changelog v1.3.0:
    - Associação entre frames: rejeição de trilhas com saltos incoerentes
      (ângulo de direção, razão de passo, aceleração excessiva)
    - Rejeição de falsos positivos: borda rejeita (não apenas flag),
      hot pixel heurístico, SNR mínimo, brilho caótico, cinemática implausível
    - Score recalibrado: limiares centralizados em T.*, dois novos componentes
      (elongation e consistência morfológica), score máximo agora 12
    - Morfologia: elongation via PCA dos pixels, FWHM estimada, consistência
      entre frames adicionados à análise
    - Logs de decisão explícitos: cada candidato registra razões de penalização
      e rejeição rastreáveis em flags e no campo "razoes_decisao" do JSON
    - Thresholds centralizados na classe T para tuning sem caçar valores
    - Compatível com saídas anteriores (JSON, TXT, MPC, PNG)

Changelog v1.2.0:
    - Score decomposto por componente exportado no JSON e TXT
    - Trilha completa frame a frame serializada no JSON
    - Flags diagnósticas por candidato (JSON + TXT)
    - Status MPC normalizado com schema estruturado
    - Métricas globais do conjunto no JSON e TXT
    - TXT reorganizado com seção de priorização operacional
    - PNGs melhorados: rank, score, flags, vetor de trajetória
    - README revisado com posicionamento técnico
=============================================================================
"""

import os
import sys
import json
import time
import argparse
import logging
import warnings
import multiprocessing as mp
from datetime import datetime
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

VERSAO = "1.4.2"

# ─────────────────────────────────────────────
# THRESHOLDS CENTRALIZADOS
# Edite aqui para tuning futuro sem caçar valores pelo código.
# ─────────────────────────────────────────────
class T:
    # Detecção
    SIGMA_DETECCAO       = 5.5       # sigma local de detecção (sobrescrito por --sigma)
    FWHM_PX              = 3.0       # FWHM assumido para o finder
    BACKGROUND_BOX       = (50, 50)  # malha local do céu para Background2D
    BACKGROUND_FILTER    = (3, 3)    # suaviza a malha sem apagar gradientes reais
    PSF_FIT_RAIO_PX      = 8         # janela local para ajuste PSF sub-pixel
    PSF_MODELO           = "moffat"  # Moffat modela melhor asas de seeing atmosférico
    PSF_MAX_ITER         = 100
    FWHM_MIN_PX          = 1.2       # abaixo disso tende a hot pixel/raio cósmico
    FWHM_MAX_PX          = 8.0       # acima disso tende a galáxia, blend ou trail
    SHARPNESS_MIN        = 0.2
    SHARPNESS_MAX        = 1.0
    ROUNDNESS_MAX        = 0.8       # |roundness| máximo no finder

    # Associação entre frames
    RAIO_MATCH_PX        = 20.0      # janela de casamento mútuo [px]

    # Filtro de movimento (pós-deriva)
    MOVE_MIN_PX          = 2.0       # deslocamento mínimo total [px]
    MOVE_MAX_PX          = 120.0     # deslocamento máximo total [px]
    RESIDUO_MIN_PX       = 1.8       # resíduo mínimo em relação à deriva

    # Refinamento WCS via Gaia DR3
    GAIA_RAIO_QUERY_ARCMIN  = 6.0    # raio do cone search em torno do CRVAL
    GAIA_MAG_LIM            = 19.0   # G_mag máximo para evitar fontes faint demais
    GAIA_MAG_MIN            = 12.0   # evitar saturadas (Pan-STARRS satura cedo)
    GAIA_MIN_MATCHES        = 8      # estrelas casadas mínimas para refinar (passada fina)
    GAIA_RAIO_MATCH_ARCSEC  = 2.0    # tolerância de cross-match por estrela (passada fina)
    GAIA_DIFF_MIN_ARCSEC    = 0.3    # só tenta refinar se RMS_pre > este valor
    GAIA_DIFF_MAX_ARCSEC    = 0.9    # aceita resultado Gaia se RMS_pos < este (bruto ou afim)
    GAIA_TIMEOUT_S          = 30     # timeout da query
    GAIA_RAIO_BRUTO_ARCSEC  = 15.0   # raio da passada bruta (WCS Pan-STARRS tem ~5-7 arcsec erro)
    GAIA_MIN_MATCHES_BRUTO  = 5      # mínimo de matches para estimar offset grosseiro
    GAIA_OFFSET_BIN_ARCSEC  = 1.5     # bin da votação 2D de offset bruto
    GAIA_RAIO_REFINO_ARCSEC = 2.5     # cluster em torno do pico da votação bruta

    # Rejeição de falso positivo
    MARGEM_BORDA_PX      = 30.0      # distância mínima à borda do frame
    HOT_PIXEL_PONT_MAX   = 0.70      # pontualidade > threshold → hot pixel suspect
    HOT_PIXEL_SNR_MIN    = 3.0       # SNR mínimo para não ser descartado como hot pixel
    BRILHO_CV_CAOS       = 1.5       # CV acima disto → brilho caótico (rejeição forte)
    SNR_MINIMO           = 2.0       # SNR mínimo em pelo menos 3 dos 4 frames

    # Score: linearidade (0–3)
    LIN_EXCELENTE        = 0.8       # residuo < → 3 pts
    LIN_BOA              = 1.5       # residuo < → 2 pts
    LIN_MARGINAL         = 2.5       # residuo < → 1 pt; senão 0

    # Score: velocidade/consistência (0–2)
    VEL_UNIFORME         = 3.0       # σ passos < → 2 pts
    VEL_MODERADA         = 5.0       # σ passos < → 1 pt; senão 0

    # Score: fotometria / CV (0–2)
    FOT_ESTAVEL          = 0.25      # CV < → 2 pts
    FOT_MODERADA         = 0.50      # CV < → 1 pt; senão 0  (era 0.45)

    # Score: morfologia / pontualidade (0–2)
    MOR_PONTUAL          = 0.12      # pont > → 2 pts  (era 0.10)
    MOR_MARGINAL         = 0.06      # pont > → 1 pt; senão 0  (era 0.05)

    # Score: elongation (0–1, novo componente)
    ELON_BOA             = 1.6       # elong < → 1 pt; senão 0

    # Score: consistência morfológica entre frames (0–1, novo)
    MOR_CONSIST_MAX_STD  = 0.08      # std(pontual_por_frame) < → 1 pt

    # Score: faixa de velocidade típica de MBA (0–1)
    FAIXA_VEL_MIN_PX     = 4.0
    FAIXA_VEL_MAX_PX     = 50.0

    # Classificação final
    SCORE_FORTE          = 8         # score >= → FORTE
    SCORE_MODERADO       = 6         # score >= → MODERADO
    SCORE_FRACO          = 4         # score >= → FRACO; senão DESCARTA

    # Mascaramento de pixels inválidos (Pan-STARRS marca saturados com 65535)
    PIXEL_SAT_VALOR      = 65535      # valor que indica saturação
    PIXEL_SAT_TOLERANCIA = 35         # margem (mascara também 65500..65535)
    PIXEL_MIN_VALIDO     = 1          # abaixo disso é "buraco" / pixel ruim
    FRAME_SAT_AVISO_PCT  = 5.0        # >% saturados → avisa qualidade ruim
    FRAME_SAT_REJEITA_PCT = 25.0      # >% saturados → aborta o pipeline

    # Deduplicação
    DEDUP_RAIO_PX        = 25.0

    # MPC: score mínimo para entrar no relatório 80-colunas
    MPC_SCORE_MIN        = 6

    # Deriva: threshold para fonte "estável"
    DERIVA_ESTAVEL_PX    = 2.0
    DERIVA_MIN_ESTAVEIS  = 5

    # Astrometria: incerteza e rejeição de fonte estática Gaia
    GAIA_STATIC_MATCH_ARCSEC = 1.5
    GAIA_STATIC_MIN_FRAMES   = 3


# ─────────────────────────────────────────────
# Cores para terminal + logging
# ─────────────────────────────────────────────
class Cor:
    VERDE    = "\033[92m"
    AMARELO  = "\033[93m"
    VERMELHO = "\033[91m"
    AZUL     = "\033[94m"
    NEGRITO  = "\033[1m"
    RESET    = "\033[0m"

_logger = logging.getLogger("TRITON")

def log(msg, nivel="INFO"):
    cores = {"INFO": Cor.AZUL, "OK": Cor.VERDE,
             "WARN": Cor.AMARELO, "ERRO": Cor.VERMELHO}
    cor = cores.get(nivel, Cor.RESET)
    print(f"{cor}[{nivel}]{Cor.RESET} {msg}")
    nivel_py = {"INFO": logging.INFO, "OK": logging.INFO,
                "WARN": logging.WARNING, "ERRO": logging.ERROR}.get(nivel, logging.INFO)
    _logger.log(nivel_py, msg)


def configurar_log_arquivo(pasta_output: Path, nome_conjunto: str) -> Path:
    arq_log = pasta_output / f"{nome_conjunto}_pipeline.log"
    handler = logging.FileHandler(arq_log, mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    _logger.setLevel(logging.INFO)
    for h in list(_logger.handlers):
        _logger.removeHandler(h)
    _logger.addHandler(handler)
    return arq_log


# ─────────────────────────────────────────────
# 1. CARREGAMENTO DAS IMAGENS
# ─────────────────────────────────────────────
def _construir_wcs_panstarrs(header: fits.Header) -> WCS:
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
        raise ValueError("Header FITS não contém chaves WCS mínimas "
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

def _consultar_gaia_dr3(ra_centro: float, dec_centro: float,
                        raio_arcmin: float = None,
                        mag_lim: float = None,
                        mag_min: float = None):
    """
    Consulta o catálogo Gaia DR3 via TAP do ESA.

    Física/astrometria:
        Gaia DR3 (epoch 2016.0) é o catálogo astrométrico mais preciso
        disponível, com erros típicos < 1 mas em RA/Dec. Limites de
        magnitude evitam estrelas saturadas no Pan-STARRS (G < 12) e
        fontes faint cuja detecção falha no detector (G > 19).
        RUWE < 1.4 filtra fontes com astrometria degradada (binárias
        não resolvidas, extensões reais).

    Retorna astropy.Table com colunas ra, dec, pmra, pmdec,
    phot_g_mean_mag, ref_epoch. Retorna None se falhar ou < GAIA_MIN_MATCHES.
    """
    try:
        from astroquery.gaia import Gaia
        Gaia.ROW_LIMIT = 5000
        if raio_arcmin is None:
            raio_arcmin = T.GAIA_RAIO_QUERY_ARCMIN
        if mag_lim is None:
            mag_lim = T.GAIA_MAG_LIM
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


def _aplicar_movimento_proprio_gaia(tabela_gaia, data_obs_jd: float):
    """
    Propaga posições Gaia DR3 (epoch 2016.0) para a época da observação.

    Física:
        Estrelas de alto movimento próprio podem se deslocar dezenas de
        mas em poucos anos. Sem corrigir, o cross-match falha exatamente
        para as estrelas mais brilhantes e próximas — as mais úteis para
        calibrar o WCS porque têm menor incerteza astrométrica.
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


def _delta_ra_graus(ra_a, ra_b):
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
        é válida para campos de ~10 arcmin onde o erro de curvatura é
        < 0.01 arcsec — muito abaixo da precisão alvo.
    """
    from scipy.spatial import cKDTree
    if raio_arcsec is None:
        raio_arcsec = T.GAIA_RAIO_MATCH_ARCSEC

    ra_med  = float(np.median(gaia_radec[:, 0]))
    dec_med = float(np.median(gaia_radec[:, 1]))
    cosdec  = np.cos(np.radians(dec_med))

    det_xy = np.column_stack([
        _delta_ra_graus(deteccoes_radec[:, 0], ra_med) * cosdec * 3600.0,
        (deteccoes_radec[:, 1] - dec_med) * 3600.0,
    ])
    gaia_xy = np.column_stack([
        _delta_ra_graus(gaia_radec[:, 0], ra_med) * cosdec * 3600.0,
        (gaia_radec[:, 1] - dec_med) * 3600.0,
    ])

    tree = cKDTree(gaia_xy)
    dist, idx = tree.query(det_xy, k=1, distance_upper_bound=raio_arcsec)
    valido = np.isfinite(dist) & (dist < raio_arcsec)
    return np.where(valido)[0], idx[valido]


def _estimar_offset_grosseiro(
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
        O WCS Pan-STARRS pode ter erro sistemático de vários arcsec. Em campo
        denso, vizinho mais próximo direto é enviesado: a estrela Gaia errada
        pode estar mais perto que a correspondente verdadeira. Por isso esta
        etapa usa votação 2D de todos os pares dentro do raio bruto e escolhe
        o pico de translação coerente antes de calcular a mediana robusta.
    """
    from scipy.spatial import cKDTree
    if raio_arcsec is None:
        raio_arcsec = T.GAIA_RAIO_BRUTO_ARCSEC

    ra_med  = float(np.median(gaia_radec[:, 0]))
    dec_med = float(np.median(gaia_radec[:, 1]))
    cosdec  = np.cos(np.radians(dec_med))

    det_xy = np.column_stack([
        _delta_ra_graus(deteccoes_radec[:, 0], ra_med) * cosdec * 3600.0,
        (deteccoes_radec[:, 1] - dec_med) * 3600.0,
    ])
    gaia_xy = np.column_stack([
        _delta_ra_graus(gaia_radec[:, 0], ra_med) * cosdec * 3600.0,
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

    raio_refino = float(T.GAIA_RAIO_REFINO_ARCSEC)
    dist_peak = np.hypot(deltas[:, 0] - centro[0], deltas[:, 1] - centro[1])
    cluster = dist_peak <= raio_refino

    # Garante pareamento 1:1 dentro do pico: para cada detecção, mantém o par
    # mais próximo do centro do pico. Isso reduz blends e duplicatas Gaia.
    candidatos = np.where(cluster)[0]
    if candidatos.size:
        ordem = candidatos[np.argsort(dist_peak[candidatos])]
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
    if n_val < T.GAIA_MIN_MATCHES_BRUTO:
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


def _refinar_wcs_gaia(frame: dict, sigma_deteccao: float = None) -> dict:
    """
    Refina o WCS de um frame contra Gaia DR3 em duas passadas.

    Passada bruta (raio 15 arcsec):
        Estima offset puro de translação — suficiente para cobrir o erro
        sistemático inicial de 5-7 arcsec do WCS Pan-STARRS. Aplica o
        offset ao CRVAL antes da passada fina.

    Passada fina (raio 2 arcsec com WCS pré-corrigido):
        Cross-match preciso + ajuste afim (CD matrix) minimizando resíduos.

    Se a passada fina falhar, o offset bruto é preservado no WCS com status
    "gaia_offset_bruto_apenas" — melhor do que descartar a correção inteira.

    Status possíveis:
        "gaia_refinado"              — ajuste afim completo aceito
        "gaia_offset_bruto_apenas"   — só translação; passada fina falhou
        "gaia_skipped_pre_baixo"     — WCS Pan-STARRS já era bom
        "gaia_falhou_rede"           — Gaia não respondeu
        "gaia_falhou_match"          — estrelas insuficientes no campo
        "gaia_falhou_pos_alto"       — ajuste não convergiu para < limiar
        "wcs_invalido"               — frame não tem WCS válido
    """
    if not frame.get("wcs_ok"):
        frame["wcs_status"]          = "wcs_invalido"
        frame["wcs_rms_pre_arcsec"]  = None
        frame["wcs_rms_pos_arcsec"]  = None
        frame["gaia_n_matches"]      = 0
        frame["gaia_n_matches_bruto"] = 0
        frame["gaia_n_matches_fino"]  = 0
        frame["gaia_offset_bruto_arcsec"] = [0.0, 0.0]
        return frame

    wcs_original = frame["wcs"]
    if sigma_deteccao is None:
        sigma_deteccao = T.SIGMA_DETECCAO

    # ── Detectar fontes para cross-match ──────────────────────────────────
    data = np.asarray(frame["data"], dtype=np.float64)
    mascara = _mascara_pixels_invalidos(data)
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
    if tabela_det is None or len(tabela_det) < T.GAIA_MIN_MATCHES_BRUTO:
        frame["wcs_status"]               = "gaia_falhou_match"
        frame["wcs_rms_pre_arcsec"]       = None
        frame["wcs_rms_pos_arcsec"]       = None
        frame["gaia_n_matches"]           = 0
        frame["gaia_n_matches_bruto"]     = 0
        frame["gaia_n_matches_fino"]      = 0
        frame["gaia_offset_bruto_arcsec"] = [0.0, 0.0]
        return frame

    x_col = "x_centroid" if "x_centroid" in tabela_det.colnames else "xcentroid"
    y_col = "y_centroid" if "y_centroid" in tabela_det.colnames else "ycentroid"
    px_det = np.column_stack([tabela_det[x_col], tabela_det[y_col]])

    try:
        sky = wcs_original.all_pix2world(px_det, 0)
        det_radec = sky.astype(np.float64)
    except Exception:
        frame["wcs_status"]               = "gaia_falhou_match"
        frame["wcs_rms_pre_arcsec"]       = None
        frame["wcs_rms_pos_arcsec"]       = None
        frame["gaia_n_matches"]           = 0
        frame["gaia_n_matches_bruto"]     = 0
        frame["gaia_n_matches_fino"]      = 0
        frame["gaia_offset_bruto_arcsec"] = [0.0, 0.0]
        return frame

    # ── Consultar Gaia ────────────────────────────────────────────────────
    ra_c  = float(wcs_original.wcs.crval[0])
    dec_c = float(wcs_original.wcs.crval[1])
    tabela_gaia = _consultar_gaia_dr3(ra_c, dec_c)
    if tabela_gaia is None:
        frame["wcs_status"]               = "gaia_falhou_rede"
        frame["wcs_rms_pre_arcsec"]       = None
        frame["wcs_rms_pos_arcsec"]       = None
        frame["gaia_n_matches"]           = 0
        frame["gaia_n_matches_bruto"]     = 0
        frame["gaia_n_matches_fino"]      = 0
        frame["gaia_offset_bruto_arcsec"] = [0.0, 0.0]
        return frame

    tabela_gaia = _aplicar_movimento_proprio_gaia(tabela_gaia, frame["jd"])
    frame["gaia_catalogo_refinado"] = tabela_gaia
    gaia_radec  = np.column_stack([
        np.array(tabela_gaia["ra"],  dtype=np.float64),
        np.array(tabela_gaia["dec"], dtype=np.float64),
    ])

    # ── PASSADA BRUTA: offset grosseiro de translação ─────────────────────
    off_ra, off_dec, n_bruto, idx_det_bruto, idx_gaia_bruto = _estimar_offset_grosseiro(
        det_radec, gaia_radec, retornar_indices=True
    )
    frame["gaia_n_matches_bruto"]     = n_bruto
    frame["gaia_offset_bruto_arcsec"] = [round(off_ra, 3), round(off_dec, 3)]

    log(f"  Gaia bruto [{frame['arquivo']}]: N={n_bruto}, "
        f"offset=({off_ra:+.2f}, {off_dec:+.2f}) arcsec", "INFO")

    if n_bruto < T.GAIA_MIN_MATCHES_BRUTO:
        frame["wcs_status"]         = "gaia_falhou_match"
        frame["wcs_rms_pre_arcsec"] = None
        frame["wcs_rms_pos_arcsec"] = None
        frame["gaia_n_matches"]     = 0
        frame["gaia_n_matches_fino"] = 0
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
        off_ra2, off_dec2, n_iter2 = _estimar_offset_grosseiro(
            det_rd_iter, gaia_radec, raio_arcsec=T.GAIA_RAIO_MATCH_ARCSEC * 3
        )
        if n_iter2 >= T.GAIA_MIN_MATCHES_BRUTO and (abs(off_ra2) > 0.05 or abs(off_dec2) > 0.05):
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
        frame["wcs_status"]       = "gaia_offset_bruto_apenas"
        frame["wcs_rms_pre_arcsec"] = None
        frame["wcs_rms_pos_arcsec"] = None
        frame["gaia_n_matches"]   = n_bruto
        frame["gaia_n_matches_fino"] = 0
        return frame

    # ── PASSADA FINA: cross-match estreito + ajuste afim ─────────────────
    idx_det, idx_gaia = _cross_match_gaia(det_radec_corr, gaia_radec)
    n_fino = len(idx_det)
    frame["gaia_n_matches_fino"] = n_fino

    log(f"  Gaia fino  [{frame['arquivo']}]: N={n_fino}", "INFO")

    if n_fino >= T.GAIA_MIN_MATCHES:
        idx_det_ajuste = idx_det
        idx_gaia_ajuste = idx_gaia
        match_origem = "fino"
    elif n_bruto >= T.GAIA_MIN_MATCHES:
        # Se a passada fina ainda fica curta, os pares do pico bruto são mais
        # informativos que descartar tudo: eles já representam uma translação
        # coerente do campo e permitem estimar/validar a solução afim.
        idx_det_ajuste = idx_det_bruto
        idx_gaia_ajuste = idx_gaia_bruto
        match_origem = "bruto"
        log(f"  Gaia ajuste[{frame['arquivo']}]: usando {len(idx_det_ajuste)} "
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
        d_ra  = _delta_ra_graus(det_match[:, 0], gaia_match[:, 0]) * cosdec_med * 3600.0
        d_dec = (det_match[:, 1] - gaia_match[:, 1]) * 3600.0
        rms_pre = float(np.sqrt(np.mean(d_ra**2 + d_dec**2)))
        frame["wcs_rms_pre_arcsec"] = round(rms_pre, 4)
        frame["gaia_n_matches"]     = n_ajuste

        if rms_pre <= T.GAIA_DIFF_MIN_ARCSEC:
            frame["wcs"]              = wcs_bruto
            frame["wcs_status"]       = "gaia_skipped_pre_baixo"
            frame["wcs_rms_pos_arcsec"] = round(rms_pre, 4)
            frame["gaia_n_matches_usados"] = n_ajuste
            return frame

        # Se o RMS_pre já está dentro do critério de aceite, o WCS bruto
        # (translação pura) é suficiente — o ajuste afim não tem sinal de
        # rotação/escala para aprender com erros de centroide aleatórios.
        if rms_pre < T.GAIA_DIFF_MAX_ARCSEC:
            frame["wcs"]              = wcs_bruto
            frame["wcs_status"]       = "gaia_refinado"
            frame["wcs_rms_pos_arcsec"] = round(rms_pre, 4)
            frame["gaia_n_matches_usados"] = n_ajuste
            return frame
    else:
        # Passada fina insuficiente — aceita offset bruto como melhor resultado
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_offset_bruto_apenas"
        frame["wcs_rms_pre_arcsec"] = None
        frame["wcs_rms_pos_arcsec"] = None
        frame["gaia_n_matches"]   = n_bruto
        return frame

    # ── Ajuste afim via lstsq com verificação de condicionamento ─────────
    # Física: se os resíduos pós-passada-bruta são uniformes (sem gradiente
    # espacial), o lstsq não tem sinal suficiente para estimar rotação/escala
    # e produz coeficientes instáveis que pioram o WCS. Verificamos o número
    # de condição da matriz A; se > 1e6 (mal condicionada), não há informação
    # suficiente para o ajuste afim — mantemos o WCS da passada bruta.
    px_match    = px_det[idx_det_ajuste]
    ra_med_fit  = float(np.mean(gaia_match[:, 0]))
    dec_med_fit = float(np.mean(gaia_match[:, 1]))
    cosdec_fit  = np.cos(np.radians(dec_med_fit))
    target_x    = _delta_ra_graus(gaia_match[:, 0], ra_med_fit) * cosdec_fit * 3600.0
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
        frame["wcs_status"]       = "gaia_offset_bruto_apenas"
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
        frame["wcs_status"]       = "gaia_offset_bruto_apenas"
        frame["wcs_rms_pos_arcsec"] = None
        return frame

    if coef_x is None or coef_y is None or int(mascara_fit.sum()) < T.GAIA_MIN_MATCHES:
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_offset_bruto_apenas"
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
    # ajusta o deslocamento residual medido pelo termo constante do lstsq.
    # O termo coef_x[2]/coef_y[2] representa o offset residual (em arcsec
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
        d_ra_n   = _delta_ra_graus(sky_novo[:, 0], gaia_match_eval[:, 0]) * cosdec_fit * 3600.0
        d_dec_n  = (sky_novo[:, 1] - gaia_match_eval[:, 1]) * 3600.0
        rms_pos  = float(np.sqrt(np.mean(d_ra_n**2 + d_dec_n**2)))
    except Exception:
        frame["wcs"]              = wcs_bruto
        frame["wcs_status"]       = "gaia_offset_bruto_apenas"
        frame["wcs_rms_pos_arcsec"] = None
        return frame

    frame["wcs_rms_pos_arcsec"] = round(rms_pos, 4)
    log(f"  Gaia RMS   [{frame['arquivo']}]: pre={frame['wcs_rms_pre_arcsec']} "
        f"pos={rms_pos:.4f} arcsec", "INFO")

    if rms_pos < T.GAIA_DIFF_MAX_ARCSEC:
        frame["wcs"]        = wcs_novo
        frame["wcs_status"] = "gaia_refinado"
    else:
        frame["wcs"]        = wcs_bruto
        frame["wcs_status"] = "gaia_offset_bruto_apenas"

    return frame


def _encontrar_hdu_ciencia(hdul: fits.HDUList) -> int:
    for i, h in enumerate(hdul):
        if h.data is not None and h.data.ndim == 2:
            return i
    raise ValueError("Nenhuma HDU com imagem 2D encontrada no FITS.")


def _extrair_tempos_observacao(header: fits.Header, arquivo: str) -> dict:
    """Extrai início e meio da exposição em JD, aceitando DATE-OBS ou MJD-OBS."""
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
            raise ValueError("DATE-OBS/MJD-OBS ausentes")
    except Exception:
        log(
            f"Timestamp inválido ou ausente em {arquivo}: "
            f"DATE-OBS={date_obs!r}, MJD-OBS={mjd_obs!r}",
            "ERRO",
        )
        sys.exit(1)

    t_mid = t_inicio + (exptime / 2.0) * u.second
    return {
        "date_obs_inicio": date_obs_inicio,
        "date_obs_mid": t_mid.isot,
        "jd_inicio": float(t_inicio.jd),
        "jd_mid": float(t_mid.jd),
        "exptime": exptime,
    }


def carregar_fits(pasta: Path) -> list[dict]:
    """Carrega os 4 arquivos FITS da pasta, ordenados por timestamp."""
    arquivos = sorted(pasta.glob("*.fits")) + sorted(pasta.glob("*.fit"))
    if len(arquivos) != 4:
        log(f"Esperava exatamente 4 arquivos FITS, encontrou {len(arquivos)}.", "ERRO")
        sys.exit(1)

    frames = []
    for arq in arquivos:
        with fits.open(arq) as hdul:
            idx = _encontrar_hdu_ciencia(hdul)
            header = hdul[0].header.copy()
            if idx != 0:
                for k, v in hdul[idx].header.items():
                    if k and k not in ("SIMPLE", "BITPIX", "NAXIS",
                                       "EXTEND", "COMMENT", "HISTORY",
                                       "XTENSION", "PCOUNT", "GCOUNT"):
                        header[k] = v
            data = hdul[idx].data.astype(np.float64)

            try:
                wcs = _construir_wcs_panstarrs(header)
                wcs_ok = True
            except Exception as e:
                log(f"WCS inválido em {arq.name}: {e}", "ERRO")
                wcs = None
                wcs_ok = False

            tempos = _extrair_tempos_observacao(header, arq.name)

            frame_dict = {
                "arquivo" : arq.name,
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
                "exptime" : tempos["exptime"],
            }

            if wcs_ok:
                # Calcula fração de pixels inválidos para detectar frames
                # com saturação excessiva antes que contaminem a detecção.
                mascara_frame = _mascara_pixels_invalidos(data)
                sat_pct = float(mascara_frame.mean() * 100.0)
                frame_dict["sat_pct"] = sat_pct
                if sat_pct > T.FRAME_SAT_REJEITA_PCT:
                    log(
                        f"Frame {arq.name} tem {sat_pct:.1f}% de pixels "
                        f"saturados/inválidos (limite: {T.FRAME_SAT_REJEITA_PCT}%). "
                        f"Baixe um frame substituto e reexecute o pipeline.",
                        "ERRO",
                    )
                    sys.exit(1)
                elif sat_pct > T.FRAME_SAT_AVISO_PCT:
                    log(
                        f"Frame {arq.name}: {sat_pct:.1f}% de pixels saturados "
                        f"— qualidade reduzida, mascaramento ativo.",
                        "WARN",
                    )
            else:
                frame_dict["sat_pct"] = 0.0

            frames.append(frame_dict)
        log(f"Carregado: {arq.name}  ({frame_dict['date_obs']})  "
            f"[{data.shape[1]}x{data.shape[0]}]  WCS={'OK' if wcs_ok else 'FALHOU'}")

    if not any(f["wcs_ok"] for f in frames):
        log("Nenhum frame tem WCS válido. Impossível calcular RA/Dec.", "ERRO")
        sys.exit(1)

    frames.sort(key=lambda x: x["jd"])

    # Refinamento astrométrico via Gaia DR3 — etapa separada, após WCS Pan-STARRS.
    # Cada frame é processado independentemente; fallback transparente se Gaia falhar.
    log("Refinando WCS contra Gaia DR3...")
    for i, f in enumerate(frames):
        if f["wcs_ok"]:
            frames[i] = _refinar_wcs_gaia(f)
            status = frames[i].get("wcs_status", "?")
            rms_pre = frames[i].get("wcs_rms_pre_arcsec")
            rms_pos = frames[i].get("wcs_rms_pos_arcsec")
            n_match = frames[i].get("gaia_n_matches", 0)
            pre_str = f"{rms_pre:.3f}" if rms_pre is not None else "N/A"
            pos_str = f"{rms_pos:.3f}" if rms_pos is not None else "N/A"
            log(f"  WCS {frames[i]['arquivo']}: {status} "
                f"(RMS pre={pre_str} arcsec, pos={pos_str} arcsec, N={n_match})")
        else:
            frames[i]["wcs_status"]         = "wcs_invalido"
            frames[i]["wcs_rms_pre_arcsec"] = None
            frames[i]["wcs_rms_pos_arcsec"] = None
            frames[i]["gaia_n_matches"]     = 0

    return frames


# ─────────────────────────────────────────────
# 2. DETECÇÃO DE FONTES EM CADA FRAME
# ─────────────────────────────────────────────
def _mascara_pixels_invalidos(data: np.ndarray) -> np.ndarray:
    """
    Retorna máscara booleana True para pixels que NÃO devem ser usados
    em detecção, fotometria ou estimativa de background.

    Pan-STARRS marca pixels saturados, defeituosos ou off-CCD com 65535
    (limite do uint16) e ocasionalmente com 0 ou valores muito baixos.
    Esses pixels não carregam sinal físico válido; passá-los ao
    Background2D infla a RMS local, e ao DAOStarFinder gera detecções
    espúrias que poluem o casamento entre frames.
    """
    sat_min = float(T.PIXEL_SAT_VALOR - T.PIXEL_SAT_TOLERANCIA)
    return (
        (data >= sat_min)
        | (data < T.PIXEL_MIN_VALIDO)
        | ~np.isfinite(data)
    )


def _estimar_background_local(img: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Estima background e RMS locais com `Background2D`.

    Física/astrometria:
        Imagens reais raramente têm céu plano: vinheta, Lua, halos de estrelas
        brilhantes e gradientes do detector mudam o ruído em escalas de dezenas
        de pixels. Um threshold global superestima fontes em regiões limpas e
        subestima ruído em bordas. A malha 50x50 mede o céu em escala maior que
        a PSF estelar, enquanto o filtro 3x3 suaviza células ruidosas sem apagar
        gradientes lentos. Assim, o critério n-sigma passa a ser local.
    """
    data = np.asarray(img, dtype=np.float64)
    # Pan-STARRS marca saturação com 65535. Sem essa máscara, o Background2D
    # infla a RMS em torno de detritos saturados (ex.: rastro de satélite),
    # elevando o threshold local e produzindo clusters de "fontes" espúrias.
    mascara = _mascara_pixels_invalidos(data)
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


def _recorte_float(img: np.ndarray, y: float, x: float,
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


def _ajustar_psf_subpixel(img_sub: np.ndarray, x0: float, y0: float,
                          fwhm_px: float = 3.0,
                          mascara_invalida: np.ndarray | None = None) -> dict | None:
    """
    Refina uma detecção por ajuste PSF e retorna centro sub-pixel.

    Física/astrometria:
        O centroide simples é uma média ponderada dos pixels e se desloca com
        ruído, background residual, pixels saturados e blends. O ajuste PSF
        compara a fonte a um modelo contínuo, obtendo o centro do perfil óptico
        em coordenadas sub-pixel. Moffat é preferido em seeing ruim porque suas
        asas de potência representam espalhamento atmosférico melhor que uma
        Gaussiana pura; em campos bem comportados a Gaussiana fica como fallback
        estável.
    """
    cut, y_origin, x_origin = _recorte_float(img_sub, y0, x0, T.PSF_FIT_RAIO_PX)
    if cut.size < 16 or min(cut.shape) < 5:
        return None

    # Se a janela de ajuste contém >30% de pixels saturados/inválidos, o modelo
    # PSF seria dominado por artefatos — descarta sem ajustar.
    if mascara_invalida is not None:
        cut_mask, _, _ = _recorte_float(
            mascara_invalida.astype(np.float64), y0, x0, T.PSF_FIT_RAIO_PX
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
    modelo_nome = T.PSF_MODELO.lower()

    try:
        if modelo_nome == "moffat":
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
            modelo_nome = "gaussian"
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
        "modelo": modelo_nome,
    }


def detectar_fontes(img: np.ndarray, sigma: float = 5.5,
                    fwhm_px: float = 3.0) -> list[tuple]:
    """
    Detecta fontes pontuais com background local + refinamento PSF.

    Retorna lista de (y, x, flux, tamanho), preservando o contrato público.
    As coordenadas vêm do ajuste PSF em float64, sem arredondamento científico.
    """
    data = np.asarray(img, dtype=np.float64)
    # Pan-STARRS marca saturação com 65535. Sem essa máscara, clusters saturados
    # aparecem como picos acima do threshold e entram no casamento entre frames
    # gerando trilhas com cinemática absurda que elimina também trilhas reais.
    mascara = _mascara_pixels_invalidos(data)
    background, rms = _estimar_background_local(data)
    img_sub = (data - background).astype(np.float64)
    detection_img = (img_sub / rms).astype(np.float64)
    # Zerar pixels mascarados evita que (65535 - bg) / rms produza picos
    # artificiais nos buracos antes de o DAOStarFinder receber a imagem.
    detection_img[mascara] = 0.0

    # DAOStarFinder fica como gerador de sementes em uma imagem já corrigida
    # por background local. A medição final não é o centroide do DAO, e sim o
    # centro do modelo PSF ajustado no recorte da fonte.
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

    fontes = []
    for linha in tabela:
        y0 = np.float64(linha[y_col])
        x0 = np.float64(linha[x_col])
        psf = _ajustar_psf_subpixel(img_sub, x0=x0, y0=y0, fwhm_px=fwhm_px,
                                    mascara_invalida=mascara)
        if psf is None:
            continue

        # FWHM físico: raios cósmicos/hot pixels tendem a FWHM sub-pixel,
        # enquanto blends, galáxias e trilhas têm FWHM grande demais para a PSF
        # estelar do conjunto. Esse filtro protege a astrometria antes do WCS.
        fwhm_fit = float(psf["fwhm"])
        if not (T.FWHM_MIN_PX <= fwhm_fit <= T.FWHM_MAX_PX):
            continue

        y = float(np.float64(psf["y"]))
        x = float(np.float64(psf["x"]))
        flux = float(np.float64(psf["flux"]))
        tam = float(np.pi * (fwhm_fit / 2.0) ** 2)
        fontes.append((y, x, flux, tam))
    return fontes


def _detectar_fontes_worker(args: tuple[int, np.ndarray, float, float]) -> tuple[int, list[tuple]]:
    """Worker top-level para multiprocessing com `spawn` no macOS/Apple Silicon."""
    idx, data, sigma, fwhm_px = args
    return idx, detectar_fontes(data, sigma=sigma, fwhm_px=fwhm_px)


def _detectar_fontes_frames_paralelo(frames: list[dict], sigma: float) -> list[list[tuple]]:
    """
    Processa os 4 frames em paralelo.

    Em Apple Silicon M2, quatro imagens FITS independentes são uma carga
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
        return [detectar_fontes(f["data"], sigma=sigma) for f in frames]

    try:
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=n_proc) as pool:
            resultados = pool.map(_detectar_fontes_worker, tarefas)
        resultados.sort(key=lambda item: item[0])
        return [fontes for _, fontes in resultados]
    except Exception as e:
        log(f"Multiprocessing indisponível; usando processamento serial: {e}", "WARN")
        return [detectar_fontes(f["data"], sigma=sigma) for f in frames]


# ─────────────────────────────────────────────
# 3. IDENTIFICAÇÃO DE OBJETOS EM MOVIMENTO
# ─────────────────────────────────────────────
def _casar_mutuo(fontes_a: list[tuple], fontes_b: list[tuple],
                 raio_px: float = None) -> dict[int, int]:
    """
    Casamento mútuo de vizinho mais próximo entre duas listas de fontes.
    Retorna {idx_a: idx_b} onde o casamento é recíproco e dentro de raio_px.
    """
    if raio_px is None:
        raio_px = T.RAIO_MATCH_PX
    if not fontes_a or not fontes_b:
        return {}

    pa = np.array([(f[0], f[1]) for f in fontes_a])
    pb = np.array([(f[0], f[1]) for f in fontes_b])

    d2 = ((pa[:, None, :] - pb[None, :, :]) ** 2).sum(axis=2)

    melhor_a_para_b = np.argmin(d2, axis=1)
    melhor_b_para_a = np.argmin(d2, axis=0)

    pares = {}
    for i, j in enumerate(melhor_a_para_b):
        if melhor_b_para_a[j] == i and d2[i, j] < raio_px ** 2:
            pares[i] = int(j)
    return pares


def _escala_pixel_arcsec(wcs: WCS | None) -> float:
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


def _fontes_para_plano_celeste(fontes: list[tuple], frame: dict,
                               ra_ref: float, dec_ref: float) -> np.ndarray:
    """
    Projeta fontes detectadas para um plano tangente local em arcsec.
    Usa o WCS refinado do frame; se o WCS falhar, retorna coordenadas em pixel.
    """
    if not fontes:
        return np.empty((0, 2), dtype=np.float64)

    pix = np.array([(f[1], f[0]) for f in fontes], dtype=np.float64)
    if not frame.get("wcs_ok") or frame.get("wcs") is None:
        return np.array([(f[0], f[1]) for f in fontes], dtype=np.float64)

    try:
        radec = frame["wcs"].all_pix2world(pix, 0).astype(np.float64)
        cosdec = np.cos(np.radians(float(dec_ref)))
        return np.column_stack([
            _delta_ra_graus(radec[:, 0], ra_ref) * cosdec * 3600.0,
            (radec[:, 1] - dec_ref) * 3600.0,
        ])
    except Exception:
        return np.array([(f[0], f[1]) for f in fontes], dtype=np.float64)


def _media_ra_circular_graus(ras: list[float]) -> float:
    """Média circular de RA em graus, correta em campos cruzando 0/360."""
    ang = np.radians(np.asarray(ras, dtype=np.float64))
    ra = np.degrees(np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang))))
    return float(ra % 360.0)


def _casar_mutuo_celeste(fontes_a: list[tuple], fontes_b: list[tuple],
                         frame_a: dict, frame_b: dict,
                         raio_px: float = None) -> dict[int, int]:
    """
    Casamento mútuo em coordenadas celestes projetadas pelo WCS refinado.
    O raio operacional continua definido em pixels, convertido para arcsec.
    """
    if raio_px is None:
        raio_px = T.RAIO_MATCH_PX
    if not fontes_a or not fontes_b:
        return {}

    try:
        ra_ref = _media_ra_circular_graus([
            float(frame_a["wcs"].wcs.crval[0]),
            float(frame_b["wcs"].wcs.crval[0]),
        ])
        dec_ref = float(np.mean([
            float(frame_a["wcs"].wcs.crval[1]),
            float(frame_b["wcs"].wcs.crval[1]),
        ]))
        escala = float(np.mean([
            _escala_pixel_arcsec(frame_a.get("wcs")),
            _escala_pixel_arcsec(frame_b.get("wcs")),
        ]))
        pa = _fontes_para_plano_celeste(fontes_a, frame_a, ra_ref, dec_ref)
        pb = _fontes_para_plano_celeste(fontes_b, frame_b, ra_ref, dec_ref)
        raio = float(raio_px) * escala
    except Exception:
        return _casar_mutuo(fontes_a, fontes_b, raio_px=raio_px)

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



def encontrar_movedores(frames: list[dict], sigma: float = 5.5) -> tuple[list[dict], dict]:
    """
    Compara fontes entre os 4 frames e identifica objetos com movimento
    residual após remover a deriva do campo.

    Com WCS refinado por Gaia, o casamento mútuo opera sobre coordenadas
    celestes projetadas em plano tangente local, não apenas em pixels brutos.

    Retorna (movedores, metricas_globais).
    """
    todas_fontes = _detectar_fontes_frames_paralelo(frames, sigma=sigma)
    n_fontes_por_frame = []
    for i, f in enumerate(frames):
        fontes = todas_fontes[i]
        n_fontes_por_frame.append(len(fontes))
        log(f"  Frame {i+1} ({f['date_obs'][:19]}): {len(fontes)} fontes")

    pares_12 = _casar_mutuo_celeste(todas_fontes[0], todas_fontes[1],
                                    frames[0], frames[1])
    pares_23 = _casar_mutuo_celeste(todas_fontes[1], todas_fontes[2],
                                    frames[1], frames[2])
    pares_34 = _casar_mutuo_celeste(todas_fontes[2], todas_fontes[3],
                                    frames[2], frames[3])

    candidatos_brutos = []
    n_trilhas_tentadas = 0

    for i0, i1 in pares_12.items():
        if i1 not in pares_23:
            continue
        i2 = pares_23[i1]
        if i2 not in pares_34:
            continue
        i3 = pares_34[i2]
        n_trilhas_tentadas += 1

        f0 = todas_fontes[0][i0]
        f1 = todas_fontes[1][i1]
        f2 = todas_fontes[2][i2]
        f3 = todas_fontes[3][i3]
        trilha = [(f0[0], f0[1], f0[2]),
                  (f1[0], f1[1], f1[2]),
                  (f2[0], f2[1], f2[2]),
                  (f3[0], f3[1], f3[2])]

        pos  = np.array([(t[0], t[1]) for t in trilha], dtype=np.float64)
        move = float(np.hypot(pos[-1, 0] - pos[0, 0], pos[-1, 1] - pos[0, 1]))
        candidatos_brutos.append({
            "trilha"    : trilha,
            "move_total": move,
            "brilhos"   : [t[2] for t in trilha],
        })

    deriva_dx, deriva_dy = 0.0, 0.0
    n_estaveis = 0
    if candidatos_brutos:
        vetores = np.array([
            (c["trilha"][-1][1] - c["trilha"][0][1],
             c["trilha"][-1][0] - c["trilha"][0][0])
            for c in candidatos_brutos
        ])
        deriva_dx = float(np.median(vetores[:, 0]))
        deriva_dy = float(np.median(vetores[:, 1]))
        res_inicial = np.hypot(vetores[:, 0] - deriva_dx,
                               vetores[:, 1] - deriva_dy)
        mascara_estaveis = res_inicial < T.DERIVA_ESTAVEL_PX
        if mascara_estaveis.sum() >= T.DERIVA_MIN_ESTAVEIS:
            deriva_dx = float(np.median(vetores[mascara_estaveis, 0]))
            deriva_dy = float(np.median(vetores[mascara_estaveis, 1]))
        n_estaveis = int(mascara_estaveis.sum())
        log(f"Deriva instrumental estimada (dx, dy): "
            f"({deriva_dx:+.2f}, {deriva_dy:+.2f}) px  "
            f"[{n_estaveis} estrelas usadas]")

        movedores = []
        for c, v in zip(candidatos_brutos, vetores):
            res = float(np.hypot(v[0] - deriva_dx, v[1] - deriva_dy))
            if T.MOVE_MIN_PX < c["move_total"] < T.MOVE_MAX_PX and res > T.RESIDUO_MIN_PX:
                c["residuo"] = res
                movedores.append(c)
    else:
        movedores = []

    metricas = {
        "n_fontes_por_frame"        : n_fontes_por_frame,
        "n_trilhas_tentadas"        : n_trilhas_tentadas,
        "deriva_dx_px"              : round(deriva_dx, 3),
        "deriva_dy_px"              : round(deriva_dy, 3),
        "n_estaveis_deriva"         : n_estaveis,
        "sigma_deteccao"            : sigma,
        "background_model"          : "photutils.Background2D",
        "background_box_size"       : list(T.BACKGROUND_BOX),
        "background_filter_size"    : list(T.BACKGROUND_FILTER),
        "psf_model"                 : T.PSF_MODELO,
        "psf_fwhm_range_px"         : [T.FWHM_MIN_PX, T.FWHM_MAX_PX],
        "multiprocessing_processes" : min(4, len(frames), os.cpu_count() or 1),
        "saturacao_pct_por_frame"   : [round(f.get("sat_pct", 0.0), 2) for f in frames],
        "gaia_refinamento"          : [
            {
                "frame"               : f["arquivo"],
                "status"              : f.get("wcs_status", "wcs_invalido"),
                "rms_pre_arcsec"      : f.get("wcs_rms_pre_arcsec"),
                "rms_pos_arcsec"      : f.get("wcs_rms_pos_arcsec"),
                "n_matches"           : f.get("gaia_n_matches", 0),
                "n_matches_bruto"     : f.get("gaia_n_matches_bruto", 0),
                "offset_bruto_arcsec" : f.get("gaia_offset_bruto_arcsec", [0.0, 0.0]),
                "n_matches_fino"      : f.get("gaia_n_matches_fino", 0),
                "n_matches_ajuste"    : f.get("gaia_n_matches_ajuste", 0),
                "n_matches_usados"    : f.get("gaia_n_matches_usados", 0),
                "match_origem"        : f.get("gaia_match_origem"),
            }
            for f in frames
        ],
    }

    return movedores, metricas


# ─────────────────────────────────────────────
# 4. ANÁLISE E SCORE DE CADA CANDIDATO
# ─────────────────────────────────────────────
def _esta_na_borda(tx: float, ty: float, shape: tuple,
                   margem_px: float = None) -> bool:
    """Retorna True se o centroide está a menos de margem_px da borda."""
    if margem_px is None:
        margem_px = T.MARGEM_BORDA_PX
    ny, nx = shape
    return tx < margem_px or tx > (nx - margem_px) \
        or ty < margem_px or ty > (ny - margem_px)


def _morfologia_por_frame(img: np.ndarray, ty: float, tx: float,
                          r: int = 20) -> dict:
    """
    Calcula métricas morfológicas simples a partir do recorte local.

    Retorna:
        pontual     : razão pico/total (0–1; alta = ponto, baixa = estendida)
        elongation  : razão eixo maior/menor via PCA dos pixels acima do limiar
                      (1.0 = circular; >1.6 = estendida/artefato de trail)
        fwhm_est_px : estimativa de FWHM via limiar de meia potência [px]
        valido      : False se recorte degenerado (borda, fluxo nulo)
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
        return {"pontual": 0.0, "elongation": 1.0, "fwhm_est_px": None, "valido": False}

    bg = float(np.percentile(cut, 30))
    cut_sub = np.clip(cut - bg, 0.0, None)
    total = float(cut_sub.sum())

    if total <= 0.0:
        return {"pontual": 0.0, "elongation": 1.0, "fwhm_est_px": None, "valido": False}

    pontual = float(cut_sub.max() / total)

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

    psf_fit = _ajustar_psf_subpixel(cut_sub, x0=cx_cut, y0=cy_cut,
                                    fwhm_px=T.FWHM_PX)
    fwhm_psf = float(psf_fit["fwhm"]) if psf_fit is not None else None
    modelo_psf = psf_fit["modelo"] if psf_fit is not None else None

    return {
        "pontual"    : pontual,
        "elongation" : elongation,
        "fwhm_est_px": fwhm_psf if fwhm_psf is not None else fwhm_est,
        "modelo_psf" : modelo_psf,
        "valido"     : True,
    }


def _verificar_falso_positivo(c: dict, frames: list[dict],
                               snrs: list[float],
                               metricas_morfo: list[dict],
                               brilhos: list[float]) -> list[str]:
    """
    Avalia heurísticas de falso positivo e retorna lista de razões de rejeição.
    Lista vazia significa candidato sem evidência de ser artefato.

    Heurísticas:
    - Borda do frame (qualquer posição na trilha)
    - SNR insuficiente em mais de 1 frame
    - Perfil de hot pixel (pontualidade extrema + SNR fraco)
    - Brilho completamente caótico (CV >> limiar)
    - Morfologia estendida grave e inconsistente (elongation muito alta)
    """
    razoes = []
    trilha   = c["trilha"]
    forma_img = frames[0]["data"].shape

    # Borda: rejeição (não apenas flag)
    frames_na_borda = [
        fi for fi in range(4)
        if _esta_na_borda(trilha[fi][1], trilha[fi][0], forma_img)
    ]
    if frames_na_borda:
        razoes.append(
            f"posição na borda em frame(s) {[f+1 for f in frames_na_borda]} "
            f"(margem {T.MARGEM_BORDA_PX:.0f} px)"
        )

    # SNR insuficiente
    snrs_validos = [s for s in snrs if s is not None]
    n_snr_baixo = sum(1 for s in snrs_validos if s < T.SNR_MINIMO)
    if n_snr_baixo >= 2:
        razoes.append(
            f"SNR insuficiente em {n_snr_baixo}/4 frames "
            f"(mínimo {T.SNR_MINIMO:.1f}; valores={[round(s,1) for s in snrs_validos]})"
        )

    # Hot pixel: pontualidade muito alta + SNR borderline
    pontuais_frame = [m["pontual"] for m in metricas_morfo if m["valido"]]
    pont_media = float(np.mean(pontuais_frame)) if pontuais_frame else 0.0
    snr_medio  = float(np.mean(snrs_validos)) if snrs_validos else 0.0
    if pont_media > T.HOT_PIXEL_PONT_MAX and snr_medio < T.HOT_PIXEL_SNR_MIN:
        razoes.append(
            f"perfil suspeito de hot pixel: pontualidade={pont_media:.3f} "
            f"(>{T.HOT_PIXEL_PONT_MAX}) com SNR médio={snr_medio:.1f} "
            f"(<{T.HOT_PIXEL_SNR_MIN})"
        )

    # Brilho caótico
    brilhos_arr = np.array(brilhos, dtype=np.float64)
    if brilhos_arr.mean() > 0:
        cv = float(brilhos_arr.std() / brilhos_arr.mean())
        if cv > T.BRILHO_CV_CAOS:
            razoes.append(
                f"brilho caótico: CV={cv:.2f} (>{T.BRILHO_CV_CAOS}) — "
                "inconsistente com objeto real"
            )

    # Elongation grave (artefato de trail CCD, raio cósmico, etc.)
    elongs = [m["elongation"] for m in metricas_morfo if m["valido"]]
    if elongs:
        elong_max = float(max(elongs))
        elong_med = float(np.mean(elongs))
        # Só rejeita se TODOS os frames mostram elongation grave
        if elong_med > 3.0 and elong_max > 4.0:
            razoes.append(
                f"morfologia sistematicamente estendida: "
                f"elongation média={elong_med:.2f}, máx={elong_max:.2f} "
                "(possível artefato, trail CCD ou galáxia de fundo)"
            )

    # FWHM incompatível com a PSF estelar: a posição MPC precisa vir do centro
    # de uma fonte pontual. FWHM sub-pixel é assinatura clássica de raio cósmico
    # ou hot pixel; FWHM muito largo desloca o centro por blend/galáxia/trail.
    fwhms = [float(m["fwhm_est_px"]) for m in metricas_morfo
             if m["valido"] and m.get("fwhm_est_px") is not None
             and np.isfinite(m["fwhm_est_px"])]
    if fwhms:
        fwhm_med = float(np.median(fwhms))
        if fwhm_med < T.FWHM_MIN_PX or fwhm_med > T.FWHM_MAX_PX:
            razoes.append(
                f"FWHM incompatível com PSF estelar: mediana={fwhm_med:.2f} px "
                f"(faixa aceita {T.FWHM_MIN_PX:.1f}-{T.FWHM_MAX_PX:.1f} px)"
            )

    return razoes


def _calcular_snr_local(img: np.ndarray, ty: float, tx: float,
                        raio_sinal: int = 4, raio_bg: int = 12) -> float:
    """
    SNR local simples: sinal = soma de pixels no raio_sinal após subtrair
    background estimado em annulus entre raio_sinal e raio_bg.
    """
    ny, nx = img.shape
    yc = int(np.rint(float(ty)))
    xc = int(np.rint(float(tx)))
    y0s = max(0, yc - raio_sinal); y1s = min(ny, yc + raio_sinal + 1)
    x0s = max(0, xc - raio_sinal); x1s = min(nx, xc + raio_sinal + 1)
    corte_sinal = img[y0s:y1s, x0s:x1s]

    y0b = max(0, yc - raio_bg); y1b = min(ny, yc + raio_bg + 1)
    x0b = max(0, xc - raio_bg); x1b = min(nx, xc + raio_bg + 1)
    corte_bg = img[y0b:y1b, x0b:x1b]

    bg_med  = float(np.median(corte_bg))
    bg_std  = float(np.std(corte_bg))
    sinal   = float(np.sum(corte_sinal - bg_med))
    n_pix   = corte_sinal.size
    ruido   = bg_std * np.sqrt(n_pix) if n_pix > 0 else 1.0
    return float(sinal / ruido) if ruido > 0 else 0.0


def _calcular_incertezas_astrometricas(
    frames: list[dict],
    metricas_morfo: list[dict],
    snrs: list[float],
) -> dict:
    """Estimativa por frame: erro centroidal FWHM/SNR combinado ao RMS Gaia."""
    por_frame = []
    sigma_total = []
    for fi, frame in enumerate(frames):
        snr = float(snrs[fi]) if fi < len(snrs) and snrs[fi] is not None else 0.0
        morfo = metricas_morfo[fi] if fi < len(metricas_morfo) else {}
        fwhm_px = morfo.get("fwhm_est_px") or T.FWHM_PX
        escala = _escala_pixel_arcsec(frame.get("wcs"))
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

        por_frame.append({
            "frame_index": fi,
            "sigma_centroid_arcsec": round(sigma_centroid, 4)
                                      if sigma_centroid is not None else None,
            "sigma_wcs_arcsec": round(sigma_wcs, 4) if sigma_wcs is not None else None,
            "sigma_pos_arcsec": round(sigma_pos, 4) if sigma_pos is not None else None,
            "metodo": "FWHM/(2.355*SNR) combinado em quadratura com RMS Gaia",
        })

    return {
        "por_frame": por_frame,
        "sigma_pos_mediana_arcsec": round(float(np.median(sigma_total)), 4)
                                    if sigma_total else None,
    }


def _avaliar_fonte_estatica_gaia(c: dict, frames: list[dict]) -> dict:
    """Marca trilhas compatíveis com a mesma fonte Gaia em vários frames."""
    matches = []
    source_ids = []
    for fi, frame in enumerate(frames):
        tabela = frame.get("gaia_catalogo_refinado")
        if tabela is None or not frame.get("wcs_ok") or len(tabela) == 0:
            continue

        try:
            ra, dec = pixel_para_radec(c["trilha"][fi][1], c["trilha"][fi][0], frame["wcs"])
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
        "status": "fonte_estatica_gaia" if mesma_fonte else "sem_match_estatico",
        "n_matches": len(matches),
        "matches": matches,
    }


def analisar_candidato(c: dict, frames: list[dict]) -> dict:
    """
    Calcula métricas, score decomposto, morfologia robusta, flags diagnósticas
    e razões de decisão para um candidato.

    Score (0–12):
        linearidade        : 0–3  (resíduo da regressão linear na trilha)
        velocidade         : 0–2  (consistência entre passos consecutivos)
        fotometria         : 0–2  (coeficiente de variação do brilho)
        morfologia         : 0–2  (pontualidade média do perfil)
        consistencia_morfo : 0–1  (estabilidade morfológica entre frames)
        elongation         : 0–1  (compacidade; penaliza fontes estendidas)
        faixa_vel          : 0–1  (velocidade total na faixa típica de MBA)

    Candidatos com razões de rejeição recebem classe=DESCARTA independente do score.
    """
    trilha    = c["trilha"]
    pos       = np.array([(t[1], t[0]) for t in trilha], dtype=np.float64)  # (x, y)
    t_idx     = np.arange(4, dtype=np.float64)
    forma_img = frames[0]["data"].shape

    razoes_penalizacao: list[str] = []   # penaliza score mas não rejeita
    razoes_rejeicao:    list[str] = []   # rejeição direta → DESCARTA

    # ── Linearidade ──────────────────────────────────────────────────────────
    px = np.polyfit(t_idx, pos[:, 0], 1)
    py = np.polyfit(t_idx, pos[:, 1], 1)
    res_x = pos[:, 0] - np.polyval(px, t_idx)
    res_y = pos[:, 1] - np.polyval(py, t_idx)
    lin = float(np.mean(np.sqrt(res_x ** 2 + res_y ** 2)))

    # ── Velocidade / consistência de passos ──────────────────────────────────
    passos_x = np.diff(pos[:, 0])
    passos_y = np.diff(pos[:, 1])
    vel_cons = float(np.std(passos_x) + np.std(passos_y))

    # Velocidade angular em arcsec/min
    vel_arcsec_min = None
    if frames[0]["wcs_ok"] and frames[-1]["wcs_ok"]:
        try:
            ra0, dec0 = frames[0]["wcs"].all_pix2world(
                trilha[0][1], trilha[0][0], 0)
            ra3, dec3 = frames[-1]["wcs"].all_pix2world(
                trilha[-1][1], trilha[-1][0], 0)
            c0 = SkyCoord(ra=float(ra0) * u.deg, dec=float(dec0) * u.deg)
            c3 = SkyCoord(ra=float(ra3) * u.deg, dec=float(dec3) * u.deg)
            sep_arcsec = float(c0.separation(c3).to(u.arcsec).value)
            dt_min = (frames[-1]["jd"] - frames[0]["jd"]) * 1440.0
            if dt_min > 0:
                vel_arcsec_min = round(sep_arcsec / dt_min, 4)
        except Exception:
            pass

    # ── Fotometria e morfologia robusta por frame ─────────────────────────────
    brilhos       : list[float] = []
    snrs          : list[float] = []
    metricas_morfo: list[dict]  = []

    for fi, frame in enumerate(frames):
        ty, tx = trilha[fi][0], trilha[fi][1]
        morfo = _morfologia_por_frame(frame["data"], ty, tx, r=20)
        metricas_morfo.append(morfo)
        brilhos.append(morfo["pontual"] * 1.0)   # placeholder; recalcular abaixo

        # Brilho: soma do recorte subtraído de background
        r = 20
        yc = int(np.rint(float(ty)))
        xc = int(np.rint(float(tx)))
        y0_ = max(0, yc - r); y1_ = min(frame["data"].shape[0], yc + r + 1)
        x0_ = max(0, xc - r); x1_ = min(frame["data"].shape[1], xc + r + 1)
        cut = frame["data"][y0_:y1_, x0_:x1_].astype(np.float64)
        bg  = np.percentile(cut, 30)
        brilhos[-1] = float(np.sum(np.clip(cut - bg, 0, None)))

        snrs.append(_calcular_snr_local(frame["data"], ty, tx))

    brilhos_arr = np.array(brilhos, dtype=np.float64)
    brilho_cv   = float(np.std(brilhos_arr) / (np.mean(brilhos_arr) + 1e-10)) \
                  if np.all(brilhos_arr > 0) else float("inf")

    pontuais_frame = [m["pontual"] for m in metricas_morfo if m["valido"]]
    pont           = float(np.mean(pontuais_frame)) if pontuais_frame else 0.0
    pont_std       = float(np.std(pontuais_frame))  if len(pontuais_frame) > 1 else 0.0

    elongs_frame   = [m["elongation"] for m in metricas_morfo if m["valido"]]
    elong_medio    = float(np.mean(elongs_frame)) if elongs_frame else 1.0

    fwhms_frame    = [m["fwhm_est_px"] for m in metricas_morfo
                      if m["valido"] and m["fwhm_est_px"] is not None]
    fwhm_medio     = float(np.mean(fwhms_frame)) if fwhms_frame else None

    # ── Verificação de falso positivo ────────────────────────────────────────
    razoes_rejeicao = _verificar_falso_positivo(
        c, frames, snrs, metricas_morfo, brilhos
    )
    incerteza_astrometrica = _calcular_incertezas_astrometricas(
        frames, metricas_morfo, snrs
    )
    gaia_static = _avaliar_fonte_estatica_gaia(c, frames)
    if gaia_static["status"] == "fonte_estatica_gaia":
        razoes_rejeicao.append(
            f"compatível com fonte estática Gaia em {gaia_static['n_matches']}/4 frames"
        )

    # ── Score por componente ─────────────────────────────────────────────────
    move = c["move_total"]

    if   lin < T.LIN_EXCELENTE:  s_lin = 3
    elif lin < T.LIN_BOA:        s_lin = 2
    elif lin < T.LIN_MARGINAL:   s_lin = 1
    else:
        s_lin = 0
        razoes_penalizacao.append(f"linearidade ruim (resíduo={lin:.2f} px)")

    if   vel_cons < T.VEL_UNIFORME:  s_vel = 2
    elif vel_cons < T.VEL_MODERADA:  s_vel = 1
    else:
        s_vel = 0
        razoes_penalizacao.append(f"velocidade irregular (σ={vel_cons:.2f} px)")

    if   brilho_cv < T.FOT_ESTAVEL:  s_fot = 2
    elif brilho_cv < T.FOT_MODERADA: s_fot = 1
    else:
        s_fot = 0
        razoes_penalizacao.append(f"brilho instável (CV={brilho_cv:.2f})")

    if   pont > T.MOR_PONTUAL:   s_mor = 2
    elif pont > T.MOR_MARGINAL:  s_mor = 1
    else:
        s_mor = 0
        razoes_penalizacao.append(f"morfologia estendida (pontualidade={pont:.4f})")

    # Consistência morfológica entre frames (novo)
    s_morfo_consist = 1 if pont_std < T.MOR_CONSIST_MAX_STD else 0
    if s_morfo_consist == 0:
        razoes_penalizacao.append(
            f"morfologia inconsistente entre frames (std pontualidade={pont_std:.3f})"
        )

    # Elongation (novo)
    s_elon = 1 if elong_medio < T.ELON_BOA else 0
    if s_elon == 0:
        razoes_penalizacao.append(
            f"elongation elevada (média={elong_medio:.2f}, limiar={T.ELON_BOA})"
        )

    s_faixa = 1 if T.FAIXA_VEL_MIN_PX < move < T.FAIXA_VEL_MAX_PX else 0

    score_bruto = s_lin + s_vel + s_fot + s_mor + s_morfo_consist + s_elon + s_faixa
    score_max   = 3 + 2 + 2 + 2 + 1 + 1 + 1   # = 12

    # Candidatos com razão de rejeição → DESCARTA independente de score
    if razoes_rejeicao:
        score   = 0
        classe  = "DESCARTA"
        cor_cli = Cor.VERMELHO
    else:
        score = score_bruto
        if   score >= T.SCORE_FORTE:    classe, cor_cli = "FORTE",    Cor.VERDE
        elif score >= T.SCORE_MODERADO: classe, cor_cli = "MODERADO", Cor.AMARELO
        elif score >= T.SCORE_FRACO:    classe, cor_cli = "FRACO",    Cor.AMARELO
        else:                           classe, cor_cli = "DESCARTA", Cor.VERMELHO

    # Probabilidade heurística (escala 5–95 %)
    prob = min(95, max(5, int(score / score_max * 100)))

    # ── Flags diagnósticas ───────────────────────────────────────────────────
    flags: list[str] = []

    if s_lin >= 2:   flags.append("LINEARIDADE_BOA")
    elif s_lin == 0: flags.append("LINEARIDADE_RUIM")

    if s_vel >= 2:   flags.append("VELOCIDADE_CONSISTENTE")
    elif s_vel == 0: flags.append("VELOCIDADE_IRREGULAR")

    if s_fot >= 2:   flags.append("BRILHO_ESTAVEL")
    elif s_fot == 0: flags.append("BRILHO_INSTAVEL")

    if s_mor >= 2:   flags.append("MORFOLOGIA_PONTUAL")
    elif s_mor == 0: flags.append("MORFOLOGIA_ESTENDIDA")

    if s_elon == 0:         flags.append("ELONGACAO_ALTA")
    if s_morfo_consist == 0: flags.append("MORFO_INCONSISTENTE")

    borda = any(
        _esta_na_borda(trilha[fi][1], trilha[fi][0], forma_img)
        for fi in range(4)
    )
    if borda:
        flags.append("EDGE_FRAME")

    if razoes_rejeicao:
        flags.append("REJEITADO_FP")
    if gaia_static["status"] == "fonte_estatica_gaia":
        flags.append("GAIA_FONTE_ESTATICA")

    # ── Log de decisão ───────────────────────────────────────────────────────
    if razoes_rejeicao:
        for r in razoes_rejeicao:
            log(f"  REJEIÇÃO FP: {r}", "WARN")
    elif razoes_penalizacao:
        for r in razoes_penalizacao:
            log(f"  Penalização: {r}", "WARN")
    else:
        log(f"  Promovido: linearidade={lin:.2f}px, vel_cons={vel_cons:.2f}, "
            f"CV={brilho_cv:.2f}, pont={pont:.4f}", "OK")

    return {
        **c,
        "linearidade"           : lin,
        "vel_consistencia"      : vel_cons,
        "vel_arcsec_min"        : vel_arcsec_min,
        "brilhos"               : brilhos_arr.tolist(),
        "brilho_cv"             : brilho_cv,
        "pontual"               : pont,
        "pont_std"              : pont_std,
        "elong_medio"           : elong_medio,
        "fwhm_medio_px"         : fwhm_medio,
        "metricas_morfo"        : metricas_morfo,
        "snrs"                  : snrs,
        "incerteza_astrometrica": incerteza_astrometrica,
        "gaia_static"           : gaia_static,
        "score"                 : score,
        "score_bruto"           : score_bruto,
        "score_max"             : score_max,
        "score_linearidade"     : s_lin,
        "score_velocidade"      : s_vel,
        "score_fotometria"      : s_fot,
        "score_morfologia"      : s_mor,
        "score_morfo_consist"   : s_morfo_consist,
        "score_elongation"      : s_elon,
        "score_faixa_vel"       : s_faixa,
        "probabilidade"         : prob,
        "classe"                : classe,
        "cor_cli"               : cor_cli,
        "flags"                 : flags,
        "razoes_penalizacao"    : razoes_penalizacao,
        "razoes_rejeicao"       : razoes_rejeicao,
    }


# ─────────────────────────────────────────────
# 5. CONVERSÃO PIXEL → RA/DEC
# ─────────────────────────────────────────────
def pixel_para_radec(x: float, y: float, wcs: WCS) -> tuple[float, float]:
    if wcs is None:
        raise ValueError("WCS ausente — impossível converter pixel em RA/Dec.")
    ra, dec = wcs.all_pix2world(x, y, 0)
    return float(ra), float(dec)


# ─────────────────────────────────────────────
# 6. CONSULTA AO MINOR PLANET CENTER
# ─────────────────────────────────────────────
def consultar_catalogo_posicional(ra: float, dec: float, data_obs: str,
                                  raio_arcmin: float = 2.0) -> dict:
    """
    Consulta o SkyBot (IMCCE) para procurar correspondência posicional próxima.

    Retorna dict com schema normalizado:
        status: nao_consultado | consulta_falhou | sem_match |
                match_ambiguo  | match_provavel
    """
    if not (0 <= ra < 360) or not (-90 <= dec <= 90):
        return {
            "status": "consulta_falhou",
            "motivo": f"coordenadas fora do domínio (RA={ra:.3f}, Dec={dec:.3f})",
            "match" : None,
        }

    try:
        from astroquery.imcce import Skybot

        coord    = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        epoch    = Time(data_obs, format="isot", scale="utc")
        resultado = Skybot.cone_search(coord, rad=raio_arcmin * u.arcmin,
                                       epoch=epoch)

        if resultado is not None and len(resultado) > 0:
            # Mais de um objeto no cone → ambíguo
            status = "match_ambiguo" if len(resultado) > 1 else "match_provavel"
            obj    = resultado[0]

            num  = str(obj["Number"])  if "Number" in resultado.colnames else ""
            nome = str(obj["Name"])    if "Name"   in resultado.colnames \
                   else str(obj.get("Designation", "desconhecido"))
            mag  = float(obj["V"])     if "V"      in resultado.colnames else None
            dist = float(obj["centerdist"].to(u.arcmin).value) \
                   if "centerdist" in resultado.colnames else None
            tipo = str(obj["Type"])    if "Type"   in resultado.colnames else None

            return {
                "status"  : status,
                "motivo"  : None,
                "n_objetos_no_cone": len(resultado),
                "match"   : {
                    "nome_designacao"     : f"{num} {nome}".strip(),
                    "tipo"                : tipo,
                    "magnitude_v"         : round(mag, 2) if mag is not None else None,
                    "distancia_arcmin"    : round(dist, 4) if dist is not None else None,
                },
            }

        return {
            "status"  : "sem_match",
            "motivo"  : None,
            "n_objetos_no_cone": 0,
            "match"   : None,
        }

    except ImportError:
        return {
            "status": "consulta_falhou",
            "motivo": "astroquery não instalado",
            "match" : None,
        }
    except Exception as e:
        msg = str(e)
        # SkyBot retorna "No table found" quando o cone não contém objetos
        if "no table found" in msg.lower() or "no result" in msg.lower():
            return {
                "status"           : "sem_match",
                "motivo"           : None,
                "n_objetos_no_cone": 0,
                "match"            : None,
            }
        return {
            "status": "consulta_falhou",
            "motivo": msg[:200],
            "match" : None,
        }


def _adicionar_flag_mpc(c: dict) -> None:
    """Adiciona flag MPC ao candidato com base no resultado da consulta."""
    mpc = c.get("mpc", {})
    status = mpc.get("status", "nao_consultado")
    flags = c["flags"]

    if status == "sem_match":
        flags.append("MPC_SEM_MATCH")
    elif status == "match_provavel":
        flags.append("MPC_MATCH_PROVAVEL")
    elif status == "match_ambiguo":
        flags.append("MPC_MATCH_AMBIGUO")
    elif status == "consulta_falhou":
        flags.append("MPC_CONSULTA_FALHOU")


# ─────────────────────────────────────────────
# 7. GERAÇÃO DO RELATÓRIO MPC
# ─────────────────────────────────────────────
def formatar_ra_mpc(ra_deg: float) -> str:
    ra_wrap = float(ra_deg) % 360.0
    return Angle(ra_wrap * u.deg).to_string(
        unit=u.hourangle, sep=" ", precision=2, pad=True,
    )


def formatar_dec_mpc(dec_deg: float) -> str:
    dec_clip = max(-90.0, min(90.0, float(dec_deg)))
    return Angle(dec_clip * u.deg).to_string(
        unit=u.deg, sep=" ", precision=1, pad=True, alwayssign=True,
    )


def formatar_data_mpc(jd: float) -> str:
    t  = Time(jd, format="jd", scale="utc")
    dt = t.to_datetime()
    frac = (dt.hour * 3600 + dt.minute * 60
            + dt.second + dt.microsecond / 1e6) / 86400.0
    return f"{dt.year:04d} {dt.month:02d} {dt.day + frac:08.5f}"


def gerar_relatorio_mpc(candidatos: list[dict], frames: list[dict],
                        nome_conjunto: str, pasta_output: Path,
                        observador: dict) -> Path:
    """
    Gera arquivo de relatório no formato MPC 80-colunas.
    Apenas candidatos com score >= 6 (MODERADO+) são incluídos.
    Magnitude deixada em branco — requer calibração fotométrica no Astrometrica.
    """
    linhas  = []
    incluidos = [c for c in candidatos if c["score"] >= T.MPC_SCORE_MIN]

    for i, c in enumerate(incluidos):
        trilha     = c["trilha"]
        designacao = f"TMP{i+1:04d}"

        for fi, frame in enumerate(frames):
            if not frame["wcs_ok"]:
                continue
            tx, ty   = trilha[fi][1], trilha[fi][0]
            ra, dec  = pixel_para_radec(tx, ty, frame["wcs"])
            data_mpc = formatar_data_mpc(frame["jd"])
            ra_mpc   = formatar_ra_mpc(ra)
            dec_mpc  = formatar_dec_mpc(dec)

            linha = (
                "     "
                + f"{designacao:<7s}"
                + " "
                + " "
                + "C"
                + data_mpc
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
            linhas.append(linha[:80].ljust(80))

    caminho   = pasta_output / f"{nome_conjunto}_MPC_report.txt"
    nome_obs  = observador.get("nome",  "Observador IASC")
    email_obs = observador.get("email", "observer@example.com")
    with open(caminho, "w", encoding="ascii", errors="replace") as f:
        f.write("COD F51\n")
        f.write("OBS IASC Citizen Scientist\n")
        f.write(f"MEA {nome_obs}\n")
        f.write("TEL 1.8-m f/4.4 Ritchey-Chretien + CCD\n")
        f.write(f"ACK TRITON pipeline v{VERSAO}\n")
        f.write(f"AC2 {email_obs}\n")
        f.write("----- ------------------- ------------------ "
                "------------------ ----- ---\n")
        for linha in linhas:
            f.write(linha + "\n")
        f.write("-----\n")

    return caminho


# ─────────────────────────────────────────────
# 8. VISUALIZAÇÃO
# ─────────────────────────────────────────────
def gerar_visualizacao(candidatos: list[dict], frames: list[dict],
                       pasta_output: Path, nome_conjunto: str):
    """
    Gera PNG com os top-5 candidatos (não-DESCARTA).
    Cada linha: 4 recortes do candidato nos 4 frames.
    Inclui: rank, score, flags principais, vetor de trajetória.
    """
    top = [c for c in candidatos if c["classe"] != "DESCARTA"][:5]
    if not top:
        log("Nenhum candidato válido para visualização.", "WARN")
        return

    cores_classe = {
        "FORTE"   : "#00e07a",
        "MODERADO": "#ffaa00",
        "FRACO"   : "#ff7733",
        "DESCARTA": "#ff2222",
    }
    times = [f["date_obs"][11:16] + " UT" for f in frames]
    n = len(top)
    fig, axes = plt.subplots(n, 4, figsize=(20, 4.5 * n), squeeze=False)
    fig.patch.set_facecolor("#0d0d1a")

    for ri, cand in enumerate(top):
        cor    = cores_classe[cand["classe"]]
        trilha = cand["trilha"]
        flags  = cand.get("flags", [])

        # Flags mais relevantes para exibir na imagem (máximo 3)
        flags_exibir = [f for f in flags
                        if f not in ("MPC_SEM_MATCH",)][:3]
        flags_str = "  ".join(flags_exibir) if flags_exibir else ""

        for fi, frame in enumerate(frames):
            ax = axes[ri][fi]
            ty, tx = trilha[fi][0], trilha[fi][1]
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
                ppy, ppx = trilha[pfi][0], trilha[pfi][1]
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
                    "sem_match"      : "SkyBot: sem match",
                    "match_provavel" : "SkyBot: match",
                    "match_ambiguo"  : "SkyBot: ambíguo",
                    "consulta_falhou": "SkyBot: falhou",
                    "nao_consultado" : "",
                }.get(mpc_status, "")

                ra_str = ""
                if frame["wcs_ok"]:
                    ra_v, dec_v = pixel_para_radec(tx, ty, frame["wcs"])
                    ra_str = f"RA {ra_v:.4f}° Dec {dec_v:.4f}°"

                ylabel_lines = [
                    f"#{ri+1} {cand['classe']}",
                    f"Score {cand['score']}/{cand.get('score_max', 12)}  "
                    f"({cand['probabilidade']}%)",
                    f"Lin={cand['linearidade']:.2f}px  Vel={cand['vel_consistencia']:.2f}",
                ]
                if ra_str:
                    ylabel_lines.append(ra_str)
                if mpc_label:
                    ylabel_lines.append(mpc_label)

                ax.set_ylabel(
                    "\n".join(ylabel_lines),
                    color=cor, fontsize=6.5, labelpad=4,
                )

            # Título do frame (apenas na linha 0)
            if ri == 0:
                ax.set_title(f"Frame {fi+1}\n{times[fi]}",
                             color="white", fontsize=9)

            # Score e rank no canto superior direito do frame 0 de cada candidato
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
        f"TRITON v{VERSAO} | {nome_conjunto} | "
        f"{frames[0]['date_obs'][:10]}\n"
        "● posição atual   + posição anterior   → vetor trajetória",
        color="white", fontsize=11, y=1.01,
    )
    plt.tight_layout(h_pad=0.8, w_pad=0.4)
    caminho = pasta_output / f"{nome_conjunto}_candidatos.png"
    plt.savefig(caminho, dpi=130, bbox_inches="tight", facecolor="#0d0d1a")
    plt.close()
    log(f"Visualização salva: {caminho}", "OK")


# ─────────────────────────────────────────────
# 9. RELATÓRIO TEXTO FINAL
# ─────────────────────────────────────────────
def _resumo_diagnostico(c: dict) -> str:
    """
    Gera uma frase técnica curta explicando o ranking do candidato.
    Inclui razões de penalização e rejeição quando presentes.
    """
    partes = []

    # Razões de rejeição têm precedência
    if c.get("razoes_rejeicao"):
        partes.append("REJEITADO: " + c["razoes_rejeicao"][0])
        if len(c["razoes_rejeicao"]) > 1:
            partes.append(f"+{len(c['razoes_rejeicao'])-1} razão(ões) adicionais")
        return "; ".join(partes) + "."

    lin = c["linearidade"]
    if lin < T.LIN_EXCELENTE:
        partes.append(f"trilha linear (res={lin:.2f} px)")
    elif lin < T.LIN_MARGINAL:
        partes.append(f"linearidade marginal (res={lin:.2f} px)")
    else:
        partes.append(f"linearidade ruim (res={lin:.2f} px)")

    vel = c["vel_consistencia"]
    if vel < T.VEL_UNIFORME:
        partes.append("velocidade uniforme")
    else:
        partes.append(f"velocidade irregular (σ={vel:.2f})")

    cv = c["brilho_cv"]
    if cv < T.FOT_ESTAVEL:
        partes.append("brilho estável")
    elif cv < T.FOT_MODERADA:
        partes.append(f"brilho moderadamente variável (CV={cv:.2f})")
    else:
        partes.append(f"brilho instável (CV={cv:.2f})")

    elong = c.get("elong_medio", 1.0)
    if elong >= T.ELON_BOA:
        partes.append(f"elongação elevada ({elong:.2f})")

    if c.get("razoes_penalizacao"):
        n = len(c["razoes_penalizacao"])
        partes.append(f"{n} penalização(ões) de score")

    mpc_status = c.get("mpc", {}).get("status", "nao_consultado")
    if mpc_status == "match_provavel":
        nome = c["mpc"].get("match", {}).get("nome_designacao", "?")
        dist = c["mpc"].get("match", {}).get("distancia_arcmin", "?")
        partes.append(f"match SkyBot: {nome} a {dist}'")
    elif mpc_status == "sem_match":
        partes.append("sem correspondência posicional no SkyBot")
    elif mpc_status == "match_ambiguo":
        partes.append("múltiplos objetos no cone SkyBot — verificar")

    if "EDGE_FRAME" in c.get("flags", []):
        partes.append("posição na borda do frame")

    return "; ".join(partes) + "."


def gerar_relatorio_txt(candidatos: list[dict], frames: list[dict],
                        nome_conjunto: str, pasta_output: Path,
                        metricas_globais: dict) -> Path:
    """Gera relatório operacional em texto com priorização e diagnóstico."""
    caminho = pasta_output / f"{nome_conjunto}_relatorio.txt"
    linhas  = []
    sep     = "=" * 72
    sub     = "-" * 72
    sub2    = "·" * 72

    # ── Cabeçalho ──────────────────────────────────────────────────────────
    linhas += [
        sep,
        "  TRITON — Relatório de Análise",
        f"  Versão do pipeline   : {VERSAO}",
        f"  Data de processamento: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"  Conjunto de imagens  : {nome_conjunto}",
        f"  Observatório         : Pan-STARRS / IASC (código F51)",
        f"  Data de observação   : {frames[0]['date_obs'][:10]}",
        f"  Intervalo temporal   : {frames[0]['date_obs'][11:19]} → "
        f"{frames[-1]['date_obs'][11:19]} UTC",
        sep, "",
    ]

    # ── Métricas globais do conjunto ───────────────────────────────────────
    mg = metricas_globais
    linhas += [
        "  MÉTRICAS DO CONJUNTO PROCESSADO",
        sub,
        f"  Frames válidos       : {mg.get('n_frames_validos', 4)} de 4",
        f"  Sigma de detecção    : {mg.get('sigma_deteccao', '?')}",
        f"  Fontes detectadas    : " +
        "  |  ".join(
            f"F{i+1}={n}" for i, n in enumerate(mg.get("n_fontes_por_frame", []))
        ),
        f"  Trilhas tentadas     : {mg.get('n_trilhas_tentadas', '?')}",
        f"  Candidatos únicos    : {mg.get('n_candidatos_unicos', '?')}",
        f"  Deriva do campo      : dx={mg.get('deriva_dx_px', 0):+.2f} px  "
        f"dy={mg.get('deriva_dy_px', 0):+.2f} px  "
        f"[{mg.get('n_estaveis_deriva', '?')} estrelas]",
        f"  Tempo de execução    : {mg.get('tempo_execucao_s', '?')} s",
        "", sub,
        f"  Candidatos totais    : {len(candidatos)}",
        f"  FORTE    (8-12) : {sum(1 for c in candidatos if c['classe']=='FORTE')}",
        f"  MODERADO (6-7)  : {sum(1 for c in candidatos if c['classe']=='MODERADO')}",
        f"  FRACO    (4-5)  : {sum(1 for c in candidatos if c['classe']=='FRACO')}",
        f"  DESCARTA (0-3)  : {sum(1 for c in candidatos if c['classe']=='DESCARTA')}",
    ]

    # Score resumido
    scores_validos = [c["score"] for c in candidatos if c["classe"] != "DESCARTA"]
    if scores_validos:
        linhas += [
            f"  Score (não-DESCARTA) : "
            f"min={min(scores_validos)}  max={max(scores_validos)}  "
            f"med={round(float(np.median(scores_validos)), 1)}",
        ]
    linhas += ["", sep, ""]

    # ── Seção de priorização operacional ───────────────────────────────────
    linhas += [
        "  PRIORIZAÇÃO OPERACIONAL",
        "  (use como guia de inspeção no Astrometrica — não como decisão final)",
        sub,
    ]

    grupos = {
        "INSPECIONAR PRIMEIRO"      : [],
        "INSPECIONAR SE HOUVER TEMPO": [],
        "PROVÁVEL ARTEFATO"         : [],
        "IGNORAR"                   : [],
    }
    for c in candidatos:
        cl = c["classe"]
        mpc_status = c.get("mpc", {}).get("status", "nao_consultado")
        borda = "EDGE_FRAME" in c.get("flags", [])
        if cl == "FORTE" and mpc_status in ("sem_match", "nao_consultado"):
            grupos["INSPECIONAR PRIMEIRO"].append(c)
        elif cl in ("FORTE", "MODERADO"):
            grupos["INSPECIONAR SE HOUVER TEMPO"].append(c)
        elif cl == "FRACO":
            grupos["PROVÁVEL ARTEFATO"].append(c)
        else:
            grupos["IGNORAR"].append(c)

    rank_global = {c["rank"]: c for c in candidatos if "rank" in c}

    for grupo_nome, grupo_cands in grupos.items():
        if not grupo_cands:
            continue
        linhas += [f"  ┌─ {grupo_nome} ({len(grupo_cands)}) "]
        for c in grupo_cands:
            rank_label = f"#{c.get('rank', '?')}"
            flags_str  = "  ".join(c.get("flags", []))
            resumo     = _resumo_diagnostico(c)
            linhas += [
                f"  │  {rank_label}  score {c['score']}/{c.get('score_max', 12)}  {c['classe']}",
                f"  │     Lin={c['score_linearidade']}/3  "
                f"Vel={c['score_velocidade']}/2  "
                f"Fot={c['score_fotometria']}/2  "
                f"Mor={c['score_morfologia']}/2  "
                f"Faixa={c['score_faixa_vel']}/1",
                f"  │     Flags: {flags_str if flags_str else '—'}",
                f"  │     {resumo}",
                f"  │",
            ]
        linhas[-1] = linhas[-1].replace("  │", "  └─")
        linhas += [""]

    linhas += [sep, ""]

    # ── Detalhamento por candidato ─────────────────────────────────────────
    linhas += ["  DETALHAMENTO POR CANDIDATO", sub, ""]

    for c in candidatos:
        if c["classe"] == "DESCARTA":
            continue
        trilha = c["trilha"]
        tx0, ty0 = trilha[0][1], trilha[0][0]
        if frames[0]["wcs_ok"]:
            ra, dec = pixel_para_radec(tx0, ty0, frames[0]["wcs"])
        else:
            ra, dec = float("nan"), float("nan")
        mpc = c.get("mpc", {"status": "nao_consultado"})
        mpc_status = mpc.get("status", "nao_consultado")

        # Status MPC legível
        if mpc_status == "match_provavel":
            match = mpc.get("match", {})
            status_mpc_txt = (
                f"MATCH PROVÁVEL — {match.get('nome_designacao', '?')}  "
                f"mag={match.get('magnitude_v', '?')}  "
                f"dist={match.get('distancia_arcmin', '?')}'  "
                f"tipo={match.get('tipo', '?')}"
            )
        elif mpc_status == "match_ambiguo":
            status_mpc_txt = (
                f"MATCH AMBÍGUO — {mpc.get('n_objetos_no_cone', '?')} "
                "objetos no cone de busca"
            )
        elif mpc_status == "sem_match":
            status_mpc_txt = "SEM CORRESPONDÊNCIA no cone de 2'"
        elif mpc_status == "consulta_falhou":
            status_mpc_txt = f"CONSULTA FALHOU — {mpc.get('motivo', '')}"
        else:
            status_mpc_txt = "Não consultado (modo offline)"

        flags_str = "  ".join(c.get("flags", [])) or "—"

        # Razões de decisão
        razoes_rej = c.get("razoes_rejeicao", [])
        razoes_pen = c.get("razoes_penalizacao", [])
        score_max  = c.get("score_max", 12)

        linhas += [
            f"  CANDIDATO #{c.get('rank', '?'):02}",
            sub,
            f"  Classificação    : {c['classe']}  "
            f"(score {c['score']}/{score_max}  |  prob {c['probabilidade']}%)",
            f"  Flags            : {flags_str}",
            f"  Diagnóstico      : {_resumo_diagnostico(c)}",
        ]

        if razoes_rej:
            linhas.append(f"  ── Razões de rejeição ───────────────────────────────────")
            for r in razoes_rej:
                linhas.append(f"  ✗  {r}")

        if razoes_pen:
            linhas.append(f"  ── Razões de penalização ────────────────────────────────")
            for r in razoes_pen:
                linhas.append(f"  ↓  {r}")

        linhas += [
            "",
            f"  ── Decomposição do score ────────────────────────────────────",
            f"  Linearidade        : {c['score_linearidade']}/3  "
            f"(resíduo {c['linearidade']:.3f} px)",
            f"  Velocidade         : {c['score_velocidade']}/2  "
            f"(σ passos {c['vel_consistencia']:.3f})",
            f"  Fotometria         : {c['score_fotometria']}/2  "
            f"(CV brilho {c['brilho_cv']:.3f})",
            f"  Morfologia         : {c['score_morfologia']}/2  "
            f"(pontualidade média {c['pontual']:.4f})",
            f"  Consist. morfológ. : {c.get('score_morfo_consist', '?')}/1  "
            f"(std pontualidade {c.get('pont_std', 0.0):.4f})",
            f"  Elongation         : {c.get('score_elongation', '?')}/1  "
            f"(elongação média {c.get('elong_medio', 1.0):.2f})",
            f"  Faixa vel. típica  : {c['score_faixa_vel']}/1  "
            f"({c['move_total']:.1f} px total)",
            f"  Score total        : {c['score']}/{score_max}",
            "",
            f"  ── Posição (Frame 1) ────────────────────────────────────────",
            f"  Pixel            : x={tx0:.1f}  y={ty0:.1f}",
            f"  Coordenadas      : RA = {ra:.6f}°  Dec = {dec:.6f}°",
            f"  RA  (HH MM SS)   : {formatar_ra_mpc(ra)}",
            f"  Dec (±DD MM SS)  : {formatar_dec_mpc(dec)}",
            "",
            f"  ── Movimento ────────────────────────────────────────────────",
            f"  Deslocamento total : {c['move_total']:.2f} px",
            f"  Direção (dx, dy)   : ({trilha[-1][1]-trilha[0][1]:.1f}, "
            f"{trilha[-1][0]-trilha[0][0]:.1f}) px",
        ]

        if c.get("vel_arcsec_min") is not None:
            linhas.append(f"  Velocidade angular : {c['vel_arcsec_min']:.2f} arcsec/min")

        fwhm_str = (f"{c['fwhm_medio_px']:.1f} px"
                    if c.get("fwhm_medio_px") is not None else "n/d")
        linhas += [
            "",
            f"  ── Fotometria e Morfologia ──────────────────────────────────",
            f"  Brilho por frame  : {[round(b, 1) for b in c['brilhos']]}",
            f"  SNR estimado      : {[round(s, 1) for s in c.get('snrs', [])]}",
            f"  Variação (CV)     : {c['brilho_cv']:.3f}  "
            + ("✓ estável" if c["brilho_cv"] < T.FOT_ESTAVEL else "! variável"),
            f"  Perfil pontual    : {c['pontual']:.4f}  "
            + ("✓ ponto" if c["pontual"] > T.MOR_MARGINAL else "! estendido"),
            f"  Elongação média   : {c.get('elong_medio', 1.0):.2f}  "
            + ("✓" if c.get("elong_medio", 1.0) < T.ELON_BOA else "! estendida"),
            f"  FWHM estimada     : {fwhm_str}",
            "",
            f"  ── Status MPC/SkyBot ────────────────────────────────────────",
            f"  {status_mpc_txt}",
            "",
            f"  ── Trilha frame a frame ─────────────────────────────────────",
        ]
        for fi, t in enumerate(trilha):
            if frames[fi]["wcs_ok"]:
                ra_f, dec_f = pixel_para_radec(t[1], t[0], frames[fi]["wcs"])
                coord_str = f"RA={ra_f:.5f}° Dec={dec_f:.5f}°"
            else:
                coord_str = "RA=?  Dec=?  (WCS inválido)"
            snr_str = f"  SNR≈{c['snrs'][fi]:.1f}" if c.get("snrs") else ""
            linhas.append(
                f"  Frame {fi+1} ({frames[fi]['date_obs'][11:19]} UTC)  "
                f"x={t[1]:.1f} y={t[0]:.1f}  flux={t[2]:.0f}  {coord_str}{snr_str}"
            )

        linhas += ["", sub, ""]

    # ── Metodologia ────────────────────────────────────────────────────────
    linhas += [
        sep,
        "  METODOLOGIA",
        sub,
        "  1. Background local via photutils.Background2D (box 50x50, filtro 3x3)",
        "  2. Sementes pontuais em imagem SNR local + ajuste PSF Moffat/Gauss",
        "     para centro sub-pixel e FWHM compatível com PSF estelar",
        "  3. Rastreamento em plano celeste local via WCS por frame",
        "     + casamento mútuo recíproco entre frames consecutivos",
        "  4. Refinamento WCS por Gaia DR3 com fallback para Pan-STARRS",
        "  5. Rejeição de falsos positivos: borda, SNR baixo, hot pixel heurístico,",
        "     brilho caótico, elongação/FWHM incompatível com fonte pontual",
        f"  6. Score 0-{3+2+2+2+1+1+1}: linearidade(3) + velocidade(2) + fotometria(2)",
        "     + morfologia(2) + consist.morfo(1) + elongação(1) + faixa_vel(1)",
        "  7. Morfologia: pontualidade, elongação via PCA, FWHM PSF, consist. frames",
        "  8. Conversão pixel→RA/Dec via WCS Pan-STARRS refinado por Gaia DR3",
        "  9. Correspondência posicional: SkyBot/IMCCE, cone 2 arcmin",
        " 10. MPC 80-colunas com magnitude omitida (requer calibração Astrometrica)",
        "  NOTA: o score é heurístico e ajustado para pré-triagem, não probabilidade",
        "        estatisticamente calibrada. Thresholds centralizados em T.* no código.",
        "",
        "  LIMITAÇÕES",
        sub,
        "  - Sem calibração fotométrica absoluta (magnitude = branco no MPC)",
        "  - WCS usa Pan-STARRS como solução inicial e Gaia DR3 quando disponível",
        "  - Score heurístico refinado: não substitui digest2 nem CNN",
        "  - SkyBot indica vizinhança posicional — ausência de match ≠ novidade",
        "  - Validação humana no Astrometrica é obrigatória antes de enviar ao IASC",
        "  - Heurísticas de falso positivo podem rejeitar candidatos legítimos em",
        "    campos com PSF degradada — revisar DESCARTA com flag REJEITADO_FP",
        "",
        f"  Pipeline : TRITON v{VERSAO}",
        "  Autora   : Jaciana Barbosa",
        sep,
    ]

    with open(caminho, "w", encoding="utf-8") as f:
        f.write("\n".join(linhas))

    return caminho


# ─────────────────────────────────────────────
# 10. EXPORTAÇÃO JSON
# ─────────────────────────────────────────────
def _serializar_trilha(c: dict, frames: list[dict]) -> list[dict]:
    """Serializa a trilha completa frame a frame para o JSON."""
    trilha_json = []
    for fi, t in enumerate(c["trilha"]):
        ty, tx, flux = t[0], t[1], t[2]
        frame = frames[fi]

        ra_deg, dec_deg = None, None
        ra_fmt, dec_fmt = None, None
        if frame["wcs_ok"]:
            try:
                ra_deg, dec_deg = pixel_para_radec(tx, ty, frame["wcs"])
                ra_fmt  = formatar_ra_mpc(ra_deg)
                dec_fmt = formatar_dec_mpc(dec_deg)
                ra_deg  = round(ra_deg, 6)
                dec_deg = round(dec_deg, 6)
            except Exception:
                pass

        snr = c.get("snrs", [None, None, None, None])[fi]
        incertezas_frame = c.get("incerteza_astrometrica", {}).get("por_frame") or []
        incerteza_frame = incertezas_frame[fi] if fi < len(incertezas_frame) else {}

        trilha_json.append({
            "frame_index"   : fi,
            "timestamp_utc" : frame["date_obs"],
            "timestamp_inicio_utc": frame.get("date_obs_inicio", frame["date_obs"]),
            "jd"            : round(frame["jd"], 6),
            "jd_inicio"     : round(frame.get("jd_inicio", frame["jd"]), 6),
            "jd_mid"        : round(frame.get("jd_mid", frame["jd"]), 6),
            "x"             : round(tx, 2),
            "y"             : round(ty, 2),
            "ra_deg"        : ra_deg,
            "dec_deg"       : dec_deg,
            "ra_fmt"        : ra_fmt,
            "dec_fmt"       : dec_fmt,
            "flux"          : round(float(flux), 1),
            "snr"           : round(float(snr), 2) if snr is not None else None,
            "incerteza_astrometrica": incerteza_frame,
            "wcs_valido"    : frame["wcs_ok"],
        })
    return trilha_json


def exportar_json(candidatos: list[dict], frames: list[dict],
                  nome_conjunto: str, pasta_output: Path,
                  metricas_globais: dict) -> Path:
    """
    Exporta JSON estruturado com:
    - métricas globais do conjunto
    - por candidato: score decomposto, trilha completa, flags, status MPC normalizado
    """
    # Cabeçalho global
    saida = {
        "pipeline"      : f"TRITON v{VERSAO}",
        "conjunto"      : nome_conjunto,
        "processado_em" : datetime.now().isoformat(timespec="seconds"),
        "observatorio"  : "Pan-STARRS / IASC (F51)",
        "data_obs"      : frames[0]["date_obs"][:10],
        "metricas_globais": metricas_globais,
        "candidatos"    : [],
    }

    for c in candidatos:
        tx0, ty0 = c["trilha"][0][1], c["trilha"][0][0]
        ra_r, dec_r, ra_fmt, dec_fmt = None, None, None, None
        if frames[0]["wcs_ok"]:
            try:
                ra_r, dec_r = pixel_para_radec(tx0, ty0, frames[0]["wcs"])
                ra_fmt  = formatar_ra_mpc(ra_r)
                dec_fmt = formatar_dec_mpc(dec_r)
                ra_r    = round(ra_r, 6)
                dec_r   = round(dec_r, 6)
            except Exception:
                pass

        entrada = {
            # Identificação
            "rank"              : c.get("rank"),
            "classe"            : c["classe"],
            "score_total"       : c["score"],
            "probabilidade"     : c["probabilidade"],

            # Score decomposto
            "score_componentes" : {
                "linearidade"       : c["score_linearidade"],
                "velocidade"        : c["score_velocidade"],
                "fotometria"        : c["score_fotometria"],
                "morfologia"        : c["score_morfologia"],
                "consist_morfologica": c.get("score_morfo_consist", 0),
                "elongation"        : c.get("score_elongation", 0),
                "faixa_vel"         : c["score_faixa_vel"],
                "maximos"           : {
                    "linearidade": 3, "velocidade": 2, "fotometria": 2,
                    "morfologia": 2, "consist_morfologica": 1,
                    "elongation": 1, "faixa_vel": 1,
                },
                "score_max"         : c.get("score_max", 12),
            },

            # Flags diagnósticas
            "flags"             : c.get("flags", []),

            # Razões de decisão auditáveis
            "razoes_decisao"    : {
                "penalizacoes" : c.get("razoes_penalizacao", []),
                "rejeicoes"    : c.get("razoes_rejeicao", []),
            },

            # Posição resumida (frame 1)
            "posicao_frame1"    : {
                "pixel_x" : round(tx0, 2),
                "pixel_y" : round(ty0, 2),
                "ra_deg"  : ra_r,
                "dec_deg" : dec_r,
                "ra_fmt"  : ra_fmt,
                "dec_fmt" : dec_fmt,
            },

            # Métricas de movimento
            "movimento"         : {
                "total_px"        : round(c["move_total"], 3),
                "residuo_deriva"  : round(c.get("residuo", 0.0), 3),
                "linearidade_px"  : round(c["linearidade"], 4),
                "vel_consistencia": round(c["vel_consistencia"], 4),
                "vel_arcsec_min"  : c.get("vel_arcsec_min"),
            },

            # Fotometria e morfologia
            "fotometria"        : {
                "brilho_por_frame": [round(b, 1) for b in c["brilhos"]],
                "brilho_cv"       : round(c["brilho_cv"], 4)
                                    if c["brilho_cv"] != float("inf") else None,
                "pontualidade"    : round(c["pontual"], 5),
                "pont_std"        : round(c.get("pont_std", 0.0), 5),
                "snr_por_frame"   : [round(s, 2) if s is not None else None
                                     for s in c.get("snrs", [])],
            },

            "incerteza_astrometrica": c.get("incerteza_astrometrica", {}),
            "gaia_static"           : c.get("gaia_static", {
                "status": "nao_avaliado",
                "n_matches": 0,
                "matches": [],
            }),

            "morfologia"        : {
                "elongation_medio": round(c.get("elong_medio", 1.0), 4),
                "fwhm_medio_px"   : c.get("fwhm_medio_px"),
                "por_frame"       : [
                    {
                        "elongation" : round(m["elongation"], 4),
                        "pontual"    : round(m["pontual"], 5),
                        "fwhm_est_px": round(float(m["fwhm_est_px"]), 4)
                                       if m["fwhm_est_px"] is not None else None,
                        "modelo_psf" : m.get("modelo_psf"),
                        "valido"     : m["valido"],
                    }
                    for m in c.get("metricas_morfo", [])
                ],
            },

            # Trilha completa frame a frame
            "trilha"            : _serializar_trilha(c, frames),

            # Status MPC normalizado
            "mpc"               : c.get("mpc", {"status": "nao_consultado"}),
        }
        saida["candidatos"].append(entrada)

    arq_json = pasta_output / f"{nome_conjunto}_candidatos.json"
    with open(arq_json, "w", encoding="utf-8") as f:
        json.dump(saida, f, indent=2, ensure_ascii=False)
    return arq_json


# ─────────────────────────────────────────────
# 11. PIPELINE PRINCIPAL
# ─────────────────────────────────────────────
def _resolver_observador(args) -> dict:
    nome = args.observador or os.environ.get("TRITON_OBS") or "Observador IASC"
    email = args.email or os.environ.get("TRITON_EMAIL") or "observer@example.com"
    return {"nome": nome, "email": email}


def _deduplicar(candidatos: list[dict], raio_px: float = None) -> list[dict]:
    if raio_px is None:
        raio_px = T.DEDUP_RAIO_PX
    unicos = []
    for c in candidatos:
        y0, x0 = c["trilha"][0][0], c["trilha"][0][1]
        dup = any(
            np.hypot(x0 - u["trilha"][0][1], y0 - u["trilha"][0][0]) < raio_px
            for u in unicos
        )
        if not dup:
            unicos.append(c)
    return unicos


def main():
    parser = argparse.ArgumentParser(
        description="TRITON — Pré-triagem de asteroides em FITS do IASC"
    )
    parser.add_argument("--imagens",    type=str, default="imagens")
    parser.add_argument("--output",     type=str, default="resultados")
    parser.add_argument("--sigma",      type=float, default=5.5)
    parser.add_argument("--sem-mpc",    action="store_true")
    parser.add_argument("--observador", type=str, default=None)
    parser.add_argument("--email",      type=str, default=None)
    args = parser.parse_args()

    t_inicio = time.monotonic()

    pasta_img = Path(args.imagens)
    pasta_out = Path(args.output)
    pasta_out.mkdir(parents=True, exist_ok=True)
    arquivos_ini = sorted(pasta_img.glob("*.fits")) + sorted(pasta_img.glob("*.fit"))
    nome_conjunto = (arquivos_ini[0].name.split("_")[0].split(".")[0]
                     if arquivos_ini else "pipeline")

    observador = _resolver_observador(args)

    print(f"\n{Cor.NEGRITO}{'='*60}")
    print(f"  TRITON v{VERSAO}")
    print(f"  Observador: {observador['nome']}  <{observador['email']}>")
    print(f"{'='*60}{Cor.RESET}\n")

    arq_log = configurar_log_arquivo(pasta_out, nome_conjunto)
    log(f"Log do pipeline: {arq_log}")
    log("Carregando imagens FITS...")
    frames = carregar_fits(pasta_img)
    nome_conjunto = frames[0]["arquivo"].split("_")[0].split(".")[0]

    log("\nDetectando objetos em movimento...")
    movedores, metricas_globais = encontrar_movedores(frames, sigma=args.sigma)
    log(f"Candidatos brutos encontrados: {len(movedores)}", "OK")

    if not movedores:
        log("Nenhum objeto em movimento detectado.", "WARN")
        sys.exit(0)

    movedores = _deduplicar(movedores)
    log(f"Após deduplicação: {len(movedores)} candidatos únicos", "OK")

    log("\nAnalisando candidatos...")
    candidatos = [analisar_candidato(m, frames) for m in movedores]
    candidatos.sort(key=lambda x: -x["score"])

    # Atribuir rank após ordenação
    for i, c in enumerate(candidatos):
        c["rank"] = i + 1

    if not args.sem_mpc:
        log("\nConsultando SkyBot/IMCCE...")
        for c in candidatos:
            if c["classe"] == "DESCARTA":
                c["mpc"] = {"status": "nao_consultado", "motivo": None, "match": None}
                continue
            if not frames[0]["wcs_ok"]:
                c["mpc"] = {"status": "consulta_falhou",
                            "motivo": "WCS inválido", "match": None}
                continue
            tx, ty  = c["trilha"][0][1], c["trilha"][0][0]
            ra, dec = pixel_para_radec(tx, ty, frames[0]["wcs"])
            log(f"  #{c['rank']}  RA={ra:.4f} Dec={dec:.4f}...")
            c["mpc"] = consultar_catalogo_posicional(ra, dec, frames[0]["date_obs"])
            status   = c["mpc"]["status"]
            if status == "match_provavel":
                log(f"  → match: {c['mpc']['match']['nome_designacao']}", "WARN")
            elif status == "match_ambiguo":
                log(f"  → ambíguo ({c['mpc']['n_objetos_no_cone']} objetos)", "WARN")
            elif status == "sem_match":
                log("  → sem correspondência", "OK")
            else:
                log(f"  → {status}: {c['mpc'].get('motivo', '')}", "WARN")
    else:
        for c in candidatos:
            c["mpc"] = {"status": "nao_consultado", "motivo": None, "match": None}

    # Adicionar flags MPC após consulta
    for c in candidatos:
        _adicionar_flag_mpc(c)

    # Completar métricas globais
    t_fim = time.monotonic()
    metricas_globais.update({
        "n_frames_validos"   : sum(1 for f in frames if f["wcs_ok"]),
        "n_candidatos_unicos": len(movedores),
        "n_candidatos_finais": len(candidatos),
        "dist_classes"       : {
            "FORTE"   : sum(1 for c in candidatos if c["classe"] == "FORTE"),
            "MODERADO": sum(1 for c in candidatos if c["classe"] == "MODERADO"),
            "FRACO"   : sum(1 for c in candidatos if c["classe"] == "FRACO"),
            "DESCARTA": sum(1 for c in candidatos if c["classe"] == "DESCARTA"),
        },
        "tempo_execucao_s"   : round(t_fim - t_inicio, 1),
    })

    # Resumo terminal
    score_max_str = str(candidatos[0].get("score_max", 12)) if candidatos else "12"
    print(f"\n{Cor.NEGRITO}{'─'*62}")
    print("  RESUMO DOS CANDIDATOS")
    print(f"{'─'*62}{Cor.RESET}")
    print(f"  {'#':>3}  {'Score':>7}  {'Prob':>5}  {'Classe':<10}  "
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
        "REJEITADO_FP"          : "FP✗",
    }

    for c in candidatos:
        if c["classe"] == "DESCARTA":
            continue
        mpc     = c.get("mpc", {})
        mpc_str = (mpc.get("match", {}) or {}).get("nome_designacao") \
                  or mpc.get("status", "?")
        flags_short = " ".join(
            _flag_abrev.get(f, f) for f in c.get("flags", [])
        )
        sm = c.get("score_max", 12)
        print(
            f"  {Cor.NEGRITO}{c['rank']:>3}{Cor.RESET}  "
            f"{c['cor_cli']}{c['score']:>4}/{sm}  "
            f"{c['probabilidade']:>4}%  "
            f"{c['classe']:<10}{Cor.RESET}  "
            f"{mpc_str[:28]:<28}  "
            f"{flags_short}"
        )

    # Gerar outputs
    log("\nGerando relatório de texto...")
    arq_txt = gerar_relatorio_txt(candidatos, frames, nome_conjunto,
                                  pasta_out, metricas_globais)
    log(f"Relatório salvo: {arq_txt}", "OK")

    log("Gerando relatório MPC...")
    arq_mpc = gerar_relatorio_mpc(candidatos, frames, nome_conjunto,
                                  pasta_out, observador)
    log(f"Relatório MPC salvo: {arq_mpc}", "OK")

    log("Gerando visualização...")
    gerar_visualizacao(candidatos, frames, pasta_out, nome_conjunto)

    log("Exportando JSON...")
    arq_json = exportar_json(candidatos, frames, nome_conjunto,
                             pasta_out, metricas_globais)
    log(f"JSON salvo: {arq_json}", "OK")

    print(f"\n{Cor.VERDE}{Cor.NEGRITO}Análise concluída. "
          f"Resultados em: {pasta_out}/{Cor.RESET}")
    print(f"   → {arq_txt.name}")
    print(f"   → {arq_mpc.name}")
    print(f"   → {nome_conjunto}_candidatos.png")
    print(f"   → {arq_json.name}")
    print(f"   → {arq_log.name}")
    print()


if __name__ == "__main__":
    main()
