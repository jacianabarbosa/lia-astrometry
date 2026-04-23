"""
=============================================================================
asteroid-hunter — detector.py
=============================================================================
Pipeline de pré-triagem de asteroides em imagens FITS do IASC/Pan-STARRS.

Autora  : Jaciana Barbosa
Versão  : 1.2.0
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
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
import astropy.units as u

from photutils.detection import DAOStarFinder

warnings.filterwarnings("ignore")

VERSAO = "1.3.0"

# ─────────────────────────────────────────────
# THRESHOLDS CENTRALIZADOS
# Edite aqui para tuning futuro sem caçar valores pelo código.
# ─────────────────────────────────────────────
class T:
    # Detecção
    SIGMA_DETECCAO       = 5.5       # sigma do DAOStarFinder (sobrescrito por --sigma)
    FWHM_PX              = 3.0       # FWHM assumido para o finder
    SHARPNESS_MIN        = 0.2
    SHARPNESS_MAX        = 1.0
    ROUNDNESS_MAX        = 0.8       # |roundness| máximo no finder

    # Associação entre frames
    RAIO_MATCH_PX        = 20.0      # janela de casamento mútuo [px]
    MAX_ANGULO_DEG       = 45.0      # desvio máximo de direção entre passos [graus]
    MAX_STEP_RATIO       = 2.5       # razão máxima entre o maior e menor passo
    MIN_STEP_PX          = 0.5       # passo mínimo para trilha com movimento real

    # Filtro de movimento (pós-deriva)
    MOVE_MIN_PX          = 2.0       # deslocamento mínimo total [px]
    MOVE_MAX_PX          = 120.0     # deslocamento máximo total [px]
    RESIDUO_MIN_PX       = 1.8       # resíduo mínimo em relação à deriva

    # Rejeição de falso positivo
    MARGEM_BORDA_PX      = 30.0      # distância mínima à borda do frame
    HOT_PIXEL_PONT_MAX   = 0.70      # pontualidade > threshold → hot pixel suspect
    HOT_PIXEL_SNR_MIN    = 3.0       # SNR mínimo para não ser descartado como hot pixel
    BRILHO_CV_CAOS       = 1.5       # CV acima disto → brilho caótico (rejeição forte)
    STEP_ACCEL_MAX       = 3.0       # aceleração máxima entre passos [px/frame²]
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

    # Deduplicação
    DEDUP_RAIO_PX        = 25.0

    # MPC: score mínimo para entrar no relatório 80-colunas
    MPC_SCORE_MIN        = 6

    # Deriva: threshold para fonte "estável"
    DERIVA_ESTAVEL_PX    = 2.0
    DERIVA_MIN_ESTAVEIS  = 5


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

_logger = logging.getLogger("asteroid-hunter")

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
    else:
        raise ValueError("Header FITS sem CD matrix nem CDELT — "
                         "escala do pixel desconhecida.")

    if "EQUINOX" in header:
        w.wcs.equinox = float(header["EQUINOX"])
    if "RADESYS" in header:
        w.wcs.radesys = str(header["RADESYS"]).strip()
    return w


def _encontrar_hdu_ciencia(hdul: fits.HDUList) -> int:
    for i, h in enumerate(hdul):
        if h.data is not None and h.data.ndim == 2:
            return i
    raise ValueError("Nenhuma HDU com imagem 2D encontrada no FITS.")


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

            date_obs = str(header.get("DATE-OBS", "")).strip()
            try:
                t  = Time(date_obs, format="isot", scale="utc")
                jd = float(t.jd)
            except Exception:
                log(f"DATE-OBS inválido ou ausente em {arq.name}: {date_obs!r}", "ERRO")
                sys.exit(1)

            frames.append({
                "arquivo" : arq.name,
                "data"    : data,
                "header"  : header,
                "wcs"     : wcs,
                "wcs_ok"  : wcs_ok,
                "date_obs": date_obs,
                "jd"      : jd,
                "exptime" : float(header.get("EXPTIME", 45.0)),
            })
        log(f"Carregado: {arq.name}  ({date_obs})  "
            f"[{data.shape[1]}x{data.shape[0]}]  WCS={'OK' if wcs_ok else 'FALHOU'}")

    if not any(f["wcs_ok"] for f in frames):
        log("Nenhum frame tem WCS válido. Impossível calcular RA/Dec.", "ERRO")
        sys.exit(1)

    frames.sort(key=lambda x: x["jd"])
    return frames


# ─────────────────────────────────────────────
# 2. DETECÇÃO DE FONTES EM CADA FRAME
# ─────────────────────────────────────────────
def detectar_fontes(img: np.ndarray, sigma: float = 5.5,
                    fwhm_px: float = 3.0) -> list[tuple]:
    """
    Detecta fontes pontuais usando photutils.DAOStarFinder.
    Retorna lista de (y, x, flux, tamanho).
    """
    mean, median, std = sigma_clipped_stats(img, sigma=3.0, maxiters=5)

    finder = DAOStarFinder(
        fwhm=fwhm_px,
        threshold=sigma * std,
        sharpness_range=(T.SHARPNESS_MIN, T.SHARPNESS_MAX),
        roundness_range=(-T.ROUNDNESS_MAX, T.ROUNDNESS_MAX),
        exclude_border=True,
    )
    tabela = finder(img - median)
    if tabela is None or len(tabela) == 0:
        return []

    x_col = "x_centroid" if "x_centroid" in tabela.colnames else "xcentroid"
    y_col = "y_centroid" if "y_centroid" in tabela.colnames else "ycentroid"

    fontes = []
    for linha in tabela:
        y    = float(linha[y_col])
        x    = float(linha[x_col])
        flux = float(linha["flux"])
        tam  = float(np.pi * (fwhm_px / 2.0) ** 2)
        fontes.append((y, x, flux, tam))
    return fontes


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


def _validar_coerencia_trilha(trilha: list[tuple]) -> tuple[bool, list[str]]:
    """
    Verifica se uma trilha de 4 posições é cinematicamente coerente.

    Retorna (valida, razoes_rejeicao).
    Heurísticas físicas aplicadas:
    - Passos individuais não devem mudar de direção abruptamente (> MAX_ANGULO_DEG)
    - Razão entre o passo maior e o menor não deve exceder MAX_STEP_RATIO
    - Aceleração entre passos consecutivos não deve exceder STEP_ACCEL_MAX px/frame²
    """
    pos = np.array([(t[1], t[0]) for t in trilha])   # (x, y)
    passos = np.diff(pos, axis=0)                      # 3 vetores de passo
    normas = np.linalg.norm(passos, axis=1)            # magnitude de cada passo

    razoes = []

    # Rejeita se algum passo for quase nulo enquanto outros são grandes
    # (evita associação com fonte estática num frame)
    if normas.max() > T.MIN_STEP_PX and normas.min() < T.MIN_STEP_PX:
        razoes.append(
            f"passo quase nulo em trilha com movimento real "
            f"(min={normas.min():.2f} px, max={normas.max():.2f} px)"
        )

    # Consistência de direção: ângulo entre passos consecutivos
    for i in range(len(passos) - 1):
        n1, n2 = normas[i], normas[i + 1]
        if n1 < 0.1 or n2 < 0.1:
            continue
        cos_ang = np.clip(
            np.dot(passos[i], passos[i + 1]) / (n1 * n2), -1.0, 1.0
        )
        angulo = float(np.degrees(np.arccos(cos_ang)))
        if angulo > T.MAX_ANGULO_DEG:
            razoes.append(
                f"mudança de direção excessiva entre passos {i+1}-{i+2}: "
                f"{angulo:.1f}° (máx {T.MAX_ANGULO_DEG}°)"
            )

    # Consistência de distância: razão max/min de passos
    if normas.min() > 0.1:
        ratio = float(normas.max() / normas.min())
        if ratio > T.MAX_STEP_RATIO:
            razoes.append(
                f"razão de passo inconsistente: {ratio:.2f}x "
                f"(máx {T.MAX_STEP_RATIO}x)"
            )

    # Aceleração: variação de passo entre frames consecutivos
    for i in range(len(normas) - 1):
        accel = abs(normas[i + 1] - normas[i])
        if accel > T.STEP_ACCEL_MAX:
            razoes.append(
                f"aceleração excessiva entre passos {i+1}-{i+2}: "
                f"{accel:.2f} px/frame (máx {T.STEP_ACCEL_MAX})"
            )

    return (len(razoes) == 0), razoes


def encontrar_movedores(frames: list[dict], sigma: float = 5.5) -> tuple[list[dict], dict]:
    """
    Compara fontes entre os 4 frames e identifica objetos com movimento
    residual (após remover a deriva do campo).

    Inclui validação cinemática de coerência de trilha.
    Retorna (movedores, metricas_globais).
    """
    todas_fontes = []
    n_fontes_por_frame = []
    for i, f in enumerate(frames):
        fontes = detectar_fontes(f["data"], sigma=sigma)
        todas_fontes.append(fontes)
        n_fontes_por_frame.append(len(fontes))
        log(f"  Frame {i+1} ({f['date_obs'][:19]}): {len(fontes)} fontes")

    pares_12 = _casar_mutuo(todas_fontes[0], todas_fontes[1])
    pares_23 = _casar_mutuo(todas_fontes[1], todas_fontes[2])
    pares_34 = _casar_mutuo(todas_fontes[2], todas_fontes[3])

    candidatos_brutos = []
    n_rejeitados_coerencia = 0

    for i0, i1 in pares_12.items():
        if i1 not in pares_23:
            continue
        i2 = pares_23[i1]
        if i2 not in pares_34:
            continue
        i3 = pares_34[i2]

        f0 = todas_fontes[0][i0]
        f1 = todas_fontes[1][i1]
        f2 = todas_fontes[2][i2]
        f3 = todas_fontes[3][i3]
        trilha = [(f0[0], f0[1], f0[2]),
                  (f1[0], f1[1], f1[2]),
                  (f2[0], f2[1], f2[2]),
                  (f3[0], f3[1], f3[2])]

        # Validação cinemática antes de aceitar a trilha
        valida, razoes_cinem = _validar_coerencia_trilha(trilha)
        if not valida:
            n_rejeitados_coerencia += 1
            log(f"  Trilha rejeitada por incoerência: {'; '.join(razoes_cinem)}", "WARN")
            continue

        pos  = np.array([(t[0], t[1]) for t in trilha])
        move = float(np.hypot(pos[-1, 0] - pos[0, 0], pos[-1, 1] - pos[0, 1]))
        candidatos_brutos.append({
            "trilha"    : trilha,
            "move_total": move,
            "brilhos"   : [t[2] for t in trilha],
        })

    n_trilhas_tentadas = len(candidatos_brutos) + n_rejeitados_coerencia

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

    if n_rejeitados_coerencia:
        log(f"Trilhas rejeitadas por incoerência cinemática: {n_rejeitados_coerencia}", "WARN")

    metricas = {
        "n_fontes_por_frame"          : n_fontes_por_frame,
        "n_trilhas_tentadas"          : n_trilhas_tentadas,
        "n_rejeitados_coerencia"      : n_rejeitados_coerencia,
        "deriva_dx_px"                : round(deriva_dx, 3),
        "deriva_dy_px"                : round(deriva_dy, 3),
        "n_estaveis_deriva"           : n_estaveis,
        "sigma_deteccao"              : sigma,
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
    y0 = max(0, int(ty) - r); y1 = min(ny, int(ty) + r)
    x0 = max(0, int(tx) - r); x1 = min(nx, int(tx) + r)
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

    return {
        "pontual"    : pontual,
        "elongation" : elongation,
        "fwhm_est_px": round(fwhm_est, 2) if fwhm_est is not None else None,
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

    return razoes


def _calcular_snr_local(img: np.ndarray, ty: float, tx: float,
                        raio_sinal: int = 8, raio_bg: int = 20) -> float:
    """
    SNR local simples: sinal = soma de pixels no raio_sinal após subtrair
    background estimado em annulus entre raio_sinal e raio_bg.
    """
    ny, nx = img.shape
    y0s = max(0, int(ty) - raio_sinal); y1s = min(ny, int(ty) + raio_sinal)
    x0s = max(0, int(tx) - raio_sinal); x1s = min(nx, int(tx) + raio_sinal)
    corte_sinal = img[y0s:y1s, x0s:x1s]

    y0b = max(0, int(ty) - raio_bg); y1b = min(ny, int(ty) + raio_bg)
    x0b = max(0, int(tx) - raio_bg); x1b = min(nx, int(tx) + raio_bg)
    corte_bg = img[y0b:y1b, x0b:x1b]

    bg_med  = float(np.median(corte_bg))
    bg_std  = float(np.std(corte_bg))
    sinal   = float(np.sum(corte_sinal - bg_med))
    n_pix   = corte_sinal.size
    ruido   = bg_std * np.sqrt(n_pix) if n_pix > 0 else 1.0
    return float(sinal / ruido) if ruido > 0 else 0.0


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
    pos       = np.array([(t[1], t[0]) for t in trilha])   # (x, y)
    t_idx     = np.arange(4)
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
        y0_ = max(0, int(ty) - r); y1_ = min(frame["data"].shape[0], int(ty) + r)
        x0_ = max(0, int(tx) - r); x1_ = min(frame["data"].shape[1], int(tx) + r)
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
            data_mpc = formatar_data_mpc(frame["jd"]
                                         + frame["exptime"] / 2.0 / 86400.0)
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
        f.write(f"ACK asteroid-hunter pipeline v{VERSAO}\n")
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
                    f"Score {cand['score']}/10  ({cand['probabilidade']}%)",
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
                    f"#{ri+1}  {cand['score']}/10",
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
        f"asteroid-hunter v{VERSAO} | {nome_conjunto} | "
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
        "  ASTEROID-HUNTER — Relatório de Análise",
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
        f"  FORTE    (8-10) : {sum(1 for c in candidatos if c['classe']=='FORTE')}",
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
                f"  │  {rank_label}  score {c['score']}/10  {c['classe']}",
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
        "  1. Detecção via photutils.DAOStarFinder + sigma-clipping (3σ background)",
        "  2. Rastreamento com casamento mútuo recíproco de vizinho mais próximo",
        "     + validação cinemática: ângulo de direção, razão de passo, aceleração",
        "  3. Deriva instrumental estimada na mediana de fontes estáveis",
        "  4. Rejeição de falsos positivos: borda, SNR baixo, hot pixel heurístico,",
        "     brilho caótico, elongação grave sistemática",
        f"  5. Score 0-{3+2+2+2+1+1+1}: linearidade(3) + velocidade(2) + fotometria(2)",
        "     + morfologia(2) + consist.morfo(1) + elongação(1) + faixa_vel(1)",
        "  6. Morfologia: pontualidade, elongação via PCA, FWHM estimada, consist. frames",
        "  7. Conversão pixel→RA/Dec via WCS primário Pan-STARRS",
        "  8. Correspondência posicional: SkyBot/IMCCE, cone 2 arcmin",
        "  9. MPC 80-colunas com magnitude omitida (requer calibração Astrometrica)",
        "  NOTA: o score é heurístico e ajustado para pré-triagem, não probabilidade",
        "        estatisticamente calibrada. Thresholds centralizados em T.* no código.",
        "",
        "  LIMITAÇÕES",
        sub,
        "  - Sem calibração fotométrica absoluta (magnitude = branco no MPC)",
        "  - WCS do Pan-STARRS ignora correções polinomiais proprietárias",
        "  - Score heurístico refinado: não substitui digest2 nem CNN",
        "  - SkyBot indica vizinhança posicional — ausência de match ≠ novidade",
        "  - Validação humana no Astrometrica é obrigatória antes de enviar ao IASC",
        "  - Heurísticas de falso positivo podem rejeitar candidatos legítimos em",
        "    campos com PSF degradada — revisar DESCARTA com flag REJEITADO_FP",
        "",
        f"  Pipeline : asteroid-hunter v{VERSAO}",
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

        trilha_json.append({
            "frame_index"   : fi,
            "timestamp_utc" : frame["date_obs"],
            "jd"            : round(frame["jd"], 6),
            "x"             : round(tx, 2),
            "y"             : round(ty, 2),
            "ra_deg"        : ra_deg,
            "dec_deg"       : dec_deg,
            "ra_fmt"        : ra_fmt,
            "dec_fmt"       : dec_fmt,
            "flux"          : round(float(flux), 1),
            "snr"           : round(float(snr), 2) if snr is not None else None,
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
        "pipeline"      : f"asteroid-hunter v{VERSAO}",
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

            "morfologia"        : {
                "elongation_medio": round(c.get("elong_medio", 1.0), 4),
                "fwhm_medio_px"   : c.get("fwhm_medio_px"),
                "por_frame"       : [
                    {
                        "elongation" : round(m["elongation"], 4),
                        "pontual"    : round(m["pontual"], 5),
                        "fwhm_est_px": m["fwhm_est_px"],
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
    nome  = args.observador or os.environ.get("ASTEROID_HUNTER_OBS",  "Observador IASC")
    email = args.email      or os.environ.get("ASTEROID_HUNTER_EMAIL", "observer@example.com")
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
        description="asteroid-hunter — Pré-triagem de asteroides em FITS do IASC"
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

    observador = _resolver_observador(args)

    print(f"\n{Cor.NEGRITO}{'='*60}")
    print(f"  ASTEROID-HUNTER v{VERSAO}")
    print(f"  Observador: {observador['nome']}  <{observador['email']}>")
    print(f"{'='*60}{Cor.RESET}\n")

    log("Carregando imagens FITS...")
    frames = carregar_fits(pasta_img)
    nome_conjunto = frames[0]["arquivo"].split("_")[0].split(".")[0]

    arq_log = configurar_log_arquivo(pasta_out, nome_conjunto)
    log(f"Log do pipeline: {arq_log}")

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
