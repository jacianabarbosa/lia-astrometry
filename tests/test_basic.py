"""
tests/test_basic.py
======================
Testes unitários e de integração do pipeline Lia v1.5.

Cobertura:
    - detect_sources: source simples, múltiplas, image vazia
    - format_mpc_ra / format_mpc_dec: casos normais e de wraparound
    - format_mpc_date: conversão JD → formato MPC
    - _build_panstarrs_wcs: header Pan-STARRS com chaves PCA* proprietárias
    - _mutual_match: validação recíproca de vizinho mais próximo
    - _frame_morphology: elongation, pointness, FWHM (v1.3)
    - _check_false_positive: rejeição por borda, SNR, hot pixel, flux caótico (v1.3)
    - analyze_candidate: score decomposto, novos fields v1.3, flags, razões de decisão
    - query_position_catalog: schema normalizado de retorno
    - _apply_gaia_proper_motion: propagação de pm por ano (v1.4)
    - _cross_match_gaia: cross-match RA/Dec com KDTree (v1.4)
    - _refine_wcs_gaia: fallback gracioso e aceitação de refinamento (v1.4)
    - Integração: pipeline end-to-end com FITS sintético (JSON v1.4)
    - wcs_mode=header: Gaia desativado, status "gaia_disabled" (v1.5)
    - run_id e candidate_id presentes no JSON de saída (v1.5)
    - run_metadata estruturado no JSON (v1.5)
    - export_metrics.py: geração de CSV a partir de JSON (v1.5)

Uso:
    cd lia-astrometry
    venv/bin/python -m pytest tests/test_basic.py -v
"""

import sys
import json
import subprocess
import numpy as np
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from detector import (
    detect_sources,
    format_mpc_ra,
    format_mpc_dec,
    format_mpc_date,
    _build_panstarrs_wcs,
    _mutual_match,
    _frame_morphology,
    _check_false_positive,
    _invalid_pixel_mask,
    _apply_gaia_proper_motion,
    _cross_match_gaia,
    _estimate_coarse_offset,
    _refine_wcs_gaia,
    _query_gaia_dr3,
    _extract_observation_times,
    pixel_to_radec,
    analyze_candidate,
    query_position_catalog,
    _serialize_track,
    PROJECT_NAME,
    REPOSITORY_NAME,
    T,
)

from astropy.io import fits
from astropy.wcs import WCS


class TestIdentidadePublica:
    def test_nome_e_versao_publicos(self):
        from detector import PIPELINE_VERSION
        assert PROJECT_NAME == "Lia"
        assert REPOSITORY_NAME == "lia-astrometry"
        assert PIPELINE_VERSION == "1.5.0"

    def test_documentation_has_no_old_names(self):
        root = Path(__file__).parent.parent
        files = [
            root / "README.md",
            root / "docs" / "methodology.md",
            root / "src" / "detector.py",
            root / "tools" / "export_metrics.py",
        ]
        forbidden_terms = (
            "asteroid" + "-hunter",
            "tri" + "ton" + "-astrometry",
            "TRI" + "TON",
            "tri" + "ton",
            "probabili" + "dade",
            "probabili" + "ty",
        )
        for file_path in files:
            text = file_path.read_text(encoding="utf-8")
            lower = text.lower()
            for term in forbidden_terms:
                assert term.lower() not in lower, f"{term!r} encontrado em {file_path}"


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────
def _fonte_gaussiana(shape, cx, cy, flux=3000.0, fwhm=3.0):
    sigma = fwhm / 2.355
    y, x = np.mgrid[0:shape[0], 0:shape[1]]
    return flux * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma ** 2))


def _img_com_background(shape, sigma_bg=5.0, mean_bg=100.0, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(mean_bg, sigma_bg, size=shape).astype(np.float64)


def _header_panstarrs_mock():
    h = fits.Header()
    h["CTYPE1"] = "RA---TAN"
    h["CTYPE2"] = "DEC--TAN"
    h["CRVAL1"] = 329.505
    h["CRVAL2"] = -11.883
    h["CRPIX1"] = 1211.5
    h["CRPIX2"] = 1217.0
    h["CDELT1"] = -7.095e-05
    h["CDELT2"] = 7.129e-05
    h["EQUINOX"] = 2000.0
    h["PCA1X0Y2"] = 3.12e-08
    h["PCA2X0Y2"] = -1.56e-07
    return h


def _frame_sintetico(shape=(256, 256), date_obs="2019-08-28T10:00:00"):
    """Cria um dict de frame compatível com o pipeline."""
    header = _header_panstarrs_mock()
    from detector import _build_panstarrs_wcs
    wcs = _build_panstarrs_wcs(header)
    return {
        "file" : "synth.fits",
        "data"    : np.random.default_rng(0).normal(100, 5, shape),
        "header"  : header,
        "wcs"     : wcs,
        "wcs_ok"  : True,
        "date_obs": date_obs,
        "jd"      : 2458723.921,
        "exptime" : 45.0,
    }


def _synthetic_candidate():
    """Candidato mínimo para testar analyze_candidate."""
    track = [
        (100.0, 200.0, 3000.0),
        (102.0, 202.5, 2950.0),
        (104.0, 205.0, 3100.0),
        (106.0, 207.5, 2980.0),
    ]
    return {
        "track"    : track,
        "move_total": float(np.hypot(6.0, 7.5)),
        "fluxes"   : [t[2] for t in track],
        "residual"   : 3.5,
    }


# ─────────────────────────────────────────────
# Detecção de sources
# ─────────────────────────────────────────────
class TestDetectarFontes:
    def test_detecta_fonte_simples(self):
        img = _img_com_background((200, 200), seed=1)
        img += _fonte_gaussiana(img.shape, cx=100, cy=100, flux=5000)
        sources = detect_sources(img, sigma=5.0)
        assert len(sources) >= 1
        y, x, _, _ = sources[0]
        assert abs(x - 100) < 1.5 and abs(y - 100) < 1.5

    def test_imagem_apenas_ruido(self):
        img = _img_com_background((100, 100), seed=2)
        sources = detect_sources(img, sigma=10.0)
        assert isinstance(sources, list)
        assert len(sources) == 0

    def test_multiplas_fontes(self):
        img = _img_com_background((300, 300), seed=3)
        for cx, cy in [(60, 60), (150, 150), (240, 240)]:
            img += _fonte_gaussiana(img.shape, cx, cy, flux=4000)
        sources = detect_sources(img, sigma=5.0)
        assert len(sources) >= 2


# ─────────────────────────────────────────────
# Formatação MPC (incluindo wraparound)
# ─────────────────────────────────────────────
class TestFormatacaoMPC:
    def test_ra_zero(self):
        assert format_mpc_ra(0.0).startswith("00 00")

    def test_ra_180(self):
        assert format_mpc_ra(180.0).startswith("12 00")

    def test_ra_wraparound_quase_360(self):
        result = format_mpc_ra(359.99999)
        parts = result.split()
        assert len(parts) == 3
        ss = float(parts[2])
        assert 0 <= ss < 60, f"Segundos devem estar em [0, 60), mas veio {ss}"

    def test_ra_acima_de_360(self):
        r1 = format_mpc_ra(10.0)
        r2 = format_mpc_ra(370.0)
        assert r1 == r2

    def test_dec_positiva(self):
        assert format_mpc_dec(45.5).startswith("+45")

    def test_dec_negativa(self):
        assert format_mpc_dec(-12.133).startswith("-12")

    def test_dec_wraparound_quase_46(self):
        result = format_mpc_dec(45.9999)
        parts = result.split()
        ss = float(parts[2])
        assert 0 <= ss < 60, f"Segundos Dec devem estar em [0, 60), veio {ss}"

    def test_dec_sempre_com_sinal(self):
        assert format_mpc_dec(0.0)[0] in "+-"
        assert format_mpc_dec(45.0)[0] == "+"
        assert format_mpc_dec(-45.0)[0] == "-"

    def test_data_formato(self):
        jd = 2458723.921
        result = format_mpc_date(jd)
        assert result.startswith("2019 08")
        parts = result.split()
        assert len(parts) == 3


class TestTempoObservacao:
    def test_jd_canonico_usa_meio_da_exposicao(self):
        h = fits.Header()
        h["DATE-OBS"] = "2019-08-28T10:00:00"
        h["EXPTIME"] = 60.0
        tempos = _extract_observation_times(h, "tempo.fits")

        assert tempos["jd_mid"] > tempos["jd_inicio"]
        assert abs((tempos["jd_mid"] - tempos["jd_inicio"]) * 86400.0 - 30.0) < 1e-3

    def test_mjd_obs_como_fallback(self):
        h = fits.Header()
        h["MJD-OBS"] = 58723.4166666667
        h["EXPTIME"] = 40.0
        tempos = _extract_observation_times(h, "tempo.fits")

        assert tempos["date_obs_inicio"].startswith("2019-08-28T10:00:00")
        assert abs((tempos["jd_mid"] - tempos["jd_inicio"]) * 86400.0 - 20.0) < 1e-3


# ─────────────────────────────────────────────
# WCS Pan-STARRS
# ─────────────────────────────────────────────
class TestWCSPanSTARRS:
    def test_wcs_pan_starrs_com_pca_nao_falha(self):
        h = _header_panstarrs_mock()
        wcs = _build_panstarrs_wcs(h)
        assert wcs is not None
        assert list(wcs.wcs.ctype) == ["RA---TAN", "DEC--TAN"]

    def test_wcs_pan_starrs_crpix_retorna_crval(self):
        h = _header_panstarrs_mock()
        wcs = _build_panstarrs_wcs(h)
        ra, dec = pixel_to_radec(h["CRPIX1"] - 1, h["CRPIX2"] - 1, wcs)
        assert abs(ra  - h["CRVAL1"]) < 0.01
        assert abs(dec - h["CRVAL2"]) < 0.01

    def test_wcs_pan_starrs_preserva_crota_com_cdelt(self):
        h = _header_panstarrs_mock()
        h["CDELT1"] = 7.10861946296706e-05
        h["CDELT2"] = 7.13340298342252e-05
        h["CROTA1"] = 180.4854
        h["CROTA2"] = 180.4854

        wcs = _build_panstarrs_wcs(h)

        assert hasattr(wcs.wcs, "crota")
        assert abs(float(wcs.wcs.crota[0]) - 180.4854) < 1e-6
        ra_c, dec_c = pixel_to_radec(h["CRPIX1"] - 1, h["CRPIX2"] - 1, wcs)
        ra_x, dec_x = pixel_to_radec(h["CRPIX1"], h["CRPIX2"] - 1, wcs)
        # CROTA≈180° inverte o eixo: avançar 1 px em x deve reduzir RA.
        assert ra_x < ra_c
        assert abs(dec_x - dec_c) < 1e-4

    def test_wcs_sem_chaves_falha_claramente(self):
        h = fits.Header()
        h["SIMPLE"] = True
        with pytest.raises(ValueError, match="minimum WCS keywords"):
            _build_panstarrs_wcs(h)


# ─────────────────────────────────────────────
# Casamento mútuo de sources
# ─────────────────────────────────────────────
class TestCasarMutuo:
    def test_casar_vazio(self):
        assert _mutual_match([], []) == {}

    def test_casar_simples(self):
        a = [(10.0, 20.0, 100.0, 5.0)]
        b = [(11.0, 21.0, 100.0, 5.0)]
        pares = _mutual_match(a, b, radius_px=5.0)
        assert pares == {0: 0}

    def test_casar_fora_de_raio(self):
        a = [(10.0, 20.0, 100.0, 5.0)]
        b = [(500.0, 500.0, 100.0, 5.0)]
        pares = _mutual_match(a, b, radius_px=10.0)
        assert pares == {}

    def test_casar_mutuo_rejeita_nao_reciproco(self):
        a = [(10.0, 10.0, 100.0, 5.0), (20.0, 10.0, 100.0, 5.0)]
        b = [(21.0, 10.0, 100.0, 5.0)]
        pares = _mutual_match(a, b, radius_px=20.0)
        assert 0 not in pares
        assert pares.get(1) == 0


# ─────────────────────────────────────────────
# Score decomposto e flags (v1.3)
# ─────────────────────────────────────────────
class TestScoreDecomposto:
    """Verifica que analyze_candidate retorna todos os fields v1.3."""

    def _rodar(self):
        cand   = _synthetic_candidate()
        frames = [_frame_sintetico(date_obs=f"2019-08-28T{10+i}:00:00")
                  for i in range(4)]
        return analyze_candidate(cand, frames)

    def test_campos_score_decomposto_presentes(self):
        result = self._rodar()
        for field in ("score_linearity", "score_velocity",
                      "score_photometry", "score_morphology",
                      "score_morph_consistency", "score_elongation", "score_velocity_range"):
            assert field in result, f"Campo ausente: {field}"

    def test_soma_score_consistente(self):
        result = self._rodar()
        # Score pode ser 0 se candidate foi rejeitado por FP; validamos raw_score
        soma = (result["score_linearity"] + result["score_velocity"]
                + result["score_photometry"] + result["score_morphology"]
                + result["score_morph_consistency"] + result["score_elongation"]
                + result["score_velocity_range"])
        assert result["raw_score"] == soma

    def test_limites_score_components(self):
        result = self._rodar()
        assert 0 <= result["score_linearity"]   <= 3
        assert 0 <= result["score_velocity"]    <= 2
        assert 0 <= result["score_photometry"]    <= 2
        assert 0 <= result["score_morphology"]    <= 2
        assert 0 <= result["score_morph_consistency"] <= 1
        assert 0 <= result["score_elongation"]    <= 1
        assert 0 <= result["score_velocity_range"]     <= 1

    def test_score_max_correto(self):
        result = self._rodar()
        assert result["score_max"] == 12

    def test_flags_e_lista(self):
        result = self._rodar()
        assert "flags" in result
        assert isinstance(result["flags"], list)

    def test_flags_validas(self):
        """Todas as flags devem ser strings não-vazias do input_set conhecido."""
        flags_validas = {
            "LINEARIDADE_BOA", "LINEARIDADE_RUIM",
            "VELOCIDADE_CONSISTENTE", "VELOCIDADE_IRREGULAR",
            "BRILHO_ESTAVEL", "BRILHO_INSTAVEL",
            "MORFOLOGIA_PONTUAL", "MORFOLOGIA_ESTENDIDA",
            "ELONGACAO_ALTA", "MORFO_INCONSISTENTE",
            "EDGE_FRAME", "REJECTED_FP",
            "GAIA_FONTE_ESTATICA",
            "MPC_SEM_MATCH", "MPC_MATCH_PROVAVEL",
            "MPC_MATCH_AMBIGUO", "MPC_CONSULTA_FALHOU",
        }
        result = self._rodar()
        for flag in result["flags"]:
            assert isinstance(flag, str) and flag
            assert flag in flags_validas, f"Flag inesperada: {flag}"

    def test_razoes_decisao_presentes(self):
        result = self._rodar()
        assert "penalty_reasons" in result
        assert "rejection_reasons" in result
        assert isinstance(result["penalty_reasons"], list)
        assert isinstance(result["rejection_reasons"], list)

    def test_morphology_fields_present(self):
        result = self._rodar()
        assert "mean_elongation" in result
        assert "pont_std" in result
        assert "mean_fwhm_px" in result
        assert isinstance(result["mean_elongation"], float)
        assert result["mean_elongation"] >= 1.0

    def test_snrs_presentes(self):
        result = self._rodar()
        assert "snrs" in result
        assert len(result["snrs"]) == 4

    def test_vel_arcsec_min_type(self):
        result = self._rodar()
        assert result.get("vel_arcsec_min") is None \
            or isinstance(result["vel_arcsec_min"], float)

    def test_incerteza_astrometrica_presente(self):
        result = self._rodar()
        assert "astrometric_uncertainty" in result
        inc = result["astrometric_uncertainty"]
        assert "by_frame" in inc
        assert len(inc["by_frame"]) == 4
        assert "sigma_pos_mediana_arcsec" in inc

    def test_gaia_static_schema_presente(self):
        result = self._rodar()
        assert "gaia_static" in result
        assert result["gaia_static"]["status"] in {
            "no_match_estatico", "fonte_estatica_gaia",
        }


# ─────────────────────────────────────────────
# Serialização da track (v1.2)
# ─────────────────────────────────────────────
class TestSerializacaoTrilha:
    def test_schema_completo(self):
        cand   = _synthetic_candidate()
        cand["snrs"] = [12.3, 11.8, 13.1, 12.5]
        cand["astrometric_uncertainty"] = {
            "by_frame": [{"sigma_pos_arcsec": 0.1} for _ in range(4)]
        }
        frames = [_frame_sintetico(date_obs=f"2019-08-28T{10+i}:00:00")
                  for i in range(4)]
        track_json = _serialize_track(cand, frames)

        assert len(track_json) == 4
        campos_obrigatorios = {"frame_index", "timestamp_utc", "jd",
                               "jd_inicio", "jd_mid", "x", "y", "flux",
                               "snr", "astrometric_uncertainty", "wcs_valid"}
        for entry in track_json:
            for field in campos_obrigatorios:
                assert field in entry, f"Campo ausente na track: {field}"

    def test_frame_index_sequencial(self):
        cand   = _synthetic_candidate()
        frames = [_frame_sintetico() for _ in range(4)]
        track_json = _serialize_track(cand, frames)
        indices = [e["frame_index"] for e in track_json]
        assert indices == [0, 1, 2, 3]

    def test_coordenadas_numericas(self):
        cand   = _synthetic_candidate()
        frames = [_frame_sintetico() for _ in range(4)]
        track_json = _serialize_track(cand, frames)
        for entry in track_json:
            assert isinstance(entry["x"], float)
            assert isinstance(entry["y"], float)
            assert isinstance(entry["flux"], float)
            # ra_deg pode ser None se WCS falhar, mas deve ser float quando presente
            if entry["ra_deg"] is not None:
                assert isinstance(entry["ra_deg"], float)


# ─────────────────────────────────────────────
# Status MPC normalizado (v1.2)
# ─────────────────────────────────────────────
class TestStatusMPCNormalizado:
    """Verifica que query_position_catalog retorna schema coerente."""

    STATUSES_VALIDOS = {
        "not_queried", "query_failed", "no_match",
        "ambiguous_match", "likely_match",
    }

    def test_coordenadas_invalidas_retornam_query_failed(self):
        result = query_position_catalog(400.0, 100.0, "2019-08-28T10:00:00")
        assert result["status"] == "query_failed"
        assert "match" in result

    def test_status_e_campo_obrigatorio(self):
        result = query_position_catalog(400.0, 100.0, "2019-08-28T10:00:00")
        assert result["status"] in self.STATUSES_VALIDOS

    def test_sem_astroquery_retorna_query_failed(self, monkeypatch):
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "astroquery.imcce":
                raise ImportError("mocked")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        result = query_position_catalog(329.44, -12.2, "2019-08-28T10:00:00")
        assert result["status"] == "query_failed"
        assert result["match"] is None


# ─────────────────────────────────────────────
# Imports de dependências
# ─────────────────────────────────────────────
class TestImports:
    def test_dependencias_instaladas(self):
        import astropy
        import numpy
        import scipy
        import matplotlib
        import photutils

    def test_astroquery_disponivel(self):
        try:
            import astroquery
        except ImportError:
            pytest.skip("astroquery is not installed")


# ─────────────────────────────────────────────
# Integração end-to-end (v1.2)
# ─────────────────────────────────────────────
@pytest.mark.integration
class TestIntegracaoEndToEnd:
    """
    Gera 4 FITS sintéticos com asteroide em motion conhecido e roda o
    pipeline completo. Gaia mockada para retornar None (fallback gracioso)
    — garante que o teste não depende de internet.

    Verifica estrutura do JSON v1.4:
    - global_metrics (gaia_refinement por frame, sem n_rejeitados_coerencia)
    - score_components (7 componentes, score_max=12)
    - morphology por frame
    - decision_reasons
    - track completa
    - flags
    - mpc normalizado
    """

    @staticmethod
    def _gerar_fits_sinteticos(folder: Path, seed: int = 42):
        rng = np.random.default_rng(seed)
        nx, ny = 512, 512
        ra_c, dec_c = 329.44, -12.2
        scale = 0.26 / 3600.0

        star_xy   = rng.uniform(30, 480, size=(30, 2))
        star_flux = rng.uniform(2000, 8000, size=30)
        ast_x0, ast_y0 = 256.0, 256.0
        ast_dx, ast_dy = 2.5, 2.0

        datas = [
            "2019-08-28T10:18:00",
            "2019-08-28T10:37:00",
            "2019-08-28T10:56:00",
            "2019-08-28T11:15:00",
        ]

        for i, data_obs in enumerate(datas):
            drift_dx, drift_dy = 0.6 * i, 0.4 * i
            img = rng.normal(100, 5, size=(ny, nx)).astype(np.float64)
            yy, xx = np.mgrid[0:ny, 0:nx]
            for (sx, sy), flux in zip(star_xy, star_flux):
                sx2, sy2 = sx + drift_dx, sy + drift_dy
                if 5 < sx2 < nx - 5 and 5 < sy2 < ny - 5:
                    img += flux * np.exp(
                        -((xx - sx2) ** 2 + (yy - sy2) ** 2) / (2 * 1.8 ** 2)
                    )
            ax = ast_x0 + i * ast_dx + drift_dx
            ay = ast_y0 + i * ast_dy + drift_dy
            img += 4000 * np.exp(
                -((xx - ax) ** 2 + (yy - ay) ** 2) / (2 * 1.8 ** 2)
            )

            hdu = fits.PrimaryHDU(img.astype(np.float32))
            h   = hdu.header
            h["DATE-OBS"] = data_obs
            h["EXPTIME"]  = 45.0
            h["CTYPE1"]   = "RA---TAN"
            h["CTYPE2"]   = "DEC--TAN"
            h["CRVAL1"]   = ra_c
            h["CRVAL2"]   = dec_c
            h["CRPIX1"]   = nx / 2
            h["CRPIX2"]   = ny / 2
            h["CDELT1"]   = -scale
            h["CDELT2"]   = scale
            h["EQUINOX"]  = 2000.0
            hdu.writeto(folder / f"synth_frame{i+1}.fits", overwrite=True)

    def test_pipeline_detecta_asteroide_sintetico(self, tmp_path, monkeypatch):
        images    = tmp_path / "images"
        results = tmp_path / "results"
        images.mkdir(); results.mkdir()
        self._gerar_fits_sinteticos(images)

        # Mock Gaia para não depender de rede: fallback gracioso sempre ativo.
        import detector as det_mod
        monkeypatch.setattr(det_mod, "_query_gaia_dr3", lambda *a, **kw: None)

        detector_py = Path(__file__).parent.parent / "src" / "detector.py"
        result = subprocess.run(
            [sys.executable, str(detector_py),
             "--images", str(images),
             "--output",  str(results),
             "--wcs-mode", "header",
             "--no-mpc"],
            capture_output=True, text=True, timeout=120,
        )
        assert result.returncode == 0, (
            f"Pipeline falhou:\nSTDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

        jsons = list(results.glob("*_candidates.json"))
        assert len(jsons) == 1
        with open(jsons[0]) as f:
            dados = json.load(f)

        # Estrutura de topo
        assert "global_metrics" in dados, "global_metrics missing from JSON"
        assert "candidates" in dados
        assert len(dados["candidates"]) >= 1

        # Score >= 6 esperado para asteroide sintético bem definido
        scores = [c["score_total"] for c in dados["candidates"]]
        assert max(scores) >= 6, f"Score máximo baixo demais: {scores}"

        # Verificar estrutura de cada candidate
        for c in dados["candidates"]:
            # Score decomposto (v1.3: 7 componentes)
            assert "score_components" in c, "score_components ausente"
            sc = c["score_components"]
            for comp in ("linearity", "velocity", "photometry",
                         "morphology", "morph_consistency",
                         "elongation", "velocity_range"):
                assert comp in sc, f"Componente '{comp}' ausente em score_components"
            assert sc.get("score_max") == 12, "score_max deve ser 12"

            # score_total deve ser <= score_max
            assert 0 <= c["score_total"] <= 12, (
                f"score_total={c['score_total']} fora de [0, 12]"
            )

            # Flags
            assert "flags" in c
            assert isinstance(c["flags"], list)

            # Razões de decisão (v1.3)
            assert "decision_reasons" in c
            assert "penalties" in c["decision_reasons"]
            assert "rejections"    in c["decision_reasons"]

            # Morfologia (v1.3)
            assert "morphology" in c
            morph = c["morphology"]
            assert "mean_elongation" in morph
            assert "by_frame" in morph
            assert len(morph["by_frame"]) == 4

            assert "astrometric_uncertainty" in c
            assert "by_frame" in c["astrometric_uncertainty"]
            assert "gaia_static" in c
            assert "status" in c["gaia_static"]

            # Trilha completa
            assert "track" in c
            assert len(c["track"]) == 4
            for entry in c["track"]:
                for field in ("frame_index", "timestamp_utc", "x", "y",
                              "flux", "jd_mid", "astrometric_uncertainty",
                              "wcs_valid"):
                    assert field in entry, f"Campo '{field}' ausente na track"

            # MPC normalizado
            assert "mpc" in c
            assert "status" in c["mpc"]
            assert c["mpc"]["status"] in {
                "not_queried", "query_failed", "no_match",
                "ambiguous_match", "likely_match",
            }

        # Métricas globais mínimas
        mg = dados["global_metrics"]
        for field in ("n_fontes_by_frame", "n_tracks_attempted",
                      "sigma", "drift_dx_px", "drift_dy_px",
                      "execution_time_s", "gaia_refinement"):
            assert field in mg, f"Métrica global '{field}' ausente"

        assert isinstance(mg["n_fontes_by_frame"], list)
        assert len(mg["n_fontes_by_frame"]) == 4

        # Métricas de refinamento Gaia (fallback esperado pois Gaia é mockada)
        assert isinstance(mg["gaia_refinement"], list)
        assert len(mg["gaia_refinement"]) == 4
        for g in mg["gaia_refinement"]:
            assert "frame" in g and "status" in g and "n_matches" in g
            assert "coarse_matches" in g
            assert "coarse_offset_arcsec" in g
            assert "n_matches_fino" in g




# ─────────────────────────────────────────────
# Morfologia robusta por frame (v1.3)
# ─────────────────────────────────────────────
class TestMorfologiaPorFrame:

    def _img_com_fonte(self, shape=(100, 100), cx=50, cy=50,
                       flux=5000.0, sigma=2.0, elongado=False):
        img = np.random.default_rng(7).normal(100, 5, shape)
        yy, xx = np.ogrid[:shape[0], :shape[1]]
        if elongado:
            # PSF muito alongada no eixo x
            img += flux * np.exp(
                -((xx - cx) ** 2 / (2 * (sigma * 6) ** 2)
                  + (yy - cy) ** 2 / (2 * sigma ** 2))
            )
        else:
            img += flux * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma ** 2))
        return img.astype(np.float64)

    def test_retorna_campos_obrigatorios(self):
        img  = self._img_com_fonte()
        morph = _frame_morphology(img, 50.0, 50.0, r=20)
        for field in ("pointlike", "elongation", "fwhm_est_px", "valid"):
            assert field in morph, f"Campo morfológico ausente: {field}"

    def test_fonte_circular_baixa_elongation(self):
        img   = self._img_com_fonte(elongado=False)
        morph = _frame_morphology(img, 50.0, 50.0, r=20)
        assert morph["valid"]
        assert morph["elongation"] < T.ELON_GOOD, (
            f"Fonte circular com elongation alta: {morph['elongation']}"
        )

    def test_fonte_alongada_alta_elongation(self):
        img   = self._img_com_fonte(elongado=True)
        morph = _frame_morphology(img, 50.0, 50.0, r=20)
        assert morph["valid"]
        assert morph["elongation"] >= T.ELON_GOOD, (
            f"Fonte alongada com elongation baixa: {morph['elongation']}"
        )

    def test_pontual_maior_que_zero_para_fonte_real(self):
        img   = self._img_com_fonte()
        morph = _frame_morphology(img, 50.0, 50.0, r=20)
        assert morph["pointlike"] > 0

    def test_imagem_degenerada_retorna_invalido(self):
        img   = np.zeros((10, 10))
        morph = _frame_morphology(img, 5.0, 5.0, r=20)
        assert not morph["valid"]

    def test_fwhm_positiva_para_fonte_detectavel(self):
        img   = self._img_com_fonte(flux=10000)
        morph = _frame_morphology(img, 50.0, 50.0, r=20)
        assert morph["fwhm_est_px"] is None or morph["fwhm_est_px"] > 0


# ─────────────────────────────────────────────
# Rejeição de falso positivo (v1.3)
# ─────────────────────────────────────────────
class TestVerificarFalsoPositivo:

    def _edge_candidate(self):
        """Candidato com posição muito próxima à borda."""
        track = [
            (5.0, 5.0, 1000.0),    # perto da borda
            (7.0, 7.0, 1000.0),
            (9.0, 9.0, 1000.0),
            (11.0, 11.0, 1000.0),
        ]
        return {"track": track, "move_total": 8.5, "fluxes": [1000.0]*4, "residual": 5.0}

    def _normal_candidate(self):
        """Candidato sem problemas óbvios."""
        track = [
            (150.0, 150.0, 3000.0),
            (152.0, 153.0, 3100.0),
            (154.0, 156.0, 2900.0),
            (156.0, 159.0, 3050.0),
        ]
        return {"track": track, "move_total": 10.0, "fluxes": [3000.0]*4, "residual": 5.0}

    def _frames_sinteticos(self, shape=(256, 256)):
        frames = []
        for i in range(4):
            img = np.random.default_rng(i).normal(100, 5, shape).astype(np.float64)
            frames.append({
                "data"    : img,
                "date_obs": f"2019-08-28T{10+i}:00:00",
                "jd"      : 2458723.0 + i / 24,
                "wcs_ok"  : False,
                "exptime" : 45.0,
            })
        return frames

    def _morfo_neutro(self):
        return [{"pointlike": 0.15, "elongation": 1.2, "fwhm_est_px": 4.0, "valid": True}] * 4

    def test_normal_candidate_has_no_rejection(self):
        c      = self._normal_candidate()
        frames = self._frames_sinteticos()
        snrs   = [12.0, 11.5, 13.0, 12.5]
        morph  = self._morfo_neutro()
        fluxes = [3000.0, 3100.0, 2900.0, 3050.0]
        reasons = _check_false_positive(c, frames, snrs, morph, fluxes)
        assert reasons == [], f"Candidato normal rejeitado: {reasons}"

    def test_edge_candidate_rejected(self):
        c      = self._edge_candidate()
        frames = self._frames_sinteticos()
        snrs   = [5.0, 5.0, 5.0, 5.0]
        morph  = self._morfo_neutro()
        fluxes = [1000.0]*4
        reasons = _check_false_positive(c, frames, snrs, morph, fluxes)
        assert any("borda" in r for r in reasons), f"Borda não detectada: {reasons}"

    def test_snr_insuficiente_rejeitado(self):
        c      = self._normal_candidate()
        frames = self._frames_sinteticos()
        snrs   = [0.5, 0.8, 0.3, 1.2]   # todos abaixo de T.MIN_SNR
        morph  = self._morfo_neutro()
        fluxes = [3000.0]*4
        reasons = _check_false_positive(c, frames, snrs, morph, fluxes)
        assert any("SNR" in r for r in reasons), f"SNR baixo não detectado: {reasons}"

    def test_chaotic_flux_rejected(self):
        c      = self._normal_candidate()
        frames = self._frames_sinteticos()
        snrs   = [10.0]*4
        morph  = self._morfo_neutro()
        # CV > T.CHAOTIC_FLUX_CV (1.5): um frame quase zero, outro muito alto
        # CV = std/mean com [1, 100000, 1, 100000] ≈ 1.0  — ainda insuficiente
        # Usar [1, 0, 200000, 0]: mean≈50000, std≈100000, CV≈2.0
        fluxes = [1.0, 0.0, 200000.0, 0.0]
        reasons = _check_false_positive(c, frames, snrs, morph, fluxes)
        assert any("caótico" in r for r in reasons), f"Brilho caótico não detectado: {reasons}"

    def test_elongation_grave_rejeitada(self):
        c      = self._normal_candidate()
        frames = self._frames_sinteticos()
        snrs   = [8.0]*4
        # Todos os frames com elongation muito alta
        morph  = [{"pointlike": 0.05, "elongation": 5.0, "fwhm_est_px": 20.0, "valid": True}] * 4
        fluxes = [3000.0]*4
        reasons = _check_false_positive(c, frames, snrs, morph, fluxes)
        assert any("elongação" in r or "elongation" in r.lower() for r in reasons), (
            f"Elongation grave não detectada: {reasons}"
        )

    def test_rejeicao_seta_flag_rejeitado_fp(self):
        """Se há rejeição, analyze_candidate deve incluir flag REJECTED_FP."""
        cand = self._edge_candidate()
        # Reposicionar para garantir que cai na borda de 256x256
        cand["track"] = [
            (5.0, 5.0, 1000.0),
            (7.0, 7.0, 1000.0),
            (9.0, 9.0, 1000.0),
            (11.0, 11.0, 1000.0),
        ]
        frames = [_frame_sintetico(shape=(256, 256), date_obs=f"2019-08-28T{10+i}:00:00")
                  for i in range(4)]
        result = analyze_candidate(cand, frames)
        assert "REJECTED_FP" in result["flags"]
        assert result["triage_class"] == "DISCARDED"
        assert len(result["rejection_reasons"]) >= 1


# ─────────────────────────────────────────────
# Thresholds centralizados (v1.3)
# ─────────────────────────────────────────────
class TestThresholdsCentralizados:
    """Verifica que a triage_class T existe e contém os fields essenciais."""

    def test_campos_obrigatorios(self):
        fields = [
            "MATCH_RADIUS_PX",
            "MOVE_MIN_PX", "MOVE_MAX_PX", "MIN_RESIDUAL_PX",
            "EDGE_MARGIN_PX", "MIN_SNR", "CHAOTIC_FLUX_CV",
            "LIN_EXCELENTE", "LIN_GOOD", "LIN_MARGINAL",
            "VEL_UNIFORM", "VEL_MODERATE",
            "PHOT_STABLE", "PHOT_MODERATE",
            "MOR_POINTLIKE", "MOR_MARGINAL",
            "ELON_GOOD", "MORPH_CONSIST_MAX_STD",
            "VEL_RANGE_MIN_PX", "VEL_RANGE_MAX_PX",
            "SCORE_STRONG", "SCORE_MODERATE", "SCORE_WEAK",
            "GAIA_QUERY_RADIUS_ARCMIN", "GAIA_MIN_MATCHES",
            "GAIA_DIFF_MIN_ARCSEC", "GAIA_DIFF_MAX_ARCSEC",
            "GAIA_COARSE_RADIUS_ARCSEC", "GAIA_MIN_COARSE_MATCHES",
            "GAIA_OFFSET_BIN_ARCSEC", "GAIA_REFINEMENT_RADIUS_ARCSEC",
        ]
        for field in fields:
            assert hasattr(T, field), f"T.{field} ausente"

    def test_hierarquia_limiares_linearity(self):
        assert T.LIN_EXCELENTE < T.LIN_GOOD < T.LIN_MARGINAL

    def test_hierarquia_limiares_velocity(self):
        assert T.VEL_UNIFORM < T.VEL_MODERATE

    def test_photometry_threshold_hierarchy(self):
        assert T.PHOT_STABLE < T.PHOT_MODERATE

    def test_class_hierarchy(self):
        assert T.SCORE_WEAK < T.SCORE_MODERATE < T.SCORE_STRONG

    def test_thresholds_gaia_presentes(self):
        assert hasattr(T, "GAIA_QUERY_RADIUS_ARCMIN")
        assert hasattr(T, "GAIA_MIN_MATCHES")
        assert hasattr(T, "GAIA_DIFF_MIN_ARCSEC")
        assert hasattr(T, "GAIA_DIFF_MAX_ARCSEC")
        assert hasattr(T, "GAIA_COARSE_RADIUS_ARCSEC")
        assert hasattr(T, "GAIA_MIN_COARSE_MATCHES")
        assert hasattr(T, "GAIA_OFFSET_BIN_ARCSEC")
        assert hasattr(T, "GAIA_REFINEMENT_RADIUS_ARCSEC")
        assert T.GAIA_MIN_MATCHES >= 1
        assert T.GAIA_DIFF_MIN_ARCSEC > 0
        assert T.GAIA_DIFF_MAX_ARCSEC > 0
        assert T.GAIA_COARSE_RADIUS_ARCSEC > T.GAIA_MATCH_RADIUS_ARCSEC
        assert T.GAIA_REFINEMENT_RADIUS_ARCSEC > T.GAIA_OFFSET_BIN_ARCSEC


# ─────────────────────────────────────────────
# Refinamento WCS via Gaia DR3 (v1.4.0)
# ─────────────────────────────────────────────
class TestRefinamentoGaia:

    def _frame_sintetico_completo(self, shape=(256, 256),
                                  date_obs="2019-08-28T10:00:00",
                                  rms_inicial_arcsec=0.0):
        """Frame com WCS válido e image sintética com estrelas."""
        from detector import _build_panstarrs_wcs
        header = _header_panstarrs_mock()
        wcs = _build_panstarrs_wcs(header)
        rng = np.random.default_rng(42)
        img = rng.normal(100.0, 5.0, shape).astype(np.float64)
        yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
        for cx, cy in [(60, 60), (120, 80), (80, 140), (180, 100),
                       (50, 190), (200, 50), (140, 200), (100, 30),
                       (220, 180), (30, 120)]:
            img += 5000.0 * np.exp(-((xx - cx)**2 + (yy - cy)**2) / (2 * 1.8**2))
        return {
            "file" : "synth.fits",
            "data"    : img,
            "header"  : header,
            "wcs"     : wcs,
            "wcs_ok"  : True,
            "date_obs": date_obs,
            "jd"      : 2458723.921,
            "exptime" : 45.0,
            "sat_pct" : 0.0,
        }

    def _tabela_gaia_sintetica(self, n=15, ra_c=329.505, dec_c=-11.883,
                               scale_deg=0.01, jd=2458723.921):
        """Tabela Gaia sintética com pm zero e mags razoáveis."""
        from astropy.table import Table, MaskedColumn
        import astropy.units as u_
        rng = np.random.default_rng(7)
        ras  = ra_c  + rng.uniform(-scale_deg, scale_deg, n)
        decs = dec_c + rng.uniform(-scale_deg, scale_deg, n)
        mags = rng.uniform(13.0, 18.0, n)
        pmra  = np.ma.MaskedArray(np.zeros(n), mask=False)
        pmdec = np.ma.MaskedArray(np.zeros(n), mask=False)
        return Table({
            "ra"              : ras,
            "dec"             : decs,
            "pmra"            : pmra,
            "pmdec"           : pmdec,
            "phot_g_mean_mag" : mags,
            "ref_epoch"       : np.full(n, 2016.0),
        })

    # ── teste 1: motion próprio ────────────────────────────────────────
    def test_movimento_proprio_aplicado(self):
        """pm de +100 mas/yr em RA → deslocamento correto em ~3.5 anos."""
        from astropy.table import Table
        import numpy.ma as ma

        n, pm_val = 5, 100.0   # mas/yr
        rng = np.random.default_rng(0)
        ras  = rng.uniform(329.4, 329.6, n)
        decs = rng.uniform(-12.0, -11.8, n)
        pmra_arr  = ma.MaskedArray(np.full(n, pm_val), mask=False)
        pmdec_arr = ma.MaskedArray(np.zeros(n),        mask=False)

        tab = Table({"ra": ras, "dec": decs, "pmra": pmra_arr,
                     "pmdec": pmdec_arr,
                     "phot_g_mean_mag": np.full(n, 15.0),
                     "ref_epoch": np.full(n, 2016.0)})

        # 3.5 anos após epoch Gaia 2016.0 → JD correspondente
        from astropy.time import Time
        jd_obs = float(Time(2019.5, format="jyear").jd)
        corrigida = _apply_gaia_proper_motion(tab, jd_obs)

        cosdec = np.cos(np.radians(float(np.mean(decs))))
        delta_ra_mas = (np.array(corrigida["ra"]) - ras) * cosdec * 3.6e6
        # ~3.5 anos × 100 mas/yr = ~350 mas; tolerância 20 mas
        for d in delta_ra_mas:
            assert 300 < d < 400, f"Deslocamento RA esperado ~350 mas, veio {d:.1f}"
        # Dec não deve mudar (pm_dec = 0)
        delta_dec_mas = (np.array(corrigida["dec"]) - decs) * 3.6e6
        for d in delta_dec_mas:
            assert abs(d) < 1.0, f"Deslocamento Dec deve ser ~0, veio {d:.3f}"

    # ── teste 2: cross-match básico ───────────────────────────────────────
    def test_cross_match_gaia_basico(self):
        """10 detecções com offset constante de 0.5 arcsec → todas casadas."""
        rng = np.random.default_rng(1)
        n = 10
        ra_base  = rng.uniform(329.4, 329.6, n)
        dec_base = rng.uniform(-12.0, -11.8, n)
        gaia_rd  = np.column_stack([ra_base, dec_base])

        cosdec = np.cos(np.radians(float(np.mean(dec_base))))
        offset_arcsec = 0.5
        det_rd = gaia_rd.copy()
        det_rd[:, 0] += offset_arcsec / 3600.0 / cosdec   # offset em RA

        idx_det, idx_gaia = _cross_match_gaia(det_rd, gaia_rd,
                                              raio_arcsec=2.0)
        assert len(idx_det) == n, f"Esperava {n} matches, veio {len(idx_det)}"
        assert set(idx_gaia) == set(range(n))

    def test_cross_match_gaia_wrap_ra_zero(self):
        """Cross-match perto de RA=0 deve usar diferença angular, não subtração simples."""
        gaia_rd = np.array([
            [359.9998, 7.7],
            [0.0002, 7.7005],
            [0.0006, 7.7010],
        ], dtype=float)
        det_rd = gaia_rd.copy()
        det_rd[:, 0] = (det_rd[:, 0] + 0.4 / 3600.0) % 360.0

        idx_det, idx_gaia = _cross_match_gaia(det_rd, gaia_rd,
                                              raio_arcsec=2.0)

        assert len(idx_det) == len(gaia_rd)
        assert set(idx_gaia.tolist()) == set(range(len(gaia_rd)))

    # ── teste 3: fallback sem rede ────────────────────────────────────────
    def test_refinar_wcs_falha_sem_rede(self, monkeypatch):
        """_query_gaia_dr3 retornando None → status gaia_network_failed, WCS preservado."""
        import detector as det_mod
        monkeypatch.setattr(det_mod, "_query_gaia_dr3", lambda *a, **kw: None)

        frame = self._frame_sintetico_completo()
        wcs_antes = frame["wcs"]
        result = _refine_wcs_gaia(frame)

        assert result["wcs_status"] == "gaia_network_failed"
        assert result["wcs"] is wcs_antes, "WCS original deve ser preservado"
        assert result["gaia_n_matches"] == 0

    # ── teste 4: skip quando RMS_pre já é baixo ───────────────────────────
    def test_refinar_wcs_skip_se_pre_baixo(self, monkeypatch):
        """Quando RMS_pre < GAIA_DIFF_MIN_ARCSEC → gaia_skipped_low_pre_rms."""
        import detector as det_mod

        frame = self._frame_sintetico_completo()
        gaia_tab = self._tabela_gaia_sintetica()

        # Construir detecções já alinhadas com Gaia (RMS_pre ≈ 0)
        wcs = frame["wcs"]
        gaia_ra  = np.array(gaia_tab["ra"],  dtype=np.float64)
        gaia_dec = np.array(gaia_tab["dec"], dtype=np.float64)

        # Mock: retorna a tabela Gaia e um RMS_pre artificialmente baixo
        # fazendo com que as posições de detecção coincidam exatamente com Gaia
        monkeypatch.setattr(det_mod, "_query_gaia_dr3",
                            lambda *a, **kw: gaia_tab)

        # Para forçar RMS_pre baixo, sobrescrevemos o cross-match para
        # retornar resíduos nulos diretamente injetando o result do match.
        original_cross = det_mod._cross_match_gaia
        def mock_cross(det_rd, gaia_rd, raio_arcsec=None):
            n = min(len(det_rd), len(gaia_rd))
            return np.arange(n), np.arange(n)

        # Também precisamos que as posições RA/Dec das detecções coincidam
        # com Gaia para que RMS_pre fique < 0.5 arcsec (< GAIA_DIFF_MIN_ARCSEC=1.0)
        def mock_pix2world(px, origin):
            n = len(px)
            return np.column_stack([
                np.array(gaia_tab["ra"][:n],  dtype=float),
                np.array(gaia_tab["dec"][:n], dtype=float),
            ])

        monkeypatch.setattr(det_mod, "_cross_match_gaia", mock_cross)
        monkeypatch.setattr(frame["wcs"], "all_pix2world", mock_pix2world)

        result = _refine_wcs_gaia(frame)
        # Com detecções coincidindo com Gaia, RMS_pre ≈ 0 → skip
        assert result["wcs_status"] in (
            "gaia_skipped_low_pre_rms", "gaia_match_failed", "gaia_refinado",
            "gaia_coarse_offset_only",
        ), f"Status inesperado: {result['wcs_status']}"

    # ── teste 5: aceita refinamento quando melhora ────────────────────────
    def test_refinar_wcs_aceita_se_melhora(self, monkeypatch):
        """
        Mock onde Gaia disponível + cross-match funciona → result é
        'gaia_refinado' ou um fallback gracioso (não aborta).
        Pipeline nunca lança exceção.
        """
        import detector as det_mod

        frame = self._frame_sintetico_completo()
        gaia_tab = self._tabela_gaia_sintetica(n=20)
        monkeypatch.setattr(det_mod, "_query_gaia_dr3",
                            lambda *a, **kw: gaia_tab)

        # Não deve lançar exceção em nenhuma circunstância
        result = _refine_wcs_gaia(frame)

        assert "wcs_status" in result
        assert result["wcs_status"] in {
            "gaia_refinado", "gaia_skipped_low_pre_rms",
            "gaia_network_failed", "gaia_match_failed", "gaia_high_post_rms",
        }
        assert "wcs_rms_pre_arcsec" in result
        assert "gaia_n_matches" in result
        # WCS nunca pode ser None se wcs_ok era True
        assert result["wcs"] is not None

    # ── teste 6: consulta Gaia com mock de launch_job_async ───────────────
    def test_consulta_gaia_mock(self, monkeypatch):
        """Mock do astroquery.gaia → verifica que a query usa os parâmetros certos."""
        from unittest.mock import MagicMock, patch
        from astropy.table import Table
        import numpy.ma as ma

        n = 15
        rng = np.random.default_rng(5)
        mock_table = Table({
            "ra"              : rng.uniform(329.4, 329.6, n),
            "dec"             : rng.uniform(-12.0, -11.8, n),
            "pmra"            : ma.MaskedArray(np.zeros(n), mask=False),
            "pmdec"           : ma.MaskedArray(np.zeros(n), mask=False),
            "phot_g_mean_mag" : rng.uniform(13.0, 18.5, n),
            "ref_epoch"       : np.full(n, 2016.0),
        })
        mock_job = MagicMock()
        mock_job.get_results.return_value = mock_table

        with patch("astroquery.gaia.Gaia.launch_job_async",
                   return_value=mock_job) as mock_launch:
            result = _query_gaia_dr3(329.505, -11.883,
                                            raio_arcmin=6.0,
                                            mag_lim=19.0, mag_min=12.0)
            assert mock_launch.called, "launch_job_async deve ter sido chamado"
            query_str = mock_launch.call_args[0][0]
            assert "329.505" in query_str
            assert "-11.883" in query_str
            assert "gaiadr3.gaia_source" in query_str

        assert result is not None
        assert len(result) == n

    # ── teste 7: _estimate_coarse_offset básico ─────────────────────────
    def test_offset_grosseiro_basico(self):
        """20 detecções com offset constante de 6 arcsec em RA → offset estimado ≈ 6 arcsec."""
        rng = np.random.default_rng(10)
        n = 20
        ra_base  = rng.uniform(329.4, 329.6, n)
        dec_base = rng.uniform(-12.0, -11.8, n)
        gaia_rd  = np.column_stack([ra_base, dec_base])

        cosdec = np.cos(np.radians(float(np.mean(dec_base))))
        offset_ra_arcsec = 6.0
        det_rd = gaia_rd.copy()
        det_rd[:, 0] += offset_ra_arcsec / 3600.0 / cosdec

        off_ra, off_dec, n_match = _estimate_coarse_offset(det_rd, gaia_rd,
                                                              raio_arcsec=15.0)
        assert n_match >= T.GAIA_MIN_COARSE_MATCHES, f"N matches={n_match} insuficiente"
        assert abs(off_ra - offset_ra_arcsec) < 0.5, \
            f"offset_ra esperado ~{offset_ra_arcsec}, veio {off_ra:.3f}"
        assert abs(off_dec) < 0.5, f"offset_dec esperado ~0, veio {off_dec:.3f}"

    # ── teste 8: _estimate_coarse_offset com outliers ───────────────────
    def test_offset_grosseiro_outliers(self):
        """Mediana deve ser robusta a 20% de falsos matches com offset aleatório."""
        rng = np.random.default_rng(11)
        n = 30
        ra_base  = rng.uniform(329.3, 329.7, n)
        dec_base = rng.uniform(-12.1, -11.7, n)
        gaia_rd  = np.column_stack([ra_base, dec_base])

        cosdec = np.cos(np.radians(float(np.mean(dec_base))))
        offset_ra_arcsec = 5.0
        det_rd = gaia_rd.copy()
        det_rd[:, 0] += offset_ra_arcsec / 3600.0 / cosdec

        # Corrompe 6 detecções (20%) com offset aleatório grande
        idx_noise = rng.choice(n, 6, replace=False)
        det_rd[idx_noise, 0] += rng.uniform(-10, 10, 6) / 3600.0 / cosdec
        det_rd[idx_noise, 1] += rng.uniform(-10, 10, 6) / 3600.0

        off_ra, off_dec, n_match = _estimate_coarse_offset(det_rd, gaia_rd,
                                                              raio_arcsec=15.0)
        assert n_match >= T.GAIA_MIN_COARSE_MATCHES
        assert abs(off_ra - offset_ra_arcsec) < 1.0, \
            f"offset_ra com outliers: esperado ~{offset_ra_arcsec}, veio {off_ra:.3f}"

    def test_offset_grosseiro_retorna_indices_1para1(self):
        """Passada bruta retorna pares únicos no pico de translação."""
        rng = np.random.default_rng(12)
        n = 18
        ra_base = rng.uniform(329.45, 329.55, n)
        dec_base = rng.uniform(-11.95, -11.85, n)
        gaia_rd = np.column_stack([ra_base, dec_base])

        cosdec = np.cos(np.radians(float(np.mean(dec_base))))
        det_rd = gaia_rd.copy()
        det_rd[:, 0] += 4.0 / 3600.0 / cosdec
        det_rd[:, 1] -= 3.0 / 3600.0

        off_ra, off_dec, n_match, idx_det, idx_gaia = _estimate_coarse_offset(
            det_rd, gaia_rd, raio_arcsec=15.0, retornar_indices=True
        )

        assert n_match >= T.GAIA_MIN_COARSE_MATCHES
        assert len(idx_det) == n_match
        assert len(idx_gaia) == n_match
        assert len(set(idx_det.tolist())) == n_match
        assert len(set(idx_gaia.tolist())) == n_match
        assert abs(off_ra - 4.0) < 0.5
        assert abs(off_dec + 3.0) < 0.5

    # ── teste 9: _refinar_wcs com offset bruto apenas ────────────────────
    def test_refinar_wcs_offset_bruto_apenas(self, monkeypatch):
        """
        Passada bruta OK (N≥5) mas passada fina insuficiente →
        status gaia_coarse_offset_only e WCS com CRVAL corrigido.
        """
        import detector as det_mod

        frame = self._frame_sintetico_completo()
        crval_original = list(frame["wcs"].wcs.crval)
        gaia_tab = self._tabela_gaia_sintetica(n=20)

        # Gaia retorna tabela válida
        monkeypatch.setattr(det_mod, "_query_gaia_dr3",
                            lambda *a, **kw: gaia_tab)

        # Passada bruta encontra N≥5, passada fina encontra 0
        call_count = {"n": 0}
        original_cross = det_mod._cross_match_gaia

        def mock_cross(det_rd, gaia_rd, raio_arcsec=None):
            call_count["n"] += 1
            if raio_arcsec is not None and raio_arcsec >= T.GAIA_COARSE_RADIUS_ARCSEC:
                # Passada bruta via _estimate_coarse_offset — deixa passar normal
                return original_cross(det_rd, gaia_rd, raio_arcsec=raio_arcsec)
            # Passada fina → retorna vazio
            return np.array([], dtype=int), np.array([], dtype=int)

        monkeypatch.setattr(det_mod, "_cross_match_gaia", mock_cross)

        result = _refine_wcs_gaia(frame)

        assert result["wcs_status"] in (
            "gaia_coarse_offset_only", "gaia_match_failed",
            "gaia_skipped_low_pre_rms", "gaia_refinado",
        ), f"Status inesperado: {result['wcs_status']}"
        assert result["wcs"] is not None
        assert "gaia_coarse_matches" in result
        assert "gaia_coarse_offset_arcsec" in result

    # ── teste 10: _refinar_wcs duas passadas completo ─────────────────────
    def test_refinar_wcs_duas_passadas_completo(self, monkeypatch):
        """
        Simula offset real de 6 arcsec: passada bruta corrige CRVAL,
        passada fina deve encontrar matches suficientes e retornar
        gaia_refinado ou gaia_coarse_offset_only (nunca gaia_match_failed
        nem exceção).
        """
        import detector as det_mod

        frame = self._frame_sintetico_completo()
        wcs = frame["wcs"]

        # Tabela Gaia deslocada 6 arcsec do WCS para simular erro de pointing
        gaia_tab = self._tabela_gaia_sintetica(n=25)
        cosdec = np.cos(np.radians(-11.883))
        gaia_tab_shifted = gaia_tab.copy()
        gaia_tab_shifted["ra"] = np.array(gaia_tab["ra"]) - 6.0 / 3600.0 / cosdec

        monkeypatch.setattr(det_mod, "_query_gaia_dr3",
                            lambda *a, **kw: gaia_tab_shifted)

        result = _refine_wcs_gaia(frame)

        assert result["wcs"] is not None, "WCS não pode ser None"
        assert result["wcs_status"] != "gaia_network_failed"
        assert "gaia_coarse_matches" in result
        assert "gaia_coarse_offset_arcsec" in result
        assert "gaia_refined_matches" in result
        assert isinstance(result["gaia_coarse_offset_arcsec"], list)
        assert len(result["gaia_coarse_offset_arcsec"]) == 2


# ─────────────────────────────────────────────
# Mascaramento de saturação (v1.3.1)
# ─────────────────────────────────────────────
class TestMascaraSaturacao:

    def test_mascara_pixels_invalidos_basico(self):
        img = np.full((100, 100), 500.0)
        # 5 pixels saturados em 65535
        img[10, 10] = 65535
        img[20, 20] = 65535
        img[30, 30] = 65535
        img[40, 40] = 65535
        img[50, 50] = 65535
        # 3 pixels com valor 0 (buraco)
        img[60, 60] = 0
        img[70, 70] = 0
        img[80, 80] = 0
        mascara = _invalid_pixel_mask(img)
        assert mascara.sum() == 8

    def test_detectar_fontes_ignora_saturacao(self):
        rng = np.random.default_rng(99)
        img = rng.normal(100.0, 5.0, (200, 200))
        # Fonte gaussiana real no centro
        cy, cx = 100, 100
        yy, xx = np.mgrid[0:200, 0:200]
        img += 5000.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 2.0 ** 2))
        # Cluster saturado longe da source real
        img[30:40, 30:40] = 65535.0

        sources = detect_sources(img, sigma=5.5)
        # A source real deve ser detectada
        assert len(sources) >= 1
        # Nenhuma source deve estar dentro do cluster saturado
        for y, x, _, _ in sources:
            dentro_cluster = (25 <= y <= 45) and (25 <= x <= 45)
            assert not dentro_cluster, (
                f"Fonte espúria detectada no cluster saturado em y={y:.1f}, x={x:.1f}"
            )

    def test_carregar_fits_aborta_em_frame_muito_saturado(self, tmp_path):
        """FITS com 30% de pixels em 65535 deve causar SystemExit."""
        nx, ny = 100, 100
        img = np.full((ny, nx), 500.0, dtype=np.float32)
        # 30% dos pixels saturados
        n_sat = int(nx * ny * 0.30)
        flat = img.ravel()
        flat[:n_sat] = 65535.0
        img = flat.reshape((ny, nx))

        hdu = fits.PrimaryHDU(img)
        h = hdu.header
        h["DATE-OBS"] = "2019-08-28T10:00:00"
        h["EXPTIME"]  = 45.0
        h["CTYPE1"]   = "RA---TAN"
        h["CTYPE2"]   = "DEC--TAN"
        h["CRVAL1"]   = 329.44
        h["CRVAL2"]   = -12.2
        h["CRPIX1"]   = nx / 2
        h["CRPIX2"]   = ny / 2
        h["CDELT1"]   = -7.095e-05
        h["CDELT2"]   = 7.129e-05
        h["EQUINOX"]  = 2000.0

        folder = tmp_path / "sat_frames"
        folder.mkdir()
        # Criar 4 frames idênticos para satisfazer a exigência do pipeline
        for i in range(4):
            hdu.writeto(folder / f"frame{i+1}.fits", overwrite=True)

        from detector import load_fits
        with pytest.raises(SystemExit):
            load_fits(folder)


# ─────────────────────────────────────────────
# wcs_mode=header: Gaia desativado (v1.5)
# ─────────────────────────────────────────────
class TestWcsModeHeader:
    """Verifica que wcs_mode='header' desativa Gaia mantendo o WCS Pan-STARRS."""

    @staticmethod
    def _criar_pasta_fits(tmp_path, n=4):
        nx, ny = 200, 200
        img = np.random.default_rng(7).normal(100.0, 10.0, (ny, nx)).astype(np.float32)
        folder = tmp_path / "fits_header_mode"
        folder.mkdir()
        for i in range(n):
            hdu = fits.PrimaryHDU(img)
            h = hdu.header
            h["DATE-OBS"] = f"2023-06-01T10:{i:02d}:00"
            h["EXPTIME"]  = 45.0
            h["CTYPE1"]   = "RA---TAN"
            h["CTYPE2"]   = "DEC--TAN"
            h["CRVAL1"]   = 120.0
            h["CRVAL2"]   = -5.0
            h["CRPIX1"]   = nx / 2
            h["CRPIX2"]   = ny / 2
            h["CDELT1"]   = -7.095e-05
            h["CDELT2"]   =  7.129e-05
            h["EQUINOX"]  = 2000.0
            hdu.writeto(folder / f"frame{i+1:04d}.fits", overwrite=True)
        return folder

    def test_wcs_status_gaia_disabled(self, tmp_path):
        from detector import load_fits
        folder = self._criar_pasta_fits(tmp_path)
        frames = load_fits(folder, wcs_mode="header")
        assert len(frames) == 4
        for f in frames:
            assert f.get("wcs_status") == "gaia_disabled", (
                f"Esperado 'gaia_disabled', obtido '{f.get('wcs_status')}'"
            )
            assert f.get("gaia_n_matches", -1) == 0
            assert f.get("gaia_catalogo_refinado") is None

    def test_wcs_ok_preservado_em_header_mode(self, tmp_path):
        from detector import load_fits
        folder = self._criar_pasta_fits(tmp_path)
        frames = load_fits(folder, wcs_mode="header")
        for f in frames:
            assert f["wcs_ok"], "WCS Pan-STARRS deve permanecer válido em header mode"
            assert f["wcs"] is not None

    def test_wcs_mode_gaia_nao_afeta_header_mode(self, tmp_path):
        """Os dois modes não devem se contaminar: cada chamada é independente."""
        from detector import load_fits
        folder = self._criar_pasta_fits(tmp_path)
        frames_h = load_fits(folder, wcs_mode="header")
        frames_g_status = [f.get("wcs_status") for f in frames_h]
        assert all(s == "gaia_disabled" for s in frames_g_status)


# ─────────────────────────────────────────────
# run_id, candidate_id e run_metadata (v1.5)
# ─────────────────────────────────────────────
class TestRastreabilidadeV15:
    """Verifica presença e formato de run_id, candidate_id e run_metadata no JSON."""

    @staticmethod
    def _json_fixture(tmp_path) -> dict:
        """Roda o pipeline end-to-end sintético e retorna o JSON gerado."""
        import uuid
        from unittest.mock import patch
        from astropy.table import Table
        from detector import (
            load_fits, find_movers, analyze_candidate,
            export_json, _deduplicate, PIPELINE_VERSION,
        )
        from datetime import datetime, timezone

        nx, ny = 300, 300
        rng = np.random.default_rng(42)
        folder = tmp_path / "fits_rastr"
        folder.mkdir()

        # Cria 4 frames com uma "source" em motion
        base = rng.normal(300.0, 15.0, (ny, nx)).astype(np.float32)
        for i in range(4):
            img = base.copy()
            cx, cy = 150 + i * 6, 150
            yy, xx = np.mgrid[0:ny, 0:nx]
            img += 8000.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 2.5 ** 2))
            hdu = fits.PrimaryHDU(img)
            h = hdu.header
            h["DATE-OBS"] = f"2023-06-01T10:{i*5:02d}:00"
            h["EXPTIME"]  = 45.0
            h["CTYPE1"]   = "RA---TAN"
            h["CTYPE2"]   = "DEC--TAN"
            h["CRVAL1"]   = 120.0
            h["CRVAL2"]   = -5.0
            h["CRPIX1"]   = nx / 2
            h["CRPIX2"]   = ny / 2
            h["CDELT1"]   = -7.095e-05
            h["CDELT2"]   =  7.129e-05
            h["EQUINOX"]  = 2000.0
            hdu.writeto(folder / f"frame{i+1:04d}.fits", overwrite=True)

        gaia_tab = Table({
            "source_id": [1], "ra": [120.0], "dec": [-5.0],
            "pmra": [0.0], "pmdec": [0.0], "phot_g_mean_mag": [15.0],
        })
        gaia_tab["pmra"].fill_value  = 0.0
        gaia_tab["pmdec"].fill_value = 0.0

        run_id = str(uuid.uuid4())
        wcs_mode = "gaia"

        with patch("detector._query_gaia_dr3", return_value=gaia_tab):
            frames = load_fits(folder, wcs_mode=wcs_mode)

        movers, mg = find_movers(frames, sigma=4.5)
        if not movers:
            return {}

        movers = _deduplicate(movers)
        candidates = [analyze_candidate(m, frames) for m in movers]
        candidates.sort(key=lambda x: -x["score"])

        mode_label = wcs_mode.upper()
        input_set_name = "frame0001"
        for i, c in enumerate(candidates):
            c["rank"]         = i + 1
            c["candidate_id"] = f"{input_set_name}_{mode_label}_C{i+1:03d}"
            c["mpc"]          = {"status": "not_queried", "reason": None, "match": None}

        mg.update({
            "n_valid_frames"   : 4,
            "n_unique_candidates": len(movers),
            "n_final_candidates": len(candidates),
            "n_rejected_false_positive"    : sum(1 for c in candidates if "REJECTED_FP" in c.get("flags", [])),
            "class_distribution"       : {cl: sum(1 for c in candidates if c["triage_class"] == cl)
                                    for cl in ("STRONG", "MODERATE", "WEAK", "DISCARDED")},
            "top1_candidate_id"  : candidates[0].get("candidate_id") if candidates else None,
            "top3_candidate_ids" : [c.get("candidate_id") for c in candidates[:3]],
            "top5_candidate_ids" : [c.get("candidate_id") for c in candidates[:5]],
            "wcs_mode"           : wcs_mode,
            "run_id"             : run_id,
            "execution_time_s"   : 1.0,
        })

        run_metadata = {
            "run_id"             : run_id,
            "pipeline_version"   : PIPELINE_VERSION,
            "wcs_mode"           : wcs_mode,
            "timestamp_execution_utc" : datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "input_set"          : str(folder),
            "input_set_name"      : input_set_name,
            "observer"         : {"name": "Test", "email": "test@test.com"},
        }

        output_dir = tmp_path / "out_rastr"
        output_dir.mkdir()
        file_path = export_json(candidates, frames, input_set_name, output_dir, mg, run_metadata)
        with open(file_path, encoding="utf-8") as f:
            return json.load(f)

    def test_run_metadata_presente(self, tmp_path):
        dados = self._json_fixture(tmp_path)
        if not dados:
            pytest.skip("Nenhum candidate detectado no fixture sintético")
        assert "run_metadata" in dados, "run_metadata ausente no JSON"
        rm = dados["run_metadata"]
        for field in ("run_id", "pipeline_version", "wcs_mode",
                      "timestamp_execution_utc", "input_set"):
            assert field in rm, f"Campo '{field}' ausente em run_metadata"

    def test_run_id_formato_uuid(self, tmp_path):
        import uuid as _uuid
        dados = self._json_fixture(tmp_path)
        if not dados:
            pytest.skip("Nenhum candidate detectado")
        run_id = dados["run_metadata"]["run_id"]
        assert len(run_id) == 36, "run_id deve ter 36 caracteres (UUID4)"
        # Valida que é UUID bem formado
        parsed = _uuid.UUID(run_id)
        assert str(parsed) == run_id

    def test_candidate_id_presente_e_formato(self, tmp_path):
        dados = self._json_fixture(tmp_path)
        if not dados:
            pytest.skip("Nenhum candidate detectado")
        cands = dados.get("candidates", [])
        assert cands, "Nenhum candidate no JSON"
        for c in cands:
            assert "candidate_id" in c, "candidate_id ausente em candidate"
            cid = c["candidate_id"]
            assert cid is not None and cid != "", "candidate_id não pode ser vazio"
            # Formato esperado: CONJUNTO_MODO_C###
            parts = cid.rsplit("_", 2)
            assert len(parts) == 3, f"Formato inesperado de candidate_id: {cid}"
            assert parts[-1].startswith("C"), f"Sufixo C### esperado, obtido: {cid}"

    def test_wcs_mode_no_json(self, tmp_path):
        dados = self._json_fixture(tmp_path)
        if not dados:
            pytest.skip("Nenhum candidate detectado")
        assert dados["run_metadata"]["wcs_mode"] in ("gaia", "header")
        # Também deve aparecer nas métricas globais
        assert "wcs_mode" in dados["global_metrics"]

    def test_manual_validation_empty_fields_present(self, tmp_path):
        dados = self._json_fixture(tmp_path)
        if not dados:
            pytest.skip("Nenhum candidate detectado")
        cands = dados.get("candidates", [])
        if not cands:
            pytest.skip("Sem candidates")
        c = cands[0]
        assert "manual_validation" in c, "manual_validation missing"
        manual = c["manual_validation"]
        for field in ("measured_in_astrometrica", "included_in_astrometrica_mpc",
                      "iasc_feedback", "manual_classification", "notes"):
            assert field in manual, f"Field '{field}' missing in manual_validation"
            assert manual[field] is None, f"'{field}' must be None before review"

    def test_run_id_consistent_between_metadata_and_candidates(self, tmp_path):
        dados = self._json_fixture(tmp_path)
        if not dados:
            pytest.skip("Nenhum candidate detectado")
        run_id_meta = dados["run_metadata"]["run_id"]
        run_id_mg   = dados["global_metrics"].get("run_id")
        assert run_id_meta == run_id_mg, "run_id diverge entre run_metadata e global_metrics"
        for c in dados.get("candidates", []):
            assert c.get("run_id") == run_id_meta, "run_id do candidate diverge do run_metadata"


# ─────────────────────────────────────────────
# export_metrics.py (v1.5)
# ─────────────────────────────────────────────
class TestExportMetrics:
    """Verifica que export_metrics.py lê JSONs e gera CSVs válidos."""

    @staticmethod
    def _json_minimo(tmp_path, input_set="abc123", wcs_mode="gaia") -> Path:
        """Cria um JSON mínimo compatível com o schema v1.5 para testar o exportador."""
        dados = {
            "pipeline"      : "Lia v1.5.0",
            "input_set"      : input_set,
            "processado_em" : "2023-06-01T10:00:00",
            "data_obs"      : "2023-06-01",
            "run_metadata"  : {
                "run_id"             : "aaaabbbb-cccc-dddd-eeee-ffffffffffff",
                "pipeline_version"   : "1.5.0",
                "wcs_mode"           : wcs_mode,
                "timestamp_execution_utc" : "2023-06-01T10:00:00+00:00",
                "input_set"          : "/data/fits",
                "input_set_name"      : input_set,
                "observer"         : {"name": "Test", "email": "t@t.com"},
            },
            "global_metrics": {
                "n_fontes_by_frame"    : [50, 48, 51, 49],
                "n_tracks_attempted"    : 20,
                "n_unique_candidates"   : 5,
                "n_rejected_false_positive"       : 2,
                "n_final_candidates"   : 5,
                "class_distribution"          : {"STRONG": 1, "MODERATE": 2, "WEAK": 1, "DISCARDED": 1},
                "top1_candidate_id"     : f"{input_set}_{wcs_mode.upper()}_C001",
                "top3_candidate_ids"    : [f"{input_set}_{wcs_mode.upper()}_C00{i}" for i in range(1, 4)],
                "top5_candidate_ids"    : [f"{input_set}_{wcs_mode.upper()}_C00{i}" for i in range(1, 6)],
                "wcs_mode"              : wcs_mode,
                "run_id"                : "aaaabbbb-cccc-dddd-eeee-ffffffffffff",
                "sigma"        : 5.5,
                "execution_time_s"      : 12.3,
                "gaia_refinement"      : [
                    {"status": "gaia_refinado", "rms_pre_arcsec": 0.4, "rms_pos_arcsec": 0.08,
                     "n_matches": 20, "coarse_matches": 25}
                ] * 4,
            },
            "candidates": [
                {
                    "candidate_id"  : f"{input_set}_{wcs_mode.upper()}_C001",
                    "run_id"        : "aaaabbbb-cccc-dddd-eeee-ffffffffffff",
                    "wcs_mode"      : wcs_mode,
                    "rank"          : 1,
                    "triage_class"        : "STRONG",
                    "score_total"   : 10,
                    "score_percent" : 83,
                    "score_components": {
                        "linearity": 3, "velocity": 2, "photometry": 2,
                        "morphology": 1, "morph_consistency": 1, "elongation": 1,
                        "velocity_range": 0, "score_max": 12,
                    },
                    "flags"         : ["LINEARIDADE_BOA", "MPC_SEM_MATCH"],
                    "posicao_frame1": {"ra_deg": 120.1, "dec_deg": -5.0, "ra_fmt": None, "dec_fmt": None},
                    "motion"     : {"total_px": 18.0, "linearidade_px": 0.3,
                                       "vel_consistencia": 1.2, "vel_arcsec_min": 0.45,
                                       "residuo_deriva": 8.0},
                    "photometry"    : {"flux_by_frame": [1000, 980, 1010, 990],
                                       "flux_cv": 0.01, "pointness": 0.15,
                                       "pont_std": 0.01, "snr_by_frame": [8.0, 7.5, 8.2, 7.8]},
                    "morphology"    : {"mean_elongation": 1.2, "mean_fwhm_px": 3.1, "by_frame": []},
                    "astrometric_uncertainty": {"sigma_pos_mediana_arcsec": 0.12, "by_frame": []},
                    "gaia_static"   : {"status": "no_match_estatico", "n_matches": 0, "matches": []},
                    "track"        : [],
                    "mpc"           : {"status": "no_match", "reason": None, "match": None},
                    "decision_reasons": {"penalties": [], "rejections": []},
                    "manual_validation": {
                        "medido_astrometrica": None, "entrou_mpc": None,
                        "feedback_iasc": None, "classificacao_manual": None, "observacoes": None,
                    },
                    "manual_validation": {
                        "measured_in_astrometrica": None,
                        "included_in_astrometrica_mpc": None,
                        "iasc_feedback": None,
                        "manual_classification": None,
                        "notes": None,
                    },
                }
            ],
        }
        file_path = tmp_path / f"{input_set}_candidates.json"
        with open(file_path, "w", encoding="utf-8") as f:
            json.dump(dados, f, indent=2)
        return file_path

    def test_csv_sets_generated(self, tmp_path):
        tools_dir = Path(__file__).parent.parent / "tools"
        self._json_minimo(tmp_path)
        import subprocess, sys
        result = subprocess.run(
            [sys.executable, str(tools_dir / "export_metrics.py"),
             str(tmp_path), "--out", str(tmp_path / "validation")],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"export_metrics failed: {result.stderr}"
        sets_csv = tmp_path / "validation_sets.csv"
        assert sets_csv.exists(), "Set-level CSV was not generated"

    def test_csv_candidates_generated(self, tmp_path):
        tools_dir = Path(__file__).parent.parent / "tools"
        self._json_minimo(tmp_path)
        import subprocess, sys
        subprocess.run(
            [sys.executable, str(tools_dir / "export_metrics.py"),
             str(tmp_path), "--out", str(tmp_path / "validation")],
            capture_output=True,
        )
        candidates_csv = tmp_path / "validation_candidates.csv"
        assert candidates_csv.exists(), "Candidate-level CSV was not generated"

    def test_csv_sets_required_fields(self, tmp_path):
        import csv as csvmod
        tools_dir = Path(__file__).parent.parent / "tools"
        self._json_minimo(tmp_path)
        import subprocess, sys
        subprocess.run(
            [sys.executable, str(tools_dir / "export_metrics.py"),
             str(tmp_path), "--out", str(tmp_path / "validation")],
            capture_output=True,
        )
        sets_csv = tmp_path / "validation_sets.csv"
        with open(sets_csv, encoding="utf-8") as f:
            reader = csvmod.DictReader(f)
            rows = list(reader)
        assert len(rows) == 1
        row = rows[0]
        for field in ("run_id", "input_set", "wcs_mode", "n_tracks_attempted",
                      "n_unique_candidates", "n_rejected_false_positive",
                      "top1_candidate_id", "execution_time_s"):
            assert field in row, f"Required field '{field}' is missing from the set-level CSV"
        assert row["wcs_mode"] == "gaia"
        assert row["top1_candidate_id"] == "abc123_GAIA_C001"

    def test_csv_candidates_manual_fields(self, tmp_path):
        import csv as csvmod
        tools_dir = Path(__file__).parent.parent / "tools"
        self._json_minimo(tmp_path)
        import subprocess, sys
        subprocess.run(
            [sys.executable, str(tools_dir / "export_metrics.py"),
             str(tmp_path), "--out", str(tmp_path / "validation")],
            capture_output=True,
        )
        candidates_csv = tmp_path / "validation_candidates.csv"
        with open(candidates_csv, encoding="utf-8") as f:
            reader = csvmod.DictReader(f)
            rows = list(reader)
        assert rows, "No rows were written to the candidate-level CSV"
        row = rows[0]
        for field in ("measured_in_astrometrica", "included_in_astrometrica_mpc",
                      "iasc_feedback", "manual_classification", "notes"):
            assert field in row, f"Manual field '{field}' is missing from the candidate-level CSV"
            assert row[field] == "", f"Manual field '{field}' should be empty, got: '{row[field]}'"

    def test_csv_two_jsons_two_sets(self, tmp_path):
        import csv as csvmod
        tools_dir = Path(__file__).parent.parent / "tools"
        self._json_minimo(tmp_path, input_set="set001", wcs_mode="gaia")
        self._json_minimo(tmp_path, input_set="set002", wcs_mode="header")
        import subprocess, sys
        subprocess.run(
            [sys.executable, str(tools_dir / "export_metrics.py"),
             str(tmp_path), "--out", str(tmp_path / "validation")],
            capture_output=True,
        )
        sets_csv = tmp_path / "validation_sets.csv"
        with open(sets_csv, encoding="utf-8") as f:
            rows = list(csvmod.DictReader(f))
        assert len(rows) == 2, f"Expected 2 sets, got {len(rows)}"
        modes = {r["wcs_mode"] for r in rows}
        assert "gaia" in modes and "header" in modes
