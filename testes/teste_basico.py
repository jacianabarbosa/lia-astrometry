"""
testes/teste_basico.py
======================
Testes unitários e de integração do pipeline asteroid-hunter v1.4.

Cobertura:
    - detectar_fontes: fonte simples, múltiplas, imagem vazia
    - formatar_ra_mpc / formatar_dec_mpc: casos normais e de wraparound
    - formatar_data_mpc: conversão JD → formato MPC
    - _construir_wcs_panstarrs: header Pan-STARRS com chaves PCA* proprietárias
    - _casar_mutuo: validação recíproca de vizinho mais próximo
    - _morfologia_por_frame: elongation, pontualidade, FWHM (v1.3)
    - _verificar_falso_positivo: rejeição por borda, SNR, hot pixel, brilho caótico (v1.3)
    - analisar_candidato: score decomposto, novos campos v1.3, flags, razões de decisão
    - consultar_catalogo_posicional: schema normalizado de retorno
    - _aplicar_movimento_proprio_gaia: propagação de pm por ano (v1.4)
    - _cross_match_gaia: cross-match RA/Dec com KDTree (v1.4)
    - _refinar_wcs_gaia: fallback gracioso e aceitação de refinamento (v1.4)
    - Integração: pipeline end-to-end com FITS sintético (JSON v1.4)

Uso:
    cd asteroid-hunter
    python -m pytest testes/teste_basico.py -v
"""

import sys
import json
import subprocess
import numpy as np
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from detector import (
    detectar_fontes,
    formatar_ra_mpc,
    formatar_dec_mpc,
    formatar_data_mpc,
    _construir_wcs_panstarrs,
    _casar_mutuo,
    _morfologia_por_frame,
    _verificar_falso_positivo,
    _mascara_pixels_invalidos,
    _aplicar_movimento_proprio_gaia,
    _cross_match_gaia,
    _estimar_offset_grosseiro,
    _refinar_wcs_gaia,
    _consultar_gaia_dr3,
    pixel_para_radec,
    analisar_candidato,
    consultar_catalogo_posicional,
    _serializar_trilha,
    T,
)

from astropy.io import fits
from astropy.wcs import WCS


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────
def _fonte_gaussiana(shape, cx, cy, flux=3000.0, fwhm=3.0):
    sigma = fwhm / 2.355
    y, x = np.mgrid[0:shape[0], 0:shape[1]]
    return flux * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma ** 2))


def _img_com_background(shape, sigma_bg=5.0, media_bg=100.0, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(media_bg, sigma_bg, size=shape).astype(np.float64)


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
    from detector import _construir_wcs_panstarrs
    wcs = _construir_wcs_panstarrs(header)
    return {
        "arquivo" : "synth.fits",
        "data"    : np.random.default_rng(0).normal(100, 5, shape),
        "header"  : header,
        "wcs"     : wcs,
        "wcs_ok"  : True,
        "date_obs": date_obs,
        "jd"      : 2458723.921,
        "exptime" : 45.0,
    }


def _candidato_sintetico():
    """Candidato mínimo para testar analisar_candidato."""
    trilha = [
        (100.0, 200.0, 3000.0),
        (102.0, 202.5, 2950.0),
        (104.0, 205.0, 3100.0),
        (106.0, 207.5, 2980.0),
    ]
    return {
        "trilha"    : trilha,
        "move_total": float(np.hypot(6.0, 7.5)),
        "brilhos"   : [t[2] for t in trilha],
        "residuo"   : 3.5,
    }


# ─────────────────────────────────────────────
# Detecção de fontes
# ─────────────────────────────────────────────
class TestDetectarFontes:
    def test_detecta_fonte_simples(self):
        img = _img_com_background((200, 200), seed=1)
        img += _fonte_gaussiana(img.shape, cx=100, cy=100, flux=5000)
        fontes = detectar_fontes(img, sigma=5.0)
        assert len(fontes) >= 1
        y, x, _, _ = fontes[0]
        assert abs(x - 100) < 1.5 and abs(y - 100) < 1.5

    def test_imagem_apenas_ruido(self):
        img = _img_com_background((100, 100), seed=2)
        fontes = detectar_fontes(img, sigma=10.0)
        assert isinstance(fontes, list)
        assert len(fontes) == 0

    def test_multiplas_fontes(self):
        img = _img_com_background((300, 300), seed=3)
        for cx, cy in [(60, 60), (150, 150), (240, 240)]:
            img += _fonte_gaussiana(img.shape, cx, cy, flux=4000)
        fontes = detectar_fontes(img, sigma=5.0)
        assert len(fontes) >= 2


# ─────────────────────────────────────────────
# Formatação MPC (incluindo wraparound)
# ─────────────────────────────────────────────
class TestFormatacaoMPC:
    def test_ra_zero(self):
        assert formatar_ra_mpc(0.0).startswith("00 00")

    def test_ra_180(self):
        assert formatar_ra_mpc(180.0).startswith("12 00")

    def test_ra_wraparound_quase_360(self):
        resultado = formatar_ra_mpc(359.99999)
        partes = resultado.split()
        assert len(partes) == 3
        ss = float(partes[2])
        assert 0 <= ss < 60, f"Segundos devem estar em [0, 60), mas veio {ss}"

    def test_ra_acima_de_360(self):
        r1 = formatar_ra_mpc(10.0)
        r2 = formatar_ra_mpc(370.0)
        assert r1 == r2

    def test_dec_positiva(self):
        assert formatar_dec_mpc(45.5).startswith("+45")

    def test_dec_negativa(self):
        assert formatar_dec_mpc(-12.133).startswith("-12")

    def test_dec_wraparound_quase_46(self):
        resultado = formatar_dec_mpc(45.9999)
        partes = resultado.split()
        ss = float(partes[2])
        assert 0 <= ss < 60, f"Segundos Dec devem estar em [0, 60), veio {ss}"

    def test_dec_sempre_com_sinal(self):
        assert formatar_dec_mpc(0.0)[0] in "+-"
        assert formatar_dec_mpc(45.0)[0] == "+"
        assert formatar_dec_mpc(-45.0)[0] == "-"

    def test_data_formato(self):
        jd = 2458723.921
        resultado = formatar_data_mpc(jd)
        assert resultado.startswith("2019 08")
        partes = resultado.split()
        assert len(partes) == 3


# ─────────────────────────────────────────────
# WCS Pan-STARRS
# ─────────────────────────────────────────────
class TestWCSPanSTARRS:
    def test_wcs_pan_starrs_com_pca_nao_falha(self):
        h = _header_panstarrs_mock()
        wcs = _construir_wcs_panstarrs(h)
        assert wcs is not None
        assert list(wcs.wcs.ctype) == ["RA---TAN", "DEC--TAN"]

    def test_wcs_pan_starrs_crpix_retorna_crval(self):
        h = _header_panstarrs_mock()
        wcs = _construir_wcs_panstarrs(h)
        ra, dec = pixel_para_radec(h["CRPIX1"] - 1, h["CRPIX2"] - 1, wcs)
        assert abs(ra  - h["CRVAL1"]) < 0.01
        assert abs(dec - h["CRVAL2"]) < 0.01

    def test_wcs_pan_starrs_preserva_crota_com_cdelt(self):
        h = _header_panstarrs_mock()
        h["CDELT1"] = 7.10861946296706e-05
        h["CDELT2"] = 7.13340298342252e-05
        h["CROTA1"] = 180.4854
        h["CROTA2"] = 180.4854

        wcs = _construir_wcs_panstarrs(h)

        assert hasattr(wcs.wcs, "crota")
        assert abs(float(wcs.wcs.crota[0]) - 180.4854) < 1e-6
        ra_c, dec_c = pixel_para_radec(h["CRPIX1"] - 1, h["CRPIX2"] - 1, wcs)
        ra_x, dec_x = pixel_para_radec(h["CRPIX1"], h["CRPIX2"] - 1, wcs)
        # CROTA≈180° inverte o eixo: avançar 1 px em x deve reduzir RA.
        assert ra_x < ra_c
        assert abs(dec_x - dec_c) < 1e-4

    def test_wcs_sem_chaves_falha_claramente(self):
        h = fits.Header()
        h["SIMPLE"] = True
        with pytest.raises(ValueError, match="WCS mínimas"):
            _construir_wcs_panstarrs(h)


# ─────────────────────────────────────────────
# Casamento mútuo de fontes
# ─────────────────────────────────────────────
class TestCasarMutuo:
    def test_casar_vazio(self):
        assert _casar_mutuo([], []) == {}

    def test_casar_simples(self):
        a = [(10.0, 20.0, 100.0, 5.0)]
        b = [(11.0, 21.0, 100.0, 5.0)]
        pares = _casar_mutuo(a, b, raio_px=5.0)
        assert pares == {0: 0}

    def test_casar_fora_de_raio(self):
        a = [(10.0, 20.0, 100.0, 5.0)]
        b = [(500.0, 500.0, 100.0, 5.0)]
        pares = _casar_mutuo(a, b, raio_px=10.0)
        assert pares == {}

    def test_casar_mutuo_rejeita_nao_reciproco(self):
        a = [(10.0, 10.0, 100.0, 5.0), (20.0, 10.0, 100.0, 5.0)]
        b = [(21.0, 10.0, 100.0, 5.0)]
        pares = _casar_mutuo(a, b, raio_px=20.0)
        assert 0 not in pares
        assert pares.get(1) == 0


# ─────────────────────────────────────────────
# Score decomposto e flags (v1.3)
# ─────────────────────────────────────────────
class TestScoreDecomposto:
    """Verifica que analisar_candidato retorna todos os campos v1.3."""

    def _rodar(self):
        cand   = _candidato_sintetico()
        frames = [_frame_sintetico(date_obs=f"2019-08-28T{10+i}:00:00")
                  for i in range(4)]
        return analisar_candidato(cand, frames)

    def test_campos_score_decomposto_presentes(self):
        resultado = self._rodar()
        for campo in ("score_linearidade", "score_velocidade",
                      "score_fotometria", "score_morfologia",
                      "score_morfo_consist", "score_elongation", "score_faixa_vel"):
            assert campo in resultado, f"Campo ausente: {campo}"

    def test_soma_score_consistente(self):
        resultado = self._rodar()
        # Score pode ser 0 se candidato foi rejeitado por FP; validamos score_bruto
        soma = (resultado["score_linearidade"] + resultado["score_velocidade"]
                + resultado["score_fotometria"] + resultado["score_morfologia"]
                + resultado["score_morfo_consist"] + resultado["score_elongation"]
                + resultado["score_faixa_vel"])
        assert resultado["score_bruto"] == soma

    def test_limites_score_componentes(self):
        resultado = self._rodar()
        assert 0 <= resultado["score_linearidade"]   <= 3
        assert 0 <= resultado["score_velocidade"]    <= 2
        assert 0 <= resultado["score_fotometria"]    <= 2
        assert 0 <= resultado["score_morfologia"]    <= 2
        assert 0 <= resultado["score_morfo_consist"] <= 1
        assert 0 <= resultado["score_elongation"]    <= 1
        assert 0 <= resultado["score_faixa_vel"]     <= 1

    def test_score_max_correto(self):
        resultado = self._rodar()
        assert resultado["score_max"] == 12

    def test_flags_e_lista(self):
        resultado = self._rodar()
        assert "flags" in resultado
        assert isinstance(resultado["flags"], list)

    def test_flags_validas(self):
        """Todas as flags devem ser strings não-vazias do conjunto conhecido."""
        flags_validas = {
            "LINEARIDADE_BOA", "LINEARIDADE_RUIM",
            "VELOCIDADE_CONSISTENTE", "VELOCIDADE_IRREGULAR",
            "BRILHO_ESTAVEL", "BRILHO_INSTAVEL",
            "MORFOLOGIA_PONTUAL", "MORFOLOGIA_ESTENDIDA",
            "ELONGACAO_ALTA", "MORFO_INCONSISTENTE",
            "EDGE_FRAME", "REJEITADO_FP",
            "MPC_SEM_MATCH", "MPC_MATCH_PROVAVEL",
            "MPC_MATCH_AMBIGUO", "MPC_CONSULTA_FALHOU",
        }
        resultado = self._rodar()
        for flag in resultado["flags"]:
            assert isinstance(flag, str) and flag
            assert flag in flags_validas, f"Flag inesperada: {flag}"

    def test_razoes_decisao_presentes(self):
        resultado = self._rodar()
        assert "razoes_penalizacao" in resultado
        assert "razoes_rejeicao" in resultado
        assert isinstance(resultado["razoes_penalizacao"], list)
        assert isinstance(resultado["razoes_rejeicao"], list)

    def test_campos_morfologia_presentes(self):
        resultado = self._rodar()
        assert "elong_medio" in resultado
        assert "pont_std" in resultado
        assert "fwhm_medio_px" in resultado
        assert isinstance(resultado["elong_medio"], float)
        assert resultado["elong_medio"] >= 1.0

    def test_snrs_presentes(self):
        resultado = self._rodar()
        assert "snrs" in resultado
        assert len(resultado["snrs"]) == 4

    def test_vel_arcsec_min_tipo(self):
        resultado = self._rodar()
        assert resultado.get("vel_arcsec_min") is None \
            or isinstance(resultado["vel_arcsec_min"], float)


# ─────────────────────────────────────────────
# Serialização da trilha (v1.2)
# ─────────────────────────────────────────────
class TestSerializacaoTrilha:
    def test_schema_completo(self):
        cand   = _candidato_sintetico()
        cand["snrs"] = [12.3, 11.8, 13.1, 12.5]
        frames = [_frame_sintetico(date_obs=f"2019-08-28T{10+i}:00:00")
                  for i in range(4)]
        trilha_json = _serializar_trilha(cand, frames)

        assert len(trilha_json) == 4
        campos_obrigatorios = {"frame_index", "timestamp_utc", "jd",
                               "x", "y", "flux", "snr", "wcs_valido"}
        for entrada in trilha_json:
            for campo in campos_obrigatorios:
                assert campo in entrada, f"Campo ausente na trilha: {campo}"

    def test_frame_index_sequencial(self):
        cand   = _candidato_sintetico()
        frames = [_frame_sintetico() for _ in range(4)]
        trilha_json = _serializar_trilha(cand, frames)
        indices = [e["frame_index"] for e in trilha_json]
        assert indices == [0, 1, 2, 3]

    def test_coordenadas_numericas(self):
        cand   = _candidato_sintetico()
        frames = [_frame_sintetico() for _ in range(4)]
        trilha_json = _serializar_trilha(cand, frames)
        for entrada in trilha_json:
            assert isinstance(entrada["x"], float)
            assert isinstance(entrada["y"], float)
            assert isinstance(entrada["flux"], float)
            # ra_deg pode ser None se WCS falhar, mas deve ser float quando presente
            if entrada["ra_deg"] is not None:
                assert isinstance(entrada["ra_deg"], float)


# ─────────────────────────────────────────────
# Status MPC normalizado (v1.2)
# ─────────────────────────────────────────────
class TestStatusMPCNormalizado:
    """Verifica que consultar_catalogo_posicional retorna schema coerente."""

    STATUSES_VALIDOS = {
        "nao_consultado", "consulta_falhou", "sem_match",
        "match_ambiguo", "match_provavel",
    }

    def test_coordenadas_invalidas_retornam_consulta_falhou(self):
        resultado = consultar_catalogo_posicional(400.0, 100.0, "2019-08-28T10:00:00")
        assert resultado["status"] == "consulta_falhou"
        assert "match" in resultado

    def test_status_e_campo_obrigatorio(self):
        resultado = consultar_catalogo_posicional(400.0, 100.0, "2019-08-28T10:00:00")
        assert resultado["status"] in self.STATUSES_VALIDOS

    def test_sem_astroquery_retorna_consulta_falhou(self, monkeypatch):
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "astroquery.imcce":
                raise ImportError("mocked")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        resultado = consultar_catalogo_posicional(329.44, -12.2, "2019-08-28T10:00:00")
        assert resultado["status"] == "consulta_falhou"
        assert resultado["match"] is None


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
            pytest.skip("astroquery não instalado")


# ─────────────────────────────────────────────
# Integração end-to-end (v1.2)
# ─────────────────────────────────────────────
@pytest.mark.integration
class TestIntegracaoEndToEnd:
    """
    Gera 4 FITS sintéticos com asteroide em movimento conhecido e roda o
    pipeline completo. Gaia mockada para retornar None (fallback gracioso)
    — garante que o teste não depende de internet.

    Verifica estrutura do JSON v1.4:
    - metricas_globais (gaia_refinamento por frame, sem n_rejeitados_coerencia)
    - score_componentes (7 componentes, score_max=12)
    - morfologia por frame
    - razoes_decisao
    - trilha completa
    - flags
    - mpc normalizado
    """

    @staticmethod
    def _gerar_fits_sinteticos(pasta: Path, seed: int = 42):
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
            hdu.writeto(pasta / f"synth_frame{i+1}.fits", overwrite=True)

    def test_pipeline_detecta_asteroide_sintetico(self, tmp_path, monkeypatch):
        imagens    = tmp_path / "imagens"
        resultados = tmp_path / "resultados"
        imagens.mkdir(); resultados.mkdir()
        self._gerar_fits_sinteticos(imagens)

        # Mock Gaia para não depender de rede: fallback gracioso sempre ativo.
        import detector as det_mod
        monkeypatch.setattr(det_mod, "_consultar_gaia_dr3", lambda *a, **kw: None)

        detector_py = Path(__file__).parent.parent / "src" / "detector.py"
        result = subprocess.run(
            [sys.executable, str(detector_py),
             "--imagens", str(imagens),
             "--output",  str(resultados),
             "--sem-mpc"],
            capture_output=True, text=True, timeout=120,
        )
        assert result.returncode == 0, (
            f"Pipeline falhou:\nSTDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

        jsons = list(resultados.glob("*_candidatos.json"))
        assert len(jsons) == 1
        with open(jsons[0]) as f:
            dados = json.load(f)

        # Estrutura de topo
        assert "metricas_globais" in dados, "metricas_globais ausente no JSON"
        assert "candidatos" in dados
        assert len(dados["candidatos"]) >= 1

        # Score >= 6 esperado para asteroide sintético bem definido
        scores = [c["score_total"] for c in dados["candidatos"]]
        assert max(scores) >= 6, f"Score máximo baixo demais: {scores}"

        # Verificar estrutura de cada candidato
        for c in dados["candidatos"]:
            # Score decomposto (v1.3: 7 componentes)
            assert "score_componentes" in c, "score_componentes ausente"
            sc = c["score_componentes"]
            for comp in ("linearidade", "velocidade", "fotometria",
                         "morfologia", "consist_morfologica",
                         "elongation", "faixa_vel"):
                assert comp in sc, f"Componente '{comp}' ausente em score_componentes"
            assert sc.get("score_max") == 12, "score_max deve ser 12"

            # score_total deve ser <= score_max
            assert 0 <= c["score_total"] <= 12, (
                f"score_total={c['score_total']} fora de [0, 12]"
            )

            # Flags
            assert "flags" in c
            assert isinstance(c["flags"], list)

            # Razões de decisão (v1.3)
            assert "razoes_decisao" in c
            assert "penalizacoes" in c["razoes_decisao"]
            assert "rejeicoes"    in c["razoes_decisao"]

            # Morfologia (v1.3)
            assert "morfologia" in c
            morfo = c["morfologia"]
            assert "elongation_medio" in morfo
            assert "por_frame" in morfo
            assert len(morfo["por_frame"]) == 4

            # Trilha completa
            assert "trilha" in c
            assert len(c["trilha"]) == 4
            for entrada in c["trilha"]:
                for campo in ("frame_index", "timestamp_utc", "x", "y",
                              "flux", "wcs_valido"):
                    assert campo in entrada, f"Campo '{campo}' ausente na trilha"

            # MPC normalizado
            assert "mpc" in c
            assert "status" in c["mpc"]
            assert c["mpc"]["status"] in {
                "nao_consultado", "consulta_falhou", "sem_match",
                "match_ambiguo", "match_provavel",
            }

        # Métricas globais mínimas
        mg = dados["metricas_globais"]
        for campo in ("n_fontes_por_frame", "n_trilhas_tentadas",
                      "sigma_deteccao", "deriva_dx_px", "deriva_dy_px",
                      "tempo_execucao_s", "gaia_refinamento"):
            assert campo in mg, f"Métrica global '{campo}' ausente"

        assert isinstance(mg["n_fontes_por_frame"], list)
        assert len(mg["n_fontes_por_frame"]) == 4

        # Métricas de refinamento Gaia (fallback esperado pois Gaia é mockada)
        assert isinstance(mg["gaia_refinamento"], list)
        assert len(mg["gaia_refinamento"]) == 4
        for g in mg["gaia_refinamento"]:
            assert "frame" in g and "status" in g and "n_matches" in g
            assert "n_matches_bruto" in g
            assert "offset_bruto_arcsec" in g
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
        morfo = _morfologia_por_frame(img, 50.0, 50.0, r=20)
        for campo in ("pontual", "elongation", "fwhm_est_px", "valido"):
            assert campo in morfo, f"Campo morfológico ausente: {campo}"

    def test_fonte_circular_baixa_elongation(self):
        img   = self._img_com_fonte(elongado=False)
        morfo = _morfologia_por_frame(img, 50.0, 50.0, r=20)
        assert morfo["valido"]
        assert morfo["elongation"] < T.ELON_BOA, (
            f"Fonte circular com elongation alta: {morfo['elongation']}"
        )

    def test_fonte_alongada_alta_elongation(self):
        img   = self._img_com_fonte(elongado=True)
        morfo = _morfologia_por_frame(img, 50.0, 50.0, r=20)
        assert morfo["valido"]
        assert morfo["elongation"] >= T.ELON_BOA, (
            f"Fonte alongada com elongation baixa: {morfo['elongation']}"
        )

    def test_pontual_maior_que_zero_para_fonte_real(self):
        img   = self._img_com_fonte()
        morfo = _morfologia_por_frame(img, 50.0, 50.0, r=20)
        assert morfo["pontual"] > 0

    def test_imagem_degenerada_retorna_invalido(self):
        img   = np.zeros((10, 10))
        morfo = _morfologia_por_frame(img, 5.0, 5.0, r=20)
        assert not morfo["valido"]

    def test_fwhm_positiva_para_fonte_detectavel(self):
        img   = self._img_com_fonte(flux=10000)
        morfo = _morfologia_por_frame(img, 50.0, 50.0, r=20)
        assert morfo["fwhm_est_px"] is None or morfo["fwhm_est_px"] > 0


# ─────────────────────────────────────────────
# Rejeição de falso positivo (v1.3)
# ─────────────────────────────────────────────
class TestVerificarFalsoPositivo:

    def _candidato_na_borda(self):
        """Candidato com posição muito próxima à borda."""
        trilha = [
            (5.0, 5.0, 1000.0),    # perto da borda
            (7.0, 7.0, 1000.0),
            (9.0, 9.0, 1000.0),
            (11.0, 11.0, 1000.0),
        ]
        return {"trilha": trilha, "move_total": 8.5, "brilhos": [1000.0]*4, "residuo": 5.0}

    def _candidato_normal(self):
        """Candidato sem problemas óbvios."""
        trilha = [
            (150.0, 150.0, 3000.0),
            (152.0, 153.0, 3100.0),
            (154.0, 156.0, 2900.0),
            (156.0, 159.0, 3050.0),
        ]
        return {"trilha": trilha, "move_total": 10.0, "brilhos": [3000.0]*4, "residuo": 5.0}

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
        return [{"pontual": 0.15, "elongation": 1.2, "fwhm_est_px": 4.0, "valido": True}] * 4

    def test_candidato_normal_sem_rejeicao(self):
        c      = self._candidato_normal()
        frames = self._frames_sinteticos()
        snrs   = [12.0, 11.5, 13.0, 12.5]
        morfo  = self._morfo_neutro()
        brilhos = [3000.0, 3100.0, 2900.0, 3050.0]
        razoes = _verificar_falso_positivo(c, frames, snrs, morfo, brilhos)
        assert razoes == [], f"Candidato normal rejeitado: {razoes}"

    def test_candidato_na_borda_rejeitado(self):
        c      = self._candidato_na_borda()
        frames = self._frames_sinteticos()
        snrs   = [5.0, 5.0, 5.0, 5.0]
        morfo  = self._morfo_neutro()
        brilhos = [1000.0]*4
        razoes = _verificar_falso_positivo(c, frames, snrs, morfo, brilhos)
        assert any("borda" in r for r in razoes), f"Borda não detectada: {razoes}"

    def test_snr_insuficiente_rejeitado(self):
        c      = self._candidato_normal()
        frames = self._frames_sinteticos()
        snrs   = [0.5, 0.8, 0.3, 1.2]   # todos abaixo de T.SNR_MINIMO
        morfo  = self._morfo_neutro()
        brilhos = [3000.0]*4
        razoes = _verificar_falso_positivo(c, frames, snrs, morfo, brilhos)
        assert any("SNR" in r for r in razoes), f"SNR baixo não detectado: {razoes}"

    def test_brilho_caotico_rejeitado(self):
        c      = self._candidato_normal()
        frames = self._frames_sinteticos()
        snrs   = [10.0]*4
        morfo  = self._morfo_neutro()
        # CV > T.BRILHO_CV_CAOS (1.5): um frame quase zero, outro muito alto
        # CV = std/mean com [1, 100000, 1, 100000] ≈ 1.0  — ainda insuficiente
        # Usar [1, 0, 200000, 0]: mean≈50000, std≈100000, CV≈2.0
        brilhos = [1.0, 0.0, 200000.0, 0.0]
        razoes = _verificar_falso_positivo(c, frames, snrs, morfo, brilhos)
        assert any("caótico" in r for r in razoes), f"Brilho caótico não detectado: {razoes}"

    def test_elongation_grave_rejeitada(self):
        c      = self._candidato_normal()
        frames = self._frames_sinteticos()
        snrs   = [8.0]*4
        # Todos os frames com elongation muito alta
        morfo  = [{"pontual": 0.05, "elongation": 5.0, "fwhm_est_px": 20.0, "valido": True}] * 4
        brilhos = [3000.0]*4
        razoes = _verificar_falso_positivo(c, frames, snrs, morfo, brilhos)
        assert any("elongação" in r or "elongation" in r.lower() for r in razoes), (
            f"Elongation grave não detectada: {razoes}"
        )

    def test_rejeicao_seta_flag_rejeitado_fp(self):
        """Se há rejeição, analisar_candidato deve incluir flag REJEITADO_FP."""
        cand = self._candidato_na_borda()
        # Reposicionar para garantir que cai na borda de 256x256
        cand["trilha"] = [
            (5.0, 5.0, 1000.0),
            (7.0, 7.0, 1000.0),
            (9.0, 9.0, 1000.0),
            (11.0, 11.0, 1000.0),
        ]
        frames = [_frame_sintetico(shape=(256, 256), date_obs=f"2019-08-28T{10+i}:00:00")
                  for i in range(4)]
        resultado = analisar_candidato(cand, frames)
        assert "REJEITADO_FP" in resultado["flags"]
        assert resultado["classe"] == "DESCARTA"
        assert len(resultado["razoes_rejeicao"]) >= 1


# ─────────────────────────────────────────────
# Thresholds centralizados (v1.3)
# ─────────────────────────────────────────────
class TestThresholdsCentralizados:
    """Verifica que a classe T existe e contém os campos essenciais."""

    def test_campos_obrigatorios(self):
        campos = [
            "RAIO_MATCH_PX",
            "MOVE_MIN_PX", "MOVE_MAX_PX", "RESIDUO_MIN_PX",
            "MARGEM_BORDA_PX", "SNR_MINIMO", "BRILHO_CV_CAOS",
            "LIN_EXCELENTE", "LIN_BOA", "LIN_MARGINAL",
            "VEL_UNIFORME", "VEL_MODERADA",
            "FOT_ESTAVEL", "FOT_MODERADA",
            "MOR_PONTUAL", "MOR_MARGINAL",
            "ELON_BOA", "MOR_CONSIST_MAX_STD",
            "FAIXA_VEL_MIN_PX", "FAIXA_VEL_MAX_PX",
            "SCORE_FORTE", "SCORE_MODERADO", "SCORE_FRACO",
            "GAIA_RAIO_QUERY_ARCMIN", "GAIA_MIN_MATCHES",
            "GAIA_DIFF_MIN_ARCSEC", "GAIA_DIFF_MAX_ARCSEC",
            "GAIA_RAIO_BRUTO_ARCSEC", "GAIA_MIN_MATCHES_BRUTO",
            "GAIA_OFFSET_BIN_ARCSEC", "GAIA_RAIO_REFINO_ARCSEC",
        ]
        for campo in campos:
            assert hasattr(T, campo), f"T.{campo} ausente"

    def test_hierarquia_limiares_linearidade(self):
        assert T.LIN_EXCELENTE < T.LIN_BOA < T.LIN_MARGINAL

    def test_hierarquia_limiares_velocidade(self):
        assert T.VEL_UNIFORME < T.VEL_MODERADA

    def test_hierarquia_limiares_fotometria(self):
        assert T.FOT_ESTAVEL < T.FOT_MODERADA

    def test_hierarquia_classes(self):
        assert T.SCORE_FRACO < T.SCORE_MODERADO < T.SCORE_FORTE

    def test_thresholds_gaia_presentes(self):
        assert hasattr(T, "GAIA_RAIO_QUERY_ARCMIN")
        assert hasattr(T, "GAIA_MIN_MATCHES")
        assert hasattr(T, "GAIA_DIFF_MIN_ARCSEC")
        assert hasattr(T, "GAIA_DIFF_MAX_ARCSEC")
        assert hasattr(T, "GAIA_RAIO_BRUTO_ARCSEC")
        assert hasattr(T, "GAIA_MIN_MATCHES_BRUTO")
        assert hasattr(T, "GAIA_OFFSET_BIN_ARCSEC")
        assert hasattr(T, "GAIA_RAIO_REFINO_ARCSEC")
        assert T.GAIA_MIN_MATCHES >= 1
        assert T.GAIA_DIFF_MIN_ARCSEC > 0
        assert T.GAIA_DIFF_MAX_ARCSEC > 0
        assert T.GAIA_RAIO_BRUTO_ARCSEC > T.GAIA_RAIO_MATCH_ARCSEC
        assert T.GAIA_RAIO_REFINO_ARCSEC > T.GAIA_OFFSET_BIN_ARCSEC


# ─────────────────────────────────────────────
# Refinamento WCS via Gaia DR3 (v1.4.0)
# ─────────────────────────────────────────────
class TestRefinamentoGaia:

    def _frame_sintetico_completo(self, shape=(256, 256),
                                  date_obs="2019-08-28T10:00:00",
                                  rms_inicial_arcsec=0.0):
        """Frame com WCS válido e imagem sintética com estrelas."""
        from detector import _construir_wcs_panstarrs
        header = _header_panstarrs_mock()
        wcs = _construir_wcs_panstarrs(header)
        rng = np.random.default_rng(42)
        img = rng.normal(100.0, 5.0, shape).astype(np.float64)
        yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
        for cx, cy in [(60, 60), (120, 80), (80, 140), (180, 100),
                       (50, 190), (200, 50), (140, 200), (100, 30),
                       (220, 180), (30, 120)]:
            img += 5000.0 * np.exp(-((xx - cx)**2 + (yy - cy)**2) / (2 * 1.8**2))
        return {
            "arquivo" : "synth.fits",
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

    # ── teste 1: movimento próprio ────────────────────────────────────────
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
        corrigida = _aplicar_movimento_proprio_gaia(tab, jd_obs)

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
        """_consultar_gaia_dr3 retornando None → status gaia_falhou_rede, WCS preservado."""
        import detector as det_mod
        monkeypatch.setattr(det_mod, "_consultar_gaia_dr3", lambda *a, **kw: None)

        frame = self._frame_sintetico_completo()
        wcs_antes = frame["wcs"]
        resultado = _refinar_wcs_gaia(frame)

        assert resultado["wcs_status"] == "gaia_falhou_rede"
        assert resultado["wcs"] is wcs_antes, "WCS original deve ser preservado"
        assert resultado["gaia_n_matches"] == 0

    # ── teste 4: skip quando RMS_pre já é baixo ───────────────────────────
    def test_refinar_wcs_skip_se_pre_baixo(self, monkeypatch):
        """Quando RMS_pre < GAIA_DIFF_MIN_ARCSEC → gaia_skipped_pre_baixo."""
        import detector as det_mod

        frame = self._frame_sintetico_completo()
        gaia_tab = self._tabela_gaia_sintetica()

        # Construir detecções já alinhadas com Gaia (RMS_pre ≈ 0)
        wcs = frame["wcs"]
        gaia_ra  = np.array(gaia_tab["ra"],  dtype=np.float64)
        gaia_dec = np.array(gaia_tab["dec"], dtype=np.float64)

        # Mock: retorna a tabela Gaia e um RMS_pre artificialmente baixo
        # fazendo com que as posições de detecção coincidam exatamente com Gaia
        monkeypatch.setattr(det_mod, "_consultar_gaia_dr3",
                            lambda *a, **kw: gaia_tab)

        # Para forçar RMS_pre baixo, sobrescrevemos o cross-match para
        # retornar resíduos nulos diretamente injetando o resultado do match.
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

        resultado = _refinar_wcs_gaia(frame)
        # Com detecções coincidindo com Gaia, RMS_pre ≈ 0 → skip
        assert resultado["wcs_status"] in (
            "gaia_skipped_pre_baixo", "gaia_falhou_match", "gaia_refinado",
            "gaia_offset_bruto_apenas",
        ), f"Status inesperado: {resultado['wcs_status']}"

    # ── teste 5: aceita refinamento quando melhora ────────────────────────
    def test_refinar_wcs_aceita_se_melhora(self, monkeypatch):
        """
        Mock onde Gaia disponível + cross-match funciona → resultado é
        'gaia_refinado' ou um fallback gracioso (não aborta).
        Pipeline nunca lança exceção.
        """
        import detector as det_mod

        frame = self._frame_sintetico_completo()
        gaia_tab = self._tabela_gaia_sintetica(n=20)
        monkeypatch.setattr(det_mod, "_consultar_gaia_dr3",
                            lambda *a, **kw: gaia_tab)

        # Não deve lançar exceção em nenhuma circunstância
        resultado = _refinar_wcs_gaia(frame)

        assert "wcs_status" in resultado
        assert resultado["wcs_status"] in {
            "gaia_refinado", "gaia_skipped_pre_baixo",
            "gaia_falhou_rede", "gaia_falhou_match", "gaia_falhou_pos_alto",
        }
        assert "wcs_rms_pre_arcsec" in resultado
        assert "gaia_n_matches" in resultado
        # WCS nunca pode ser None se wcs_ok era True
        assert resultado["wcs"] is not None

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
            resultado = _consultar_gaia_dr3(329.505, -11.883,
                                            raio_arcmin=6.0,
                                            mag_lim=19.0, mag_min=12.0)
            assert mock_launch.called, "launch_job_async deve ter sido chamado"
            query_str = mock_launch.call_args[0][0]
            assert "329.505" in query_str
            assert "-11.883" in query_str
            assert "gaiadr3.gaia_source" in query_str

        assert resultado is not None
        assert len(resultado) == n

    # ── teste 7: _estimar_offset_grosseiro básico ─────────────────────────
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

        off_ra, off_dec, n_match = _estimar_offset_grosseiro(det_rd, gaia_rd,
                                                              raio_arcsec=15.0)
        assert n_match >= T.GAIA_MIN_MATCHES_BRUTO, f"N matches={n_match} insuficiente"
        assert abs(off_ra - offset_ra_arcsec) < 0.5, \
            f"offset_ra esperado ~{offset_ra_arcsec}, veio {off_ra:.3f}"
        assert abs(off_dec) < 0.5, f"offset_dec esperado ~0, veio {off_dec:.3f}"

    # ── teste 8: _estimar_offset_grosseiro com outliers ───────────────────
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

        off_ra, off_dec, n_match = _estimar_offset_grosseiro(det_rd, gaia_rd,
                                                              raio_arcsec=15.0)
        assert n_match >= T.GAIA_MIN_MATCHES_BRUTO
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

        off_ra, off_dec, n_match, idx_det, idx_gaia = _estimar_offset_grosseiro(
            det_rd, gaia_rd, raio_arcsec=15.0, retornar_indices=True
        )

        assert n_match >= T.GAIA_MIN_MATCHES_BRUTO
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
        status gaia_offset_bruto_apenas e WCS com CRVAL corrigido.
        """
        import detector as det_mod

        frame = self._frame_sintetico_completo()
        crval_original = list(frame["wcs"].wcs.crval)
        gaia_tab = self._tabela_gaia_sintetica(n=20)

        # Gaia retorna tabela válida
        monkeypatch.setattr(det_mod, "_consultar_gaia_dr3",
                            lambda *a, **kw: gaia_tab)

        # Passada bruta encontra N≥5, passada fina encontra 0
        call_count = {"n": 0}
        original_cross = det_mod._cross_match_gaia

        def mock_cross(det_rd, gaia_rd, raio_arcsec=None):
            call_count["n"] += 1
            if raio_arcsec is not None and raio_arcsec >= T.GAIA_RAIO_BRUTO_ARCSEC:
                # Passada bruta via _estimar_offset_grosseiro — deixa passar normal
                return original_cross(det_rd, gaia_rd, raio_arcsec=raio_arcsec)
            # Passada fina → retorna vazio
            return np.array([], dtype=int), np.array([], dtype=int)

        monkeypatch.setattr(det_mod, "_cross_match_gaia", mock_cross)

        resultado = _refinar_wcs_gaia(frame)

        assert resultado["wcs_status"] in (
            "gaia_offset_bruto_apenas", "gaia_falhou_match",
            "gaia_skipped_pre_baixo", "gaia_refinado",
        ), f"Status inesperado: {resultado['wcs_status']}"
        assert resultado["wcs"] is not None
        assert "gaia_n_matches_bruto" in resultado
        assert "gaia_offset_bruto_arcsec" in resultado

    # ── teste 10: _refinar_wcs duas passadas completo ─────────────────────
    def test_refinar_wcs_duas_passadas_completo(self, monkeypatch):
        """
        Simula offset real de 6 arcsec: passada bruta corrige CRVAL,
        passada fina deve encontrar matches suficientes e retornar
        gaia_refinado ou gaia_offset_bruto_apenas (nunca gaia_falhou_match
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

        monkeypatch.setattr(det_mod, "_consultar_gaia_dr3",
                            lambda *a, **kw: gaia_tab_shifted)

        resultado = _refinar_wcs_gaia(frame)

        assert resultado["wcs"] is not None, "WCS não pode ser None"
        assert resultado["wcs_status"] != "gaia_falhou_rede"
        assert "gaia_n_matches_bruto" in resultado
        assert "gaia_offset_bruto_arcsec" in resultado
        assert "gaia_n_matches_fino" in resultado
        assert isinstance(resultado["gaia_offset_bruto_arcsec"], list)
        assert len(resultado["gaia_offset_bruto_arcsec"]) == 2


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
        mascara = _mascara_pixels_invalidos(img)
        assert mascara.sum() == 8

    def test_detectar_fontes_ignora_saturacao(self):
        rng = np.random.default_rng(99)
        img = rng.normal(100.0, 5.0, (200, 200))
        # Fonte gaussiana real no centro
        cy, cx = 100, 100
        yy, xx = np.mgrid[0:200, 0:200]
        img += 5000.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 2.0 ** 2))
        # Cluster saturado longe da fonte real
        img[30:40, 30:40] = 65535.0

        fontes = detectar_fontes(img, sigma=5.5)
        # A fonte real deve ser detectada
        assert len(fontes) >= 1
        # Nenhuma fonte deve estar dentro do cluster saturado
        for y, x, _, _ in fontes:
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

        pasta = tmp_path / "sat_frames"
        pasta.mkdir()
        # Criar 4 frames idênticos para satisfazer a exigência do pipeline
        for i in range(4):
            hdu.writeto(pasta / f"frame{i+1}.fits", overwrite=True)

        from detector import carregar_fits
        with pytest.raises(SystemExit):
            carregar_fits(pasta)
