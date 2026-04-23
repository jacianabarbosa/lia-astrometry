"""
testes/teste_basico.py
======================
Testes unitários e de integração do pipeline asteroid-hunter v1.3.

Cobertura:
    - detectar_fontes: fonte simples, múltiplas, imagem vazia
    - formatar_ra_mpc / formatar_dec_mpc: casos normais e de wraparound
    - formatar_data_mpc: conversão JD → formato MPC
    - _construir_wcs_panstarrs: header Pan-STARRS com chaves PCA* proprietárias
    - _casar_mutuo: validação recíproca de vizinho mais próximo
    - _validar_coerencia_trilha: coerência cinemática de trilhas (v1.3)
    - _morfologia_por_frame: elongation, pontualidade, FWHM (v1.3)
    - _verificar_falso_positivo: rejeição por borda, SNR, hot pixel, brilho caótico (v1.3)
    - analisar_candidato: score decomposto, novos campos v1.3, flags, razões de decisão
    - consultar_catalogo_posicional: schema normalizado de retorno
    - Integração: pipeline end-to-end com FITS sintético (JSON v1.3)

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
    _validar_coerencia_trilha,
    _morfologia_por_frame,
    _verificar_falso_positivo,
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
    pipeline completo. Verifica estrutura do JSON v1.3:
    - metricas_globais (incluindo n_rejeitados_coerencia)
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

    def test_pipeline_detecta_asteroide_sintetico(self, tmp_path):
        imagens    = tmp_path / "imagens"
        resultados = tmp_path / "resultados"
        imagens.mkdir(); resultados.mkdir()
        self._gerar_fits_sinteticos(imagens)

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
                      "n_rejeitados_coerencia",
                      "sigma_deteccao", "deriva_dx_px", "deriva_dy_px",
                      "tempo_execucao_s"):
            assert campo in mg, f"Métrica global '{campo}' ausente"

        assert isinstance(mg["n_fontes_por_frame"], list)
        assert len(mg["n_fontes_por_frame"]) == 4
        assert isinstance(mg["n_rejeitados_coerencia"], int)


# ─────────────────────────────────────────────
# Coerência cinemática de trilhas (v1.3)
# ─────────────────────────────────────────────
class TestValidarCoerenciaTrilha:

    def _trilha(self, passos):
        """Constrói trilha com passos acumulados a partir de (100, 100)."""
        pontos = [(100.0, 100.0, 1000.0)]
        for dy, dx in passos:
            y0, x0 = pontos[-1][0], pontos[-1][1]
            pontos.append((y0 + dy, x0 + dx, 1000.0))
        return pontos

    def test_trilha_linear_aceita(self):
        trilha = self._trilha([(2.0, 3.0), (2.0, 3.0), (2.0, 3.0)])
        valida, razoes = _validar_coerencia_trilha(trilha)
        assert valida, f"Trilha linear rejeitada: {razoes}"

    def test_trilha_com_inversao_de_direcao_rejeitada(self):
        # Passo 1: vai para (+2,+3); passo 2: inverte completamente (-2,-3)
        trilha = self._trilha([(2.0, 3.0), (-2.0, -3.0), (2.0, 3.0)])
        valida, razoes = _validar_coerencia_trilha(trilha)
        assert not valida
        assert any("direção" in r for r in razoes)

    def test_trilha_com_razao_passo_excessiva_rejeitada(self):
        # Passo 1 pequeno, passo 2 muito grande (ratio >> T.MAX_STEP_RATIO)
        trilha = self._trilha([(1.0, 1.0), (20.0, 20.0), (1.0, 1.0)])
        valida, razoes = _validar_coerencia_trilha(trilha)
        assert not valida
        assert any("razão" in r or "aceleração" in r for r in razoes)

    def test_trilha_com_aceleracao_excessiva_rejeitada(self):
        # Passos crescendo demais: 1, 1, 10 px/frame
        trilha = self._trilha([(1.0, 1.0), (1.0, 1.0), (10.0, 10.0)])
        valida, razoes = _validar_coerencia_trilha(trilha)
        assert not valida

    def test_trilha_quase_linear_com_pequeno_desvio_aceita(self):
        # Desvio angular pequeno deve ser aceito (< T.MAX_ANGULO_DEG)
        trilha = self._trilha([(2.0, 3.0), (2.1, 2.9), (1.9, 3.1)])
        valida, razoes = _validar_coerencia_trilha(trilha)
        assert valida, f"Trilha quase linear rejeitada: {razoes}"

    def test_retorna_lista_razoes(self):
        trilha = self._trilha([(2.0, 3.0), (2.0, 3.0), (2.0, 3.0)])
        valida, razoes = _validar_coerencia_trilha(trilha)
        assert isinstance(razoes, list)

    def test_thresholds_centralizados(self):
        """T.MAX_ANGULO_DEG e T.MAX_STEP_RATIO devem existir e ser positivos."""
        assert T.MAX_ANGULO_DEG > 0
        assert T.MAX_STEP_RATIO > 1


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
            "RAIO_MATCH_PX", "MAX_ANGULO_DEG", "MAX_STEP_RATIO",
            "MOVE_MIN_PX", "MOVE_MAX_PX", "RESIDUO_MIN_PX",
            "MARGEM_BORDA_PX", "SNR_MINIMO", "BRILHO_CV_CAOS",
            "LIN_EXCELENTE", "LIN_BOA", "LIN_MARGINAL",
            "VEL_UNIFORME", "VEL_MODERADA",
            "FOT_ESTAVEL", "FOT_MODERADA",
            "MOR_PONTUAL", "MOR_MARGINAL",
            "ELON_BOA", "MOR_CONSIST_MAX_STD",
            "FAIXA_VEL_MIN_PX", "FAIXA_VEL_MAX_PX",
            "SCORE_FORTE", "SCORE_MODERADO", "SCORE_FRACO",
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
