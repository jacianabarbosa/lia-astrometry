[Read in English](README.md)

# Lia

**Lia** é um pipeline em Python para pré-triagem de candidatos a objetos em movimento em sequências de imagens FITS do IASC, com refinamento astrométrico via Gaia DR3 e suporte à validação manual no Astrometrica.

A Lia **não substitui** o Astrometrica, a validação pelo MPC nem o fluxo oficial do IASC. Ela é uma camada de pré-triagem e priorização que reduz o esforço de inspeção manual e gera saídas auditáveis para a validação humana posterior.

> O nome **Lia** é utilizado como nome do projeto e dedicação pessoal da autora. Não é uma sigla.

**Autora:** Jaciana Barbosa  
**Repositório:** `lia-astrometry`  
**Versão:** `v1.5.0`  
**Licença:** MIT

## Visão Geral

Os conjuntos de prática e de campanha do IASC/Pan-STARRS contêm quatro frames FITS do mesmo campo. A tarefa científica é encontrar fontes que se movem de forma coerente enquanto a maioria das estrelas e galáxias permanece fixa. A Lia automatiza o primeiro passo de triagem: detecta fontes pontuais, constrói trilhas candidatas, rejeita artefatos óbvios, classifica os candidatos restantes e produz saídas que podem ser revisadas manualmente no Astrometrica.

O resultado é intencionalmente conservador. Um candidato bem classificado não é uma descoberta confirmada de asteroide — é um candidato que vale a pena medir ou descartar no fluxo normal do Astrometrica.

## Exemplo de Saída

![Candidatos de exemplo](images/example_candidates.png)

Visualização produzida pela Lia para uma sequência do IASC. Cada linha é um candidato; as colunas são os quatro frames em ordem cronológica. Círculos coloridos marcam o centroide medido e setas mostram a direção do movimento.

## Posicionamento

A Lia é um pipeline de pré-triagem. Ela prioriza candidatos; não confirma descobertas, não envia observações e não substitui o processo final de medição.

As medidas astrométricas finais devem ser feitas ou revisadas no Astrometrica. Os relatórios no formato MPC gerados pela Lia são rascunhos auxiliares apenas. O envio oficial da campanha deve continuar seguindo o fluxo padrão Astrometrica/IASC.

O score de triagem é heurístico. É útil para classificação e auditoria, mas não é uma estimativa de confiança estatisticamente calibrada.

## Funcionalidades

- Carregamento de FITS com tratamento de imagens em `float64` para estabilidade numérica.
- Extração de JD/MJD de meio da exposição a partir de `DATE-OBS` ou `MJD-OBS` mais o tempo de exposição.
- Estimativa local do céu com `photutils.Background2D`.
- Mascaramento de pixels Pan-STARRS marcados próximos de `65535`, típico de regiões saturadas, defeituosas ou fora do CCD, evitando detecções espúrias e estimativas infladas de RMS local do fundo.
- Detecção de fontes com `DAOStarFinder`.
- Refinamento de centroide sub-pixel usando modelo PSF Moffat com fallback Gaussiano.
- Verificações de FWHM e morfologia para rejeitar hot pixels, raios cósmicos, blends e fontes estendidas.
- Refinamento WCS opcional com Gaia DR3 usando correspondência cruzada em dois estágios e rastreamento de resíduos.
- Modo WCS apenas pelo header para comparações internas entre Gaia e header.
- Construção de trilhas de quatro frames em representação de plano tangente local.
- Validação cinemática linear com resíduos e `R^2`.
- Verificação de vizinhança posicional via SkyBot/IMCCE para objetos conhecidos do Sistema Solar.
- Exportação em JSON, relatório de texto, rascunho MPC, recortes PNG e planilha CSV para validação.

## Fluxo de Trabalho

```text
conjunto FITS do IASC
-> pré-triagem com a Lia
-> lista priorizada de candidatos
-> inspeção manual e medição no Astrometrica
-> relatório MPC gerado pelo Astrometrica
-> envio ao IASC quando apropriado
-> retorno operacional do IASC, quando disponível
```

## Instalação

```bash
git clone https://github.com/jacianabarbosa/lia-astrometry.git
cd lia-astrometry

python3 -m venv venv
source venv/bin/activate

venv/bin/python -m pip install --upgrade pip
venv/bin/python -m pip install -r requirements.txt
venv/bin/python -m pytest tests/test_basic.py -v
```

> **macOS com pyenv:** use `pyenv local 3.13.1` antes de criar o venv se `python3` não encontrar a versão correta.

A suíte atual tem 92 testes cobrindo identidade pública, detecção de fontes, tratamento WCS, refinamento Gaia, score cinemático, mascaramento de saturação, modo header, rastreabilidade JSON e exportação CSV de métricas.

## Uso

Execução básica:

```bash
venv/bin/python src/detector.py --images images/XY14_p10 --output results/XY14_p10
```

Ajustar o threshold de detecção:

```bash
venv/bin/python src/detector.py --sigma 5.0   # mais sensível — mais falsos positivos
venv/bin/python src/detector.py --sigma 6.5   # menos sensível — mais conservador
```

Pular a consulta de vizinhança SkyBot/MPC:

```bash
venv/bin/python src/detector.py --no-mpc
```

Metadados do observador podem ser passados diretamente ou por variáveis de ambiente:

```bash
venv/bin/python src/detector.py --observer "Jaciana Barbosa" --email "nome@exemplo.com"
export LIA_OBS="Jaciana Barbosa"
export LIA_EMAIL="nome@exemplo.com"
```

## Opções de Linha de Comando

```text
usage: detector.py [-h] [--images IMAGES] [--output OUTPUT] [--sigma SIGMA]
                   [--no-mpc] [--observer OBSERVER] [--email EMAIL]
                   [--wcs-mode {gaia,header}]

Lia - pre-screen moving-object candidates in IASC FITS sequences

options:
  -h, --help            mostra a ajuda e encerra
  --images IMAGES       pasta de entrada com exatamente quatro frames FITS
  --output OUTPUT       pasta de saída para relatórios, JSON, logs e figuras
  --sigma SIGMA         threshold de detecção em SNR local, padrão 5.5
  --no-mpc              pula verificações de vizinhança SkyBot/IMCCE
  --observer OBSERVER   nome do observador para rascunhos auxiliares
  --email EMAIL         email do observador para rascunhos auxiliares
  --wcs-mode {gaia,header}
                        modo WCS: refinamento Gaia DR3 ou WCS reconstruído do header
```

## Modos WCS: Gaia e Apenas Header

O refinamento Gaia DR3 é o modo padrão:

```bash
venv/bin/python src/detector.py --images fits/conjunto01 --output results/conjunto01_gaia --wcs-mode gaia
```

O modo apenas header desativa somente o refinamento Gaia. A estimativa de fundo, o ajuste PSF, o rastreamento, o score, os filtros, os thresholds e a geração de saídas permanecem idênticos:

```bash
venv/bin/python src/detector.py --images fits/conjunto01 --output results/conjunto01_header --wcs-mode header
```

Nota técnica: a Lia não chama `astropy.wcs.WCS(header)` diretamente para frames Pan-STARRS. Esses headers FITS podem incluir chaves proprietárias de distorção como `PCA1X0Y2` e `PCA2X0Y2`, que o parser padrão pode interpretar como entradas de matriz PC e falhar com matrizes singulares. A Lia reconstrói o WCS do header a partir das chaves astrométricas primárias (`CTYPE`, `CRVAL`, `CRPIX`, `CD`/`CDELT` e `CROTA`). Os dois modos partem dessa reconstrução; apenas a etapa de refinamento Gaia DR3 muda.

Isso permite uma comparação interna usando o mesmo código:

```text
recuperação de candidatos
mudanças no ranking
diferenças de RA/Dec
separação angular em relação às medições do Astrometrica
RMS e número de matches do Gaia
casos em que o Gaia melhora, não altera ou piora as coordenadas preliminares
```

Apenas o relatório final gerado pelo Astrometrica deve ser enviado ao IASC. A comparação Gaia/header é metodológica e interna; não envie relatórios duplicados.

## Saídas

Para cada conjunto processado, a Lia gera:

- `*_candidates.json`: metadados da execução, candidatos, trilhas, scores, morfologia, status WCS, status Gaia, status SkyBot e campos de validação manual (vazios para preenchimento posterior).
- `*_report.txt`: relatório operacional legível para revisão dos candidatos.
- `*_MPC_report.txt`: rascunho auxiliar no formato MPC, apenas para revisão.
- `*_candidates.png`: recortes visuais dos principais candidatos.
- `*_pipeline.log`: trilha de auditoria científica da execução.

O JSON contém `run_metadata` com `run_id`, `pipeline_name`, `pipeline_version`, `repository`, `input_set`, `timestamp_execution_utc`, `wcs_mode` e `sigma`.

Cada candidato inclui `candidate_id`, `triage_class`, `heuristic_score`, `score_components`, `decision_reasons`, `track`, `motion`, `photometry`, `morphology` e campos de revisão manual para avaliação posterior no Astrometrica/IASC.

## Exportação de Métricas

A Lia inclui um exportador CSV para estudos de validação:

```bash
venv/bin/python tools/export_metrics.py results/ --out validation_metrics --recursive
```

Isso gera:

- `validation_metrics_sets.csv`: uma linha por conjunto processado.
- `validation_metrics_candidates.csv`: uma linha por candidato.

Os arquivos CSV incluem colunas vazias como `measured_in_astrometrica`, `included_in_astrometrica_mpc`, `iasc_feedback`, `manual_classification` e `notes`, onde os resultados da campanha podem ser inseridos após a revisão humana.

## Score e Flags

O score de triagem varia de 0 a 12 e combina:

| Componente | Máx. | Justificativa |
|---|---:|---|
| Linearidade | 3 | resíduo da trilha de quatro frames em torno de uma trajetória linear |
| Consistência de velocidade | 2 | uniformidade dos passos entre frames consecutivos |
| Estabilidade fotométrica | 2 | estabilidade relativa do fluxo ao longo dos frames |
| Morfologia pontual | 2 | rejeição de fontes difusas, blends ou artefatos |
| Consistência morfológica | 1 | estabilidade da forma da fonte entre frames |
| Elongação | 1 | penalização de detecções com trail ou extensão |
| Faixa típica de movimento | 1 | deslocamento total compatível com candidatos do cinturão principal em sequências curtas do IASC |

Total: 12. O score é uma heurística de ranking, não uma confiança calibrada.

Classes:

| Score | Classe | Significado |
|---:|---|---|
| 8–12 | `STRONG` | inspecionar primeiro |
| 6–7 | `MODERATE` | candidato plausível |
| 4–5 | `WEAK` | candidato de baixa prioridade |
| 0–3 | `DISCARDED` | provável artefato ou fonte rejeitada |

O score é uma heurística de triagem reproduzível, não um valor de confiança estatisticamente calibrado.

## Plano de Validação

### Avaliação Operacional Principal

Execute a Lia no modo Gaia em conjuntos do IASC:

```text
conjunto FITS do IASC
-> Lia com Gaia DR3
-> lista priorizada de candidatos
-> inspeção manual no Astrometrica
-> relatório MPC gerado pelo Astrometrica
-> envio ao IASC
-> retorno operacional do IASC, quando disponível
```

Métricas principais a coletar:

- redução do espaço de busca de candidatos;
- utilidade dos candidatos top-1/top-3/top-5;
- candidatos mensuráveis no Astrometrica;
- candidatos incluídos em um relatório MPC gerado pelo Astrometrica;
- retorno do IASC, quando disponível;
- falsos positivos e modos de falha.

### Comparação Interna Gaia vs Header

Execute os mesmos conjuntos nos dois modos:

```bash
venv/bin/python src/detector.py --wcs-mode gaia
venv/bin/python src/detector.py --wcs-mode header
```

Compare a recuperação de candidatos, mudanças de ranking, RMS do Gaia, contagem de matches e diferenças angulares em relação às medições do Astrometrica. Essa comparação é interna; apenas o relatório final gerado pelo Astrometrica deve entrar no fluxo oficial do IASC.

### Calibração Estatística Futura

O score heurístico é uma ferramenta de ranking. Um trabalho futuro pode calibrá-lo contra resultados validados pelo IASC usando regressão logística, análise ROC ou calibração beta para estimar confiança calibrada de candidato real. Isso exige volume suficiente de resultados positivos e negativos no fluxo Astrometrica/IASC.

## Limitações

- A Lia não substitui o Astrometrica.
- O score de triagem é heurístico, não é um valor de confiança estatisticamente calibrado.
- O refinamento Gaia DR3 pode falhar ou não melhorar todos os campos.
- O WCS do header pode já ser suficiente em alguns frames do IASC.
- O ajuste PSF pode falhar para fontes com baixo SNR, em blend, saturadas ou na borda.
- Os parâmetros do `Background2D` precisam de validação empírica com dados de prática e campanha do IASC.
- A rejeição de falsos positivos pode descartar candidatos reais em condições de seeing degradado ou com artefatos severos.
- O retorno do IASC é validação operacional, não verdade absoluta universal.
- A Lia não realiza determinação de órbita, score digest2, shift-and-stack, classificação por CNN, fluxo com interface gráfica nem envio automático ao MPC/IASC.

## Estrutura do Projeto

```text
lia-astrometry/
├── README.md
├── README.pt-BR.md
├── docs/
│   ├── methodology.md
│   └── methodology.pt-BR.md
├── src/
│   └── detector.py
├── tests/
│   └── test_basic.py
├── tools/
│   └── export_metrics.py
├── requirements.txt
└── pytest.ini
```

## Metodologia

A metodologia técnica completa está em [docs/methodology.pt-BR.md](docs/methodology.pt-BR.md). Ela cobre o tratamento de metadados FITS, tempo de meio da exposição, Background2D, ajuste PSF, correspondência cruzada com Gaia DR3, resíduos WCS, validação cinemática, rejeição de falsos positivos, esquema de saídas e validação planejada.

## Referências

- Gaia Collaboration et al. (2023), *Gaia Data Release 3*, Astronomy & Astrophysics, 674, A1.
- Stetson, P. B. (1987), *DAOPHOT: A Computer Program for Crowded-Field Stellar Photometry*, Publications of the Astronomical Society of the Pacific, 99, 191.
- Bradley et al. (2024), *astropy/photutils: source detection and photometry tools*.
- Astropy Collaboration et al. (2022), *The Astropy Project: sustaining and growing a community-oriented open-source project*.
- Documentação de cone search do SkyBot/IMCCE, usada via `astroquery.imcce.Skybot`.
- Smullen et al. (2025), *TRIPP: TRansient Image Processing Pipeline*, arXiv:2501.18142.
- Materiais do International Astronomical Search Collaboration (IASC) e do fluxo legado AIsteroid.
- Documentação de formato de observação e astrometria do Minor Planet Center.
- Materiais públicos de campanha do IASC e orientações do fluxo com Astrometrica.

## Licença

Licença MIT. Consulte [LICENSE](LICENSE).

## Histórico de Versões

### v1.5.0

- Projeto renomeado para **Lia**.
- Identidade do repositório atualizada para `lia-astrometry`.
- Documentação atualizada para publicação internacional.
- Campos de rastreabilidade científica adicionados ou revisados para validação com os fluxos IASC/Astrometrica.
- Esclarecimento de que a Lia é um pipeline de pré-triagem e não substitui o Astrometrica nem a validação oficial MPC/IASC.
- Documentação preparada para comparação interna entre o refinamento Gaia DR3 e o modo WCS apenas pelo header.

### Versões Anteriores

Versões internas anteriores desenvolveram o núcleo de detecção, Background2D, ajuste PSF, refinamento WCS com Gaia DR3, verificações SkyBot, logs de auditoria e geração de rascunho MPC. Lia é a identidade pública atual desse trabalho.
