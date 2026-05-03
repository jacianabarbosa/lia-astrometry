[Read in English](README.md)

# Lia

**Lia** é um pipeline Python de pré-triagem astrométrica para candidatos a objetos em movimento em sequências FITS do IASC/Pan-STARRS.

A ferramenta combina estimativa local de fundo, detecção de fontes, refinamento de centróides, ajuste astrométrico com Gaia DR3, rastreamento entre frames e ranqueamento heurístico para priorizar candidatos antes da revisão manual no Astrometrica.

A Lia é uma ferramenta de apoio. Ela **não** confirma descobertas de asteroides, não substitui o Astrometrica, não substitui o Minor Planet Center (MPC) e não substitui o fluxo oficial de campanhas do IASC.

> O nome **Lia** é usado como nome do projeto e dedicação pessoal da autora. Não é uma sigla.

- **Autora:** Jaciana Barbosa
- **Repositório:** `lia-astrometry`
- **Versão:** `v1.5.0`
- **Licença:** MIT

---

## Sumário

- [Visão geral](#visão-geral)
- [Por que isso importa](#por-que-isso-importa)
- [Objetivo científico](#objetivo-científico)
- [O que a Lia faz](#o-que-a-lia-faz)
- [Metodologia resumida](#metodologia-resumida)
- [Exemplo de saída](#exemplo-de-saída)
- [Instalação rápida](#instalação-rápida)
- [Uso básico](#uso-básico)
- [Modos astrométricos: Gaia e header](#modos-astrométricos-gaia-e-header)
- [Saídas geradas](#saídas-geradas)
- [Entendendo o score](#entendendo-o-score)
- [Exportação de métricas](#exportação-de-métricas)
- [Plano de validação](#plano-de-validação)
- [Limitações](#limitações)
- [Mais documentação](#mais-documentação)
- [Estrutura do projeto](#estrutura-do-projeto)
- [Referências](#referências)
- [Licença](#licença)

---

## Visão geral

A Lia foi criada para auxiliar a análise de conjuntos FITS usados em campanhas de busca de asteroides, especialmente no contexto do **International Astronomical Search Collaboration (IASC)**.

Em um conjunto típico do IASC, quatro imagens FITS do mesmo campo do céu são obtidas em sequência. Estrelas e galáxias permanecem praticamente fixas nesse intervalo curto, enquanto objetos do Sistema Solar, como candidatos a asteroides, podem aparecer como fontes pontuais que se deslocam levemente entre os frames.

A Lia processa esses quatro frames, identifica fontes, procura trilhas com movimento coerente e gera uma lista ranqueada de candidatos para inspeção humana. O objetivo não é substituir o observador, mas reduzir o espaço inicial de busca e orientar a validação manual no Astrometrica.

---

## Por que isso importa

Asteroides são pequenos corpos rochosos remanescentes da formação do Sistema Solar. Encontrar e medir asteroides ajuda astrônomos a melhorar órbitas, identificar novos objetos do cinturão principal e apoiar o trabalho de defesa planetária relacionado a objetos próximos da Terra.

O **International Astronomical Search Collaboration (IASC)** é um projeto de ciência cidadã da NASA Science em que equipes inspecionam imagens de telescópios profissionais e submetem medições validadas de asteroides por meio de um processo oficial de campanha.

Nas campanhas gerais de caça a asteroides, os conjuntos de imagens do IASC são enviados às equipes participantes pelos organizadores da campanha. Segundo o IASC, essas imagens são fornecidas pelo Institute for Astronomy da University of Hawaii e obtidas com o telescópio Pan-STARRS de 1,8 m em Haleakalā, ao longo da eclíptica, onde muitos asteroides são encontrados.

Essas são observações reais de telescópio, não imagens geradas por IA nem imagens sintéticas.

A Lia foi criada para facilitar a primeira triagem desses dados, especialmente para estudantes, professores, equipes de ciência cidadã e observadores que desejam organizar melhor a inspeção antes da validação manual.

Fontes:

- [NASA Science: International Astronomical Search Collaboration](https://science.nasa.gov/citizen-science/international-astronomical-search-collaboration/)
- [Registro de campanhas do IASC](https://iasc.cosmosearch.org/Home/Registration)

---

## Objetivo científico

A Lia foi desenvolvida para responder uma pergunta prática no fluxo de campanhas do IASC:

> É possível reduzir a quantidade de candidatos a serem inspecionados manualmente e priorizar objetos que possam ser recuperados, medidos e reportados no Astrometrica?

O pipeline não busca confirmar asteroides de forma automática. Seu papel é:

- reduzir o espaço inicial de busca;
- priorizar candidatos com comportamento cinemático e morfológico plausível;
- gerar saídas auditáveis;
- apoiar a inspeção humana;
- registrar métricas para comparação posterior com Astrometrica e, quando disponível, retorno operacional do IASC.

A avaliação científica da Lia será baseada na proporção de candidatos priorizados que podem ser recuperados, revisados, medidos e eventualmente incluídos em relatórios MPC gerados pelo Astrometrica.

---

## O que a Lia faz

Conjuntos de prática e campanha do IASC/Pan-STARRS geralmente contêm quatro frames FITS reais do mesmo campo do céu.

A Lia:

- carrega as quatro imagens FITS;
- mascara pixels saturados, inválidos ou não finitos;
- estima o fundo local do céu;
- detecta fontes pontuais;
- refina centróides com precisão sub-pixel;
- reconstrói o WCS a partir do header FITS;
- opcionalmente refina a solução astrométrica com Gaia DR3;
- associa detecções entre os quatro frames;
- valida a coerência cinemática das trilhas;
- rejeita ou penaliza falsos positivos prováveis;
- classifica candidatos de `STRONG` a `DISCARDED`;
- consulta objetos conhecidos próximos via SkyBot/IMCCE, quando habilitado;
- gera saídas em JSON, texto, PNG, CSV e rascunho auxiliar em estilo MPC.

O resultado é uma lista priorizada de candidatos. Um score alto significa:

> “inspecione isto primeiro”

e não:

> “isto é um asteroide confirmado”.

---

## Metodologia resumida

A Lia processa cada conjunto como uma sequência temporal curta de quatro frames FITS.

Primeiro, o pipeline lê os dados da imagem, preservando precisão numérica em `float64`, extrai metadados de tempo e reconstrói o WCS a partir do header FITS. Quando possível, o tempo de meio da exposição é calculado e usado como referência temporal para as posições dos candidatos.

Em seguida, a Lia estima o fundo local do céu com `Background2D`, o que torna a detecção mais robusta em imagens com gradientes, halos, regiões mascaradas ou variação local de ruído. As fontes são detectadas com `DAOStarFinder` e podem ter seus centróides refinados por ajuste PSF/sub-pixel.

No modo `gaia`, a Lia consulta o Gaia DR3, cruza estrelas de referência com fontes detectadas e tenta refinar a solução astrométrica. No modo `header`, o pipeline usa apenas o WCS reconstruído a partir do header FITS. Essa separação permite comparar internamente o impacto do refinamento Gaia sem alterar o restante do pipeline.

As detecções são associadas entre os quatro frames para formar trilhas candidatas. Cada trilha é avaliada por critérios cinemáticos, fotométricos e morfológicos, incluindo linearidade do movimento, consistência da velocidade entre frames, estabilidade de brilho, pontualidade, elongação, FWHM, proximidade de borda, SNR e possíveis artefatos ou fontes estáticas.

A saída final é uma lista ranqueada de candidatos para inspeção manual no Astrometrica.

A metodologia completa está em:

- [Metodologia em Português do Brasil](docs/methodology.pt-BR.md)
- [English methodology](docs/methodology.md)

---

## Exemplo de saída

![Candidatos de exemplo](images/example_candidates.png)

Cada linha representa um candidato. As colunas mostram os quatro frames em ordem cronológica. Círculos marcam o centroide medido, e setas indicam a direção do movimento estimado.

---

## Instalação rápida

```bash
git clone https://github.com/jacianabarbosa/lia-astrometry.git
cd lia-astrometry

python3 -m venv venv
source venv/bin/activate

venv/bin/python -m pip install --upgrade pip
venv/bin/python -m pip install -r requirements.txt
venv/bin/python -m pytest tests/ -v
```

> **macOS com pyenv:** use `pyenv local 3.13.1` antes de criar o ambiente virtual se `python3` não encontrar a versão correta.

A suíte atual tem **92 testes passando**, cobrindo detecção de fontes, WCS, refinamento Gaia, score, máscara de saturação, rastreabilidade JSON e exportação CSV.

---

## Uso básico

Execute a Lia em uma pasta com exatamente quatro arquivos `.fits` ou `.fit`:

```bash
venv/bin/python src/detector.py --images images/XY14_p10 --output results/XY14_p10
```

Opções úteis:

```bash
venv/bin/python src/detector.py --sigma 5.0       # mais sensível, mais falsos positivos
venv/bin/python src/detector.py --sigma 6.5       # mais conservador
venv/bin/python src/detector.py --no-mpc          # pula consulta SkyBot de objetos conhecidos
venv/bin/python src/detector.py --wcs-mode header # usa apenas o WCS reconstruído do header
venv/bin/python src/detector.py --wcs-mode gaia   # padrão: refinamento WCS com Gaia DR3
```

Metadados do observador podem ser passados diretamente ou por variáveis de ambiente:

```bash
venv/bin/python src/detector.py --observer "Seu Nome" --email "nome@exemplo.com"

export LIA_OBS="Seu Nome"
export LIA_EMAIL="nome@exemplo.com"
```

---

## Modos astrométricos: Gaia e header

A Lia possui dois modos astrométricos principais.

### `--wcs-mode gaia`

Modo padrão. A Lia reconstrói o WCS inicial a partir do header FITS e tenta refinar a solução astrométrica usando estrelas de referência do Gaia DR3.

Esse modo é recomendado para o uso operacional do pipeline antes da inspeção no Astrometrica.

```bash
venv/bin/python src/detector.py --images images/XY14_p10 --output results/XY14_p10_gaia --wcs-mode gaia
```

### `--wcs-mode header`

Modo de comparação interna. A Lia usa apenas o WCS reconstruído a partir do header FITS, sem refinamento Gaia DR3.

Todo o restante do pipeline permanece idêntico: estimativa de fundo, detecção, PSF, tracking, score, filtros e geração de saídas.

```bash
venv/bin/python src/detector.py --images images/XY14_p10 --output results/XY14_p10_header --wcs-mode header
```

Esse modo existe para avaliar o impacto do refinamento Gaia nas coordenadas preliminares e no ranking dos candidatos. Ele não representa um fluxo separado de envio ao IASC.

---

## Saídas geradas

Para cada conjunto processado, a Lia gera:

| Arquivo | Função |
|---|---|
| `*_candidates.json` | Dados estruturados dos candidatos, metadados da execução e métricas globais |
| `*_report.txt` | Relatório legível para inspeção humana |
| `*_MPC_report.txt` | Rascunho auxiliar em estilo MPC, apenas para revisão |
| `*_candidates.png` | Recortes visuais dos principais candidatos |
| `*_pipeline.log` | Log de auditoria da execução |

O envio oficial da campanha deve seguir o processo normal Astrometrica/IASC. Use a Lia para decidir onde olhar primeiro; depois, inspecione e meça os candidatos manualmente.

---

## Entendendo o score

O score é um índice heurístico de prioridade operacional, variando de 0 a 12.

Ele indica quais candidatos devem ser inspecionados primeiro, mas **não** representa uma probabilidade estatística de o objeto ser um asteroide.

| Score | Classe | Significado operacional |
|---:|---|---|
| 8–12 | `STRONG` | inspecionar primeiro |
| 6–7 | `MODERATE` | candidato plausível |
| 4–5 | `WEAK` | candidato de baixa prioridade |
| 0–3 | `DISCARDED` | provável artefato, fonte rejeitada ou candidato fraco |

O score combina:

- linearidade do movimento;
- consistência da velocidade entre frames;
- estabilidade fotométrica;
- morfologia pontual;
- consistência morfológica;
- elongação;
- faixa de movimento esperada para sequências curtas do IASC.

Mesmo candidatos `STRONG` precisam ser revisados no Astrometrica. Candidatos `WEAK` ou `DISCARDED` podem ainda ser úteis em casos específicos, especialmente se o campo estiver degradado ou se a fonte tiver baixo SNR.

---

## Exportação de métricas

Depois de rodar vários conjuntos de imagens, exporte CSVs para validação:

```bash
venv/bin/python tools/export_metrics.py results/ --out validation_metrics --recursive
```

Isso cria:

- `validation_metrics_sets.csv`: uma linha por conjunto processado;
- `validation_metrics_candidates.csv`: uma linha por candidato.

O CSV inclui colunas automáticas, como:

- versão do pipeline;
- modo WCS;
- fontes por frame;
- trilhas tentadas;
- candidatos finais;
- distribuição de classes;
- top candidates;
- status Gaia;
- RMS Gaia;
- tempo de execução.

Também inclui colunas vazias para revisão manual posterior, como:

- `measured_in_astrometrica`;
- `included_in_astrometrica_mpc`;
- `iasc_feedback`;
- `manual_classification`;
- `notes`.

Esses campos permitem comparar a saída da Lia com a inspeção no Astrometrica e, quando disponível, com o retorno operacional do IASC.

---

## Plano de validação

A validação da Lia será conduzida em duas frentes: avaliação operacional e comparação astrométrica interna.

### 1. Avaliação operacional

Fluxo principal:

```text
conjunto FITS do IASC
→ Lia com Gaia DR3
→ lista priorizada de candidatos
→ inspeção manual no Astrometrica
→ relatório MPC gerado pelo Astrometrica
→ envio ao IASC
→ retorno operacional do IASC, quando disponível
```

Métricas principais:

- redução do espaço de busca de candidatos;
- utilidade dos candidatos top-1, top-3 e top-5;
- proporção de candidatos mensuráveis no Astrometrica;
- proporção de candidatos incluídos em relatório MPC gerado pelo Astrometrica;
- retorno operacional do IASC, quando disponível;
- falsos positivos e modos de falha;
- diferença entre coordenadas preliminares da Lia e medições no Astrometrica.

### 2. Comparação interna Gaia vs header

Cada conjunto pode ser executado nos dois modos:

```bash
venv/bin/python src/detector.py --wcs-mode gaia
venv/bin/python src/detector.py --wcs-mode header
```

Essa comparação permite avaliar:

- recuperação de candidatos;
- mudanças de ranking;
- diferenças de RA/Dec;
- separação angular em relação às medições no Astrometrica;
- RMS do refinamento Gaia;
- número de estrelas Gaia usadas;
- casos em que Gaia melhora, não altera ou piora a posição preliminar.

Apenas o relatório final gerado pelo Astrometrica deve ser enviado ao IASC. A comparação `gaia`/`header` é interna e metodológica.

---

## Limitações

- A Lia não confirma descobertas de asteroides.
- A Lia não substitui o Astrometrica.
- A Lia não substitui o Minor Planet Center.
- A Lia não substitui o fluxo oficial do IASC.
- O score é heurístico e não é um valor de confiança estatisticamente calibrado.
- O refinamento Gaia DR3 pode falhar ou não melhorar todos os campos.
- O WCS reconstruído a partir do header pode ser suficiente em alguns casos.
- O ajuste PSF pode falhar para fontes com baixo SNR, em blend, saturadas, com trail ou próximas à borda.
- Os parâmetros do `Background2D` precisam de validação empírica em conjuntos de prática e campanha do IASC.
- As heurísticas de rejeição de falsos positivos podem descartar candidatos reais em condições degradadas.
- O retorno do IASC é uma validação operacional, não uma verdade absoluta universal.
- O rascunho MPC gerado pela Lia é apenas auxiliar; medições finais devem ser revisadas ou produzidas no Astrometrica.

Funcionalidades não implementadas nesta versão:

- determinação de órbita;
- score digest2;
- shift-and-stack;
- classificação por CNN;
- interface gráfica;
- envio automático ao MPC/IASC;
- calibração estatística do score.

---

## Mais documentação

Para a explicação científica e técnica completa, leia:

- [Metodologia em Português do Brasil](docs/methodology.pt-BR.md)
- [English methodology](docs/methodology.md)

Esses documentos explicam tempo em FITS, estimativa de fundo, ajuste PSF, refinamento Gaia DR3, modo WCS apenas pelo header, SkyBot, score, validação, limitações e por que a revisão humana continua necessária.

---

## Estrutura do projeto

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

---

## Referências

- [NASA Science: International Astronomical Search Collaboration](https://science.nasa.gov/citizen-science/international-astronomical-search-collaboration/)
- [Site oficial do IASC](https://iasc.cosmosearch.org/)
- [Registro de campanhas do IASC](https://iasc.cosmosearch.org/Home/Registration)
- [Colaboradores do IASC](https://iasc.cosmosearch.org/Home/Collaborators)
- Gaia Collaboration et al. (2023), *Gaia Data Release 3*, Astronomy & Astrophysics, 674, A1.
- Stetson, P. B. (1987), *DAOPHOT: A Computer Program for Crowded-Field Stellar Photometry*, Publications of the Astronomical Society of the Pacific, 99, 191.
- Astropy Collaboration et al. (2022), *The Astropy Project*.
- Bradley et al. (2024), *astropy/photutils: source detection and photometry tools*.
- Documentação do SkyBot/IMCCE, usado por meio de `astroquery.imcce.Skybot`.
- Documentação de astrometria e formato de observação do Minor Planet Center.

---

## Licença

Licença MIT. Consulte [LICENSE](LICENSE).
