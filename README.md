# asteroid-hunter v1.3

Pipeline Python de pré-triagem de asteroides em imagens FITS do IASC/Pan-STARRS.

**Autora:** Jaciana Barbosa — Cientista da Computação & Astrônoma Amadora  
Icapuí, Ceará, Brasil  
Participante do IASC (International Astronomical Search Collaboration) e do Caça Asteroides MCTI

---

## Sumário

- [Sobre o Projeto](#sobre-o-projeto)
- [Posicionamento](#posicionamento)
- [Como Funciona](#como-funciona)
- [Instalação](#instalação)
- [Uso](#uso)
- [Outputs Gerados](#outputs-gerados)
- [Score e Flags](#score-e-flags)
- [Limitações](#limitações)
- [Estrutura do Projeto](#estrutura-do-projeto)
- [Metodologia](#metodologia)
- [Referências](#referências)
- [Licença](#licença)

---

## Sobre o Projeto

O **asteroid-hunter** automatiza a pré-triagem de candidatos a asteroides em conjuntos de 4 imagens FITS distribuídos pelo programa IASC.

O pipeline detecta objetos em movimento, calcula coordenadas astrométricas, consulta o catálogo SkyBot/IMCCE para verificar se existe objeto conhecido próximo, atribui um score multidimensional a cada candidato e gera relatórios estruturados para uso no Astrometrica.

A ferramenta foi construída para uso pessoal no contexto dos programas IASC e Caça Asteroides MCTI. Não é um produto comercial, não tem interface gráfica e não pretende substituir nenhuma etapa do processo científico formal.

---

## Posicionamento

**O que esta ferramenta faz:**
- Triagem rápida dos candidatos mais prováveis para inspeção no Astrometrica
- Geração de coordenadas RA/Dec para facilitar a localização no Astrometrica
- Consulta de vizinhança posicional via SkyBot (objeto conhecido próximo?)
- Exportação de dados estruturados e auditáveis para revisão posterior

**O que esta ferramenta não faz:**
- Não substitui a validação astrométrica no Astrometrica
- Não calcula magnitudes absolutas (sem calibração fotométrica)
- Não realiza determinação orbital
- Não garante que candidatos com score alto sejam asteroides reais
- Não deve ser usada para enviar relatórios ao IASC sem revisão humana

Todo candidato apontado pelo pipeline deve ser inspecionado visualmente no Astrometrica antes de qualquer envio. O score é um auxiliar de priorização, não uma decisão.

---

## Como Funciona

```
4 arquivos .fits
      │
      ▼
[1] Carregamento         Lê FITS, extrai WCS e timestamps
      │
      ▼
[2] Detecção de fontes   DAOStarFinder + sigma-clipping
      │
      ▼
[3] Rastreamento         Casamento mútuo recíproco + validação cinemática
                         (ângulo de direção, razão de passo, aceleração)
      │
      ▼
[4] Remoção de deriva    Mediana dos vetores de fontes estáveis
      │
      ▼
[5] Rejeição de FP       Borda, SNR baixo, hot pixel, brilho caótico,
                         elongação grave e sistemática
      │
      ▼
[6] Score (0–12)         Linearidade + Velocidade + Fotometria + Morfologia
                         + Consistência morfológica + Elongação + Faixa vel.
      │
      ▼
[7] Conversão RA/Dec     WCS primário do header Pan-STARRS
      │
      ▼
[8] Consulta SkyBot      Vizinhança posicional via astroquery/IMCCE
      │
      ▼
[9] Outputs              .json  .txt  _MPC_report.txt  .png  .log
```

---

## Instalação

**Pré-requisitos:** Python 3.10+, pip, conexão com internet (para SkyBot).

```bash
# Clone ou extraia o repositório
git clone https://github.com/jacianabarbosa/asteroid-hunter.git
cd asteroid-hunter

# Ambiente virtual (recomendado)
python3 -m venv venv
source venv/bin/activate        # macOS/Linux
# venv\Scripts\activate         # Windows

# Dependências
source venv/bin/activate  

# Verificação
python -m pytest testes/teste_basico.py -v
```

> **macOS com pyenv:** use `pyenv local 3.13.1` antes de criar o venv se `python3` não encontrar a versão correta.

---

## Uso

```bash
# Modo básico (4 arquivos .fits na pasta imagens/)
python src/detector.py

# Pasta personalizada
python src/detector.py --imagens /meus/fits/ --output /meus/resultados/

# Ajustar sensibilidade (padrão: 5.5σ)
python src/detector.py --sigma 5.0   # mais sensível — mais falsos positivos
python src/detector.py --sigma 6.5   # menos sensível — mais conservador

# Modo offline (sem consulta ao SkyBot)
python src/detector.py --sem-mpc

# Identificação do observador (necessário antes de enviar ao IASC)
python src/detector.py --observador "Jaciana Barbosa" --email "seu@email.com"

# Ou via variável de ambiente (evita repetir a cada execução)
export ASTEROID_HUNTER_OBS="Jaciana Barbosa"
export ASTEROID_HUNTER_EMAIL="seu@email.com"
```

**Todas as opções:**

```
--imagens     Pasta com os 4 arquivos FITS     (padrão: imagens/)
--output      Pasta para salvar resultados     (padrão: resultados/)
--sigma       Threshold de detecção em sigma   (padrão: 5.5)
--sem-mpc     Pular consulta ao catálogo
--observador  Nome do observador
--email       Email do observador
```

---

## Outputs Gerados

Para cada conjunto processado, o pipeline gera 5 arquivos em `resultados/`:

### `*_candidatos.json` — Registro estruturado e auditável

O JSON é o output principal para revisão e integração com outras ferramentas. Contém:

- **Cabeçalho global:** pipeline, conjunto, data de processamento, observatório
- **Métricas globais do conjunto:** fontes por frame, trilhas tentadas, deriva do campo, sigma usado, tempo de execução, distribuição de classes
- **Por candidato:**
  - `rank`, `classe`, `score_total`, `probabilidade`
  - `score_componentes`: decomposição completa (linearidade/3, velocidade/2, fotometria/2, morfologia/2, faixa_vel/1)
  - `flags`: lista de flags diagnósticas derivadas do pipeline
  - `posicao_frame1`: pixel x/y, RA/Dec, formatos MPC
  - `movimento`: deslocamento total, resíduo de deriva, linearidade, consistência de velocidade, velocidade em arcsec/min
  - `fotometria`: brilho por frame, CV, pontualidade, SNR estimado por frame
  - `trilha`: array com posição, coordenadas, fluxo e SNR de cada frame
  - `mpc`: status normalizado da consulta SkyBot com campos estruturados quando há match

Exemplo de estrutura:

```json
{
  "pipeline": "asteroid-hunter v1.2.0",
  "conjunto": "o8723g0213o",
  "metricas_globais": {
    "n_fontes_por_frame": [42, 39, 41, 40],
    "n_trilhas_tentadas": 35,
    "deriva_dx_px": 0.62,
    "deriva_dy_px": 0.41,
    "sigma_deteccao": 5.5,
    "tempo_execucao_s": 12.3
  },
  "candidatos": [
    {
      "rank": 1,
      "classe": "FORTE",
      "score_total": 9,
      "score_componentes": {
        "linearidade": 3,
        "velocidade": 2,
        "fotometria": 2,
        "morfologia": 2,
        "faixa_vel": 0,
        "maximos": {"linearidade": 3, "velocidade": 2, "fotometria": 2,
                    "morfologia": 2, "faixa_vel": 1}
      },
      "flags": ["LINEARIDADE_BOA", "VELOCIDADE_CONSISTENTE",
                "BRILHO_ESTAVEL", "MORFOLOGIA_PONTUAL", "MPC_SEM_MATCH"],
      "trilha": [
        {
          "frame_index": 0,
          "timestamp_utc": "2019-08-28T10:18:00",
          "x": 835.3, "y": 1472.1,
          "ra_deg": 329.441200, "dec_deg": -12.198700,
          "flux": 3241.5, "snr": 14.2, "wcs_valido": true
        }
      ],
      "mpc": {
        "status": "sem_match",
        "n_objetos_no_cone": 0,
        "match": null
      }
    }
  ]
}
```

### `*_relatorio.txt` — Relatório operacional

Organizado para uso prático no fluxo do Astrometrica:

1. **Métricas do conjunto processado** — contexto da execução
2. **Priorização operacional** — candidatos agrupados por urgência de inspeção
3. **Detalhamento por candidato** — score decomposto, flags, trilha frame a frame, status MPC

### `*_MPC_report.txt` — Rascunho no formato MPC 80-colunas

Contém apenas candidatos com score ≥ 6 (MODERADO+). Magnitude deixada em branco — requer calibração fotométrica no Astrometrica (UCAC4/Gaia). **Não enviar ao IASC sem revisão e re-medição no Astrometrica.**

### `*_candidatos.png` — Suporte visual

Recortes dos 4 frames para cada candidato (top-5 não-DESCARTA), com:
- Círculo na posição atual
- Marcadores de posição anterior
- Vetor indicando direção de movimento
- Rank, score e flags resumidas

### `*_pipeline.log` — Log estruturado da execução

Registro cronológico de todas as etapas, incluindo contagem de fontes por frame, deriva estimada e resultados das consultas SkyBot.

---

## Score e Flags

### Score heurístico (0–12)

O score é composto por 7 componentes independentes. Os thresholds estão centralizados na classe `T` em `detector.py` para facilitar tuning futuro sem caçar valores espalhados no código.

| Componente              | Máx | Critério                                               |
|-------------------------|-----|--------------------------------------------------------|
| Linearidade             |  3  | Resíduo da regressão linear na trilha (px)             |
| Velocidade              |  2  | Uniformidade dos passos entre frames (σ px)            |
| Fotometria              |  2  | Coeficiente de variação do brilho integrado (4 frames) |
| Morfologia              |  2  | Pontualidade média do perfil (pico/fluxo total)        |
| Consist. morfológica    |  1  | Estabilidade da pontualidade entre frames              |
| Elongação               |  1  | Razão eixo maior/menor via PCA dos pixels              |
| Faixa velocidade (MBA)  |  1  | Deslocamento total na faixa típica de MBA (4–50 px)    |

| Score | Classe   | Significado operacional                                   |
|-------|----------|-----------------------------------------------------------|
| 8–12  | FORTE    | Candidato com múltiplas características de asteroide real |
| 6–7   | MODERADO | Candidato plausível, verificar no Astrometrica            |
| 4–5   | FRACO    | Baixa confiança, verificar se houver tempo                |
| 0–3   | DESCARTA | Provável artefato, ruído ou fonte estendida               |

Candidatos com **razão de rejeição** (borda, SNR insuficiente, hot pixel, brilho caótico, elongação grave) recebem classe DESCARTA independente do score calculado.

**O score é heurístico e ajustado para pré-triagem.** Não é uma probabilidade estatisticamente calibrada. Um FORTE pode ser uma estrela variável; um FRACO pode ser um asteroide real com trilha irregular em campo degradado. Toda decisão final deve ser tomada no Astrometrica.

### Rejeição de falsos positivos

Antes da classificação final, o pipeline avalia heurísticas de falso positivo. Se qualquer condição abaixo for satisfeita, o candidato recebe `classe=DESCARTA` e flag `REJEITADO_FP`, com a razão explicitada no JSON (`razoes_decisao.rejeicoes`) e no TXT:

| Condição                    | Critério                                                |
|-----------------------------|---------------------------------------------------------|
| Borda do frame              | Centroide a < 30 px da borda em qualquer frame         |
| SNR insuficiente            | SNR < 2.0 em ≥ 2 dos 4 frames                          |
| Hot pixel suspeito          | Pontualidade > 0.70 e SNR médio < 3.0                  |
| Brilho caótico              | CV do brilho > 1.5                                      |
| Elongação grave sistemática | Elongação média > 3.0 e máxima > 4.0 em todos os frames|

Penalizações de score (sem rejeição completa) são registradas em `razoes_decisao.penalizacoes` e logadas no terminal durante o processamento.

### Validação cinemática de trilhas

Antes de aceitar uma trilha formada pelo casamento de fontes entre frames, o pipeline verifica coerência cinemática:

| Critério                           | Limiar padrão          |
|------------------------------------|------------------------|
| Desvio de direção entre passos     | ≤ 45°                  |
| Razão entre maior e menor passo    | ≤ 2.5×                 |
| Aceleração entre passos            | ≤ 3.0 px/frame²        |

Trilhas que violam esses critérios são rejeitadas antes da análise de score. O número de trilhas rejeitadas é registrado em `metricas_globais.n_rejeitados_coerencia` no JSON.

### Flags diagnósticas

| Flag                     | Significado                                               |
|--------------------------|-----------------------------------------------------------|
| `LINEARIDADE_BOA`        | Resíduo linear < 1.5 px                                  |
| `LINEARIDADE_RUIM`       | Resíduo linear ≥ 2.5 px                                  |
| `VELOCIDADE_CONSISTENTE` | σ dos passos < 3                                         |
| `VELOCIDADE_IRREGULAR`   | σ dos passos ≥ 5                                         |
| `BRILHO_ESTAVEL`         | CV do brilho < 0.25                                      |
| `BRILHO_INSTAVEL`        | CV do brilho ≥ 0.50                                      |
| `MORFOLOGIA_PONTUAL`     | Pontualidade média > 0.12                                |
| `MORFOLOGIA_ESTENDIDA`   | Pontualidade média ≤ 0.06                                |
| `ELONGACAO_ALTA`         | Elongação média ≥ 1.6 (fonte possivelmente não estelar)  |
| `MORFO_INCONSISTENTE`    | std(pontualidade entre frames) ≥ 0.08                    |
| `EDGE_FRAME`             | Centroide a < 30 px da borda em algum frame              |
| `REJEITADO_FP`           | Candidato descartado por heurística de falso positivo    |
| `MPC_SEM_MATCH`          | SkyBot não retornou objetos no cone de 2'                |
| `MPC_MATCH_PROVAVEL`     | Um objeto conhecido encontrado no cone                   |
| `MPC_MATCH_AMBIGUO`      | Mais de um objeto no cone — verificar manualmente        |
| `MPC_CONSULTA_FALHOU`    | Erro na consulta (sem internet, astroquery, coords inv.) |

---

## Limitações

- **Sem calibração fotométrica absoluta.** O brilho medido é relativo (ADU). A magnitude no relatório MPC fica em branco e deve ser preenchida pelo Astrometrica após calibração via UCAC4 ou Gaia DR3.

- **WCS simplificado.** O pipeline ignora as correções polinomiais proprietárias do Pan-STARRS (chaves `PCA*`). O erro astrométrico resultante é pequeno na maior parte do campo, mas pode ser mais significativo nas bordas. Use as coordenadas do Astrometrica para a medição final.

- **Score heurístico, não probabilístico.** Os pesos e limiares foram ajustados para melhorar a pré-triagem com base em critérios físicos observacionais, mas não foram calibrados estatisticamente em conjuntos rotulados. Um score 10/12 não equivale a 83% de probabilidade de ser asteroide.

- **Heurísticas de FP podem rejeitar candidatos legítimos.** Em campos com PSF degradada, seeing ruim ou imagens com artefatos globais, as heurísticas de elongação ou SNR podem descartar objetos reais. Sempre revise candidatos com flag `REJEITADO_FP` se o campo suspeitar de problemas de qualidade.

- **Validação cinemática pode filtrar objetos com movimento não-uniforme.** NEAs rápidos ou objetos com movimento medido em campos com deriva grande podem ter trilhas que violam os thresholds de ângulo ou razão de passo. Os thresholds estão centralizados em `T.*` e podem ser ajustados.

- **SkyBot indica vizinhança, não identidade.** Um `MPC_MATCH_PROVAVEL` significa que existe objeto catalogado próximo à posição e época consultadas. A confirmação de que se trata do mesmo objeto requer comparação de trilha e velocidade no Astrometrica.

- **4 frames, campos esparsos.** O pipeline foi desenvolvido e testado com sequências de 4 frames do Pan-STARRS/IASC. Não foi testado com outros telescópios, escalas de pixel ou campos densos.

- **O que deve ser validado com dados reais:** (1) os thresholds cinemáticos em `T.MAX_ANGULO_DEG`, `T.MAX_STEP_RATIO` e `T.STEP_ACCEL_MAX` podem precisar de ajuste em campos com seeing variável; (2) `T.HOT_PIXEL_PONT_MAX` e `T.SNR_MINIMO` devem ser verificados em campos com CCD de ruído diferente do Pan-STARRS; (3) o componente de elongação pode penalizar estrelas binárias não resolvidas — monitorar com practice sets do IASC.

---

## Estrutura do Projeto

```
asteroid-hunter/
├── src/
│   └── detector.py          ← Pipeline principal
├── imagens/                 ← Coloque seus arquivos .fits aqui
├── resultados/              ← Outputs gerados automaticamente
├── testes/
│   └── teste_basico.py      ← Testes unitários e de integração
├── docs/
│   └── metodologia.md       ← Documentação técnica detalhada
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Metodologia

### Detecção de fontes

Usa `photutils.DAOStarFinder` com estatísticas de background calculadas por `astropy.stats.sigma_clipped_stats` (σ=3.0, 5 iterações). O threshold de detecção é configurável via `--sigma` (padrão 5.5σ).

### Rastreamento com validação cinemática

Casamento mútuo recíproco de vizinho mais próximo (raio 20 px) entre frames consecutivos. Um objeto só é rastreado se o casamento for recíproco nas três transições (1→2, 2→3, 3→4).

Antes de aceitar uma trilha, o pipeline verifica coerência cinemática:
- **Ângulo de direção:** desvio máximo entre passos consecutivos ≤ 45° (parâmetro `T.MAX_ANGULO_DEG`)
- **Razão de passo:** razão máxima/mínima dos passos ≤ 2.5× (`T.MAX_STEP_RATIO`)
- **Aceleração:** variação de passo entre frames ≤ 3.0 px/frame² (`T.STEP_ACCEL_MAX`)

Trilhas com saltos incoerentes são descartadas antes da análise de score.

### Remoção de deriva instrumental

A deriva entre frames é estimada como a mediana dos vetores de movimento de todas as fontes rastreadas. Fontes com resíduo < 2 px são usadas para refinar a estimativa. Objetos com resíduo > 1.8 px após subtração da deriva são candidatos a mover.

### Morfologia

Para cada candidato e cada frame, o pipeline extrai:
- **Pontualidade:** razão pico/fluxo total no recorte de 40×40 px
- **Elongação:** razão dos eixos via PCA ponderado dos pixels acima de 10% do pico
- **FWHM estimada:** raio onde a intensidade cai à metade do pico (em pixels)
- **Consistência entre frames:** desvio padrão da pontualidade nos 4 frames

### Score (0–12)

Sete componentes mensuráveis, com thresholds centralizados em `T.*`:
- **Linearidade (0–3):** resíduo médio da regressão linear nas posições
- **Velocidade (0–2):** soma dos desvios padrão dos passos em x e y
- **Fotometria (0–2):** coeficiente de variação do brilho integrado (apertura 20 px)
- **Morfologia (0–2):** pontualidade média do perfil
- **Consistência morfológica (0–1):** estabilidade da pontualidade entre frames
- **Elongação (0–1):** fonte compacta (razão eixos < 1.6) ganha ponto
- **Faixa velocidade (0–1):** deslocamento total entre 4 e 50 px (típico de MBA)

### Rejeição de falsos positivos

Avaliada antes da classificação final usando heurísticas físicas sobre os dados já disponíveis (borda, SNR, morfologia, fotometria). Candidatos rejeitados recebem `DESCARTA` com razão explícita no JSON e no TXT.

### Consulta SkyBot

`astroquery.imcce.Skybot.cone_search` com raio de 2' na época exata do frame 1. A ausência de correspondência não implica que o objeto seja desconhecido.

---

## Referências

- **Miller et al. (2024)** — *The International Astronomical Search Collaboration (IASC)*. PASP 136, 024502. [ResearchGate](https://www.researchgate.net/publication/378410609)
- **Ospina et al. (2018)** — *AIsteroid*. [GitHub](https://github.com/seap-udea/AIsteroid)
- **Minor Planet Center** — MPC Submission Guidelines. [minorplanetcenter.net](https://minorplanetcenter.net/iau/info/Astrometry.html)
- **Berthier et al.** — SkyBot / IMCCE. [ssd.imcce.fr/webservices/skybot](https://ssd.imcce.fr/webservices/skybot/)

### Stack técnica

| Biblioteca   | Versão mínima | Uso                                    |
|--------------|---------------|----------------------------------------|
| astropy      | 6.0           | FITS I/O, WCS, coordenadas, tempo      |
| photutils    | 1.13          | Detecção de fontes (DAOStarFinder)     |
| astroquery   | 0.4.7         | Consultas ao SkyBot/IMCCE              |
| numpy        | 1.26          | Operações matriciais                   |
| matplotlib   | 3.9           | Visualizações                          |
| pytest       | 8.0           | Testes                                 |

---

## Changelog

### v1.3.0 (2026-04)

**Foco: melhorar a pré-triagem — candidatos plausíveis mais no topo, menos artefatos entre os primeiros ranks, decisões mais auditáveis.**

- **Associação mais robusta:** validação cinemática de trilha (ângulo, razão de passo, aceleração) antes de aceitar a associação entre frames; trilhas incoerentes são rejeitadas cedo e contadas em `metricas_globais.n_rejeitados_coerencia`
- **Rejeição de falsos positivos:** borda agora rejeita (não apenas flag), + SNR mínimo, hot pixel heurístico, brilho caótico, elongação grave sistemática; flag `REJEITADO_FP` + razão explícita no JSON e TXT
- **Score recalibrado (0–12):** dois novos componentes — consistência morfológica entre frames e elongação; limiares de fotometria e morfologia ajustados; thresholds centralizados na classe `T`
- **Morfologia robusta:** elongação via PCA ponderado dos pixels, FWHM estimada por limiar de meia potência, consistência da pontualidade entre os 4 frames — integradas ao score e ao JSON
- **Logs de decisão explícitos:** campo `razoes_decisao` no JSON com listas separadas de penalizações e rejeições; logs no terminal informam motivo de cada penalização/rejeição/promoção
- **Flags novas:** `ELONGACAO_ALTA`, `MORFO_INCONSISTENTE`, `REJEITADO_FP`
- **Thresholds centralizados:** classe `T` em `detector.py` agrupa todos os limiares para tuning sem caçar valores no código
- **Testes novos:** `TestValidarCoerenciaTrilha`, `TestMorfologiaPorFrame`, `TestVerificarFalsoPositivo`, `TestThresholdsCentralizados`

### v1.2.0 (2026-04)

- Score decomposto por componente exportado no JSON (`score_componentes`) e detalhado no TXT
- Trilha completa frame a frame serializada no JSON com timestamp, coordenadas, fluxo e SNR por frame
- Flags diagnósticas por candidato derivadas do pipeline, exportadas no JSON e resumidas no TXT
- Status MPC normalizado com schema estruturado
- Métricas globais do conjunto: fontes por frame, trilhas tentadas, deriva, sigma, tempo de execução
- Relatório TXT reorganizado com seção de priorização operacional e diagnóstico técnico por candidato
- PNGs melhorados: rank, score, flags, vetor de trajetória
- SNR estimado por frame (apertura local simples)
- Velocidade angular em arcsec/min quando WCS disponível
- Flag `EDGE_FRAME` para candidatos próximos à borda

### v1.1.0 (2026-04)

- WCS Pan-STARRS construído manualmente (ignora chaves PCA* proprietárias)
- Detecção via photutils.DAOStarFinder (robusta em campos esparsos)
- Formatação sexagesimal com wraparound correto
- Magnitude omitida no MPC (sem calibração)
- Rastreamento com validação mútua
- Deriva estimada em fontes estáveis
- HDU da ciência detectada automaticamente
- Observador via CLI ou variável de ambiente
- Logging estruturado em arquivo

### v1.0.0 (2026-04)

- Versão inicial.

---

## Licença

MIT License — veja [LICENSE](LICENSE) para detalhes.

---

**Jaciana Barbosa**  
Cientista da Computação | Astrônoma Amadora  
Icapuí, Ceará, Brasil
