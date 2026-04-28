# asteroid-hunter v1.4.2

Pipeline Python de pré-triagem de asteroides em imagens FITS do IASC/Pan-STARRS.

**Autora:** Jaciana Barbosa — Cientista da Computação & Astrônoma Amadora  
Icapuí, Ceará, Brasil  
Participante do IASC (International Astronomical Search Collaboration)

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

O **asteroid-hunter** é uma ferramenta de apoio para procurar possíveis asteroides em imagens astronômicas. Ele foi pensado para o fluxo do IASC, em que o participante recebe quatro imagens do mesmo pedaço do céu, tiradas em horários diferentes, e precisa procurar pequenos pontos que mudaram de posição.

O IASC (International Astronomical Search Collaboration) é um programa internacional de ciência cidadã, realizado em parceria com a NASA, em que estudantes, professores e astrônomos amadores analisam imagens reais de telescópios profissionais para identificar e medir asteroides. No fluxo usado aqui, as imagens vêm do Pan-STARRS, conjunto de telescópios localizado no Observatório Haleakalā, na ilha de Maui, Havaí. As observações validadas podem contribuir para o acompanhamento de objetos já conhecidos e para a triagem de novos candidatos.

Essas imagens vêm no formato FITS, que é o formato mais usado em astronomia porque guarda não só a imagem, mas também informações importantes no cabeçalho: horário da observação, escala da imagem, orientação do campo e uma primeira solução de coordenadas do céu. Cada imagem é chamada aqui de frame.

A ideia física é simples. As estrelas estão tão distantes que, dentro de uma sequência curta de observação, parecem ficar paradas umas em relação às outras. Um asteroide está muito mais perto, dentro do Sistema Solar, então aparece como um ponto fraco que muda de lugar entre um frame e outro. O trabalho do pipeline é procurar esses pontos móveis sem confundir com ruído, pixel ruim, estrela saturada, galáxia pequena ou defeito de imagem.

O programa faz uma primeira triagem automática: detecta fontes pontuais, acompanha essas fontes entre os quatro frames, mede movimento, brilho e formato, rejeita casos claramente ruins e organiza os candidatos para revisão no Astrometrica. Ele não decide sozinho o que é asteroide. Ele ajuda a pessoa a gastar menos tempo procurando no escuro.

Na versão atual, o pipeline não depende apenas da posição em pixels para acompanhar uma fonte entre os frames. Primeiro ele refina o WCS com Gaia DR3 quando possível; depois usa esse WCS para projetar as fontes em coordenadas celestes locais. Isso reduz erros quando há pequeno deslocamento de apontamento, rotação ou diferença astrométrica entre as imagens.

Além da detecção de movimento, o pipeline calcula coordenadas RA/Dec, consulta o SkyBot/IMCCE para verificar se existe objeto conhecido próximo, atribui um score multidimensional a cada candidato e gera relatórios estruturados para uso no Astrometrica.

A ferramenta foi construída para uso pessoal no contexto do IASC. Não é um produto comercial, não tem interface gráfica e não pretende substituir nenhuma etapa do processo científico formal.

---

## Posicionamento

**O que esta ferramenta faz:**
- Triagem rápida dos candidatos mais prováveis para inspeção no Astrometrica
- Geração de coordenadas RA/Dec para facilitar a localização no Astrometrica
- Rastreamento entre frames usando o WCS refinado, quando disponível
- Consulta de vizinhança posicional via SkyBot (objeto conhecido próximo?)
- Exportação de dados estruturados e auditáveis para revisão posterior

**O que esta ferramenta não faz:**
- Não substitui a validação astrométrica no software Astrometrica disponibilizado pelo IASC.
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
[1] Carregamento         Lê FITS, extrai WCS inicial, timestamps e qualidade
      │
      ▼
[2] Máscara de pixels    Saturação/inválidos: >=65500, <1 e não-finitos
      │
      ▼
[3] Detecção de fontes   Background2D + DAOStarFinder + PSF sub-pixel
      │
      ▼
[4] Refinamento WCS      Gaia DR3 + movimento próprio + ajuste afim
      │
      ▼
[5] Rastreamento         Fontes projetadas em plano celeste pelo WCS refinado
                         + casamento mútuo recíproco entre frames
      │
      ▼
[6] Rejeição de FP       Borda, SNR baixo, hot pixel, brilho caótico,
                         elongação grave e sistemática
      │
      ▼
[7] Score (0–12)         Linearidade + Velocidade + Fotometria + Morfologia
                         + Consistência morfológica + Elongação + Faixa vel.
      │
      ▼
[8] Conversão RA/Dec     WCS Pan-STARRS refinado por Gaia quando disponível
      │
      ▼
[9] Consulta SkyBot      Vizinhança posicional via astroquery/IMCCE
      │
      ▼
[10] Outputs             .json  .txt  _MPC_report.txt  .png  .log
```

---

## Instalação

**Pré-requisitos:** Python 3.10+, pip, conexão com internet para Gaia DR3 e SkyBot. Sem rede, o pipeline mantém fallback transparente para o WCS Pan-STARRS e marca o status da consulta.

```bash
# Clone ou extraia o repositório
git clone https://github.com/jacianabarbosa/asteroid-hunter.git
cd asteroid-hunter

# Ambiente virtual (recomendado)
python3 -m venv venv
source venv/bin/activate        # macOS/Linux
# venv\Scripts\activate         # Windows

# Dependências
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# Verificação
python -m pytest testes/teste_basico.py -v
```

> **macOS com pyenv:** use `pyenv local 3.13.1` antes de criar o venv se `python3` não encontrar a versão correta.

### Recriar o venv se a pasta mudou de nome

Se o projeto foi renomeado ou movido, o `pip` do ambiente virtual pode continuar apontando para o caminho antigo e falhar com erro parecido com `No such file or directory`. Nesse caso, recrie o `venv`:

```bash
deactivate 2>/dev/null
rm -rf venv
python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Use `python -m pip` em vez de `pip` para garantir que as dependências sejam instaladas no Python ativo do ambiente virtual.

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
- **Métricas globais do conjunto:** fontes por frame, trilhas tentadas, deriva diagnóstica do campo, saturação por frame, status do refinamento Gaia, sigma usado, tempo de execução, distribuição de classes
- **Por candidato:**
  - `rank`, `classe`, `score_total`, `probabilidade`
  - `score_componentes`: decomposição completa (linearidade/3, velocidade/2, fotometria/2, morfologia/2, consistência morfológica/1, elongação/1, faixa_vel/1)
  - `flags`: lista de flags diagnósticas derivadas do pipeline
  - `posicao_frame1`: pixel x/y, RA/Dec, formatos MPC
  - `movimento`: deslocamento total em pixels, resíduo de deriva, linearidade, consistência de velocidade, velocidade em arcsec/min
  - `fotometria`: brilho por frame, CV, pontualidade, SNR estimado por frame
  - `trilha`: array com posição, coordenadas, fluxo e SNR de cada frame
  - `mpc`: status normalizado da consulta SkyBot com campos estruturados quando há match

Exemplo de estrutura:

```json
{
  "pipeline": "asteroid-hunter v1.4.2",
  "conjunto": "o8723g0213o",
  "metricas_globais": {
    "n_fontes_por_frame": [42, 39, 41, 40],
    "n_trilhas_tentadas": 35,
    "deriva_dx_px": 0.62,
    "deriva_dy_px": 0.41,
    "sigma_deteccao": 5.5,
    "tempo_execucao_s": 12.3,
    "saturacao_pct_por_frame": [0.02, 0.01, 0.02, 0.01],
    "gaia_refinamento": [
      {"frame": 1, "status": "gaia_refinado", "n_matches": 18,
       "rms_pre_arcsec": 1.42, "rms_pos_arcsec": 0.31,
       "n_matches_bruto": 32, "offset_bruto_arcsec": [3.2, 0.2],
       "n_matches_fino": 30, "n_matches_ajuste": 30,
       "n_matches_usados": 30, "match_origem": "fino"}
    ]
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
        "consist_morfologica": 1,
        "elongation": 1,
        "faixa_vel": 0,
        "maximos": {"linearidade": 3, "velocidade": 2, "fotometria": 2,
                    "morfologia": 2, "consist_morfologica": 1,
                    "elongation": 1, "faixa_vel": 1},
        "score_max": 12
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

Registro cronológico de todas as etapas, incluindo contagem de fontes por frame, qualidade/saturação, status do refinamento Gaia, deriva diagnóstica e resultados das consultas SkyBot.

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

### Refinamento astrométrico Gaia

O WCS é a parte do FITS que diz qual ponto do céu corresponde a cada pixel da imagem. Se esse mapa estiver alguns segundos de arco fora do lugar, o candidato ainda pode aparecer visualmente, mas a posição RA/Dec fica imprecisa e o rastreamento entre frames fica mais vulnerável.

Gaia DR3 é a terceira liberação de dados da missão espacial Gaia, da Agência Espacial Europeia (ESA). Ela funciona como um mapa de referência muito preciso de posições e movimentos de estrelas, por isso é usada aqui para corrigir o WCS dos frames.

O refinamento opera em **duas passadas** para contornar erro inicial do WCS Pan-STARRS. A construção do WCS preserva `CROTA1/2` quando o header usa `CDELT`, e todas as diferenças de RA usam wrap 0/360 para campos próximos de RA=0°:

1. **Passada bruta (raio 15 arcsec):** estima o offset puro de translação via mediana dos resíduos. Corrige o CRVAL antes de tentar o ajuste afim.
2. **Passada fina (raio 2 arcsec):** cross-match preciso com WCS pré-corrigido + ajuste afim da CD matrix por mínimos quadrados.

Se a passada fina falhar (matches insuficientes ou RMS alto), o offset bruto é preservado no WCS com status `gaia_offset_bruto_apenas`. Quando a passada fina já deixa RMS abaixo de 1.0 arcsec, o status é `gaia_skipped_pre_baixo`; o WCS bruto corrigido é mantido porque já está dentro da tolerância operacional.

Se Gaia não estiver disponível ou não houver estrelas no campo, o WCS Pan-STARRS inicial é mantido. O status por frame fica em `metricas_globais.gaia_refinamento` com campos `n_matches_bruto`, `offset_bruto_arcsec`, `n_matches_fino`, `n_matches_ajuste`, `n_matches_usados` e `match_origem`.

### Rastreamento em coordenadas celestes

Depois que as fontes são detectadas em cada imagem, o pipeline transforma os centroides de pixel em posições no céu usando o WCS de cada frame. Essas posições são projetadas em um plano local, em segundos de arco, e só então ocorre o casamento entre frames.

Na prática, o raio continua sendo controlado por `T.RAIO_MATCH_PX`, porque é mais intuitivo ajustar em pixels. Internamente, esse raio é convertido para segundos de arco pela escala real do WCS. O casamento continua recíproco: uma fonte do frame 1 precisa escolher uma fonte do frame 2, e essa fonte do frame 2 precisa escolher a mesma fonte de volta.

Esse ponto é importante. Antes, duas fontes eram comparadas diretamente por distância em pixel. Agora, o rastreamento acompanha a posição no céu. Isso é mais compatível com o que se mede em astrometria, principalmente quando os frames têm pequenos offsets ou quando o header FITS tem rotação.

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

- **WCS com fallback.** O pipeline usa Pan-STARRS como solução inicial e refina com Gaia DR3 quando possível. Se Gaia estiver indisponível ou o ajuste não atingir o RMS exigido, o WCS inicial é preservado. Use as coordenadas medidas no Astrometrica para a submissão final.

- **Rastreamento depende da qualidade do WCS.** Quando Gaia refina bem o campo, o casamento celeste tende a ser mais estável que o casamento bruto em pixel. Se Gaia falhar, o pipeline ainda usa o WCS Pan-STARRS inicial; nesse caso, a qualidade do header FITS passa a pesar mais.

- **Score heurístico, não probabilístico.** Os pesos e limiares foram ajustados para melhorar a pré-triagem com base em critérios físicos observacionais, mas não foram calibrados estatisticamente em conjuntos rotulados. Um score 10/12 não equivale a 83% de probabilidade de ser asteroide.

- **Heurísticas de FP podem rejeitar candidatos legítimos.** Em campos com PSF degradada, seeing ruim ou imagens com artefatos globais, as heurísticas de elongação ou SNR podem descartar objetos reais. Sempre revise candidatos com flag `REJEITADO_FP` se o campo suspeitar de problemas de qualidade.

- **Dependência de rede para catálogos.** O refinamento Gaia DR3 e a consulta SkyBot dependem de serviços externos. O pipeline tem fallback gracioso, mas a execução offline reduz a informação astrométrica/catalógica disponível.

- **SkyBot indica vizinhança, não identidade.** Um `MPC_MATCH_PROVAVEL` significa que existe objeto catalogado próximo à posição e época consultadas. A confirmação de que se trata do mesmo objeto requer comparação de trilha e velocidade no Astrometrica.

- **4 frames, campos esparsos.** O pipeline foi desenvolvido e testado com sequências de 4 frames do Pan-STARRS/IASC. Não foi testado com outros telescópios, escalas de pixel ou campos densos.

- **O que deve ser validado com dados reais:** (1) thresholds de Gaia em `T.GAIA_*` podem precisar de ajuste em campos pobres ou congestionados; (2) `T.HOT_PIXEL_PONT_MAX` e `T.SNR_MINIMO` devem ser verificados em campos com CCD de ruído diferente do Pan-STARRS; (3) o componente de elongação pode penalizar estrelas binárias não resolvidas — monitorar com practice sets do IASC.

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

O detector começa estimando o fundo do céu com `photutils.Background2D`. Isso é melhor que um fundo global único, porque imagens reais têm gradientes, halos, bordas ruidosas e variação local de ruído.

Depois, o `DAOStarFinder` entra como gerador de sementes: ele encontra pontos acima do limiar local de SNR. Cada semente passa por ajuste PSF sub-pixel, com modelo Moffat por padrão e Gaussiana como fallback. O centro final usado pelo pipeline vem desse ajuste, não apenas do centroide inicial do DAO.

Pixels saturados, zerados ou não-finitos são mascarados antes do cálculo de fundo, da detecção e do ajuste PSF. O threshold de detecção é configurável via `--sigma` (padrão 5.5σ).

### Rastreamento

O rastreamento usa casamento mútuo recíproco de vizinho mais próximo entre frames consecutivos, mas a distância não é mais calculada diretamente no plano de pixels.

O fluxo atual é:

1. detecta fontes em cada frame;
2. converte o centroide de cada fonte para RA/Dec usando o WCS do próprio frame;
3. projeta RA/Dec em um plano celeste local em segundos de arco;
4. faz o casamento mútuo recíproco nesse plano;
5. exige que a trilha exista nas três transições: 1→2, 2→3 e 3→4.

O raio padrão continua sendo `T.RAIO_MATCH_PX = 20`, mas é convertido internamente para segundos de arco usando a escala média do WCS. Campos próximos de RA=0° usam média circular de RA e diferença angular com wrap 0/360.

Na v1.4, a validação cinemática antiga por ângulo/razão/aceleração foi removida. A coerência da trilha passou a ser avaliada depois, no score de linearidade, velocidade, fotometria e morfologia. A causa principal dos falsos saltos entre frames, que era astrométrica, passou a ser atacada antes, com o WCS refinado.

### Deriva diagnóstica

A deriva mediana entre frames ainda é estimada a partir das trilhas casadas e registrada nas métricas globais. Ela ajuda a separar deslocamento comum do campo de movimento real do candidato. Na prática, estrelas tendem a compartilhar a mesma deriva; um asteroide deve sobrar como resíduo.

Essa etapa não substitui o WCS. Ela é um diagnóstico e um filtro residual em pixels, útil para manter compatibilidade com o restante do score e com a visualização.

### Refinamento WCS via Gaia DR3

Para cada frame, o pipeline consulta Gaia DR3 via TAP/ESA, filtra estrelas por magnitude e RUWE, propaga movimento próprio desde a época 2016.0 e faz cross-match com as fontes detectadas. Um ajuste afim com sigma-clipping pode reconstruir a matriz CD do WCS quando o RMS inicial ainda exige refinamento. Se o offset bruto já deixa RMS abaixo de 1.0 arcsec, o WCS corrigido é mantido sem ajuste afim adicional.

Quando o refinamento Gaia não está disponível, o pipeline continua rodando com o WCS Pan-STARRS do header. O relatório deixa esse estado explícito; não é uma falha silenciosa.

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
| astropy-healpix | 1.0        | Suporte astrométrico/catálogos         |
| numpy        | 1.26          | Operações matriciais                   |
| matplotlib   | 3.9           | Visualizações                          |
| pytest       | 8.0           | Testes                                 |

---

## Changelog

### Atualização de rastreamento celeste (2026-04)

**Foco: usar o WCS refinado também no rastreamento entre frames.**

- O casamento entre fontes deixou de depender só da distância em pixels.
- As fontes detectadas são convertidas para RA/Dec pelo WCS de cada frame.
- As posições celestes são projetadas em um plano local em segundos de arco.
- O raio `T.RAIO_MATCH_PX` continua sendo configurado em pixels, mas agora é convertido pela escala média do WCS.
- Campos próximos de RA=0° usam média circular de RA e diferença angular com wrap 0/360.
- Se o WCS falhar, o pipeline preserva fallback para o casamento em pixels.
- A suíte local segue passando com 72 testes.

### v1.4.2 (2026-04)

**Foco: corrigir a astrometria Gaia em campos Pan-STARRS com `CROTA` e RA perto de 0°.**

- **WCS Pan-STARRS corrigido:** `_construir_wcs_panstarrs` agora preserva `CROTA1/2` quando o header usa `CDELT`; isso evita perder a rotação de ~180° dos FITS `o8730g0300o`
- **Offset bruto preserva rotação:** o WCS intermediário usado após a passada bruta copia `CROTA`, evitando que a passada fina volte para uma solução sem rotação
- **Wrap correto de RA:** cross-match, RMS e ajuste afim usam a menor diferença angular em RA, essencial para campos próximos de 0/360°
- **Passada bruta mais robusta:** votação 2D de translação com pareamento 1:1 no pico de offset, em vez de vizinho mais próximo simples em campo denso
- **Ajuste afim com sigma-clipping:** pares incoerentes podem ser removidos antes de avaliar o RMS pós-ajuste
- **Logs completos:** o arquivo `.log` agora é inicializado antes do carregamento/refinamento Gaia
- **Métricas expandidas:** `gaia_refinamento` inclui `n_matches_ajuste`, `n_matches_usados` e `match_origem`
- **Resultado validado em `o8730g0300o`:** 29-31 matches finos por frame, RMS Gaia de 0.659-0.744 arcsec e status `gaia_skipped_pre_baixo`
- **Testes:** cobertura para `CROTA`, RA wrap em Gaia e índices 1:1 da passada bruta; suíte com 72 testes

### v1.4.1 (2026-04)

**Foco: corrigir o cross-match Gaia que falhava porque o erro inicial do WCS Pan-STARRS (~5–7 arcsec) era maior que o raio de busca (2 arcsec).**

- **Cross-match em duas passadas:** passada bruta (raio 15 arcsec) estima offset de translação via mediana robusta dos resíduos; passada fina (raio 2 arcsec) executa o ajuste afim completo com WCS pré-corrigido
- **Preservação do offset bruto:** mesmo que a passada fina falhe por matches insuficientes ou RMS alto, o CRVAL corrigido é mantido no WCS com status `gaia_offset_bruto_apenas`
- **Nova função `_estimar_offset_grosseiro`:** isolada para testabilidade, robusta a ~20% de falsos matches por mediana
- **Novos thresholds:** `GAIA_RAIO_BRUTO_ARCSEC = 15.0` e `GAIA_MIN_MATCHES_BRUTO = 5` na classe `T`
- **Métricas expandidas:** `gaia_refinamento` por frame inclui `n_matches_bruto`, `offset_bruto_arcsec` e `n_matches_fino`
- **Log detalhado:** offset bruto e N de matches de cada passada logados por frame durante a execução
- **4 novos testes:** `test_offset_grosseiro_basico`, `test_offset_grosseiro_outliers`, `test_refinar_wcs_offset_bruto_apenas`, `test_refinar_wcs_duas_passadas_completo` — total 69 testes, todos passando

### v1.4.0 (2026-04)

**Foco: corrigir a causa raiz astrométrica com refinamento absoluto por Gaia DR3.**

- **Refinamento WCS por frame:** consulta Gaia DR3 via TAP/ESA, filtro G 12–19 e `RUWE < 1.4`, propagação de movimento próprio desde epoch 2016.0 e cross-match por KDTree
- **Ajuste afim do WCS:** reconstrução da CD matrix por mínimos quadrados; refinamento aceito apenas com RMS pós-ajuste < 0.5 arcsec
- **Fallback transparente:** falha de rede, poucos matches ou RMS alto preservam o WCS Pan-STARRS inicial e registram status por frame
- **Limpeza de workarounds:** removida a validação cinemática por ângulo/razão/aceleração e removida a lógica de deriva por par de frames
- **Métricas novas:** `metricas_globais.gaia_refinamento` com status, matches, RMS pré e RMS pós por frame
- **Testes Gaia:** cobertura para movimento próprio, cross-match, fallback sem rede, skip por RMS baixo, aceite de refinamento e mock de TAP
- **Dependências:** `astropy-healpix>=1.0` adicionado ao `requirements.txt`

### v1.3.2 (2026-04)

**Foco: workaround temporário para offsets de pointing/WCS entre frames.**

- Deriva estimada por par de frames (1→2, 2→3, 3→4) como mediana dos vetores entre estrelas casadas
- Reordenação das operações para estimar deriva antes da validação cinemática
- Validação passou a receber trilha corrigida pela deriva acumulada, preservando a trilha original para RA/Dec e MPC
- Logs explícitos dos três pares, com aviso quando `|deriva| > 3 px`
- Thresholds `DERIVA_PAR_*` e testes `TestDerivaPorPar` adicionados nesta versão

### v1.3.1 (2026-04)

**Foco: impedir que pixels saturados fossem detectados como fontes reais.**

- Nova máscara de pixels inválidos para valores `>= 65500`, `< 1` e não-finitos
- `Background2D`, imagem de detecção e `DAOStarFinder` passaram a receber/considerar a máscara
- Ajuste PSF sub-pixel descarta fonte quando a janela local tem mais de 30% de pixels mascarados
- `carregar_fits` calcula `sat_pct` por frame, avisa acima de 5% e aborta acima de 25%
- Novos thresholds `PIXEL_SAT_*` e `FRAME_SAT_*` na classe `T`
- Testes `TestMascaraSaturacao` adicionados

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
