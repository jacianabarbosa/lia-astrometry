# TRITON — Tracking and Robust Identification of Transient Objects in the Night

**Pipeline científico de pré-triagem de asteroides em imagens FITS do IASC/Pan-STARRS**

**Rastreamento e Identificação Robusta de Objetos Transientes na Noite**

**Autora:** Jaciana Barbosa — Cientista da Computação & Astrônoma Amadora  
Icapuí, Ceará, Brasil  
Participante do IASC (International Astronomical Search Collaboration)

> *TRITON foi criado para rastrear pontos fracos e transitórios no céu noturno. Em uma sequência FITS, a maior parte do campo é estática; o objeto de interesse é pequeno, discreto e aparece pelo movimento coerente entre frames.*

---

## Sumário

- [Problema e Solução](#problema-e-solução)
- [Diferenciais Técnicos](#diferenciais-técnicos)
- [Workflow](#workflow)
- [Como Funciona](#como-funciona)
- [Instalação](#instalação)
- [Uso](#uso)
- [Outputs Gerados](#outputs-gerados)
- [Score e Flags](#score-e-flags)
- [Limitações](#limitações)
- [Estrutura do Projeto](#estrutura-do-projeto)
- [Referências](#referências)
- [Licença](#licença)

---

## Problema e Solução

O IASC (International Astronomical Search Collaboration) é um programa internacional de ciência cidadã, realizado em parceria com a NASA, em que estudantes, professores e astrônomos amadores analisam imagens reais de telescópios profissionais para identificar e medir asteroides. No fluxo usado aqui, as imagens vêm do Pan-STARRS, no Observatório Haleakalā, em Maui - Havaí.

O IASC distribui conjuntos com sequências de quatro imagens FITS do mesmo campo para que participantes identifiquem asteroides visualmente. Cada imagem contém muitas fontes pontuais. Um asteroide é um ponto que se move coerentemente entre os quatro frames — mas boa parte do que parece se mover pode ser artefato: pixels saturados, raios cósmicos, hot pixels, halos de estrelas brilhantes, ruído de leitura ou associação errada entre fontes.

**O problema central é de fadiga cognitiva e escala.** Inspecionar manualmente cada ponto suspeito num conjunto de quatro imagens de 2400×2400 pixels não é viável. É preciso uma triagem automatizada que apresente ao observador apenas os candidatos com maior coerência física — e faça isso de forma auditável, registrando o motivo de cada decisão.

O TRITON resolve isso com um pipeline de dez etapas que parte da imagem crua e entrega uma lista ranqueada de candidatos com coordenadas RA/Dec prontas para revisão no software Astrometrica, disponibilizado pelo IASC para fazer os relatórios. Em execuções com conjuntos reais do IASC, o pipeline reduz dezenas de fontes detectadas por frame a uma lista pequena de trilhas candidatas, registrando por que cada uma foi promovida, penalizada ou descartada.

---

## Diferenciais Técnicos

### Astrometria referenciada a Gaia DR3

Usar o WCS do header FITS diretamente é simples, mas pode deixar erro sistemático de alguns segundos de arco. Em imagens Pan-STARRS, isso equivale a vários pixels e pode comprometer tanto a posição reportada quanto o rastreamento entre frames.

O TRITON refina o WCS de cada frame independentemente contra o catálogo Gaia DR3 da ESA, usando um cross-match em duas passadas (bruta 15 arcsec + fina 2 arcsec) com propagação de movimento próprio das estrelas de referência. Em testes com o conjunto `o8730g0300o`, a passada bruta identificou offset sistemático de +3.2 arcsec em RA em todos os quatro frames, corrigindo a astrometria antes do rastreamento.

### Rastreamento em coordenadas celestes

Após o refinamento WCS, o casamento entre fontes é feito no plano celeste (arcsec), não em pixels. Isso é fundamental quando frames consecutivos têm pequenos offsets de apontamento ou rotação — situação comum no Pan-STARRS.

### Detecção robusta a saturação

Frames Pan-STARRS podem trazer regiões saturadas em 65535 (limite uint16). Sem mascaramento, esses pixels contaminam a estimativa local de fundo e podem gerar fontes espúrias nas bordas das regiões saturadas. O TRITON mascara pixels ≥ 65500, < 1 e não-finitos antes de qualquer etapa de detecção.

### Score multidimensional auditável

Cada candidato recebe um score de 0 a 12 com sete componentes independentes. Toda decisão — penalização, rejeição ou promoção — é registrada no JSON e no log com a razão explícita. Nenhuma rejeição é silenciosa.

### Fallback gracioso em cadeia

Gaia indisponível → mantém WCS Pan-STARRS e registra status. WCS inválido → fallback para casamento em pixels. SkyBot offline → marca campo sem consulta. Em frames com saturação excessiva, o pipeline aborta cedo para evitar resultado enganoso.

---

## Workflow

O TRITON se encaixa entre o recebimento dos FITS do IASC e a revisão final no Astrometrica:

```
┌─────────────────────────────────────────────────────────────────┐
│  IASC distribui conjunto de 4 FITS  (e-mail / portal)           │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌───────────────────────────────────────────────────────────────────────┐
│  TRITON  (este pipeline)                                              │
│                                                                       │
│  venv/bin/python src/detector.py --imagens fits/ --output resultados/ │
│                                                                       │
│  Saída:                                                               │
│    candidatos.json   ← lista ranqueada, scores, coordenadas           │
│    relatorio.txt     ← priorização operacional                        │
│    MPC_report.txt    ← rascunho 80 colunas para candidatos ≥ 6        │
│    candidatos.png    ← recortes visuais dos top candidatos            │
│    pipeline.log      ← auditoria completa das decisões                │
└─────────────────────────────┬─────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Observador inspeciona candidatos FORTE e MODERADO no           │
│  Astrometrica, re-mede posições e preenche magnitudes           │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│  Submissão ao IASC/MPC após validação humana                    │
└─────────────────────────────────────────────────────────────────┘
```

**Tempo típico por conjunto:** 15–45 segundos (depende da latência da consulta Gaia DR3 e SkyBot).

---

## Como Funciona

```
4 arquivos .fits
      │
      ▼
[1] Carregamento         Lê FITS, extrai WCS (com CROTA preservado),
                         timestamps e fração de pixels saturados
      │
      ▼
[2] Máscara de pixels    Saturação ≥ 65500, valores < 1, não-finitos
      │
      ▼
[3] Detecção de fontes   Background2D (malha local 50×50 px) +
                         DAOStarFinder (threshold --sigma) +
                         ajuste PSF sub-pixel (Moffat / Gaussiana)
      │
      ▼
[4] Refinamento WCS      Gaia DR3: passada bruta (15 arcsec) →
                         corrige CRVAL → passada fina (2 arcsec) →
                         aceita se RMS < 0.9 arcsec
      │
      ▼
[5] Rastreamento         Fontes convertidas para RA/Dec pelo WCS
                         de cada frame → casamento mútuo recíproco
                         no plano celeste (arcsec)
      │
      ▼
[6] Rejeição de FP       Borda, SNR baixo, hot pixel suspeito,
                         brilho caótico, elongação grave
      │
      ▼
[7] Score (0–12)         Linearidade · Velocidade · Fotometria ·
                         Morfologia · Consistência · Elongação ·
                         Faixa de velocidade MBA
      │
      ▼
[8] Conversão RA/Dec     WCS refinado por Gaia quando disponível
      │
      ▼
[9] Consulta SkyBot      Vizinhança posicional via IMCCE (raio 2')
      │
      ▼
[10] Outputs             .json  .txt  _MPC_report.txt  .png  .log
```

---

## Instalação

**Pré-requisitos:** Python 3.10+, pip, conexão com internet para Gaia DR3 e SkyBot.

Sem rede, o pipeline mantém fallback transparente para o WCS Pan-STARRS e marca o status da consulta.

```bash
# Clone ou extraia o repositório
git clone https://github.com/jacianabarbosa/triton.git
cd triton

# Ambiente virtual (recomendado)
python3 -m venv venv
source venv/bin/activate        # macOS/Linux

# Dependências
venv/bin/python -m pip install --upgrade pip
venv/bin/python -m pip install -r requirements.txt

# Verificação
venv/bin/python -m pytest testes/teste_basico.py -v
```

> **macOS com pyenv:** use `pyenv local 3.13.1` antes de criar o venv se `python3` não encontrar a versão correta.

### Recriar o venv se a pasta mudou de nome

Se o projeto foi renomeado ou movido, recrie o venv:

```bash
deactivate 2>/dev/null
rm -rf venv
python3 -m venv venv
source venv/bin/activate
venv/bin/python -m pip install --upgrade pip
venv/bin/python -m pip install -r requirements.txt
```

Use `venv/bin/python -m pip` em vez de `pip` direto no macOS para garantir que as dependências sejam instaladas no Python do projeto.

---

## Uso

```bash
# Modo básico (4 arquivos .fits na pasta imagens/)
venv/bin/python src/detector.py

# Pasta personalizada
venv/bin/python src/detector.py --imagens /meus/fits/ --output /meus/resultados/

# Ajustar sensibilidade (padrão: 5.5σ)
venv/bin/python src/detector.py --sigma 5.0   # mais sensível — mais falsos positivos
venv/bin/python src/detector.py --sigma 6.5   # menos sensível — mais conservador

# Modo offline (sem consulta ao SkyBot)
venv/bin/python src/detector.py --sem-mpc

# Identificação do observador (necessário antes de enviar ao IASC)
venv/bin/python src/detector.py --observador "Jaciana Barbosa" --email "seu@email.com"

# Ou via variável de ambiente (evita repetir a cada execução)
export TRITON_OBS="Jaciana Barbosa"
export TRITON_EMAIL="seu@email.com"
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

O JSON é o output principal. Contém:

- **Cabeçalho global:** pipeline, conjunto, data de processamento, observatório
- **Métricas globais:** fontes por frame, trilhas tentadas, deriva diagnóstica, saturação por frame, status Gaia por frame, sigma usado, tempo de execução, distribuição de classes
- **Por candidato:**
  - `rank`, `classe`, `score_total`, `probabilidade`
  - `score_componentes`: decomposição completa com máximos por componente
  - `flags`: lista de flags diagnósticas
  - `posicao_frame1`: pixel x/y, RA/Dec, formatos MPC
  - `movimento`: deslocamento total, resíduo de deriva, linearidade, velocidade em arcsec/min
  - `fotometria`: brilho por frame, CV, pontualidade, SNR por frame
  - `trilha`: posição, coordenadas, fluxo e SNR frame a frame
  - `mpc`: status da consulta SkyBot com campos estruturados
  - `razoes_decisao`: penalizações e rejeições com texto explícito

Exemplo de estrutura:

```json
{
  "pipeline": "TRITON v1.4.2",
  "conjunto": "o8723g0213o",
  "metricas_globais": {
    "n_fontes_por_frame": [42, 39, 41, 40],
    "n_trilhas_tentadas": 35,
    "sigma_deteccao": 5.5,
    "tempo_execucao_s": 18.4,
    "saturacao_pct_por_frame": [0.02, 0.01, 0.02, 0.01],
    "gaia_refinamento": [
      {
        "frame": "o8723g0213o.fits",
        "status": "gaia_refinado",
        "n_matches_bruto": 29,
        "offset_bruto_arcsec": [3.2, 0.4],
        "n_matches_fino": 29,
        "rms_pre_arcsec": 0.664,
        "rms_pos_arcsec": 0.664
      }
    ]
  },
  "candidatos": [
    {
      "rank": 1,
      "classe": "FORTE",
      "score_total": 9,
      "score_componentes": {
        "linearidade": 3, "velocidade": 2, "fotometria": 2,
        "morfologia": 2, "consist_morfologica": 0,
        "elongation": 0, "faixa_vel": 0,
        "score_max": 12
      },
      "flags": ["LINEARIDADE_BOA", "VELOCIDADE_CONSISTENTE",
                "BRILHO_ESTAVEL", "MORFOLOGIA_PONTUAL", "MPC_SEM_MATCH"],
      "trilha": [
        {
          "frame_index": 0,
          "timestamp_utc": "2019-08-28T10:18:00",
          "x": 835.3, "y": 1472.1,
          "ra_deg": 329.4412, "dec_deg": -12.1987,
          "flux": 3241.5, "snr": 8.2, "wcs_valido": true
        }
      ],
      "mpc": { "status": "sem_match", "n_objetos_no_cone": 0, "match": null }
    }
  ]
}
```

### `*_relatorio.txt` — Relatório operacional

Organizado para uso prático no fluxo do Astrometrica:

1. **Métricas do conjunto** — contexto da execução, status Gaia por frame
2. **Priorização operacional** — candidatos agrupados por urgência de inspeção
3. **Detalhamento por candidato** — score decomposto, flags, trilha frame a frame, status MPC, razões de decisão

### `*_MPC_report.txt` — Rascunho no formato MPC 80-colunas

Contém apenas candidatos com score ≥ 6 (MODERADO+). Magnitude deixada em branco — requer calibração fotométrica no Astrometrica via UCAC4 ou Gaia DR3. **Não enviar ao IASC sem revisão e re-medição no Astrometrica.**

### `*_candidatos.png` — Suporte visual

Recortes dos 4 frames para cada candidato (top-5 não-DESCARTA), com:
- Círculo na posição atual
- Marcadores de posição anterior
- Vetor indicando direção de movimento
- Rank, score e flags resumidas

### `*_pipeline.log` — Auditoria completa

Registro cronológico de todas as etapas: contagem de fontes, saturação por frame, status e RMS do refinamento Gaia, offset bruto por frame, deriva diagnóstica, resultados SkyBot, razões de cada rejeição.

---

## Score e Flags

### Heurística de Coerência Cinemática (0–12)

O score quantifica em que medida o candidato se comporta como um corpo do Sistema Solar em movimento retilíneo uniforme. É composto por 7 componentes independentes com thresholds centralizados na classe `T` em `detector.py`.

| Componente              | Máx | O que mede                                             |
|-------------------------|-----|--------------------------------------------------------|
| Linearidade             |  3  | Resíduo da regressão linear na trilha (px)             |
| Velocidade              |  2  | Uniformidade dos passos entre frames (σ px)            |
| Fotometria              |  2  | Coeficiente de variação do brilho integrado            |
| Morfologia              |  2  | Pontualidade do perfil (pico / fluxo total)            |
| Consist. morfológica    |  1  | Estabilidade da pontualidade entre os 4 frames         |
| Elongação               |  1  | Razão eixo maior/menor via PCA dos pixels              |
| Faixa velocidade (MBA)  |  1  | Deslocamento total na faixa típica de Main Belt (4–50 px) |

| Score | Classe   | Significado operacional                                   |
|-------|----------|-----------------------------------------------------------|
| 8–12  | FORTE    | Candidato com múltiplas características de asteroide real |
| 6–7   | MODERADO | Candidato plausível, verificar no Astrometrica            |
| 4–5   | FRACO    | Baixa confiança, verificar se houver tempo                |
| 0–3   | DESCARTA | Provável artefato, ruído ou fonte estendida               |

Candidatos com **razão de rejeição** recebem `DESCARTA` independente do score calculado.

**O score é heurístico, não uma probabilidade calibrada.** Um FORTE pode ser uma estrela variável; um FRACO pode ser um asteroide real em campo degradado. Toda decisão final deve ser tomada no Astrometrica.

### Rejeição de falsos positivos

| Condição                    | Critério                                                |
|-----------------------------|---------------------------------------------------------|
| Borda do frame              | Centroide a < 30 px da borda em qualquer frame         |
| SNR insuficiente            | SNR < 2.0 em ≥ 2 dos 4 frames                          |
| Hot pixel suspeito          | Pontualidade > 0.70 e SNR médio < 3.0                  |
| Brilho caótico              | CV do brilho > 1.5                                      |
| Elongação grave sistemática | Elongação média > 3.0 e máxima > 4.0 em todos os frames|

### Refinamento WCS Gaia DR3

| Status                      | Significado                                            |
|-----------------------------|--------------------------------------------------------|
| `gaia_refinado`             | Translação bruta aceita (RMS < 0.9 arcsec)            |
| `gaia_offset_bruto_apenas`  | Translação aplicada; passada fina insuficiente         |
| `gaia_skipped_pre_baixo`    | WCS Pan-STARRS já dentro da tolerância                |
| `gaia_falhou_rede`          | Gaia não respondeu                                     |
| `gaia_falhou_match`         | Estrelas insuficientes no campo                        |
| `wcs_invalido`              | Frame sem WCS válido                                   |

### Flags diagnósticas

| Flag                     | Significado                                               |
|--------------------------|-----------------------------------------------------------|
| `LINEARIDADE_BOA`        | Resíduo linear < 1.5 px                                  |
| `LINEARIDADE_RUIM`       | Resíduo linear ≥ 2.5 px                                  |
| `VELOCIDADE_CONSISTENTE` | σ dos passos < 3 px                                      |
| `VELOCIDADE_IRREGULAR`   | σ dos passos ≥ 5 px                                      |
| `BRILHO_ESTAVEL`         | CV do brilho < 0.25                                      |
| `BRILHO_INSTAVEL`        | CV do brilho ≥ 0.50                                      |
| `MORFOLOGIA_PONTUAL`     | Pontualidade média > 0.12                                |
| `MORFOLOGIA_ESTENDIDA`   | Pontualidade média ≤ 0.06                                |
| `ELONGACAO_ALTA`         | Elongação média ≥ 1.6                                    |
| `MORFO_INCONSISTENTE`    | std(pontualidade entre frames) ≥ 0.08                    |
| `EDGE_FRAME`             | Centroide a < 30 px da borda em algum frame              |
| `REJEITADO_FP`           | Candidato descartado por heurística de falso positivo    |
| `MPC_SEM_MATCH`          | SkyBot não retornou objetos no cone de 2'                |
| `MPC_MATCH_PROVAVEL`     | Um objeto conhecido encontrado no cone                   |
| `MPC_MATCH_AMBIGUO`      | Mais de um objeto no cone — verificar manualmente        |
| `MPC_CONSULTA_FALHOU`    | Erro na consulta (sem rede, timeout, coordenadas inválidas) |

---

## Limitações

- **Sem calibração fotométrica absoluta.** Magnitude no relatório MPC fica em branco. Preencher no Astrometrica via UCAC4 ou Gaia DR3 antes de qualquer envio.

- **WCS Pan-STARRS com fallback.** O pipeline refina com Gaia DR3 quando possível. Se Gaia estiver indisponível ou o RMS não atingir a tolerância, o WCS do header é preservado. Use sempre as coordenadas re-medidas no Astrometrica para submissão.

- **Score heurístico, não probabilístico.** Os thresholds foram ajustados com base em critérios físicos observacionais mas não foram calibrados estatisticamente em conjuntos rotulados.

- **Heurísticas de FP podem rejeitar candidatos legítimos.** Em campos com PSF degradada, seeing ruim ou artefatos globais, elongação ou SNR podem descartar objetos reais. Revise sempre candidatos `REJEITADO_FP` em campos suspeitos.

- **Dependência de rede.** Refinamento Gaia DR3 e consulta SkyBot dependem de serviços externos (ESA TAP e IMCCE). O pipeline tem fallback gracioso, mas a execução offline reduz informação astrométrica e catalógica.

- **4 frames, campos Pan-STARRS.** Desenvolvido e testado para sequências IASC/Pan-STARRS. Outros telescópios, escalas de pixel ou cadências podem exigir ajuste de thresholds na classe `T`.

---

## Estrutura do Projeto

```
triton/
├── src/
│   └── detector.py              ← Pipeline principal (TRITON)
├── imagens/                     ← Coloque seus arquivos .fits aqui
├── resultados/                  ← Outputs gerados automaticamente
├── testes/
│   └── teste_basico.py          ← 72 testes unitários e de integração
├── docs/
│   └── metodologia.md           ← Documentação técnica detalhada
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Referências

- **Gaia Collaboration (2023)** — *Gaia Data Release 3*. A&A 674, A1. [doi:10.1051/0004-6361/202243940](https://doi.org/10.1051/0004-6361/202243940)
- **Miller et al. (2024)** — *The International Astronomical Search Collaboration (IASC)*. PASP 136, 024502. [ResearchGate](https://www.researchgate.net/publication/378410609)
- **Berthier et al.** — SkyBot / IMCCE. [ssd.imcce.fr/webservices/skybot](https://ssd.imcce.fr/webservices/skybot/)
- **Minor Planet Center** — MPC Submission Guidelines. [minorplanetcenter.net](https://minorplanetcenter.net/iau/info/Astrometry.html)
- **Bradley et al. (2024)** — *photutils: Astronomical source detection and photometry*. Zenodo. [doi:10.5281/zenodo.596036](https://doi.org/10.5281/zenodo.596036)
- **Astropy Collaboration (2022)** — *The Astropy Project*. ApJ 935, 167. [doi:10.3847/1538-4357/ac7c74](https://doi.org/10.3847/1538-4357/ac7c74)

### Stack técnica

| Biblioteca      | Versão mínima | Uso                                    |
|-----------------|---------------|----------------------------------------|
| astropy         | 6.0           | FITS I/O, WCS, coordenadas, tempo      |
| photutils       | 1.13          | Detecção de fontes (DAOStarFinder)     |
| astroquery      | 0.4.7         | Gaia TAP, SkyBot/IMCCE                 |
| astropy-healpix | 1.0           | Suporte astrométrico/catálogos         |
| numpy           | 1.26          | Operações matriciais                   |
| scipy           | 1.13          | KDTree para cross-match                |
| matplotlib      | 3.9           | Visualizações                          |
| pytest          | 8.0           | Testes                                 |

---

## Changelog

### v1.4.2 (2026-04)

- Correção do WCS Pan-STARRS com `CROTA` preservado
- Passada bruta com votação 2D de translação; wrap correto de RA
- Ajuste afim com sigma-clipping; aceite direto do WCS bruto quando RMS < 0.9 arcsec
- 72 testes; status `gaia_refinado` validado em `o8730g0300o` com RMS 0.66–0.85 arcsec

### v1.4.1 (2026-04)

- Cross-match Gaia em duas passadas para contornar erro inicial de 5–7 arcsec do WCS Pan-STARRS
- Nova função `_estimar_offset_grosseiro` com mediana robusta a outliers

### v1.4.0 (2026-04)

- Refinamento astrométrico por Gaia DR3 por frame
- Remoção da heurística cinemática em pixels (substituída por astrometria absoluta)

### v1.3.x (2026-04)

- v1.3.2: estimativa de deriva por par de frames (workaround, removido em v1.4)
- v1.3.1: mascaramento de pixels saturados antes de Background2D e DAOStarFinder
- v1.3.0: score multidimensional, morfologia por PCA, rejeição de falsos positivos

---

## Licença

MIT License — veja [LICENSE](LICENSE) para detalhes.

---

**Jaciana Barbosa**  
Cientista da Computação | Astrônoma Amadora  
Icapuí, Ceará, Brasil
