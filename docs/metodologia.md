# Metodologia Técnica — asteroid-hunter

**Autora: Jaciana Barbosa**  
**Versão: 1.4.2**

---

## 1. Visão Geral do Pipeline

O asteroid-hunter implementa um pipeline de pré-triagem para objetos em movimento em sequências FITS do IASC/Pan-STARRS. O IASC (International Astronomical Search Collaboration) é um programa internacional de ciência cidadã, realizado em parceria com a NASA, em que estudantes, professores e astrônomos amadores analisam imagens reais de telescópios profissionais para identificar e medir asteroides. No fluxo tratado por este projeto, as imagens vêm do Pan-STARRS, conjunto de telescópios localizado no Observatório Haleakalā, na ilha de Maui, Havaí.

A entrada típica são quatro imagens do mesmo campo, feitas em horários diferentes. Estrelas e galáxias de fundo devem permanecer praticamente fixas; um objeto do Sistema Solar aparece como uma fonte pontual que se desloca entre os frames.

A versão 1.4.2 usa detecção fotométrica local, mascaramento de pixels inválidos, WCS Pan-STARRS com `CROTA` preservado, refinamento astrométrico em duas passadas por Gaia DR3 quando possível, rastreamento em coordenadas celestes projetadas pelo WCS de cada frame, rejeição de falsos positivos e exportação de relatórios auditáveis.

O pipeline é auxiliar: ele prioriza candidatos para inspeção humana no Astrometrica, mas não substitui a medição astrométrica final nem a submissão revisada ao IASC/MPC.

---

## 2. Dados de Entrada

### Formato das imagens

- **Formato**: FITS (Flexible Image Transport System)
- **Origem típica**: Pan-STARRS / IASC
- **Local do telescópio**: Observatório Haleakalā, Maui, Havaí
- **Código MPC**: F51
- **Cadência esperada**: 4 frames do mesmo campo
- **Metadados essenciais**: timestamp, WCS inicial e dimensões da imagem

### Metadados usados do header FITS

```text
DATE-OBS       Timestamp ISO 8601 da observação
EXPTIME        Tempo de exposição
CRVAL1/2       Coordenadas RA/Dec do centro
CRPIX1/2       Pixel de referência
CROTA1/2       Rotação do WCS quando o header usa CDELT
CD* ou CDELT*  Matriz/escala da transformação pixel-céu
NAXIS1/2       Dimensões da imagem
```

O WCS Pan-STARRS é usado como solução inicial. Na v1.4.2, cada frame pode ser refinado de forma independente contra Gaia DR3 em duas passadas. Quando o FITS traz `CDELT` + `CROTA`, a rotação é preservada; ignorar `CROTA` gera uma solução espelhada/rotacionada e inviabiliza o cross-match Gaia.

---

## 3. Qualidade de Frame e Máscara de Saturação

Antes da estimativa de background e da detecção, o pipeline constrói uma máscara de pixels inválidos:

- pixels `>= 65500`, incluindo a saturação Pan-STARRS em `65535`;
- pixels `< 1`;
- valores não finitos (`NaN`, `inf`).

A máscara é usada em três pontos:

- `Background2D`, para impedir que pixels saturados contaminem a estimativa local de céu;
- imagem de detecção, zerando pixels mascarados antes do `DAOStarFinder`;
- ajuste PSF sub-pixel, descartando fontes cuja janela local tenha mais de 30% de pixels mascarados.

O carregamento calcula `sat_pct` por frame. Frames com mais de 5% de pixels inválidos geram aviso de qualidade; acima de 25%, o pipeline aborta para evitar que uma imagem degradada contamine todos os resultados.

---

## 4. Detecção de Fontes

### 4.1 Background local

O background é estimado por `photutils.Background2D` com malha local e filtro mediano:

```python
Background2D(data, box_size=(50, 50), filter_size=(3, 3), mask=mascara)
```

Essa abordagem é mais robusta que um percentil global porque preserva variações suaves de fundo e evita que regiões saturadas inflem a estatística local.

### 4.2 DAOStarFinder e PSF

As sementes pontuais são detectadas com `photutils.DAOStarFinder`, usando threshold configurável por `--sigma` e FWHM inicial de 3 px. Em seguida, o pipeline tenta refinar o centroide por ajuste PSF sub-pixel com modelo Moffat ou Gaussiano.

Fontes com FWHM incompatível com PSF estelar, sharpness/roundness fora da faixa, baixa SNR ou janela contaminada por pixels inválidos são rejeitadas cedo.

---

## 5. Refinamento Astrométrico Gaia DR3

### 5.1 Consulta e movimento próprio

Gaia DR3 é a terceira liberação de dados da missão espacial Gaia, da Agência Espacial Europeia (ESA). O catálogo fornece posições, movimentos próprios e magnitudes de estrelas com alta precisão, funcionando como uma régua astrométrica para corrigir a relação entre pixels da imagem e coordenadas reais no céu.

Para cada frame, o pipeline consulta Gaia DR3 por TAP/ESA em torno do centro WCS inicial. O filtro padrão usa estrelas com magnitude G entre 12 e 19 e `RUWE < 1.4`, evitando estrelas saturadas e soluções astrométricas ruins.

Como Gaia DR3 está referenciado à época 2016.0, as posições são propagadas para a data da observação com movimento próprio:

```python
ra_corr  = ra  + pmra  * dt_anos / cos(dec)
dec_corr = dec + pmdec * dt_anos
```

### 5.2 Cross-match em duas passadas

O WCS Pan-STARRS pode ter erro sistemático inicial maior que a tolerância de 2 arcsec necessária para o ajuste afim. Além disso, campos próximos de RA=0° exigem cálculo de diferença angular com wrap 0/360. A solução é um cross-match em duas etapas:

**Passada bruta (raio 15 arcsec):**

```python
off_ra, off_dec, n = _estimar_offset_grosseiro(det_radec, gaia_radec)
```

Estima o offset de translação pura por votação 2D dos pares detecção-Gaia dentro do raio bruto, com pareamento 1:1 no pico de translação. O CRVAL é corrigido por esse offset antes da passada fina, preservando `CROTA` no WCS intermediário.

**Passada fina (raio 2 arcsec):**

Com o WCS pré-corrigido, as detecções são reprojetadas e casadas com Gaia dentro do raio apertado. Com matches suficientes, o pipeline pode resolver uma transformação afim por mínimos quadrados com sigma-clipping e reconstruir a CD matrix.

O refinamento afim é aceito apenas se o RMS pós-ajuste ficar abaixo de `0.5 arcsec`. Se a passada fina já deixa RMS abaixo de `1.0 arcsec`, o pipeline usa `gaia_skipped_pre_baixo` e mantém o WCS bruto corrigido. Se a passada fina falhar, o offset bruto é preservado com status `gaia_offset_bruto_apenas`.

### 5.3 Status e fallback

O refinamento é sempre gracioso. O status por frame é registrado em `metricas_globais.gaia_refinamento`:

| Status                    | Significado                                              |
|---------------------------|----------------------------------------------------------|
| `gaia_refinado`           | Ajuste afim completo aceito (RMS < 0.5 arcsec)          |
| `gaia_offset_bruto_apenas`| Só translação aplicada; passada fina insuficiente        |
| `gaia_skipped_pre_baixo`  | WCS Pan-STARRS já era bom (RMS < 1.0 arcsec)            |
| `gaia_falhou_rede`        | Gaia não respondeu                                       |
| `gaia_falhou_match`       | Estrelas insuficientes mesmo na passada bruta            |
| `gaia_falhou_pos_alto`    | Ajuste afim não convergiu                                |
| `wcs_invalido`            | Frame sem WCS válido                                     |

Campos por frame nas métricas: `n_matches_bruto`, `offset_bruto_arcsec [dRA, dDec]`, `n_matches_fino`, `n_matches_ajuste`, `n_matches_usados`, `match_origem`, `rms_pre_arcsec`, `rms_pos_arcsec`.

---

## 6. Rastreamento entre Frames

O rastreamento usa casamento mútuo recíproco de vizinho mais próximo entre frames consecutivos, mas a distância principal não é mais calculada diretamente no plano de pixels.

O fluxo atual é:

1. fontes pontuais são detectadas em cada frame;
2. cada centroide `(x, y)` é convertido para RA/Dec usando o WCS daquele frame;
3. as coordenadas RA/Dec são projetadas para um plano tangente local em segundos de arco;
4. o casamento mútuo recíproco é feito nesse plano celeste;
5. uma trilha só é aceita quando a associação existe nas transições 1→2, 2→3 e 3→4.

O raio operacional continua sendo `T.RAIO_MATCH_PX = 20`, para manter o ajuste intuitivo em pixels. Internamente, esse raio é convertido para segundos de arco pela escala média do WCS de cada par de frames. Para campos perto de RA=0°, a referência de RA usa média circular e as diferenças angulares usam wrap 0/360.

Esse refinamento é importante porque dois frames podem ter pequenos deslocamentos de apontamento, rotação ou erro inicial de WCS. Em pixel bruto, esses efeitos podem parecer movimento falso. Em coordenadas celestes, a associação fica mais próxima da medida astrométrica real.

A versão 1.4.0 removeu a validação cinemática antiga por ângulo, razão de passo e aceleração, porque ela rejeitava trilhas reais quando havia offsets sistemáticos de WCS/pointing entre frames. A correção agora ocorre no domínio astrométrico, por refinamento absoluto Gaia DR3, em vez de tentar compensar o problema apenas em pixels.

Se o WCS de um frame falhar, o pipeline preserva fallback para o casamento em pixels. Esse fallback mantém a execução possível, mas a qualidade astrométrica do header FITS passa a pesar mais.

O pipeline ainda estima a deriva mediana de fontes estáveis para separar movimento residual de campo e manter métricas diagnósticas, mas a lógica de "deriva por par de frames" da v1.3.2 foi removida como workaround obsoleto.

---

## 7. Rejeição de Falsos Positivos

Antes da classificação final, o pipeline aplica filtros físicos e diagnósticos:

| Condição                    | Critério padrão                                      |
|-----------------------------|------------------------------------------------------|
| Borda do frame              | Centroide a menos de 30 px da borda                  |
| SNR insuficiente            | SNR < 2.0 em pelo menos 2 dos 4 frames               |
| Hot pixel suspeito          | Pontualidade muito alta com SNR médio baixo          |
| Brilho caótico              | Coeficiente de variação do brilho > 1.5              |
| Elongação grave sistemática | Elongação média > 3.0 e máxima > 4.0                 |
| FWHM incompatível           | Perfil muito estreito ou muito largo para PSF estelar |

Candidatos rejeitados recebem `classe=DESCARTA`, flag `REJEITADO_FP` e razão explícita em `razoes_decisao.rejeicoes`.

---

## 8. Sistema de Score

O score é heurístico e vai de 0 a 12 pontos.

| Componente              | Máx | Critério                                               |
|-------------------------|-----|--------------------------------------------------------|
| Linearidade             | 3   | Resíduo da regressão linear da trilha                  |
| Velocidade              | 2   | Uniformidade dos passos entre frames                   |
| Fotometria              | 2   | Coeficiente de variação do brilho integrado            |
| Morfologia              | 2   | Pontualidade média do perfil                           |
| Consistência morfológica| 1   | Estabilidade da pontualidade entre frames              |
| Elongação               | 1   | Fonte compacta por razão eixo maior/menor              |
| Faixa de velocidade MBA | 1   | Deslocamento total entre 4 e 50 px                     |

| Score | Classe   | Interpretação operacional                     |
|-------|----------|-----------------------------------------------|
| 8–12  | FORTE    | Prioridade alta para inspeção no Astrometrica |
| 6–7   | MODERADO | Candidato plausível                           |
| 4–5   | FRACO    | Baixa confiança                               |
| 0–3   | DESCARTA | Provável artefato ou ruído                    |

Os thresholds ficam centralizados na classe `T` em `src/detector.py`.

---

## 9. Conversão Astrométrica e SkyBot

As posições de pixel são convertidas para RA/Dec usando o WCS disponível no frame: refinado por Gaia DR3 quando o ajuste foi aceito, ou Pan-STARRS inicial quando o fallback foi acionado.

O pipeline consulta o SkyBot/IMCCE com raio de 2 arcmin na época do primeiro frame:

```python
Skybot.cone_search(coordenadas, rad=2 * u.arcmin, epoch=Time(data_obs))
```

Interpretação:

- **Match provável**: existe objeto conhecido próximo; confirmar no Astrometrica.
- **Sem match**: não implica descoberta; exige revisão humana e comparação no Astrometrica.
- **Erro/timeout**: consulta indisponível; verificar manualmente.

---

## 10. Outputs

Para cada conjunto processado, o pipeline gera:

- `*_candidatos.json`: candidatos, trilhas, flags, scores, métricas globais e status Gaia/SkyBot;
- `*_relatorio.txt`: relatório operacional para inspeção;
- `*_MPC_report.txt`: rascunho MPC 80-colunas para candidatos MODERADO+;
- `*_candidatos.png`: recortes visuais dos principais candidatos;
- `*_pipeline.log`: log cronológico da execução.

O relatório MPC deixa magnitude em branco porque o pipeline não faz calibração fotométrica absoluta.

---

## 11. Validação

A suíte de testes da v1.4.2 cobre 72 casos, incluindo:

- máscara de saturação e pixels inválidos;
- propagação de movimento próprio Gaia;
- cross-match Gaia básico, com outliers e com RA perto de 0/360°;
- preservação de `CROTA` no WCS Pan-STARRS;
- estimativa de offset grosseiro (passada bruta);
- fallback sem rede;
- comportamento quando passada fina tem matches insuficientes (`gaia_offset_bruto_apenas`);
- aceitação/rejeição do refinamento afim por RMS;
- integração end-to-end com Gaia mockada.

Execução:

```bash
python -m pytest testes/teste_basico.py -v
```

---

## 12. Limitações

- A conexão com Gaia DR3 e SkyBot depende de rede; ambos têm fallback, mas a precisão/identificação pode ficar limitada.
- Não há calibração fotométrica absoluta; magnitude final deve ser medida no Astrometrica.
- O score não é uma probabilidade estatística calibrada.
- O pipeline foi desenvolvido para sequências IASC/Pan-STARRS de 4 frames; outros telescópios ou cadências podem exigir novos thresholds.
- A confirmação científica continua dependendo de revisão humana e medição astrométrica final.

---

*Documento atualizado para: asteroid-hunter v1.4.2*  
*Autora: Jaciana Barbosa — 2026*
