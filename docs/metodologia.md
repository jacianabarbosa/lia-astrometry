# Metodologia Técnica — asteroid-hunter

**Autora: Jaciana Barbosa**  
**Versão: 1.0.0**

---

## 1. Visão Geral do Pipeline

O asteroid-hunter implementa um pipeline de 8 estágios para detecção de objetos em movimento em imagens astronômicas do programa IASC/Pan-STARRS. Este documento descreve em detalhe cada etapa, adequado para uso como base metodológica em publicações científicas.

---

## 2. Dados de Entrada

### Formato das imagens
- **Formato**: FITS (Flexible Image Transport System)
- **Origem**: Telescópio Pan-STARRS 1 (1.8-m, Haleakalā, Havaí)
- **Código MPC**: F51
- **Cadência**: 4 frames, separados por ~19 minutos cada
- **Exposição**: 45 segundos por frame
- **Campo**: ~0.3° × 0.3°
- **Escala de placa**: ~0.26 arcsec/pixel

### Metadados utilizados do header FITS
```
DATE-OBS  →  Timestamp ISO 8601 da observação
EXPTIME   →  Tempo de exposição (segundos)
CRVAL1/2  →  Coordenadas RA/Dec do centro da imagem
CRPIX1/2  →  Pixel de referência
CD1_1 etc →  Matriz de transformação pixel→céu
NAXIS1/2  →  Dimensões da imagem
```

---

## 3. Detecção de Fontes

### 3.1 Remoção de Background

O background é estimado usando o **percentil 35** dos valores de pixel:

```python
background = np.percentile(imagem, 35)
imagem_limpa = imagem - background
```

O percentil 35 é escolhido empiricamente para imagens astronômicas tipicamente esparsas (maioria dos pixels = céu).

### 3.2 Estimativa de Ruído

```python
sigma = np.std(imagem_limpa[imagem_limpa > 0])
threshold = N × sigma    # N = 5.5 por padrão
```

### 3.3 Extração de Fontes

Pixeis acima do threshold são binarizados e agrupados em regiões conectadas (8-conectividade via `scipy.ndimage.label`). Filtros aplicados:

- **Tamanho mínimo**: 2 pixels (remove ruído de pixel único)
- **Tamanho máximo**: 200 pixels (remove galáxias estendidas e saturações)

Para cada região, extrai-se o centroide via `scipy.ndimage.center_of_mass`.

---

## 4. Rastreamento entre Frames

### 4.1 Associação de Vizinho Mais Próximo

Para cada fonte no Frame 1, busca-se a fonte mais próxima em cada frame subsequente com raio de busca de 50 pixels. Objetos encontrados em todos os 4 frames formam uma trilha candidata.

Filtro de movimento: apenas trilhas com deslocamento total entre 2 e 120 pixels são mantidas (corresponde a velocidades típicas de MBAs a essa escala).

### 4.2 Remoção de Deriva Instrumental

O campo pode sofrer pequenas translações entre frames devido a variações de apontamento. Para separar essa deriva sistemática do movimento real:

```
deriva_dx = mediana(dx de todos os candidatos)
deriva_dy = mediana(dy de todos os candidatos)

movimento_residual = sqrt((dx - deriva_dx)² + (dy - deriva_dy)²)
```

Candidatos com movimento residual > 1.8 pixels são considerados movedores reais.

---

## 5. Sistema de Score Multidimensional

### 5.1 Linearidade (0–3 pontos)

Ajuste de linha reta pelo método dos mínimos quadrados:

```
residuo_medio = media(distância_perpendicular_à_reta_ajustada)
```

| Resíduo         | Pontos |
|-----------------|--------|
| < 0.8 pixels    | 3      |
| 0.8 – 1.5 px    | 2      |
| 1.5 – 2.5 px    | 1      |
| > 2.5 pixels    | 0      |

### 5.2 Consistência de Velocidade (0–2 pontos)

```
consistencia = std(passos_x) + std(passos_y)
```

Passos são as diferenças de posição entre frames consecutivos. Baixo desvio padrão indica velocidade uniforme.

### 5.3 Estabilidade de Brilho (0–2 pontos)

Coeficiente de Variação (CV) do brilho de pico nos 4 frames:

```
CV = std(brilhos) / media(brilhos)
```

Objetos reais mantêm magnitude aproximadamente constante. Raios cósmicos e artefatos costumam aparecer em apenas 1–2 frames.

### 5.4 Morfologia Pontual (0–2 pontos)

Razão pico/total em recorte de 40×40 pixels ao redor do objeto:

```
pontual = max(recorte) / sum(recorte)
```

Fontes pontuais (estrelas, asteroides) têm razão alta. Galáxias e objetos estendidos têm razão baixa.

### 5.5 Velocidade Típica (0–1 ponto)

Asteroides do Cinturão Principal (MBA) a esta distância do telescópio movem tipicamente 4–50 pixels em 4 frames de 19 min:

```
movimento_total = sqrt(dx² + dy²)
```

---

## 6. Conversão Astrométrica

### 6.1 Transformação Pixel → RA/Dec

Usando o **World Coordinate System (WCS)** armazenado no header FITS:

```python
from astropy.wcs import WCS
wcs = WCS(header)
ra, dec = wcs.pixel_to_world(x, y)
```

As imagens Pan-STARRS já vêm com solução WCS pré-calculada pelo pipeline do telescópio, com precisão de ~0.1 arcsec.

### 6.2 Limitações

A precisão astrométrica deste pipeline depende da qualidade do WCS do header FITS. Para comparação: o Astrometrica realiza refinamento adicional usando estrelas de catálogo (UCAC4, PPMXL, Gaia), podendo atingir ~0.05 arcsec. A diferença típica é de 0.1–0.3 arcsec, o que é aceitável para observações preliminares do IASC.

---

## 7. Consulta ao Minor Planet Center

### 7.1 SkyBot (IMCCE)

O pipeline usa o serviço **SkyBot** via `astroquery.imcce.Skybot`:

```python
resultado = Skybot.cone_search(
    coordenadas,
    rad=2 * u.arcmin,
    epoch=Time(data_obs)
)
```

O SkyBot calcula as efemérides de todos os objetos conhecidos do sistema solar para a data/hora da observação e verifica quais estão dentro do raio de busca. Este é o mesmo princípio usado pelo Astrometrica internamente.

### 7.2 Interpretação dos resultados

- **Objeto encontrado**: provavelmente é o asteroide catalogado. Reportar como observação de confirmação.
- **Objeto não encontrado**: candidato a nova descoberta. Verificar no Astrometrica e reportar ao IASC.
- **Erro/timeout**: verificar manualmente. O SkyBot pode estar indisponível.

---

## 8. Geração do Relatório MPC

### 8.1 Formato MPC 80-colunas

```
Colunas  1– 5: Número do objeto (ou espaços)
Colunas  6–12: Designação provisional
Coluna    13 : Flag de descoberta
Coluna    14 : Nota 1
Coluna    15 : Nota 2
Colunas 16–32: Data (YYYY MM DD.ddddd)
Colunas 33–44: RA (HH MM SS.ss)
Colunas 45–56: Dec (±DD MM SS.s)
Colunas 57–65: Espaços
Colunas 66–70: Magnitude
Coluna    71 : Banda fotométrica
Colunas 72–77: Espaços
Colunas 78–80: Código do observatório
```

### 8.2 Conversões necessárias

**RA em graus → HH MM SS.ss**:
```
RA_horas = RA_graus / 15
HH = int(RA_horas)
MM = int((RA_horas - HH) × 60)
SS = ((RA_horas - HH) × 60 - MM) × 60
```

**Dec em graus → ±DD MM SS.s**:
```
sinal = '+' se Dec >= 0, '-' caso contrário
DD = int(|Dec|)
MM = int((|Dec| - DD) × 60)
SS = ((|Dec| - DD) × 60 - MM) × 60
```

**JD → Data MPC**:
```
Converter JD para datetime
frac_dia = (hora × 3600 + min × 60 + seg) / 86400
data_mpc = "YYYY MM DD.ddddd"
```

---

## 9. Limitações e Trabalho Futuro

### 9.1 Limitações atuais

1. **Precisão astrométrica**: ~0.1–0.3 arcsec (vs ~0.05 arcsec do Astrometrica com Gaia)
2. **Objetos muito fracos**: abaixo de 5σ podem ser perdidos
3. **Campos congestionados**: muitas fontes podem causar associações incorretas
4. **Raios cósmicos**: podem passar pelo filtro de brilho se aparecerem em múltiplos frames

### 9.2 Melhorias planejadas

- [ ] Calibração astrométrica refinada via catálogo Gaia DR3
- [ ] PSF fitting via `photutils` para medições mais precisas
- [ ] Interface gráfica simples para confirmação visual
- [ ] Suporte a conjuntos de mais de 4 frames
- [ ] Cálculo do digest2 score (índice NEO do MPC)

---

## 10. Validação

Para validar o pipeline, recomenda-se:

1. Rodar nas **imagens de prática do IASC** (que têm asteroides conhecidos)
2. Comparar as coordenadas geradas com as do Astrometrica
3. Verificar se os objetos FORTE/MODERADO aparecem no Astrometrica na mesma posição

---

*Documento gerado para: asteroid-hunter v1.0.0*  
*Autora: Jaciana Barbosa — 2026*
