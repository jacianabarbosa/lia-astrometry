[Read in English](methodology.md)

# Metodologia da Lia

**Versão:** 1.5.0  
**Repositório:** `lia-astrometry`  
**Autora:** Jaciana Barbosa

A Lia é um pipeline de pré-triagem de candidatos a objetos em movimento em sequências de imagens FITS do IASC/Pan-STARRS. Ela foi desenvolvida para reduzir o esforço de inspeção manual antes da revisão no Astrometrica. **Não substitui** o Astrometrica, a validação pelo MPC nem o fluxo oficial do IASC.

## Etapas do Pipeline

1. **Carregamento dos FITS**  
   A Lia lê o HDU de ciência de cada arquivo FITS e converte os dados da imagem para `float64`. Isso evita perda prematura de precisão numérica durante a subtração do fundo, o ajuste PSF, o refinamento de centroide e a conversão pelo WCS.

2. **Extração de metadados**  
   O pipeline lê as palavras-chave do WCS, as dimensões da imagem, o tempo de exposição e os metadados de tempo. `DATE-OBS` é preferido quando disponível; `MJD-OBS` é aceito como fallback.

3. **Timestamp e tempo de meio da exposição**  
   Para astrometria de objetos em movimento, o timestamp cientificamente relevante é o instante médio da exposição. A Lia armazena tanto o início quanto o meio da exposição em JD e MJD. O tempo de meio de exposição (MJD) é usado nos relatórios astrométricos.

4. **Estimativa local do fundo do céu**  
   O fundo é estimado com `photutils.Background2D` usando uma malha local e filtragem por mediana. Essa abordagem é mais robusta do que um único threshold global em campos com gradientes de iluminação, halos de estrelas brilhantes, estrutura do detector ou pixels mascarados.

5. **Detecção de fontes**  
   As fontes são detectadas na imagem subtraída do fundo e normalizada pelo SNR local com `DAOStarFinder`. Pixels saturados, inválidos e não finitos são mascarados antes da detecção.

6. **Refinamento PSF e centroide sub-pixel**  
   Cada semente de fonte é refinada usando um modelo PSF Moffat2D. Um modelo Gaussiano2D é usado como fallback quando o ajuste Moffat falha. O FWHM ajustado é usado para rejeitar fontes muito estreitas para o seeing atmosférico, muito largas para uma fonte pontual, ou prováveis blends, artefatos ou detecções do tipo raio cósmico.

7. **Refinamento astrométrico Gaia DR3**  
   No modo `--wcs-mode gaia`, a Lia consulta o Gaia DR3 ao redor do centro WCS inicial, propaga os movimentos próprios das estrelas de referência até a época da observação, realiza uma correspondência cruzada grosseira de translação e em seguida uma correspondência fina com atualização WCS afim opcional. O JSON registra o status do Gaia, os resíduos RMS, os offsets e a contagem de matches para cada frame.

8. **Fallback WCS apenas pelo header**  
   No modo `--wcs-mode header`, o refinamento Gaia é ignorado e o WCS construído a partir do header é usado diretamente. O restante do pipeline permanece idêntico. Esse modo existe para comparações internas controladas, não como um fluxo separado para envio ao IASC.

9. **Construção de trilhas**  
   As fontes detectadas são projetadas em uma representação de plano tangente local e vinculadas nos quatro frames usando correspondência recíproca de vizinho mais próximo. Uma trilha candidata deve estar presente em todos os quatro frames.

10. **Validação cinemática**  
    A Lia ajusta uma trajetória linear às quatro posições medidas e registra o resíduo médio e o `R^2`. Sequências curtas do IASC devem apresentar movimento aproximadamente linear; resíduos irregulares são evidência de artefatos, blends ou vinculação incorreta de fontes.

11. **Rejeição de falsos positivos**  
    Os candidatos são rejeitados ou penalizados com base em distância à borda, SNR, FWHM, elongação, morfologia de fonte pontual, variação de fluxo, heurísticas de hot pixel e verificações de fonte estática no Gaia.

12. **Score heurístico**  
    O score combina linearidade, consistência de velocidade, estabilidade fotométrica, morfologia, consistência morfológica, elongação e faixa de velocidade. É um score de triagem, não um valor de confiança estatisticamente calibrado.

13. **Verificação de vizinhança SkyBot/MPC**  
    Quando habilitada, a Lia consulta o SkyBot/IMCCE ao redor das coordenadas do candidato e da época da observação para identificar objetos conhecidos do Sistema Solar próximos. Um match no SkyBot é informação de contexto, não validação final.

14. **Geração de saídas**  
    A Lia grava JSON, texto, rascunho MPC, recortes PNG e logs. O JSON inclui `run_metadata`, `global_metrics`, rastreabilidade por candidato, campos de validação manual e proveniência tanto do modo Gaia quanto do modo header.

## Glossário

| Termo | Significado |
|---|---|
| FITS | Flexible Image Transport System — o formato padrão de imagens astronômicas usado pelos conjuntos IASC/Pan-STARRS. |
| IASC | International Astronomical Search Collaboration — o fluxo de campanha em que participantes inspecionam conjuntos de imagens FITS e enviam medições validadas pelo processo oficial. |
| Astrometrica | O software usado no fluxo do IASC para inspeção visual, medição astrométrica, calibração fotométrica e preparação de relatórios. |
| Relatório MPC | Um relatório de observação para o Minor Planet Center. A Lia pode gerar texto auxiliar de rascunho, mas o relatório oficial deve ser revisado ou produzido no Astrometrica. |
| WCS | World Coordinate System — a transformação entre coordenadas de pixel e coordenadas celestes (RA/Dec). |
| Gaia DR3 | O terceiro lançamento de dados do Gaia, usado aqui como catálogo de referência astrométrica. |
| Resíduos WCS / RMS | Diferenças entre as posições das estrelas de referência do Gaia e as posições preditas pelo WCS ajustado, geralmente resumidas como RMS em segundos de arco. |
| Background2D | Um modelo local de fundo do céu do `photutils` que estima o fundo sobre uma malha em vez de assumir um campo uniforme. |
| Ajuste PSF | Ajuste de uma função de espalhamento pontual analítica, como Moffat ou Gaussiana, para estimar um centroide sub-pixel e a largura da fonte. |
| Centroide | O centro estimado de uma fonte detectada em coordenadas de pixel. A Lia preserva valores sub-pixel internamente. |
| MJD / tempo de meio da exposição | Modified Julian Date no ponto médio de uma exposição. As posições de objetos em movimento devem estar associadas ao tempo de meio da exposição. |
| SkyBot | Serviço do IMCCE para consulta de objetos conhecidos do Sistema Solar em uma dada posição e instante. |
| Candidato a objeto em movimento | Uma trilha de fonte que parece se mover de forma coerente nos quatro frames e requer validação humana. |
| Score heurístico de triagem | Um score de classificação transparente baseado em critérios físicos e instrumentais. Não é um valor de confiança estatisticamente calibrado. |

## Justificativa Metodológica

A estimativa local do fundo reduz falsos positivos em imagens não uniformes porque uma fonte é avaliada em relação ao ruído e ao nível do céu próximo à sua própria posição, e não em relação a uma estatística global. Isso é importante em campos do Pan-STARRS onde halos, estrutura do detector e pixels mascarados podem alterar o threshold local de detecção.

O centroide sub-pixel baseado em PSF melhora as posições preliminares ao ajustar um modelo contínuo de fonte à distribuição local de brilho. O mesmo ajuste fornece uma estimativa de FWHM, que ajuda a distinguir fontes astronômicas pontuais de hot pixels, raios cósmicos, artefatos saturados, trails e blends.

O refinamento Gaia DR3 reduz offsets sistemáticos do WCS quando há estrelas de referência limpas suficientes disponíveis. A Lia registra o número de matches e os resíduos RMS para que cada frame refinado possa ser auditado. Se o Gaia falhar ou não melhorar o campo, a Lia preserva e documenta o fallback para o WCS do header.

A validação cinemática é necessária porque as sequências do IASC são curtas. Um objeto real do Sistema Solar deve mostrar movimento aproximadamente linear em quatro frames, enquanto picos de ruído, estrelas associadas incorretamente e artefatos de borda frequentemente produzem passos inconsistentes ou resíduos lineares ruins.

A validação humana continua necessária. A Lia não pode determinar uma órbita, não pode calibrar a fotometria final e não pode substituir o julgamento necessário para inspecionar blends, fontes com baixo SNR e modos de falha conhecidos no Astrometrica.

## Integridade dos Dados

A Lia converte arrays de imagens FITS para `float64` no carregamento e mantém as coordenadas sub-pixel das fontes como valores de ponto flutuante durante a subtração do fundo, o ajuste PSF, a construção de trilhas, a conversão WCS e a exportação JSON. O arredondamento é reservado para logs legíveis, relatórios e campos de exibição serializados finais.

Pixels inválidos são mascarados antes da estimativa do fundo e da detecção de fontes. A máscara cobre pixels saturados, pixels inválidos próximos de zero e valores não finitos. Isso impede que regiões saturadas inflacionem o RMS local do fundo e evita que pixels inválidos se tornem detecções artificiais.

## Avaliação Gaia vs Apenas Header

O modo operacional principal é:

```bash
python src/detector.py --wcs-mode gaia
```

O modo de comparação interna é:

```bash
python src/detector.py --wcs-mode header
```

O modo header desativa apenas o refinamento Gaia DR3. A estimativa do fundo, o ajuste PSF, o threshold de detecção, o rastreamento, o score, os filtros e as saídas permanecem idênticos. Isso permite uma comparação controlada da recuperação de candidatos, mudanças de ranking, diferenças de RA/Dec, separação angular em relação às medições do Astrometrica, RMS do Gaia e contagem de matches.

Apenas o relatório final gerado pelo Astrometrica deve ser enviado ao IASC. A comparação Gaia/header é um estudo de validação interno.

## Plano de Validação

### Avaliação Operacional

```text
conjunto FITS do IASC
-> Lia com Gaia DR3
-> lista priorizada de candidatos
-> inspeção manual no Astrometrica
-> relatório MPC gerado pelo Astrometrica
-> envio ao IASC
-> retorno operacional do IASC, quando disponível
```

Métricas a coletar:

- redução do espaço de busca de candidatos;
- utilidade dos candidatos top-1/top-3/top-5;
- candidatos mensuráveis no Astrometrica;
- candidatos incluídos em um relatório MPC gerado pelo Astrometrica;
- retorno do IASC quando disponível;
- falsos positivos e modos de falha.

### Comparação Interna Gaia vs Header

Execute cada conjunto com `--wcs-mode gaia` e `--wcs-mode header`, depois exporte ambas as saídas JSON através do `tools/export_metrics.py`. Compare mudanças de ranking, RMS do Gaia, contagem de matches, diferenças de coordenadas e separação angular em relação às medições do Astrometrica.

## Limitações

- A Lia não substitui o Astrometrica.
- O score é heurístico, não é um valor de confiança estatisticamente calibrado.
- O Gaia DR3 pode falhar ou não melhorar todos os campos.
- O WCS do header pode já ser suficiente em alguns casos.
- O ajuste PSF pode falhar para fontes com baixo SNR, em blend, saturadas, com trail ou na borda do frame.
- Os parâmetros do `Background2D` precisam de validação empírica em conjuntos de prática e campanha do IASC.
- A rejeição de falsos positivos pode descartar candidatos reais em condições degradadas.
- O retorno do IASC é validação operacional, não verdade absoluta universal.
- A Lia não implementa determinação de órbita, score digest2, shift-and-stack, classificação por CNN, fluxos com interface gráfica nem envio automático ao MPC/IASC.
