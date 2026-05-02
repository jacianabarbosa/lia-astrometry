[Read in English](README.md)

# Lia

**Lia** é uma ferramenta em Python que ajuda estudantes, professores e equipes de ciência cidadã a fazer uma pré-triagem de candidatos a objetos em movimento em conjuntos FITS do IASC antes da revisão manual no Astrometrica.

Ela analisa uma sequência de quatro imagens, detecta fontes pontuais, procura objetos que se movem de forma coerente contra o fundo de estrelas fixas, ranqueia os candidatos e gera relatórios para orientar a inspeção humana.

A Lia é uma ferramenta de apoio. Ela **não** confirma descobertas de asteroides, não substitui o Astrometrica, não substitui o Minor Planet Center (MPC) e não substitui o fluxo oficial de campanhas do IASC.

> O nome **Lia** é usado como nome do projeto e dedicação pessoal da autora. Não é uma sigla.

- **Autora:** Jaciana Barbosa
- **Repositório:** `lia-astrometry`
- **Versão:** `v1.5.0`
- **Licença:** MIT

## Por Que Isso Importa

Asteroides são pequenos corpos rochosos que sobraram da formação do Sistema Solar. Encontrar e medir asteroides ajuda astrônomos a melhorar órbitas, identificar novos objetos do cinturão principal e apoiar o trabalho de defesa planetária relacionado a objetos próximos da Terra.

O **International Astronomical Search Collaboration (IASC)** é um programa de ciência cidadã em que equipes inspecionam imagens de telescópios profissionais e submetem medições validadas por meio de um processo oficial de campanha. A NASA lista o IASC como um projeto de ciência cidadã da NASA Science, e a página de colaboradores do IASC lista a NASA entre os colaboradores dos Estados Unidos.

As campanhas do IASC são importantes porque permitem que estudantes e pessoas sem formação técnica avancada participem de trabalho astronômico real. A Lia foi criada para facilitar a primeira triagem sem esconder a necessidade de validação humana cuidadosa.

Fontes: [página do IASC na NASA Science](https://science.nasa.gov/citizen-science/international-astronomical-search-collaboration/) e [página de colaboradores do IASC](https://iasc.cosmosearch.org/Home/Collaborators).

## O Que A Lia Faz

Conjuntos de prática e campanha do IASC/Pan-STARRS geralmente contêm quatro imagens FITS do mesmo campo do céu. A maioria das estrelas fica parada de um frame para outro; possíveis asteroides se deslocam um pouco.

A Lia:

- carrega as quatro imagens FITS;
- mascara pixels saturados ou inválidos do Pan-STARRS;
- estima o fundo local do céu;
- detecta fontes pontuais;
- liga detecções entre os quatro frames;
- rejeita artefatos evidentes;
- classifica os candidatos de `STRONG` a `DISCARDED`;
- gera saídas em JSON, texto, PNG, CSV e rascunho auxiliar no estilo MPC.

O resultado é uma lista priorizada de candidatos. Um score alto significa “inspecione isto primeiro”, não “isto é um asteroide confirmado”.

## Exemplo de Saída

![Candidatos de exemplo](images/example_candidates.png)

Cada linha é um candidato. As colunas são os quatro frames em ordem cronológica. Círculos coloridos marcam o centroide medido, e setas mostram a direção do movimento.

## Instalação Rápida

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

A suíte atual tem 92 testes passando, cobrindo detecção de fontes, WCS, refinamento Gaia, score, máscara de saturação, rastreabilidade JSON e exportação CSV.

## Uso Básico

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

## Saídas

Para cada conjunto processado, a Lia gera:

- `*_candidates.json`: dados estruturados dos candidatos e metadados da execução;
- `*_report.txt`: relatório legível para inspeção;
- `*_MPC_report.txt`: rascunho auxiliar no estilo MPC, apenas para revisão;
- `*_candidates.png`: recortes visuais dos principais candidatos;
- `*_pipeline.log`: log de auditoria da execução.

O envio oficial da campanha ainda deve seguir o processo normal Astrometrica/IASC. Use a Lia para decidir onde olhar primeiro, depois inspecione e meça os candidatos manualmente.

## Entendendo O Score

O score é uma heurística transparente de ranking de 0 a 12:

| Score | Classe | Significado |
|---:|---|---|
| 8-12 | `STRONG` | inspecionar primeiro |
| 6-7 | `MODERATE` | candidato plausível |
| 4-5 | `WEAK` | candidato de baixa prioridade |
| 0-3 | `DISCARDED` | provável artefato ou fonte rejeitada |

O score combina movimento linear, consistência de velocidade entre frames, estabilidade de brilho, morfologia pontual, elongação e faixa de movimento esperada para sequências curtas do IASC. Ele não é um valor de confiança calibrado.

## Exportação de Métricas

Depois de rodar vários conjuntos de imagens, exporte CSVs para validação:

```bash
venv/bin/python tools/export_metrics.py results/ --out validation_metrics --recursive
```

Isso cria:

- `validation_metrics_sets.csv`: uma linha por conjunto processado;
- `validation_metrics_candidates.csv`: uma linha por candidato.

O CSV inclui colunas vazias de revisão manual, como `measured_in_astrometrica`, `included_in_astrometrica_mpc`, `iasc_feedback`, `manual_classification` e `notes`.

## Mais Documentação

Para a explicação científica e técnica completa, leia:

- [Metodologia em Portugues do Brasil](docs/methodology.pt-BR.md)
- [English methodology](docs/methodology.md)

Esses documentos explicam tempo em FITS, estimativa de fundo, ajuste PSF, refinamento Gaia DR3, modo WCS apenas pelo header, SkyBot, score, validação, limitações e por que a revisão humana continua necessária.

## Estrutura Do Projeto

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

## Referências

- [NASA Science: International Astronomical Search Collaboration](https://science.nasa.gov/citizen-science/international-astronomical-search-collaboration/)
- [Site oficial do IASC](https://iasc.cosmosearch.org/)
- [Colaboradores do IASC](https://iasc.cosmosearch.org/Home/Collaborators)
- Gaia Collaboration et al. (2023), *Gaia Data Release 3*, Astronomy & Astrophysics, 674, A1.
- Stetson, P. B. (1987), *DAOPHOT: A Computer Program for Crowded-Field Stellar Photometry*, PASP, 99, 191.
- Astropy Collaboration et al. (2022), *The Astropy Project*.
- Documentação do SkyBot/IMCCE, usado por meio de `astroquery.imcce.Skybot`.
- Documentação de astrometria e formato de observação do Minor Planet Center.

## Licença

Licença MIT. Consulte [LICENSE](LICENSE).
