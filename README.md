# Change Detection in SAR Images using Hypothesis Testing and Shannon Entropy based on the Rayleigh Distribution

<p align="center">
  <img src="images/capa_change_rayleigh_paper.png" alt="IEEE TGRS Paper Cover" width="750"/>
</p>

<p align="center">
  <a href="https://doi.org/10.1109/TGRS.2025.3601857"><img src="https://img.shields.io/badge/DOI-10.1109%2FTGRS.2025.3601857-00629B?style=for-the-badge&logo=ieee&logoColor=white" alt="DOI"/></a>
  <a href="https://quarto.org/"><img src="https://img.shields.io/badge/Quarto-Reveal.js-4B72C2?style=for-the-badge&logo=quarto&logoColor=white" alt="Quarto"/></a>
  <a href="https://www.ufpe.br/ppge"><img src="https://img.shields.io/badge/UFPE-PPGE-red?style=for-the-badge" alt="UFPE PPGE"/></a>
  <a href="https://jodavid.github.io/"><img src="https://img.shields.io/badge/Author-Jodavid%20Ferreira-blue?style=for-the-badge&logo=github" alt="Author"/></a>
</p>

---

## 📌 Visão Geral / Overview

Este repositório contém os slides da apresentação realizada no **Programa de Pós-Graduação em Estatística (PPGE)** da **Universidade Federal de Pernambuco (UFPE)** em 2026.

A apresentação detalha a metodologia, formulação matemática e resultados do artigo científico publicado na conceituada revista internacional **IEEE Transactions on Geoscience and Remote Sensing (TGRS)**:

> **Ferreira, J. A., & Palm, B. G. (2025).**  
> *Change Detection in SAR Images using Hypothesis Testing and Shannon Entropy based on the Rayleigh Distribution.*  
> **IEEE Transactions on Geoscience and Remote Sensing**, 2025.  
> 🔗 **DOI:** [10.1109/TGRS.2025.3601857](https://doi.org/10.1109/TGRS.2025.3601857)

Os slides foram desenvolvidos utilizando o ecossistema [Quarto](https://quarto.org/) com o formato [Reveal.js](https://revealjs.com/), combinando tipografia matemática rica ($\LaTeX$), visualizações dinâmicas e design responsivo.

---

## 📄 Resumo (Abstract)

### Português
Esta apresentação introduz um algoritmo estatístico para detecção de mudanças em imagens de Radar de Abertura Sintética (SAR) baseado em teoria da informação (Entropia de Shannon) e na distribuição Rayleigh, contemplando dados tanto no formato *single-look* quanto *multi-look*. 

Embora testes de hipóteses sejam frequentemente empregados em sensoriamento remoto por radar, poucos métodos fornecem uma estrutura unificada e analiticamente tratável para ambos os regimes de processamento. A formulação proposta baseia-se no cálculo da variância assintótica da entropia estocástica e constrói uma estatística de teste assintoticamente qui-quadrado. A abordagem foi validada por meio de extensivas simulações de Monte Carlo e aplicada a conjuntos de dados reais: imagens *single-look* do dataset **CARABAS II** (Suécia) e dados *multi-look* polarimétricos do **UAVSAR** (NASA/JPL, Los Angeles). Os resultados comprovam a precisão do método, superando detectores clássicos multivariados mesmo operando em canais individuais de amplitude.

### English
This presentation is about the article that presents a change detection algorithm based on Shannon entropy and the Rayleigh distribution for both single-look and multilook synthetic aperture radar (SAR) images. While hypothesis testing has been widely used in SAR change detection, few studies have applied this approach to both types of data. To address this, we propose a method that utilizes Shannon entropy to detect changes between two samples. The performance of the algorithm was evaluated through Monte Carlo simulations using synthetic SAR data and further validated on real-world datasets, including single-look images from the CARABAS II dataset and multilook data from the uninhabited aerial vehicle synthetic aperture radar (UAVSAR). The results demonstrate that the proposed method is effective in detecting changes, highlighting its versatility and reinforcing the way for future research into alternative entropy measures and probability distributions in change detection tasks.

---

## 🔬 Fundamentação Teórica e Metodologia

```
 ┌──────────────────────┐        ┌──────────────────────┐
 │ SAR Image (Time t₁)  │        │ SAR Image (Time t₂)  │
 └──────────┬───────────┘        └──────────┬───────────┘
            │                               │
            ▼                               ▼
 ┌──────────────────────────────────────────────────────┐
 │   Sliding Window Extraction (w × w) & Rayleigh MLE   │
 │            β̂₁ = √[ (1/2N) ∑ y₁² ]                    │
 │            β̂₂ = √[ (1/2N) ∑ y₂² ]                    │
 └──────────────────────────┬───────────────────────────┘
                            │
                            ▼
 ┌──────────────────────────────────────────────────────┐
 │       Shannon Stochastic Entropy Calculation         │
 │          H_S(β̂) = 1 + ln(β̂/√2) + γ/2                 │
 │       Key Result:  Var[H_S(β̂)] = σ_S² = 1/4          │
 └──────────────────────────┬───────────────────────────┘
                            │
                            ▼
 ┌──────────────────────────────────────────────────────┐
 │           Hypothesis Testing: H₀ vs H₁               │
 │    S_S = [4·N₁·N₂ / (N₁ + N₂)] · [H_S(β̂₁) - H_S(β̂₂)]²  │
 │            Under H₀: S_S ~ χ²(1)                     │
 └──────────────────────────┬───────────────────────────┘
                            │
                            ▼
 ┌──────────────────────────────────────────────────────┐
 │           P-value Thresholding (P(χ² > S_S) ≤ α)     │
 │               Binary Change Map                      │
 └──────────────────────────────────────────────────────┘
```

### 1. Modelagem Estatística via Distribuição Rayleigh
Em áreas homogêneas de imagens SAR de amplitude, a reflexão do sinal resultante da superposição de múltiplos espalhadores elementares independentes faz com que a amplitude $Y$ siga uma distribuição Rayleigh com parâmetro de escala $\beta > 0$:

$$f_Y(y; \beta) = \frac{y}{\beta^2} e^{-\frac{y^2}{2\beta^2}}, \quad y \ge 0$$

- **Estimador de Máxima Verossimilhança (EMV):**
  $$\hat{\beta} = \sqrt{\frac{1}{2N} \sum_{n=1}^{N} y_n^2}$$
- **Informação de Fisher:** $\mathcal{K}(\beta) = \frac{4N}{\beta^2}$
- **Normalidade Assintótica:** $\hat{\beta} \xrightarrow{\mathcal{D}} \mathcal{N}\left(\beta, \frac{\beta^2}{4N}\right)$

### 2. Entropia de Shannon e Classe $(h, \phi)$
Para a variável aleatória contínua $Y \sim \mathcal{R}(\beta)$, a entropia diferencial de Shannon é:

$$H_S(\beta) = 1 + \ln\left(\frac{\beta}{\sqrt{2}}\right) + \frac{\gamma}{2}$$

onde $\gamma \approx 0.5772$ denota a constante de Euler-Mascheroni.

### 3. Resultado Chave: Variância Invariante da Entropia
Aplicando o Método Delta assintótico:

$$\sigma_S^2(\beta) = \left(\frac{\partial H_S(\beta)}{\partial \beta}\right)^2 \cdot \mathcal{K}(\beta)^{-1} = \left(\frac{1}{\beta}\right)^2 \left(\frac{\beta^2}{4}\right) = \frac{1}{4}$$

> **Teorema Prático:** A variância da entropia estocástica de Shannon estimada sob a lei de Rayleigh é **estritamente constante** ($\sigma_S^2 = 1/4$), independente do parâmetro de escala $\beta$.

### 4. Estatística do Teste de Hipóteses
Para avaliar se duas amostras independentes $N_1$ e $N_2$ (antes e depois) provêm da mesma distribuição ($H_0: H_S(\beta_1) = H_S(\beta_2)$ vs. $H_1: H_S(\beta_1) \neq H_S(\beta_2)$):

$$S_S(\hat{\beta}_1, \hat{\beta}_2) = \frac{4 N_1 N_2}{N_1 + N_2} \left[ H_S(\hat{\beta}_1) - H_S(\hat{\beta}_2) \right]^2 \xrightarrow{\mathcal{D}} \chi^2_1$$

Para janelas idênticas ($N_1 = N_2 = N$):

$$S_S = 2 N \left[ H_S(\hat{\beta}_1) - H_S(\hat{\beta}_2) \right]^2$$

Rejeita-se $H_0$ ao nível nominal $\alpha$ se $P(\chi^2_1 > S_S) \le \alpha$.

---

## 📊 Principais Resultados

### Simulações de Monte Carlo
- **Erro Tipo I:** Verificado com 5.000 réplicas em janelas de $3\times3$ até $11\times11$, demonstrando excelente aderência aos níveis nominais ($1\%$, $5\%$ e $10\%$).
- **Curva ROC e Poder:** Área sob a curva **AUC = 0.9019** para pequenas variações paramétricas ($\Delta\beta = 0.03$), com convergência do poder do teste para $1$ à medida que o tamanho amostral aumenta.

<p align="center">
  <img src="images/ROC_curve.png" alt="Curva ROC" width="500"/>
</p>

### Dados Reais: CARABAS II (*Single-Look*)
- Sistema UWB VHF (Suécia) com 25 alvos militares ocultos sob densa cobertura florestal (*clutter* natural).
- **Desempenho:** Detecção de **25/25 alvos (100%)** com **zero alarmes falsos** registrados na área analisada.

<p align="center">
  <img src="images/th/m1_p1_entropia_editado.png" alt="Mapa de Mudanças CARABAS II" width="550"/>
</p>

### Dados Reais: UAVSAR Los Angeles (*Multi-Look*)
- Sensor banda L (NASA/JPL) em cenário urbano complexo (Abril/2009 vs. Maio/2015).
- Comparação com métodos multivariados da literatura baseados na distribuição de Wishart (Nascimento et al., 2018; Bouhlel et al., 2020):

| Detector | Canais Utilizados | FP (%) | FN (%) | FA (%) | Taxa Detecção (DR %) | Kappa ($\kappa$ %) | MCC |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| $S_{LR}$ (Razão de Verossimilhança) | Todos (PolSAR) | 0.06 | 13.43 | 13.49 | 20.41 | 24.56 | - |
| $S_{KL}$ (Kullback-Leibler) | Todos (PolSAR) | 0.05 | 14.06 | 14.11 | 19.77 | 23.02 | - |
| $S_S$ (Shannon Wishart) | Todos (PolSAR) | 0.34 | 5.47 | 5.82 | 35.39 | 52.48 | - |
| $S_R$ (Rényi) | Todos (PolSAR) | 0.43 | 3.71 | 4.14 | 42.99 | 62.27 | - |
| DRT (Determinant Ratio Test) | Todos (PolSAR) | - | - | - | 63.38 | - | - |
| HLT (Hotelling-Lawley Trace) | Todos (PolSAR) | - | - | - | 52.83 | - | - |
| **$S_S(\mathcal{R})$ (Proposto)** | **HH** | 4.13 | 2.82 | 8.08 | **79.86%** | 72.22% | 0.72 |
| **$S_S(\mathcal{R})$ (Proposto)** | **HV** | 1.72 | 5.82 | 9.35 | 69.99% | 73.79% | 0.75 |
| **$S_S(\mathcal{R})$ (Proposto)** | **VV** | 2.12 | 4.60 | 8.17 | 74.13% | **75.69%** | **0.76** |

---

## 📁 Estrutura do Repositório

```plaintext
.
├── SEMINARIO_UFPE_Rayleigh_and_Information.qmd   # Código-fonte da apresentação em Quarto Reveal.js
├── SEMINARIO_UFPE_Rayleigh_and_Information.html  # Apresentação compilada em HTML pronta para visualização
├── css/                                          # Folhas de estilo personalizadas
│   ├── jodavid.css                               # Estilização visual dos slides
│   ├── extra.css
│   └── mmp.css
├── images/                                       # Ilustrações, diagramas e mapas de detecção
│   ├── capa_change_rayleigh_paper.png
│   ├── carabas/                                  # Imagens do sensor CARABAS II
│   ├── resultados/out1/                          # Resultados UAVSAR (Los Angeles)
│   ├── power_test/                               # Gráficos de poder do teste de hipóteses
│   ├── ROC_curve.png                             # Curva ROC Monte Carlo
│   └── th/                                       # Ground truth e máscaras de mudança
├── Slide/                                        # Versões e templates auxiliares em LaTeX/Beamer
└── README.md                                     # Documentação principal do repositório
```

---

## 💻 Como Visualizar e Compilar

### 1. Visualização Direta
Você pode abrir diretamente o arquivo [`SEMINARIO_UFPE_Rayleigh_and_Information.html`](file:///media/jodavid/JODA_HD2/00_ESTATISTICA_HD2/11_APRESENTACOES/2026/SEMINARIO_UFPE_PPGE_Rayleigh_and_Information_Theory/SEMINARIO_UFPE_Rayleigh_and_Information.html) em qualquer navegador moderno (Chrome, Firefox, Safari, Edge):

- Pressione `F` para entrar em modo tela cheia.
- Pressione `O` para alternar para visão panorâmica dos slides (*Overview*).
- Pressione `S` para abrir a janela de notas do apresentador (*Speaker view*).
- Pressione `B` para abrir o quadro-negro digital (*Chalkboard*).

### 2. Compilar via Quarto CLI
Certifique-se de possuir o [Quarto](https://quarto.org/) instalado em seu sistema:

```bash
# Para visualizar localmente com recarregamento automático (live-reload):
quarto preview SEMINARIO_UFPE_Rayleigh_and_Information.qmd

# Para compilar e gerar o arquivo HTML estático final:
quarto render SEMINARIO_UFPE_Rayleigh_and_Information.qmd
```

---

## 📖 Como Citar (Citation)

Se este trabalho ou os conceitos desta apresentação forem úteis em sua pesquisa, por favor cite o artigo oficial da IEEE TGRS:

```bibtex
@article{ferreira2025change,
  author    = {Ferreira, Jodavid A. and Palm, Bruna G.},
  journal   = {IEEE Transactions on Geoscience and Remote Sensing}, 
  title     = {Change Detection in SAR Images using Hypothesis Testing and Shannon Entropy based on the Rayleigh Distribution}, 
  year      = {2025},
  volume    = {},
  number    = {},
  pages     = {1--12},
  doi       = {10.1109/TGRS.2025.3601857}
}
```

---

## 👤 Autor e Contato

**Prof. Dr. Jodavid Ferreira**  
*Departamento de Estatística - Universidade Federal de Pernambuco (UFPE)*  
*Programa de Pós-Graduação em Estatística (PPGE)*

- 🌐 **Website Pessoal:** [jodavid.github.io](https://jodavid.github.io/)
- 📄 **Currículo Lattes:** [lattes.cnpq.br/4617170601890026](http://lattes.cnpq.br/4617170601890026)
- 💼 **LinkedIn:** [jodavidferreira](https://www.linkedin.com/in/jodavidferreira/)
- ✉️ **E-mail:** [jodavid.ferreira@ufpe.br](mailto:jodavid.ferreira@ufpe.br)

---
<p align="center">
  <sub>Apresentado no Seminário do Programa de Pós-Graduação em Estatística (PPGE / UFPE) - 2026</sub>
</p>
