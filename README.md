# Single-Cell Multi-Omics And Bulk RNA-seq Analysis of Human Dendritic Cell Homeostatic Maturation

## Abstract
Homeostatic maturation of dendritic cells (DCs) is crucial for maintaining peripheral tolerance. However, the existence and underlying molecular mechanisms of this DC maturation program in humans remain poorly understood. Here, we identified a homeostatic maturation program of type 2 conventional DCs (cDC2s) within human spleens. The homeostatically maturated cDC2s (DC2hm) upregulate common DC maturation genes, preferentially express thymic stromal lymphopoietin (TSLP) receptor, and are enriched in TSLP-STAT5 signaling. DC2hm also showed higher expression of the autoimmune regulator (AIRE) gene and tissue-restricted antigens. Functionally, DC2hm are potent in inducing T regulatory cell differentiation. TSLP administration stimulates a DC2hm-like program both ex vivo and in vivo in humanized mice. Our study thus unravels a human DC homeostatic maturation program and underscores the role of TSLP in inducing DC2hm important for maintaining peripheral tolerance.

---

## Overview
This repository contains the analysis pipelines, data preprocessing scripts, and visualization tools used in our study of dendritic cell (DC) homeostatic maturation. We employ single-cell multi-omics approaches to unravel the molecular mechanisms underpinning DC maturation and its role in maintaining immune tolerance.

---

## Installation
Clone this repository and install dependencies:

```bash
git clone https://github.com/yourusername/dendritic-cell-analysis.git
cd dendritic-cell-analysis
pip install -r requirements.txt
```

Refer to the [prerequisites](#prerequisites) for additional setup information.

---

## Prerequisites
- **Python:** >= 3.8
- **R:** >= 4.0
- **Required Packages:**
  - Scanpy
  - Seurat
  - Harmony
  - scRNA-seq-specific libraries

---


## Data and Results
- **Raw Data:** Input datasets are stored in the `data/` folder.
- **Intermediate Results:** Workflow outputs, including processed matrices, are stored in `results/`.
- **Figures:** Publication-quality plots and graphs are generated in `figures/`.

---


## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contact
For questions or collaborations, please contact Harris Han(HarrisHan@hotmail.com) or the corresponding author Dr. Liang Cheng (liangcheng@whu.edu.cn)

---
