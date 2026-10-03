# Cluster → Cell Type Annotations

PBMC flow/mass cytometry panel (Leiden clusters 0–39). Values are **relative scaled expression per marker** (0 = lowest-expressing cluster, 1 = highest), so a "high" value means that cluster is among the top expressers of that marker relative to the other clusters, not an absolute fluorescence intensity.

## Legend of lineage-defining markers

| Marker | Meaning |
|---|---|
| CD3 | T lineage |
| CD4 / CD8 | helper / cytotoxic T |
| CD56, CD16, CD57, GPR56 | NK / cytotoxic markers |
| CD19, CD20 | B lineage |
| CD14, CCR2, CD86, HLA-DR | monocytes / myeloid |
| CCR7, CD27, CD28, CD45RO | naive (CCR7+CD27+CD45RO−) vs memory (CD45RO+) vs terminally differentiated (CCR7−CD27−) |
| CXCR3 / CCR4 / CCR6 | Th1 / Th2 / Th17 homing |
| CD25 / CD127 | Treg (CD25+ CD127−) and activation (CD25+) |
| PD1, TIGIT, BTLA | exhaustion / regulation |
| CD38 | activation; memory B (low) vs plasmablast (high) |

## Full annotations

| Cluster | Cell type | Key markers | Confidence |
|---|---|---|---|
| 0 | CD4+ naive T cell | CD3+ CD4+ CD28+ CCR7+ CD27+ CD127+, CD45RO low | High |
| 1 | Naive B cell | CD19+ CD20+ HLA-DR+ BTLA+ CXCR5+, CD27− CD38 low | High |
| 2 | Classical monocyte | CD14+ CCR2+ HLA-DR+ CD86+ | High |
| 3 | Classical monocyte | CD14+ CCR2+ CD86+ HLA-DR+ | High |
| 4 | CD4+ memory T cell (EM, CCR4+) | CD3+ CD4+ CD28+ CD127+ CD45RO+, CCR7 mid | High |
| 5 | CD8+ TEMRA | CD3+ CD8+ CD57+ GPR56+ TIGIT+, CCR7− CD28− CD27− CD45RO− | High |
| 6 | CD4+ central memory T (TIGIT+ PD1+) | CD45RO+ CD27+ CCR7+ CD28+, TIGIT+ PD1+ | High |
| 7 | Classical monocyte | CD14+ CCR2+ CD86+ HLA-DR+ | High |
| 8 | CD4+ naive T cell | CD3+ CD4+ CD28+ CCR7+ CD27+ CD127+, CD45RO− | High |
| 9 | CD8+ memory T cell | CD3+ CD8+ CD45RO+ CD127+ CD28+ CD27+ | High |
| 10 | Naive / transitional B cell | CD19+ CD20+ CXCR5+ BTLA+ CCR6+ HLA-DR+, CD27− | High |
| 11 | Treg-like (CD25+ CD127− CCR4+) | CD25+ CD127 low CCR4+; CD3/CD4 relatively dim | Medium |
| 12 | Classical monocyte | CD14+ CCR2+ CD86+ HLA-DR+ | High |
| 13 | CD4+ central memory T (CCR4+/Th2-skewed) | CD45RO+ CCR7+ CD27+ CD28+ CD127+, CCR4+ | High |
| 14 | CD8+ effector memory T (TIGIT+) | CD3+ CD8+ CD45RO+ GPR56+, CCR7− CD28− CD27−, TIGIT+ | High |
| 15 | Treg-like (CD25+ CD127− CCR4+) | CD25+ CD127 low CCR4+; CD3/CD4 relatively dim | Medium |
| 16 | Intermediate monocyte | CD14+ CD16+ CD86+ HLA-DR+ | High |
| 17 | NK cell (CD16 low) | CD56+ GPR56+ CD38+ TIGIT+, CD3− | High |
| 18 | Ambiguous (monocyte/T-cell mixed) | CCR4+ CCR2+ CD14+ CD86+ CD45RO+ | Low |
| 19 | CD4+ effector memory T (CD57+ GPR56+) | CD3+ CD4+ CD45RO+, CCR7− CD27 low, CD57+ GPR56+ CXCR3+ | Medium |
| 20 | CD4+ central memory T cell | CD45RO+ CCR7+ CD27+ CD28+ CD127+ | High |
| 21 | NK cell (CD56dim CD16+) | CD56+ CD16+ GPR56+ TIGIT+, CD3− | High |
| 22 | NK cell (adaptive, CD57+) | CD56+ CD57+ GPR56+, CD3− | High |
| 23 | CD8+ naive T cell | CD3+ CD8+ CCR7+ CD27+ CD28+ CD127+, CD45RO− | High |
| 24 | Memory B cell | CD19+ CD20+ CD27+ CD38 low, CXCR5+ BTLA+ HLA-DR+ | High |
| 25 | CD4+ central memory T (Th1, CXCR3+) | CD45RO+ CCR7+ CD27+ CD28+, CXCR3+ | High |
| 26 | CD56bright NK cell | CD56 high, CD16−, CD38+ | High |
| 27 | Monocyte / DC-like (CD14 dim) | HLA-DR+ CD86+ CCR2+ CD38+, CD14 dim | Medium |
| 28 | Unclassified | PD1+ CXCR3+ CCR4+ CD19+, no clean lineage | Low |
| 29 | Activated Treg-like | CD25+ CD127− CD38+ CD45RO+ CCR4+; CCR2 also high (myeloid contamination?) | Medium |
| 30 | CD3+ CD4− CD8− T cells (γδ / DN, tentative) | CD3+ CD45RO+ CD127+, CD4− CD8− | Low |
| 31 | CD3+ naive-like T (lineage dim, tentative) | CCR7+ CD127+ CD28+ CD27+, CD3/CD4 dim | Low |
| 32 | Naive B cell | CD19+ CD20+ BTLA+ CXCR5+, CD27− CD38− | High |
| 33 | Unclassified (activated HLA-DR+ CD38+ GPR56+) | no clean lineage | Low |
| 34 | Non-classical monocyte | CD16+ CD14 dim CCR2− | Medium |
| 35 | Exhausted CD8+ T (CD8 dim) | CD3+ PD1+ TIGIT+ CD57+ GPR56+ CD45RO+, CD8 dim | Medium |
| 36 | Activated CD4+ T cell (CD25hi) | CD3+ CD4+ CD25 max, CD127+ CD28+ CCR7+ CD27+ | High |
| 37 | Doublets (T cell + monocyte) | CD3+ CD4+ and CD14+ CCR2+ CD86+ simultaneously | Medium |
| 38 | CD4+ naive T cell (activated subset) | CD3+ CD4+ CD28+ CCR7+ CD27+ CD127+, CD25+ BTLA+, CD45RO int | High |
| 39 | Multiplets / artifact | CD3+ CD19+ CD14+ all high with max/zero marker pattern | Low |

## Summary by lineage

- **CD4+ T cells**
  - Naive: 0, 8, 38
  - Central memory: 13 (CCR4+), 20, 25 (CXCR3+/Th1)
  - Memory / effector memory: 4 (CCR4+), 6 (TIGIT+ PD1+), 19 (CD57+ GPR56+)
  - Activated (CD25hi): 36
  - Treg-like: 11, 15, 29 (medium/low confidence)
- **CD8+ T cells**
  - Naive: 23
  - Memory: 9
  - Effector memory: 14 (TIGIT+)
  - TEMRA: 5 (CD57+ GPR56+)
  - Exhausted: 35 (PD1+ TIGIT+, CD8 dim)
- **CD3+ CD4− CD8− T cells (tentative):** 30 (γδ/DN), 31
- **NK cells:** 17 (CD16lo), 21 (CD16+), 22 (CD57+), 26 (CD56bright)
- **B cells**
  - Naive: 1, 10, 32
  - Memory: 24
- **Monocytes**
  - Classical: 2, 3, 7, 12
  - Intermediate: 16
  - Non-classical: 34
  - CD14 dim / DC-like: 27
- **Doublets / multiplets / artifact:** 37, 39
- **Unclassified / ambiguous:** 18, 28, 33

## Caveats

1. Values are relative (scaled per marker), so lineage-negative populations can still show moderate scaled values for markers they weakly express; use the top 3–5 markers per cluster as the phenotype.
2. Clusters 11/15/29 fit the Treg pattern (CD25+ CD127− CCR4+) but have relatively dim CD3/CD4 and low CD28, which is atypical for canonical Tregs — verify against FoxP3/Helios if available.
3. Cluster 29 also has high CCR2, suggesting possible monocyte contamination or an activated myeloid component.
4. Clusters 37 and 39 look like doublets/multiplets (co-expression of T, B, and myeloid markers; 39 has a suspicious max/zero pattern).
5. Cluster 35's exhausted phenotype (PD1+ TIGIT+) is classic for exhausted CD8+ T cells, but its CD8 is relatively dim (CD8 is known to be downregulated on exhausted cells).
6. Clusters 30/31 lack clean CD4/CD8 — could be γδ T, MAIT, or DN T cells; confirm with TCR Vδ2 / CD161 staining.
