# CCC Themes — Manual Curation

Filtered: `interaction_padj < 0.05` → 130 interactions, 63 unique LR pairs, 92 source-target pairs.

Scores are differential (PPD vs. control): positive = upregulated in PPD, negative = downregulated.

---

## Theme 1 — Axon Guidance

| LR pair | n interactions | mean score |
|---|---|---|
| Sema6a^Plxna2 | 6 | +0.40 |
| Sema6d^Kdr_Plxna1 | 1 | −0.85 |
| Sema3d^Nrp2_Plxna2 | 1 | +0.96 |
| Sema3c^Nrp2_Plxna2 | 1 | +0.97 |
| Slit2^Robo1 | 1 | +1.00 |
| Slit2^Robo2 | 1 | −0.85 |
| Slit3^Robo1 | 1 | **+2.49** |
| Nell2^Robo3 | 2 | **+1.45** |
| Ncam1^Robo1 | 1 | **+1.73** |
| Ncam1^Robo3 | 1 | +0.76 |
| Efna5^Ephb1 | 1 | +0.04 |
| App^Dcc | 1 | −0.18 |
| Slit2^Dcc | 1 | −0.06 |

---

## Theme 2 — Cell Adhesion / CAMs

| LR pair | n interactions | mean score |
|---|---|---|
| Cntn4^Ptprg | 8 | +0.75 |
| Cntn1^Ptprz1 | 2 | +1.08 |
| Ncam1^Ptprz1 | 2 | +0.71 |
| L1cam^Ptprz1 | 1 | +0.14 |
| L1cam^Cntn1 | 1 | +0.50 |
| Nfasc^Cntn1_Cntnap1 | 1 | +1.06 |
| Alcam^Chl1 | 6 | −0.08 |
| Fstl5^Chl1 | 3 | +0.12 |
| Ntng1^Lrrc4c | 1 | −1.08 |
| Ntng2^Lrrc4c | 1 | −1.19 |

---

## Theme 3 — Neurotrophic / Growth Factors

| LR pair | n interactions | mean score |
|---|---|---|
| Bdnf^Ntrk2 | 3 | −0.06 |
| Bdnf^Sort1 | 3 | +0.23 |
| Ptn^Ncl | 8 | **−0.82** |
| Ptn^Ptprz1 | 5 | +0.24 |
| Ptn^Itgav_Itgb3 | 5 | −0.29 |
| Ptn^Alk | 3 | +0.28 |
| Ptn^Cdh10 | 2 | −0.30 |
| Mdk^Sorl1 | 3 | **−0.83** |
| Mdk^Ncl | 2 | **−1.74** |
| Mdk^Alk | 1 | −0.21 |
| Mdk^Ptprz1 | 1 | −0.25 |
| Fstl1^Dip2a | 4 | +0.84 |

---

## Theme 4 — ECM-Integrin Remodeling

| LR pair | n interactions | mean score |
|---|---|---|
| Fn1^Itgav_Itgb3 | 2 | +0.36 |
| Fn1^Itgav_Itgb6 | 2 | +0.36 |
| Fn1^Itga9 | 2 | +0.94 |
| Fn1^Itgav_Itgb8 | 1 | +0.77 |
| Fn1^Nt5e | 1 | +0.98 |
| Tnr^Itgav_Itgb3 | 2 | +0.82 |
| Tnr^Itgav_Itgb6 | 2 | +0.82 |
| Col4a3^Itgav_Itgb8 | 2 | +0.27 |
| Col6a1^Itgav_Itgb8 | 1 | +0.68 |
| Col18a1^Kdr | 1 | −0.12 |
| Lama1^Nt5e | 1 | +1.05 |
| Lama4^Itgav_Itgb8 | 1 | −0.12 |
| Lamb1^Itgav_Itgb8 | 1 | +0.67 |
| Vwf^Itgav_Itgb3 | 2 | +0.36 |
| Vwf^Itga9 | 2 | +0.95 |
| Gpc3^Igf1r | 3 | −0.48 |
| Gpc3^Unc5c | 1 | **−3.20** |
| Gpc3^Unc5d | 3 | −0.70 |

---

## Theme 5 — APP / Amyloid-Related

| LR pair | n interactions | mean score |
|---|---|---|
| App^Aplp2 | 1 | −0.86 |
| App^Dcc | 1 | −0.18 |
| App^Lrp6 | 1 | +0.04 |
| Spon1^App | 2 | −0.11 |
| Slit2^App | 1 | −0.86 |
| Spon1^Lrp8 | 2 | +0.83 |

---

## Theme 6 — Morphogen / Developmental Signaling

| LR pair | n interactions | mean score |
|---|---|---|
| Rspo2^Rnf43 | 2 | +0.43 |
| Bmp6^Bmpr1a_Bmpr2 | 1 | +0.68 |
| Bmp6^Bmpr1b_Bmpr2 | 1 | +0.68 |
| Timp3^Kdr | 3 | +0.07 |

---

## Top 5 Strongest Signals

| Rank | LR pair | Source | Target | Score | Theme |
|---|---|---|---|---|---|
| 1 | Gpc3^Unc5c | VLMC | VLMC | −3.20 | ECM-Integrin |
| 2 | Slit3^Robo1 | VLMC | Oligo | +2.49 | Axon Guidance |
| 3 | Nell2^Robo3 | OPC | L5 ET Glut | +2.27 | Axon Guidance |
| 4 | Ptn^Ncl | Oligo | Pvalb Gaba | −2.15 | Neurotrophic |
| 5 | Mdk^Ncl | Choroid | Pvalb Gaba | −1.84 | Neurotrophic |

---

## Notes

- `App^Dcc` and `Slit2^App` span themes 1 and 5 — assigned to theme 5 (APP ligand/receptor context takes priority).
- `Ptn^Ptprz1` and `Mdk^Ptprz1` straddle themes 2 and 3 — assigned to theme 3 (ligand drives biology).
- Theme 6 is small (4 pairs); consider merging into ECM or keeping as developmental context.
- Choroid plexus (325 CHOR NN) appears as a sender in themes 3 and 4 — may be worth highlighting on poster.