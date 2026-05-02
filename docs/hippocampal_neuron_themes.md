# Hippocampal Neuron Themes — CORT vs OIL, Postpartum Rat Brain
## Differential Gene Expression + Cell-Cell Communication

**Datasets:** Slide-tags spatial transcriptomics (female rats, postpartum day 1–23)  
**Comparison:** Corticosterone-treated (CORT) vs vehicle-treated (OIL)  
**Analyses:** EdgeR LRT (DEG) + LIANA+ (CCC), library size filter 50k  
**Date:** 2026-05-02

---

## Background

Female Sprague-Dawley rats were treated with corticosterone (CORT) or oil vehicle from postpartum day 1 through day 23, modelling sustained glucocorticoid elevation as a precipitating factor in postpartum depression (PPD). Brain tissue was profiled by Slide-tags single-nucleus spatial transcriptomics. The analysis below focuses on hippocampal neurons — primarily CA1-ProS, CA3, and DG glutamatergic populations — and their surrounding microenvironment, integrating DEG results with cell-cell communication (CCC) findings.

### Structural observation motivating this report

Non-neuronal cells (astrocytes, oligodendrocytes, OPCs, endothelium, microglia, VLMCs, pericytes) respond to CORT in a largely **homogeneous, high-magnitude** manner: near-universal upregulation of GR target genes (Fkbp5, Sgk1, Tsc22d3, Ddit4, Nfkbia) across all glial/vascular types, with DEG counts ranging from 105 (microglia) to 608 (oligodendrocytes). Hippocampal neurons, by contrast, respond in a **heterogeneous, cell-type-specific** pattern, with fewer DEGs (CA1: 84, DG: 91, CA3: 9) and directional effects that frequently oppose the glial response for the same gene. This dissociation is the central biological phenomenon described below.

---

## Theme 1 — Glucocorticoid Receptor Desensitization and HPA Axis Disruption

**Core finding:** Hippocampal neurons downregulate both glucocorticoid receptor (GR) and mineralocorticoid receptor (MR) under chronic CORT, while glial cells show the opposite — strong upregulation of GR target genes.

| Gene | Cell type | logFC | FDR | Direction vs. glia |
|---|---|---|---|---|
| Nr3c1 (GR) | CA1-ProS, DG | −0.5 to −0.9 | sig | **opposite** (glia: GR targets strongly up) |
| Nr3c2 (MR) | CA1-ProS | −0.58 | 0.034 | unique to hippocampus |
| Sgk1 | DG | −2.9 | 1.2e-13 | **opposite** (Oligo/OPC: +7.0/+0.8) |
| Fkbp5 | (not significant in CA1/DG) | — | — | absent vs. massive glial response |

The GR:MR ratio in hippocampus is a primary determinant of stress resilience. MR binds corticosterone with high affinity at basal levels (tonic signal); GR is activated at stress peak (phasic signal). Downregulation of both in CA1 under chronic CORT reflects receptor desensitization — a failure of the hippocampus to appropriately process glucocorticoid signals. This impairs negative feedback on the HPA axis, promoting hypercortisolemia, which is among the most replicated endocrine findings in PPD.

The cell-type dissociation for Sgk1 is particularly striking: oligodendrocytes show the largest single-gene upregulation in the dataset (logFC = +7.0, FDR = 1.6e-59), while DG granule cells show strong downregulation (logFC = −2.9, FDR = 1.2e-13). Sgk1 regulates ion transport, synaptic plasticity, and neuronal survival under stress — its loss in DG granule cells has direct implications for hippocampal plasticity.

**From CCC:** The Cntn4→Ptprg axis (L6 IT CTX → endothelium, astrocytes; FDR = 8.4e-05) is relevant here: PTPRG is a negative regulator of BDNF/TrkB signaling, and its upregulation in CORT could suppress the downstream pro-survival signaling that normally compensates for GR desensitization.

---

## Theme 2 — Reduced Trophic Support to Hippocampal Neurons (BDNF and PTN Axes)

**Core finding:** Hippocampal neurons experience a coordinated reduction in neurotrophic input — both from within (Bdnf in DG) and from surrounding cells (Ptn from oligodendrocytes, Mdk from choroid plexus).

### BDNF axis

| Gene | Cell type | logFC | FDR |
|---|---|---|---|
| Bdnf | DG | −0.74 | 0.017 |
| Ntrk2 (TrkB) | Astro-NT, Astro-TE, OPC | +0.58 to +0.65 | 0.021–0.048 |

Hippocampal BDNF reduction is one of the most replicated findings in rodent stress and depression models. Bdnf is downregulated specifically in DG — the site of adult neurogenesis — and not in CA1 or CA3, consistent with its role in dentate gyrus-dependent pattern separation, spatial memory, and neurogenic reserve. Antidepressants (SSRIs, ketamine) converge on BDNF-TrkB restoration.

The upregulation of TrkB (Ntrk2) in astrocytes and OPCs while BDNF is reduced in DG neurons creates an asymmetry: glial cells may be attempting to capture the diminishing BDNF supply, but the net effect on neuronal TrkB signaling — the pro-survival, pro-plasticity arm — is expected to be reduced.

**From CCC (FDR < 0.05):** The Bdnf→Sort1 interaction (DG Glut → Microglia; FDR = 8.7e-3) represents a mechanistic pivot: proBDNF (the precursor, which has opposite biological effects to mature BDNF) signals through Sortilin (Sort1) on microglia to promote apoptosis and synaptic pruning rather than survival. In CORT, Sort1 is strongly upregulated in microglia (logFC = +2.459) while Bdnf is down in DG — shifting the balance toward the proBDNF/Sort1 pro-apoptotic arm. This is a mechanistic link from hippocampal CORT exposure to microglial-driven synaptic loss.

From CCC at p < 0.05: Bdnf→Ntrk2 (DG → Astro-NT) is directionally down in CORT, confirming reduced canonical BDNF signaling from DG across the hippocampal microenvironment.

### PTN / MDK axis (oligodendrocyte and choroid plexus inputs)

| Finding | FDR | Direction |
|---|---|---|
| Ptn (DEG in oligodendrocytes) | 4.1e-06 | down (logFC = −1.9) |
| Ptn→Ncl CCC (Oligo→L6 CT; most significant Oligo interaction) | 1.6e-04 | down in CORT |
| Ptn→Itgav/Itgb3 (Oligo→L6 IT Glut) | 2.8e-05 | down in CORT |
| Mdk→Ncl / Sorl1 (Choroid Plexus → Astro + neurons) | 2.5–2.7e-02 | down in CORT |

PTN (Pleiotrophin) and MDK (Midkine) are structurally related heparin-binding growth factors that promote neurite outgrowth, oligodendrocyte differentiation, and synaptic maintenance. Their downregulation from two independent sources — mature oligodendrocytes (PTN) and the choroid plexus (MDK) — represents a convergent reduction in neurotrophic support to cortical and hippocampal target neurons. The choroid plexus is a key neuroendocrine interface for glucocorticoid entry into the CSF and has been specifically implicated in maternal brain remodeling; MDK downregulation from this structure directly reflects glucocorticoid suppression of CSF neurotrophic signaling.

**Summary:** Hippocampal neurons in CORT are receiving less Bdnf (from within), less PTN (from oligodendrocytes), and less MDK (from choroid plexus), while the proBDNF→Sortilin pro-apoptotic signal to microglia is increased.

---

## Theme 3 — Circadian Clock Disruption Originating in Glia, Converging on Hippocampal Circuitry

**Core finding:** All four canonical clock genes (Per1, Per2, Per3, Cry1) are upregulated across glial and vascular cells by CORT, representing a broad circadian phase-shift driven by GR binding to clock gene promoters. Hippocampal neurons are not the primary source of this signal, but are embedded in a brain environment where the molecular clock is dysregulated in all surrounding cell types.

| Gene | Cell types with significant upregulation | logFC range | FDR range |
|---|---|---|---|
| Per2 | Endo, Microglia, Astro-NT, Astro-TE, L6_IT, L6_CT, L4_5_IT | +0.8 to +5.2 | 1.6e-04 to 0.027 |
| Per1 | Oligo, Endo | +2.3 to +3.5 | 0.003 to 0.026 |
| Per3 | Microglia, Endo | +1.4 to +2.2 | 0.001 to 0.042 |
| Cry1 | Oligo, Endo | +1.2 to +3.5 | 2.5e-04 to 0.036 |

Per2 is the most broadly distributed, appearing in both non-neuronal cells and cortical neurons. The hippocampal neurons (CA1, CA3, DG) themselves do not show significant Per/Cry changes at FDR < 0.05 — the circadian disruption is primarily a glial/vascular signature. This matters because the hippocampal molecular clock is entrained in part by signals from surrounding astrocytes and microglia. If those cells are phase-shifted by CORT via GRE-driven Per2 upregulation, the hippocampal neuronal clock is expected to be disrupted at the circuit level even without direct transcriptional changes in hippocampal neurons.

**Biological relevance to PPD:** The postpartum period is characterized by profound circadian misalignment (disrupted sleep, altered cortisol diurnal rhythm, light-dark schedule changes from infant care). The convergent upregulation of Per2/Per3 across endothelium and glia in CORT animals models the HPA-clock coupling that underlies this misalignment. Clock gene variants (especially PER3) are associated with vulnerability to postpartum mood disorders.

**Mechanistic note:** CORT acutely phase-shifts the circadian clock by binding GREs in the Per1/Per2 promoters, directly inducing their expression. The result — upregulation of negative-arm clock components (Per genes) — paradoxically disrupts clock periodicity rather than simply phase-shifting it, contributing to circadian arrhythmia.

---

## Theme 4 — Neuroinflammatory Sensitization of Hippocampal Neurons (IL-6 Receptor Upregulation and Opposing Directionality)

**Core finding:** Hippocampal neurons upregulate Il6r (IL-6 receptor) under CORT, increasing their sensitivity to circulating or locally produced IL-6, while glial cells simultaneously activate anti-inflammatory programs. This creates opposing inflammatory trajectories in neurons vs. glia within the same tissue.

| Gene | Cell types | logFC | FDR | Note |
|---|---|---|---|---|
| Il6r | CA1-ProS, CA3, DG (+ 5 cortical types) | +0.8 to +1.7 | 7.0e-05 to 0.046 | **up in neurons** |
| Nfkbia (IκBα) | Astro-NT, Astro-TE, Oligo, Endo, Peri | +2.0 to +4.3 | 2.4e-04 to 3.8e-04 | anti-inflammatory, **up in glia** |
| Tsc22d3 (GILZ) | Astro-NT, Astro-TE, Oligo, Endo, VLMC | +0.9 to +2.6 | up to 0.019 | anti-inflammatory, **up in glia** |
| Ifngr1 | Endo, Microglia | +1.1 | sig | up in non-neuronal |
| Ifngr1 | CA1-ProS | −1.5 | sig | **down in hippocampal neurons** |
| Mif | Microglia | +3.5 | 1.3e-04 | pro-inflammatory |

The Il6r finding is the most striking: eight cell types (5 cortical layers + CA1/CA3/DG) upregulate the IL-6 receptor, sensitizing a large fraction of brain neurons to IL-6. IL-6 signaling in hippocampus is associated with stress-induced anhedonia and anxiety, and elevated IL-6 is documented in women with PPD. The receptor upregulation under CORT could prime hippocampal neurons for an exaggerated IL-6 response even at basal cytokine levels.

The Ifngr1 directionality contrast is a clean example of the neuron-glia dissociation: endothelium and microglia upregulate the IFN-γ receptor (pro-inflammatory sensitization), while CA1 neurons downregulate it (possible tolerance or desensitization). These cells are in physical proximity — the finding suggests fundamentally different inflammatory set-points are being established in neurons vs. the surrounding microenvironment.

**From CCC:** The Bdnf→Sort1 axis (DG→Microglia, FDR = 8.7e-3) links this directly to neuroinflammation: increased pro-BDNF/Sortilin signaling from DG neurons to microglia promotes microglial activation and synaptic pruning. Microglia already show upregulation of Mif (macrophage migration inhibitory factor; logFC = +3.5, FDR = 1.3e-4), a pro-inflammatory cytokine that counter-regulates glucocorticoid action — creating a positive feedback loop where CORT-activated microglia produce a cytokine that limits glucocorticoid's own anti-inflammatory effects.

**The paradox:** CORT simultaneously activates strong anti-inflammatory programs in glia (Nfkbia, Tsc22d3) and sensitizes neurons to pro-inflammatory cytokine signaling (Il6r). This duality may explain why PPD is associated with both immune suppression markers and elevated inflammatory cytokines in peripheral and CNS measurements.

---

## Theme 5 — Disrupted Excitatory-Inhibitory Balance via Hippocampal and Interneuron-Specific Signaling

**Core finding:** The hippocampal excitatory-inhibitory balance is disrupted under CORT through multiple convergent mechanisms: (a) altered interneuron→excitatory neuron CCC, (b) specific vulnerability of Sst interneurons in cortex, and (c) hippocampal CA2/CA1 neurons reducing signals that normally coordinate glutamatergic circuits.

### Interneuron-specific disruption

**Sst interneurons:** Uchl1 (Ubiquitin C-terminal hydrolase L1 / PGP9.5) is significantly downregulated in Sst interneurons specifically (logFC = −1.174, FDR = 0.0048), with a non-significant trend in Pvalb interneurons. Uchl1 is essential for maintaining the free ubiquitin pool and proteasomal clearance of misfolded proteins. Sst interneurons are among the most metabolically demanding neurons in cortex — they innervate dendritic compartments of pyramidal cells across layers, regulating gain control and coincidence detection. Their selective vulnerability under CORT, compounded by protein homeostasis failure (reduced Uchl1), likely disrupts dendritic computation in the circuits they regulate.

**From CCC — Pvalb vs CA2-FC-IG hippocampal input (Fstl5→Chl1):**

| Source | Direction of Fstl5 in CORT | Target |
|---|---|---|
| Pvalb Gaba | **up** (+1.013) | L5/L6 cortical Glut, CA1-ProS, DG |
| CA2-FC-IG Glut | **down** (−0.580) | overlapping cortical + hippocampal targets |

Pvalb interneurons increase Fstl5 (a BMP antagonist with trophic properties) in CORT, while hippocampal CA2/FC neurons decrease it, targeting an overlapping set of cortical pyramidal and hippocampal granule cells. Chl1 (receptor, a neural adhesion molecule required for inhibitory synapse maintenance) is concurrently downregulated in most targets. Together, this suggests that the balance between Pvalb-driven inhibitory signaling and hippocampal CA2 excitatory modulation is shifted in CORT, with reduced Chl1 in recipient neurons impairing the synaptic scaffolding needed for inhibitory transmission.

Pvalb interneuron disruption is among the strongest neurobiological hypotheses for the GABAergic component of PPD and depression more broadly.

### E/I balance — additional signals

From CCC at p < 0.1: Nxph1→Nrxn1/2 (Neurexophilin-Neurexin), a signaling pair that regulates inhibitory synapse formation, is enriched as a biological theme across multiple cell-type pairs, suggesting broader disruption of GABAergic synapse maintenance beyond the individual significant hits.

### Broader hippocampal circuit context

The combination of GR desensitization (Theme 1), reduced BDNF (Theme 2), circadian-disrupted surrounding glia (Theme 3), and neuroinflammatory sensitization (Theme 4) collectively impair the functional state of hippocampal neurons. The E/I disruption theme sits at the output of all these upstream changes: a hippocampus with impaired GR negative feedback, reduced neurotrophin support, disrupted clock entrainment, and increased inflammatory tone is expected to have dysregulated excitatory-inhibitory dynamics as a consequence.

---

## Cross-cutting Summary

| Theme | Primary DEG evidence | CCC evidence | Glial response direction | Hippocampal neuron direction |
|---|---|---|---|---|
| 1. GR desensitization / HPA | Nr3c1↓, Nr3c2↓, Sgk1↓ in CA1/DG | Cntn4→Ptprg (suppresses TrkB) | GR targets massively up | **opposite** — receptors down |
| 2. Neurotrophic collapse | Bdnf↓ in DG, Ptn↓ in Oligo | Bdnf→Sort1 (pro-apoptotic), Ptn↓ CCC | Ntrk2 up (compensatory) | **convergent** — both depleted |
| 3. Circadian disruption | Per2/Per1/Per3/Cry1 in glia | (no CCC-specific) | Clock genes up (phase-shifted) | Not directly changed; circuit effect |
| 4. Neuroinflammatory sensitization | Il6r up in CA1/CA3/DG | Bdnf→Sort1, Mif in microglia | Anti-inflammatory (Nfkbia, GILZ) | **opposite** — pro-inflammatory receptor up |
| 5. E/I balance disruption | Uchl1↓ (Sst), Sgk1↓ (DG) | Fstl5^Chl1 (Pvalb↑, CA2↓), Nxph1-Nrxn1/2 | Homogeneous stress response | Interneuron + hippocampal-specific |

The consistent pattern across themes 1 and 4 is **opposing directionality**: genes that are strongly upregulated in glial cells under CORT (GR target cascade, anti-inflammatory programs) are absent or downregulated in hippocampal neurons. This is not noise — it appears to reflect a fundamental difference in how glial cells (which are glucocorticoid-responsive at the transcriptional level) and hippocampal neurons (which downregulate the receptor itself) respond to chronic CORT exposure. The hippocampal neuronal response looks more like receptor desensitization and trophic withdrawal than an active stress response.

The net picture is a hippocampus that has disengaged from the glucocorticoid signaling system, lost trophic support from multiple sources, become sensitized to inflammatory cytokines, is embedded in a circadian-disrupted glial microenvironment, and shows interneuron-specific vulnerabilities that would impair local circuit computation — all consistent with the cognitive and affective phenotypes observed in PPD models and in postpartum women.
