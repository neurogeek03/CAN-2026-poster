# Additional Themes — CORT vs OIL, Postpartum Rat Brain
## Differential Gene Expression + Cell-Cell Communication (Supplement to `hippocampal_neuron_themes.md`)

**Datasets:** Slide-tags spatial transcriptomics (female rats, postpartum day 1–23)  
**Comparison:** Corticosterone-treated (CORT) vs vehicle-treated (OIL)  
**Analyses:** EdgeR LRT (DEG) + LIANA+ (CCC), library size filter 50k  
**Date:** 2026-05-02

---

## Scope

The five themes in `hippocampal_neuron_themes.md` were built from a curated gene set (the heatmap CSV). The three themes below were identified by mining the full EdgeR output (~229 k gene × cell-type rows) and the LIANA+ table (~125 k ligand-receptor pairs) for significant hits (padj < 0.05) not represented in the original five. All gene-level statistics reported here were verified directly from the source files.

---

## Theme 6 — Pan-Glial Pseudohypoxic Reprogramming (Hif3a)

**Core finding:** Hif3a (hypoxia-inducible factor 3α) is massively upregulated in every non-neuronal cell type profiled, with effect sizes and significance levels that rank among the largest in the entire dataset. Hippocampal neurons show no significant Hif3a change, producing another sharp neuron-glia dissociation.

| Cell type | logFC | padj |
|---|---|---|
| Oligo | +8.81 | 2.8e-62 |
| Astro-TE | +3.35 | 3.1e-38 |
| Astro-NT | +3.13 | 4.8e-30 |
| OPC | +6.70 | 2.1e-25 |
| Endo | +8.60 | 3.4e-16 |
| VLMC | +3.67 | 2.7e-7 |
| Microglia | +7.73 | 3.7e-7 |
| CA1-ProS / DG / CA3 | n.s. | — |

This pattern is not represented in the original five themes, which focus on GR targets (Fkbp5, Sgk1, Tsc22d3) as the dominant non-neuronal stress response. Hif3a represents a parallel and distinct arm.

**Mechanism:** Hif3a has a canonical glucocorticoid-response element (GRE) in its promoter; direct GR binding drives its transcription under elevated CORT, independent of oxygen tension. This constitutes a "pseudohypoxic" program: the cellular transcriptional machinery is responding as if oxygen-deprived, even under normoxic conditions. The effect is compounded by co-occurring metabolic reprogramming genes in oligodendrocytes: Bcat1 (+7.04, padj 1.3e-56 — branched-chain amino acid catabolism as alternative fuel), Pdk4 (+4.54, padj 2.6e-44 — pyruvate dehydrogenase kinase, shunting pyruvate away from mitochondrial oxidation), and Slc2a1 (+1.69, padj 5.2e-6 — GLUT1, increased glucose import). In endothelium, Chka (+4.42, padj 3.4e-16) and Lpl (+3.35, padj 1.1e-14) indicate phospholipid remodeling and fatty acid mobilization. Together, these constitute a coordinated metabolic reprogramming of all glial and vascular cells: less reliance on oxidative phosphorylation, more reliance on amino acid catabolism and glycolysis.

**Biological significance:** Oligodendrocytes are the most metabolically active CNS cell type by myelin volume and the most sensitive to bioenergetic failure. Hif3a↑ in the context of Pdk4↑ (which blocks acetyl-CoA entry into the TCA cycle) means oligodendrocytes are receiving less mitochondrial ATP, which could impair myelin maintenance and lipid synthesis — consistent with the Ptn↓ finding from Theme 2 (loss of PTN, an oligodendrocyte-derived neurotrophic factor, may reflect a cell in metabolic stress rather than a targeted transcriptional decision). Endothelial Hif3a↑ with Chka/Lpl↑ is consistent with BBB membrane lipid remodeling under pseudohypoxic stress, which can increase BBB permeability and alter the transcytosis of glucocorticoids and neuroactive compounds into the brain parenchyma.

---

## Theme 7 — Collapse of the Calcium–cAMP–CREB Plasticity Transcription Axis in Hippocampal Neurons

**Core finding:** Multiple components of the calcium-stimulated cAMP → CREB/Nur transcription cascade that drives activity-dependent gene expression and LTP are significantly downregulated in CA1 and DG under CORT. These are not covered in the original five themes, which focus on trophic inputs and receptor expression rather than intracellular second-messenger and transcription factor machinery.

| Gene | Cell type | logFC | padj | Function |
|---|---|---|---|---|
| Adcy8 | CA1-ProS | −1.90 | 2.6e-6 | Adenylate cyclase 8 — calcium-stimulated cAMP production; required for LTP |
| Adcy8 | CA2-FC-IG | −1.85 | 0.008 | Same |
| Nr4a3 (Nor-1) | CA1-ProS | −1.51 | 7.6e-6 | Nur-family nuclear receptor; BDNF/cAMP-inducible transcription factor for plasticity genes |
| Nr4a3 | DG | −2.23 | 0.002 | Same |
| Nr4a3 | L5 ET CTX | −1.87 | 1.7e-4 | Broader cortical involvement |
| Creb5 | DG | −3.87 | 0.004 | CREB family TF; activity-dependent transcription in hippocampus |
| Pde9a | DG | −1.36 | 4.1e-9 | PDE9A — cGMP-specific phosphodiesterase; modulates NO/cGMP → BDNF/TrkB signaling |

**Mechanistic interpretation:** Adcy8 is the calcium-calmodulin-stimulated adenylate cyclase preferentially expressed in hippocampal pyramidal neurons, where it converts synaptic calcium influx (via NMDAR activation) into cAMP that activates PKA → CREB. Its downregulation in CA1 disconnects calcium influx from the downstream plasticity transcription program. Nr4a3 (Nor-1) is a rapidly inducible nuclear receptor that lies immediately downstream of BDNF/TrkB and cAMP/PKA signaling; it directly activates BDNF, Arc, CaMKII, and other plasticity genes. Its loss in both CA1 and DG under CORT — at effect sizes of 1.5–2.2 logFC — represents a downstream amplification of the trophic withdrawal described in Theme 2 (BDNF↓ in DG): even residual BDNF signaling would be less effective because its downstream transcriptional effector is gone. Creb5↓ in DG adds a further layer: reduced CREB activity impairs the induction of many LTP-associated genes.

The Pde9a finding in DG is notable from a pharmacological angle: PDE9A inhibitors are actively investigated as cognitive enhancers, including for depression. Its downregulation in DG could reflect a compensatory attempt to preserve cGMP tone in the face of reduced NO/sGC input, or it could simply reflect the broader transcriptional depression of plasticity-related enzymes.

**Net effect:** Under chronic CORT, hippocampal neurons cannot efficiently convert synaptic activity into gene expression changes. This is a functional consequence distinct from the receptor desensitization (Theme 1) or trophic withdrawal (Theme 2) described earlier — even if GR were re-sensitized and BDNF were restored, the intracellular second-messenger-to-transcription coupling is impaired.

---

## Theme 8 — Loss of Wnt Secretion and Synaptic Adhesion Disruption in the Dentate Gyrus

**Core finding:** DG granule cells downregulate Wls (Wntless), the obligate chaperone for all Wnt ligand secretion, and co-downregulate multiple synaptic adhesion molecules. This represents a failure of autocrine/paracrine Wnt signaling from DG neurons — a mechanism directly relevant to adult hippocampal neurogenesis.

| Gene | Cell type | logFC | padj | Function |
|---|---|---|---|---|
| Wls | DG | −2.22 | 4.8e-5 | Wntless — required for all Wnt ligand lipidation and secretion |
| Wls | CA1-ProS | −2.57 | 7.1e-4 | Broader hippocampal Wnt secretion loss |
| Wls | L2/3 IT CTX | −1.16 | 1.8e-4 | Cortical involvement |
| Pcdh15 | DG | −1.61 | 1.1e-4 | Protocadherin-15 — synaptic hair-bundle adhesion molecule, expressed in DG mossy fibers |
| Pcdh7 | DG | −0.90 | 2.3e-4 | Protocadherin-7 — synaptic adhesion, dendritic arborization |
| Uqcrq | DG | −3.57 | 2.3e-4 | Ubiquinol-cytochrome c reductase subunit — mitochondrial complex III |

**From CCC:** Ntng1→Lrrc4c (DG Glut → Oligo; padj = 0.0046) is downregulated in CORT. Ntng1 (Netrin-G1) is secreted by DG granule cells and binds Lrrc4c (NGL-1) on postsynaptic partners, including oligodendrocytes and glial processes. Its downregulation reflects reduced DG→glia transsynaptic communication in addition to the cell-autonomous Wnt secretion loss.

**Biological significance of Wls loss:** Wls (also called Evi/Sprinter/GPR177) is not a Wnt ligand itself — it is required for the lipid modification (palmitoylation) and vesicular packaging of every secreted Wnt protein. A cell that has lost Wls cannot secrete any Wnt ligand, regardless of how much Wnt mRNA is present. DG granule cells are the primary local source of Wnt3a and Wnt7a in the dentate gyrus, both of which are required for adult neurogenesis (promoting proliferation and differentiation of subgranular zone progenitors via Frizzled/LRP receptors on neural stem cells). Loss of Wls in DG granule cells under CORT would eliminate this local Wnt signal, providing a mechanistic link between elevated CORT and the well-documented suppression of adult hippocampal neurogenesis in stress and depression models. This extends Theme 2 (reduced trophic support) from the BDNF/PTN axes to the Wnt axis as a third independent neurotrophic system converging on failure.

The Uqcrq↓ in DG (mitochondrial complex III, logFC −3.57) adds a metabolic layer consistent with the pseudohypoxic theme (Theme 6): DG granule cells appear to have reduced mitochondrial electron transport capacity, contrasting with the metabolic upregulation response observed in surrounding glia. This neuron-glia dissociation in metabolic gene direction echoes the pattern already described for GR targets and inflammatory receptors.

**Pcdh15 and Pcdh7 context:** Both protocadherins are downregulated specifically in DG (not CA1 or CA3). Pcdh15 is expressed in DG mossy fiber synapses and regulates the synaptic maintenance of the DG→CA3 projection. Pcdh7 organizes dendritic arborization and inhibitory synapse development. Their co-downregulation with Wls in DG — but not in other hippocampal regions — suggests the dentate gyrus is particularly vulnerable to synaptic structural destabilization under CORT.

---

## Relationship to Original Five Themes

| New theme | Extends / relates to | Distinct because |
|---|---|---|
| 6. Glial pseudohypoxia (Hif3a) | Theme 1 (non-neuronal stress response) | Different pathway — Hif3a/metabolic reprogramming vs. GR target cascade; affects every non-neuronal type at extreme effect sizes |
| 7. cAMP/CREB/Nur plasticity collapse | Theme 2 (trophic collapse), Theme 1 (GR desensitization) | Intracellular second-messenger level — even if trophic inputs and GR were restored, the transduction machinery is disrupted |
| 8. Wnt secretion and DG synaptic adhesion loss | Theme 2 (neurotrophic collapse) | Third independent trophic axis; unique to DG; mechanistically linked to adult neurogenesis via Wls |

The dominant structural observation from the original report — that hippocampal neurons respond to CORT by disengaging from glucocorticoid signaling and losing trophic support — is deepened by these findings. Themes 7 and 8 show that hippocampal neurons additionally lose intracellular plasticity transduction machinery (Adcy8, Nr4a3, Creb5) and autocrine trophic secretion capability (Wls). Theme 6 shows that the surrounding glia respond to the same CORT exposure with a metabolic reprogramming of their own — a pseudohypoxic shift that impairs their capacity to support myelination, nutrient delivery, and the neurotrophic functions documented in Theme 2.
