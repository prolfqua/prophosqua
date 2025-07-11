# Methods in Molecular Biology Article (https://link.springer.com/series/7651)

# Integrated Analysis of Post-Translational Modifications and Total Proteome: Methods for Distinguishing Expression from Usage Changes

## **Abstract**

**Background:** Post-translational modifications (PTMs), particularly phosphorylation, regulate protein function at substoichiometric levels, making their quantitative analysis technically challenging. A critical limitation in PTM studies is distinguishing between changes in modification abundance due to altered protein expression versus genuine changes in modification stoichiometry (the fraction of protein molecules that are modified).

**Methods:** This chapter presents a comprehensive workflow integrating TMT-based quantitative proteomics with specialized computational analysis to distinguish differential PTM expression (DPE) from differential PTM usage (DPU). The protocol encompasses automated sample preparation using SP3-based digestion, phosphopeptide enrichment with Ti-IMAC, high-pH offline fractionation, LC-MS/MS analysis, and integrated statistical analysis using the `prolfquapp` and `prophosqua` R packages.

**Key Features:** The workflow provides automated, scalable protocols for large sample cohorts, robust statistical frameworks for integrated PTM-protein analysis, comprehensive visualization tools including N-to-C plots for protein-centric PTM mapping, and sequence motif analysis for kinase prediction. All protocols are optimized for reproducibility and include extensive quality control measures.

**Applications:** This integrated approach enables researchers to identify true signaling changes in PTM studies, prioritize functionally relevant modification sites, understand pathway-level regulation in disease contexts, and generate robust datasets suitable for publication in high-impact journals.

# 1. Introduction

## 1.1 Biological Significance of Post-Translational Modifications

Protein phosphorylation represents a reversible post-translational modification (PTM) that modulates protein conformation, activity, localization, and interactions, functioning as a critical molecular switch in cellular signaling networks ([10.1080/14789450.2021.1976152](https://doi.org/10.1080/14789450.2021.1976152)). This dynamic modification system enables cells to respond rapidly to environmental changes and coordinate complex biological processes without requiring new protein synthesis. Dysregulated phosphorylation underlies numerous diseases, including cancer, neurodegeneration, and metabolic disorders, making large-scale phosphoproteomic profiling by mass spectrometry essential for understanding disease mechanisms and identifying therapeutic targets.

## 1.2 Technical Challenges in PTM Analysis

The analysis of post-translational modifications presents significant analytical challenges due to the low stoichiometry of phosphopeptides, which typically represent less than 1-5% of total peptide abundance ([10.1039/C5MB00024F](https://doi.org/10.1039/C5MB00024F); [10.1186/1477-5956-4-15](https://doi.org/10.1186/1477-5956-4-15)). This substoichiometric nature necessitates specialized workflows for quantification and creates substantial technical hurdles for PTM analysis ([10.1016/j.mcpro.2024.100044](https://doi.org/10.1016/j.mcpro.2024.100044)). Traditional proteomics approaches often fail to distinguish between changes in PTM abundance that result from altered protein expression versus genuine changes in modification efficiency or stoichiometry ([10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477); [10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708)). This confounding factor has led to widespread misinterpretation of PTM data, where apparent increases in phosphorylation may simply reflect higher protein abundance rather than enhanced kinase activity or altered signaling states ([10.1039/C5MB00024F](https://doi.org/10.1039/C5MB00024F)).

Furthermore, PTM datasets typically exhibit higher rates of missing values compared to total proteome data, with quantified phosphosites often representing only 60-75% of detected sites prior to imputation ([10.1101/2020.08.31.276329](https://doi.org/10.1101/2020.08.31.276329)), due to the substoichiometric nature of modifications and the requirement for enrichment strategies that can introduce technical variability ([10.1039/C5MB00024F](https://doi.org/10.1039/C5MB00024F); [10.1016/j.aca.2021.338716](https://doi.org/10.1016/j.aca.2021.338716)). Each stage of the phosphoproteomics workflow contributes to overall data variability and must be optimized thoroughly to ensure reproducible quantification ([10.1038/s41467-018-03309-6](https://doi.org/10.1038/s41467-018-03309-6); [10.1016/j.aca.2021.338716](https://doi.org/10.1016/j.aca.2021.338716)). These analytical challenges demand experimental design and computational approaches to extract biologically meaningful information from complex PTM datasets.

## 1.3 Current Methodological Approaches

### 1.3.1 Mass Spectrometry-Based Quantification Strategies

Isobaric tagging with TMT reagents has emerged as the gold standard for multiplexed PTM analysis, enabling simultaneous analysis of up to 35 samples in a single combined run and significantly reducing inter-sample technical variability ([10.1021/acs.jproteome.4c00668](https://doi.org/10.1021/acs.jproteome.4c00668)). This multiplexing approach provides enhanced statistical power while minimizing batch effects that can confound biological interpretation, creating closed systems with unique statistical properties that allow for greater confidence in comparative results ([10.1021/ac0262560](https://doi.org/10.1021/ac0262560)). Enrichment of phosphorylated peptides is typically achieved through metal affinity chromatography (IMAC) or metal oxide affinity chromatography (MOAC), with Fe³⁺, Zr⁴⁺, Ti⁴⁺ ions, or TiO₂ providing complementary selectivity profiles that enable coverage of the phosphoproteome ([10.1074/mcp.M114.045609](https://doi.org/10.1074/mcp.M114.045609); [10.1101/2020.04.13.038810](https://doi.org/10.1101/2020.04.13.038810)).

### 1.3.2 Computational Tools and Workflows

The computational landscape for quantitative protein post-translational modification analysis has evolved significantly, with numerous specialized software platforms supporting bottom-up proteomics workflows. The analysis of bottom-up proteomics data involves several critical steps that must be executed sequentially to achieve accurate PTM identification and quantification. Following mass spectrometry data acquisition, raw spectra undergo database searching to identify peptide sequences, which is followed by protein inference to determine the most likely protein identities from the identified peptides. For PTM analysis, site localization scoring represents a particularly crucial step, as it determines the confidence with which a modification can be assigned to a specific amino acid residue within a peptide sequence. This process becomes especially challenging when multiple potential modification sites exist within a single peptide, requiring algorithms to evaluate the quality of site-determining fragment ions. Specialized tools for PTM site localization scoring provide critical functionality for confident site assignment, including PTMProphet ([10.1038/s41467-020-17914-x](https://doi.org/10.1038/s41467-020-17914-x)), PhosphoRS ([10.1038/nmeth.1107](https://doi.org/10.1038/nmeth.1107)), and Ascore ([10.1038/nmeth.1107](https://doi.org/10.1038/nmeth.1107)).

The computational ecosystem includes several established DDA-TMT compatible software suites, encompassing both free and commercial options. These include Andromeda/MaxQuant (Cox et al., 2011), Proteome Discoverer (Thermo Fisher Scientific), FragPipe ([10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7)) (fragpipe.nesvilab.org), and PeptideShaker (peptide-shaker.compomics.com).

For data-dependent acquisition (DDA) analysis, platforms like FragPipe ([10.1101/2025.05.27.656447](https://doi.org/10.1101/2025.05.27.656447)) integrate the MSFragger search engine ([10.1038/nmeth.4256](https://doi.org/10.1038/nmeth.4256)) and are widely adopted for both TMT and label-free workflows. FragPipe provides PTM analysis capabilities, including site localization scoring through integrated PTMProphet functionality and generation of both site-level quantification reports and multi-site feature reports that group peptides sharing identical modification sites.

For data-independent acquisition (DIA) analysis, specialized tools such as Spectronaut (Biognosys) and DIA-NN (Demichev et al., Nat Methods 2020, [10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x)) offer PTM support with site-localization scoring algorithms and site-level quantification capabilities. Both platforms implement statistical frameworks for PTM identification and quantification in DIA workflows.

### 1.3.3 Data Analysis Frameworks

Proteomic data analysis can be conducted at multiple levels of granularity, each providing distinct analytical perspectives. A peptidoform represents a specific peptide sequence with a particular set of modifications, while site-level analysis aggregates signals for each residue position across multiple peptides. Differential analysis can therefore target either peptidoforms or individual modification sites, with each approach offering unique advantages for different research questions.

FragPipe introduces an additional analytical concept through multisite features, which refers to the set of all identified peptide-forms sharing the same set of modification sites under investigation, regardless of peptide length, sequence derivatives, cleavage state, or whether the modifications are unambiguously localized. This approach groups peptides with different sequences, which can arise from missed cleavages or alternative cleavage patterns, but identical possible modification sites into the same multisite feature, providing a balanced approach between peptidoform specificity and site-level aggregation. Discarding all peptideforms with ambiguous localization would result in the loss of a large fraction of PTM information. Multisite features preserve information about phosphorylation sites that occur together on the same peptide molecule, revealing co-modification patterns that indicate coordinated regulation by kinases or functional relationships between sites. Traditional site-level analysis treats each phosphorylation site independently, losing this valuable information about which specific combinations of sites are actually modified together in biological samples.

Site-level reports collapse signals from all peptidoforms mapping to the same PTM at an amino acid residue into a single quantitative value, ensuring each site has exactly one intensity measurement per sample. In downstream analysis, site-level intensities integrate seamlessly with site-centric enrichment tools and kinase-activity inference platforms (PhosR ([10.1016/j.cels.2021.04.007](https://doi.org/10.1016/j.cels.2021.04.007)), PTM-SEA ([10.1016/j.molcel.2019.05.030](https://doi.org/10.1016/j.molcel.2019.05.030)), RoKAI ([10.1038/s41467-021-21211-6](https://doi.org/10.1038/s41467-021-21211-6)), etc.), streamlining biological interpretation and network reconstruction analyses ([10.1186/s12014-020-09290-x](https://doi.org/10.1186/s12014-020-09290-x)).

### 1.3.4 Statistical Analysis of PTM Data

Statistical analysis of PTM data requires specialized approaches that account for the unique characteristics of modification data. Several statistical packages address the critical challenge of PTM-protein integration, including MSstatsPTM ([10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)) and msqrob2PTM ([10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708)), which model both modified and unmodified peptides to properly account for protein-level changes and distinguish genuine PTM regulation from protein abundance effects.

These two frameworks implement different analytical strategies for addressing the confounding between PTM and protein abundance changes using log-transformed intensities. MSstatsPTM employs a "model then correct" approach, separately fitting linear models to modified and unmodified peptide data, then combining the resulting statistical inferences to distinguish between absolute PTM changes and changes adjusted for protein abundance. In contrast, msqrob2PTM follows a "correct then model" strategy, first normalizing log-transformed PTM intensities by their corresponding log-transformed protein abundances at the data level, then applying differential analysis to the protein-corrected values. This methodological distinction does not lead to different statistical properties or interpretations, as both approaches ultimately estimate PTM changes adjusted for protein abundance effects.

## 1.4 Innovation of the Integrated Approach

### 1.4.1 Distinguishing PTM Expression from PTM Usage

The msqrob2PTM framework exemplifies the statistical approaches now available for PTM analysis, defining two complementary analytical strategies: Differential Peptidoform Abundance (DPA) and Differential Peptidoform Usage (DPU) ([10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708)). DPA directly models the log₂ intensities of each modified peptide (peptidoform), detecting absolute changes in PTM levels between experimental conditions. In contrast, DPU adjusts the PTM intensities by the corresponding protein-level changes, effectively testing for changes in the relative usage or stoichiometry of a modification site. This dual approach enables researchers to distinguish between PTM changes driven by protein abundance differences (DPA) versus those reflecting genuine changes in modification efficiency or regulatory activity (DPU).

### 1.4.2 Workflow Integration

This chapter introduces an analytical workflow that addresses the challenge of distinguishing Differential PTM Expression (DPE) from Differential PTM Usage (DPU), providing researchers with the tools necessary for accurate biological interpretation of PTM data.

**Differential PTM Expression (DPE)** tests raw PTM signal changes between experimental conditions, identifying any modification whose abundance changes regardless of underlying protein-level effects. This analysis captures the total effect of experimental perturbations on PTM levels, including both direct modification effects and indirect effects mediated through protein abundance changes.

**Differential PTM Usage (DPU)** tests protein-normalized PTM changes, revealing modification sites where stoichiometry genuinely changes independent of protein abundance alterations. This analysis specifically identifies sites where modification efficiency, kinase activity, or phosphatase activity has been altered by experimental conditions.

The integrated analytical approach combines optimized sample preparation protocols, mass spectrometry acquisition strategies, and computational tools (including the `prolfquapp` and `prophosqua` R packages) to provide a complete solution spanning from sample processing to biological interpretation. This framework enables researchers to extract maximum biological insight from PTM datasets while avoiding common analytical pitfalls that have historically complicated PTM data interpretation.

## 1.5 Chapter Overview and Learning Objectives

This protocol demonstrates the complete workflow using the Atg16l1 macrophage dataset (Maculins *et al.*, eLife 2021, [10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320)), which comprises TMT-11-plex measurements from six conditions (WT/KO × uninfected/early/late infection) with both total proteome and phospho-enriched samples. 

**Learning objectives include:**
1. Implementing automated, scalable sample preparation protocols
2. Performing integrated statistical analysis of PTM and protein data
3. Interpreting DPE vs. DPU results for biological insights
4. Generating publication-quality visualizations and reports
5. Conducting sequence motif analysis for kinase prediction

The workflow emphasizes reproducibility, automation, and biological interpretation, making it suitable for both experienced researchers and newcomers to integrated PTM analysis.

# 2. Materials

## 2.1 Biological Samples

As described above, we are focusing on the dataset from Maculins et al. ([10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320)) in this section. We are using the phospho-enriched and the total proteome datasets for demonstration purposes. In the original publication, the authors also performed a KGG-enrichment (ubiquitin remnant). Furthermore, the same dataset was also used in the MSstatsPTM publication (MCP, Kohler et. al, 2023, [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)). 

**Representative sample types compatible with this protocol:**
- Fresh-frozen tissue samples (≤5 mg per sample)
- Cultured cell pellets (≥1×10⁶ cells per sample)
- Primary cells isolated from tissue
- Organoids and spheroid cultures

**Sample storage requirements:**
- Store at -80°C immediately after collection
- Avoid freeze-thaw cycles
- Process within 6 months of collection for optimal results

The material and methods that were used in this publication are not identical to what we suggest here. Nevertheless, the most important parts are that chemical labeling was used and there is a PTM-enriched dataset with potentially multiple fractions (phospho-enriched here) as well as a total proteome part that can be integrated. The main differences in sample preparation between the proposed method and the method used by Maculins et al. are summarized in the table. 

| Step | Maculins et al. ([10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320))| Current Protocol (This Work) |
| ----- | ----- | ----- |
| **Lysis & Digestion** | Manual lysis, SDS buffer; trypsin digestion | SDS lysis; automated SP3 digestion on KingFisher Flex |
| **Protein Cleanup** | Acetone or methanol/chloroform precipitation | Automated SP3 bead-based cleanup |
| **TMT Labeling Timing** | Post-enrichment (for KGG peptides); pre-enrichment (flow-through) | Pre-enrichment TMTpro labeling for all samples |
| **Enrichment Targets** | KGG (ubiquitin remnant), phosphotyrosine (pY), phosphoserine/threonine (pST) | Optional pY or motif antibody enrichment \+ Fe-NTA for pST |
| **Enrichment Method** | Manual, multiple kits (e.g., CST pY1000, Thermo Fe-NTA) | Automated enrichment using KingFisher Flex with BindIt control |
| **Fractionation** | Manual StageTips or high-pH kit, 6–18 fractions | Automated high-pH fractionation on Vanquish Flex, concatenated into 36 fractions |
| **Peptide Elution & Loading** | Manual cleanup, reconstitution for LC-MS | Evotip-compatible elution and direct loading |
| **Automation** | Mostly manual steps | Fully automated digestion, enrichment, and fractionation |
| **Scalability & Reproducibility** | Moderate, sample-specific variability | High-throughput, reproducible across large sample cohorts |

## 2.2 Sample Preparation

### 2.2.1 Cell/tissue lysis, digestion, and labeling of tryptic digests with TMT reagents

1. **Lysis buffer:** 50 mM Tris-HCl, pH 7.5, 4% (w/v) SDS, optional: PhosSTOP (Merck, Germany, Cat. #4906845001; store at -20°C, use within 2 years)
2. **Nuclease (optional):** 5 U/µL Benzonase Nuclease (Merck, Germany, Cat. #70664-3; store at -20°C, stable for 2 years)
3. **Tissue homogenization:** Tissue Lyzer II (Qiagen, Germany, Cat. #85300) with 5 mm stainless steel beads (Qiagen, Cat. #69989)
4. **Sonication:** Probe-type sonicator UP200St (Hielscher, Germany, Cat. #UP200St) with 2 mm probe
5. **Mixing:** Orbital benchtop mixer for microtubes, e.g., Eppendorf Thermomixer C (Eppendorf, Germany, Cat. #5382000015; temperature range: 4-99°C, shaking speed: 300-3000 rpm)
6. **Protein quantification:** UV/Vis spectrophotometer, e.g., Lunatic (Unchained Labs, USA, Cat. #300200; measurement range: 0.03-275 mg/mL)
7. **Magnetic beads:** Speedmag carboxylated magnetic beads, hydrophilic and hydrophobic (Merck, Germany, Cat. #09-981-121 and #09-981-123; store at 4°C, stable for 2 years at 50 mg/mL stock concentration)
8. **Pipetting:** Multichannel pipette (8 or 12 channel, 10-300 µL range)
9. **Centrifugation:** Microtube centrifuge, e.g., Eppendorf 5424 (Eppendorf, Germany, Cat. #5424000010; max. speed: 21,130×g)
10. **Automation platform:** KingFisher Flex with deep-well head (Thermo, USA, Cat. #5400610; includes magnetic head and plate carriers)
11. **Plasticware:** KingFisher deep-well plate (Thermo, Cat. #95040450), microwell plate (Thermo, Cat. #95040210), tip comb (Thermo, Cat. #94460205)
12. **Automation software:** BindIT 4.0 Software (Thermo, included with KingFisher)
13. **Reducing agent:** 250 mM tris(2-carboxyethyl)phosphine (TCEP) in water (Sigma-Aldrich, Cat. #C4706; prepare fresh or store at -20°C for up to 6 months)
14. **Alkylating agent:** 500 mM chloroacetamide (ClAA) in water (Sigma-Aldrich, Cat. #C0267; prepare fresh or store at -20°C for up to 3 months, light-sensitive)
15. **Protease:** Trypsin, mass spectrometry grade (Promega, Cat. #V5280; store at -20°C, stable for 5 years)
16. **Solvents:** 80% and 100% ethanol (EtOH), LC-grade (Fisher Chemical, Cat. #E/0650DF/17)
17. **Water:** LC-grade water (Fisher Chemical, Cat. #W/0106/17)
18. **Organic solvent:** Anhydrous acetonitrile (Sigma-Aldrich, Cat. #271004; store at room temperature, water content <0.001%)
19. **Buffer:** 50 mM Tetraethylammonium bicarbonate (TEAB), pH 8.5 (Sigma-Aldrich, Cat. #T7408; prepare fresh or store at 4°C for up to 1 month)
20. **Quenching reagent:** Stop solution: 5% hydroxylamine (Thermo, Cat. #90115; prepare fresh in water)
21. **Isobaric labels:** TMTpro reagent, e.g., 18plex (Thermo, USA, Cat. #A55256; store at -20°C with desiccant, stable for 2 years)
22. **Concentration:** SpeedVac concentrator (Thermo, Cat. #SPD210-115; compatible with microtubes and deep-well plates)

### 2.2.2 Optional antibody-based enrichment

1. **Antibody:** PTMscan antibody targeting site and/or kinase substrate motif of interest, e.g., PTMScan HS Phospho-Tyrosine (Cell Signaling, USA, Cat. #8954; store at -20°C, stable for 1 year)
2. **Secondary enrichment:** High Select Fe-NTA Phosphopeptide Enrichment Kit (Thermo, USA, Cat. #A32992; store at 4°C, use within expiration date)
3. **Binding buffer:** IP binding buffer: 50 mM Tris-HCl, pH 7.5, 1% NP40 (IGEPAL CA-630) (Sigma-Aldrich, Cat. #I8896)
4. **Wash buffer:** IP wash buffer: 50 mM Tris-HCl, pH 7.5 (prepare fresh or store at 4°C for up to 1 week)
5. **Mixing:** End-over-end shaker (e.g., Thermo, Cat. #88881001; speed: 5-60 rpm)

### 2.2.3 High pH offline fractionation

1. **UHPLC system:** Vanquish Flex UHPLC with fraction collector and 100 µL sample loop (Thermo, Cat. #VF-S01-A; includes autosampler, pump, column compartment, and fraction collector)
2. **Software:** Chromeleon 7.3.1 instrument control software (Thermo, included with system)
3. **Standard:** BSA peptide standard: 500 μmol dry, diluted in 600 μL water to 83.3 μmol final concentration for 100 μL injection (Sigma-Aldrich, Cat. #A7906)
4. **Mobile phase A:** High pH buffer A: 2.5 mM ammonium hydroxide in LC-grade water (prepare fresh daily)
5. **Mobile phase B:** High pH buffer B: 95% acetonitrile/10 mM ammonium hydroxide in LC-grade water (prepare fresh daily)
6. **Separation column:** XBridge peptide BEH C18 130Å, 3.5 µm, 4.6 mm × 250 mm (Waters, Cat. #186003115; store at room temperature, avoid extreme pH)
7. **Storage solution:** 20% methanol in water (for column storage)
8. **Concentration:** SpeedVac compatible with deep-well plates, e.g., Savant SpeedVac Vacuum Concentrator SPD210 (Thermo, Cat. #SPD210-115)

### 2.2.4 Phosphopeptide enrichment and Evotip loading

1. **Automation platform:** KingFisher Flex as specified in 2.2.1
2. **Binding buffer:** Ti-IMAC buffer I: 80% acetonitrile, 5% trifluoroacetic acid, 0.1 M glycolic acid (prepare fresh daily)
3. **Wash buffer 1:** Ti-IMAC buffer II: 80% acetonitrile, 1% trifluoroacetic acid (prepare fresh daily)
4. **Wash buffer 2:** Ti-IMAC buffer III: 10% acetonitrile, 0.2% trifluoroacetic acid (prepare fresh daily)
5. **Elution buffer:** Ti-IMAC elution buffer: 1% ammonium hydroxide (prepare fresh)
6. **Neutralization:** Ti-IMAC neutralization buffer: 20% formic acid (Sigma-Aldrich, Cat. #F0507)
7. **Affinity media:** ReSyn Biosciences Ti-IMAC HP beads (ReSyn Biosciences, Cat. #R-TI.HP.1; store at 4°C as 20% slurry, stable for 1 year)
8. **Sample loading:** Evotip pure (Evosep, Denmark, Cat. #EV-1005; single-use, store at room temperature)

### 2.2.5 LC-MS/MS analysis

1. **LC system:** Evosep One (Evosep, Denmark, Cat. #EV-1000; includes autosampler and gradient pump)
2. **Mobile phase A:** Buffer A: 0.1% formic acid in LC-grade water (prepare fresh weekly)
3. **Mobile phase B:** Buffer B: 99.9% acetonitrile/0.1% formic acid (prepare fresh weekly)
4. **Analytical column:** Aurora Elite XT 15×75 C18 UHPLC column (IonOpticks, Australia, Cat. #AUR3-15075-EXT; store at room temperature)
5. **Mass spectrometer:** Orbitrap Exploris 480 with EASY-Spray Source (Thermo, USA, Cat. #IQLAAEGAAPFADBMBHQ; includes ion source and control software)
6. **Column heating:** Column heater and heater controller (IonOpticks, Australia, Cat. #AUR-COL-HEATER; temperature range: ambient to 70°C)

### 2.2.6 Data analysis

1. **Search software:** FragPipe 22.0 (free download from [fragpipe.nesvilab.org](https://fragpipe.nesvilab.org); system requirements: 16GB+ RAM, 50GB+ storage)
2. **Statistical platform:** R version ≥4.0.0 ([https://www.r-project.org/](https://www.r-project.org/)) and RStudio ([https://posit.co/](https://posit.co/))
3. **R packages:**
   - `prolfquapp` ([10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911)): `install.packages("prolfquapp")`
   - `prolfqua` ([10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441)): dependency of prolfquapp
   - `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)): `devtools::install_github("fgcz/prophosqua")`
   - Additional dependencies: `tidyverse`, `ggseqlogo`, `writexl`, `rmarkdown`

**System requirements:**
- **Minimum:** 16GB RAM, 100GB free storage, 4-core processor
- **Recommended:** 32GB+ RAM, 500GB+ SSD storage, 8+ core processor
- **Operating system:** Windows 10+, macOS 10.15+, or Linux Ubuntu 18.04+

# 3. Methods

## 3.1 Sample Preparation

### 3.1.1 Cell lysis, digestion, and labeling of tryptic digests with TMT reagents
**Expected time:** 2 days (Day 1: lysis to digestion setup, Day 2: TMT labeling)
**Expected yield:** 80-95% protein recovery, >95% TMT labeling efficiency

This protocol describes a generic lysis approach for fresh-frozen tissue or frozen cell pellets, utilizing denaturing SDS for enzymatic deactivation and efficient membrane disruption. It has been successfully used for most eukaryotic systems, including yeast and plants (see Note 4.1). Protein clean-up and digestion are based on the automated SP3 protocol, as described in Leutert et al. ([10.15252/msb.20199021](https://doi.org/10.15252/msb.20199021)). The protocol efficiently removes all components incompatible with downstream processing steps, including buffers containing primary amines.

**Day 1: Cell lysis and protein extraction**

1. Add up to 50 μL lysis buffer for small pellets (≤10 μL pellet volume) and more for larger pellets, using 3×–5× the pellet volume and completely resuspend and disrupt cells in buffer by pipetting up and down; apply 1 min Tissue Lyzer (30 Hz) cycles with glass beads if cell pellet is not completely disrupted and lysed (**Critical:** Complete cell disruption is essential for protein extraction efficiency)

2. For tissue samples, add 300 μL lysis buffer to up to 5 mg of tissue (optimal volume depends on tissue origin and protein content) and disrupt the sample using the Tissue Lyzer II with 2× 1 min cycles (30 Hz) (**Time:** 5 min total homogenization)

3. Apply a 1-minute sonication pulse to the sample (probe-type sonicator) (**Critical:** Use 50% amplitude to avoid overheating)

4. Incubate 5 min at 95°C (**Critical:** This step inactivates phosphatases and kinases)

5. Apply a 1-minute sonication pulse to the sample (probe-type sonicator)

6. If viscosity is still high due to high DNA content, dilute the sample to 1% SDS and add 5–10 U of Benzonase and incubate for 15 min at room temperature (**Quality check:** Sample should flow easily when pipetted)

7. Centrifuge 5 min at 20,000×g and proceed with supernatant (**Expected yield:** >90% of total protein in supernatant)

8. Take a small aliquot for protein concentration determination and dilute at least 1:20 for 280 nm reading on the Lunatic (**Expected concentration:** 1-5 mg/mL for cell pellets, 0.5-2 mg/mL for tissue)

9. Snap-freeze samples in liquid nitrogen and store at -20°C or proceed with the protocol for 30–50 µg of total protein with 0.5–2 µg/µL (**Stopping point:** Samples can be stored for up to 6 months)

**Reduction and alkylation:**

10. Add TCEP to a final concentration of 5 mM (**Time:** Immediate)

11. Add ClAA to a final concentration of 15 mM (**Critical:** Add in the dark, ClAA is light-sensitive)

12. Incubate 30 min at 30°C in the dark, shaking at 500 rpm (**Quality check:** No precipitation should be visible)

**SP3 bead preparation and protein binding:**

13. Prepare carboxylated magnetic beads (see Note 4.2)
    1. Use 10 μg of beads for each μg of protein (**Critical ratio for efficient binding**)
    2. Beads are available at a stock concentration of 50 μg/μL. Take the required volume of beads and mix hydrophilic and hydrophobic beads at a 1:1 ratio
    3. Wash beads three times with water and add the respective amount in 250 µL to each well

14. Add 100% EtOH to the reduced and alkylated sample to reach 60% EtOH (v/v) (1.5× of reduced and alkylated sample volume) (**Critical:** Add slowly while mixing to prevent precipitation)

15. Fill three deep-well plates with 500 µL wash buffer (80% EtOH) (**Preparation time:** 10 min)

16. Add 100 µL of trypsin in 20 mM TEAB at an enzyme:protein ratio of 1:50 to a microwell plate (**Critical:** Use MS-grade trypsin)

**Automated digestion:**

17. Process samples with KingFisher protocol as described in Leutert et al. with minor changes (**Time:** 2 hours automated processing)
    1. After transferring beads into the trypsin digestion buffer, pause the KingFisher protocol, transfer the plate to the Thermomixer heated to 37°C, and incubate overnight while shaking at 500 rpm (see Note 4.3) (**Day 1 stopping point:** 16-18 hours digestion)
    2. The next day, return the digestion plate to its original KingFisher slot, add 100 µL water to a second microwell plate, and place it in the empty KingFisher slot
    3. Continue with the KingFisher protocol for bead collection and transfer into the water elution plate (**Time:** 30 min)
    4. Pool digest, water elution, and dry samples in SpeedVac (**Time:** 2-4 hours depending on volume)

**Day 2: TMT labeling**

18. Resuspend peptides in 45 µL 50 mM TEAB and sonicate in a water bath for 10 min (**Expected recovery:** >90% of starting peptides)

19. Resuspend 100-200 µg TMTpro reagent in 15 µL anhydrous acetonitrile, vortex, and spin down (see Note 4.4) (**Critical:** Use immediately after resuspension)

20. Add TMTpro reagents to the respective sample, mix, quickly spin down, and incubate for 60 min at RT in the dark (**Expected labeling efficiency:** >95%)

**TMT labeling quality control and sample pooling:**

21. Add 3.5 µL of 5% hydroxylamine to quench the reaction (**Time:** Immediate quenching essential)

22. Pool a small fraction of each sample (e.g., 1 μL, ensure precise pipetting) and pool at equal ratio to check for sample loading variation and labeling efficiency using LC-MS/MS in QC step (**Critical QC step:** Assess before proceeding)

23. If needed, adjust the pooling ratio based on QC LC-MS/MS, otherwise, proceed with pooling all samples at an equal ratio (**Expected CV:** <20% across channels)

24. Optional: if antibody-based enrichment is required, split the sample in half

25. Dry down in SpeedVac (**Stopping point:** can be stored frozen at -80°C for up to 1 year)

26. Optional: For final clean-up, subject sample to C18 solid-phase extraction (SPE) (e.g., Sep-Pak, Waters) (**Recovery:** >85% of peptides)

**Expected outcomes:**
- **Total peptide yield:** 20-40 μg per 50 μg input protein
- **TMT labeling efficiency:** >95% (verify by LC-MS/MS)
- **Sample purity:** Suitable for downstream enrichment and fractionation
- **Storage stability:** 1 year at -80°C without significant degradation

### 3.1.2 Optional antibody-based enrichment

This protocol describes the targeted enrichment of peptides harboring specific kinase substrate motifs or phosphorylated tyrosines using PTMScan HS antibody beads from Cell Signaling Technology. Binding and wash buffers have been modified, and an additional second enrichment step using Fe-NTA spin columns is used to further increase phospho-peptide enrichment specificity (see Note 4.5).

1. Resuspend one half of the pooled labeled peptides in 400 µL cold IP binding buffer, vortex, and incubate rotating at 4 °C for 5 min, check pH (should be neutral or slightly basic, not below pH 7\)  
2. Spin down solution 5 min at 20'000 x g to clear the solution and pellet insoluble material  
3. Gently mix the antibody bead slurry to obtain a uniform bead suspension and transfer 10 µL of the antibody slurry (see Note 4.6) into an Eppendorf tube.   
4. Wash the beads with 500 µL of cold IP binding buffer by inverting the tube to resuspend the beads in the buffer, and then remove the buffer by placing the tube on a magnetic rack. Repeat for a total of four washes   
5. Transfer the soluble peptide solution to the washed beads  
6. Incubate on an end-over-end shaker for 2 hours at 4 °C (see Note 4.7)  
7. Collect buffer solution by spinning at no more than 1000 x g for 4-5 seconds  
8. Collect beads on magnetic rack and transfer supernatant into a new tube (can be saved and stored at \-80°C)  
9. Wash the beads with 400 µL of cold I**P binding buffer** by inverting the tube to resuspend the beads. Collect the beads on a magnetic rack and discard the buffer. Repeat for a total of two washes  
10. Wash the beads with 400 µL of cold **IP wash buffer** by inverting the tube to resuspend the beads. Collect the beads on a magnetic rack and discard the buffer. Repeat for a total of two washes  
11. Wash the beads with 400 µL **of cold water** by inverting the tube to resuspend the beads. Collect the beads on a magnetic rack and discard the water. Repeat for a total of two washes  
12. Add 50 µL High Select Fe-NTA binding buffer and incubate beads for 10 min at room temperature on a shaker or Thermomixer. Ensure beads stay in suspension but are not splashed to the sides of the tube  
13. Collect beads and transfer elution to High Select Fe-NTA spin column conditioned as per the manufacturer's instructions  
14. Repeat the elution step and combine with the first elution in the same spin column  
15. Follow the manufacturer's instructions for binding and washing routines  
16. Add 30 µL elution buffer and collect elution into a fresh Eppendorf tube. Repeat elution once into the same collection tube  
17. Dry peptides down to near-completeness, leaving \~5 \-10 µL of liquid in the tube if continuing with data acquisition immediately. Otherwise, dry and store peptides at \-80 °C.  
18. Acidify with 40 µL buffer A, check pH (should be at or below pH 4, see section 3.1.4)

    

### 3.1.3 High pH offline fractionation

This protocol describes the general flow of events for high-pH offline fractionation of peptides on the Thermo Vanquish Flex System, featuring automatic fraction collection and concatenation for subsequent phospho-peptide enrichment on the KingFisher Flex. 

1. For robust and sensitive performance of phosphoproteomics experiments, we typically aim for 10-25 µg of peptide per fraction. These amounts have been demonstrated to be effective for automated protocols on the KingFisher Flex ([10.1016/j.mcpro.2024.100754](https://doi.org/10.1016/j.mcpro.2024.100754), [10.1002/pmic.202100245](https://doi.org/10.1002/pmic.202100245)). We generate 36 fractions, corresponding to a total of approximately 180-540 µg of pooled peptide before fractionation.  
2. Resuspend the peptides in 100 µL of high-pH buffer A, and sonicate for 5 min to ensure complete resuspension. Pellet insoluble material and only transfer supernatant into the HPLC vial  
3. Purge and equilibrate HPLC with high pH buffer A and 2.1% high pH buffer B)  
4. Connect XBridge peptide BEH column and adjust the flow rate to 0.75 mL/min, column compartment heater to 40 °C, and equilibrate the system until pressure and UV readings stabilize (see Note 4.7)  
5. Place the KingFisher deep-well plate in the fraction collector (FC) and check the 'Use Safe Needle Height' setting in Chromeleon FC General Settings.  
6. Set up a sample queue in Chromeleon and add the custom variable specifying the fractionation starting position (needs to be specified in the Script editor)  
7. Always run QC without fractionation using BSA peptide standard and check the elution profile for proper peak separation before an actual sample.  
8. For sample fractionation, adjust the FC.FractionRange in the Script editor to 36 and the MaxTubesPerFraction to 72, collecting each fraction for 60 seconds. The needle will start again at the first position after 36 fractions, i.e. concatenating fraction 37 into fraction 1, fraction 38 into fraction 2, etc. (see Note 4.8)  
9. Set up the flow gradient to run from 2.1% B to 42.1% B in 54 minutes, followed by 5 minutes at 95% B and another 10 minutes equilibrating back to 2.1% B.  
10. After the run, equilibrate the column with 100 % acetonitrile and then return the system to 20 % methanol.  
11. Remove the deep-well plate from the fraction collector, snap-freeze, and dry down the samples in a SpeedVac. Peptides can be stored at \-80°C

    

### 3.1.4 Phosphopeptide enrichment and Evotip loading

1. Resuspend each fraction in 150 µL Ti-IMAC buffer I (see Note 4.10), sonicate for 10 min.   
2. If total proteome measurements are required as well, take a small aliquot of either the input before phospho-peptide binding or the non-bound fraction after phospho-peptide binding  
3. Prepare an appropriate volume of Ti-IMAC HP beads (see Note 4.11) from 20 % bead slurry, aiming for a peptide-to-bead ratio of 1:4   
4. Wash the beads with 200 µL Ti-IMAC buffer I, collecting the beads on a magnetic rack each time and discarding the wash buffer. Repeat for a total of three washes, then add beads to the microwell plate in a final volume of 150 µL per well   
5. Prepare three additional microwell plates, each with 150 µL of Ti-IMAC buffer I, II, and III.  
6. Start enrichment protocol from the BindIT control software as described in Leutert et al. ([10.15252/msb.20199021](https://doi.org/10.15252/msb.20199021))  
7. Add 80 µL Ti-IMAC elution buffer to each well of a microwell plate and add the plate to the indicated KingFisher slot during the pause after the last wash step, and continue with the peptide elution  
8. Immediately after elution, neutralize peptides by adding 10 μL of neutralization buffer.  
9. Follow the Evotip pure loading protocol for peptide binding to Evotips 

## 3.2 LC-MS/MS analysis

The following section describes the instrument settings for the commonly used Thermo Orbitrap Exploris 480 mass spectrometer. Other instrument platforms with sufficient mass resolving power in the reporter mass region can also be used, provided the appropriate settings are applied.

1. Peptides are separated on an Aurora Elite XT column using the 40SPD Evosep method, interfaced with the mass spectrometer through an EASY-spray source.  
2. On MS1 level, peptide masses were detected at a resolution of 120'000, with an ion target of $3\times10^6$, maximal injection time of 45 ms, and RF lens set to 40 $\%$  
3. MS2 spectra were recorded for the top 12 precursors with an intensity threshold of $10^3$, which were isolated at 0.7 Th and subjected to dynamic exclusion for 16 s. Normalized collision energy was set to 32 $\%$ and spectra were recorded with a resolution of 45'000 (phosphopeptides) or 30'000 and TMTpro on (total proteome), ion target of $1\times10^5$ and a maximal injection time of 250 ms (phosphopeptides) or Auto (total proteome)


## 3.3 Mass Spectrometry Data Processing


Raw mass spectrometry (MS) data acquired using data-dependent acquisition (DDA) combined with tandem mass tag (TMT) labeling are processed using specialized software to:
1. Database searching: Spectra are searched against a protein sequence database including species-specific sequences, common contaminants, and decoy sequences
2. Reporter ion quantification: Reporter ion intensities are extracted
3. Site localization scoring: Post-translational modification (PTM) site localization probabilities are calculated

Several free and commercial DDA-TMT compatible software suites are available, including Andromeda/MaxQuant (Cox et al., 2011), Proteome Discoverer (Thermo Fisher Scientific), FragPipe ([10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7)) (fragpipe.nesvilab.org), and PeptideShaker (peptide-shaker.compomics.com). Additional specialized tools for PTM site localization scoring include PTMProphet ([10.1038/s41467-020-17914-x](https://doi.org/10.1038/s41467-020-17914-x)), PhosphoRS ([10.1038/nmeth.1107](https://doi.org/10.1038/nmeth.1107)), and Ascore ([10.1038/nmeth.1107](https://doi.org/10.1038/nmeth.1107)).

### 3.3.1 FragPipe Method

The TMT 16-plex phospho workflow in FragPipe ([10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7)) 22.0 with MSFragger ([10.1038/nmeth.4256](https://doi.org/10.1038/nmeth.4256)) (version 4.1) was utilized.
Database searching was performed against a species-specific protein sequence database supplemented with common contaminants and reversed decoy sequences. The search parameters included:

- Fixed modifications:
  - Carbamidomethylation of cysteine (+57.0215 Da)
  - TMT labeling of lysine residues and peptide N-termini (+304.2071 Da)
- Variable modifications:
  - Phosphorylation on serine, threonine, and tyrosine residues (+79.9663 Da)
  - Oxidation of methionine (+15.9949 Da)
  - Acetylation at protein N-termini (+42.0106 Da)

Reporter ions are quantified by IonQuant (version 1.10.27), with downstream processing in TMTIntegrator (version 1.10.27). TMTIntegrator ([10.1101/2025.05.27.656447](https://doi.org/10.1101/2025.05.27.656447)) performs quantification and normalization specifically tailored for multiplexed TMT experiments by:

- Selecting the best peptide-spectrum match (PSM) based on intensity and quality metrics.
- Normalizing intensities by retention time bins and applying $\log_2$ transformations.
- Removing outliers based on interquartile range filtering.
- Aggregating quantified values to protein or PTM-site levels.

Selected TMTIntegrator parameters are:

- Best PSM selection
- Outlier removal enabled
- Allow overlabel/underlabel: OFF
- Use MS1 intensity: OFF
- Aggregation method: Median

For phospho-enriched samples, PTM sites with PTMProphet ([10.1021/acs.jproteome.9b00205](https://doi.org/10.1021/acs.jproteome.9b00205)) localization scores ≥ 0.75 are retained.

Note: FragPipe workflow parameters are stored in configuration files, which can be customized for different experimental designs. The exact configuration file used to process the data in this protocol is available for download from Zenodo [10.5281/zenodo.15850770](http://doi.org/10.5281/zenodo.15850770) or [gitlab.bfabric.org/wolski/PTM_example](https://gitlab.bfabric.org/wolski/PTM_example).

FragPipe provides two types of PTM-features reports for the phospho-enriched samples:

- Multisite features - refers to the set of all identified peptide-forms sharing the same set of modification sites under investigation, regardless of peptide length, its derivatives, cleavage state, or whether the modifications are unambiguously localized. Hence, also peptides with different sequence, ambiguous localization but identical possible modification sites are grouped into the same **multisite** feature.

- Site-level reports collapse signals from all peptidoforms mapping to the same amino acid position into a single quantitative value per PTM site, ensuring each site has exactly one intensity per sample. In downstream analysis, site-level intensities plug into site-centric enrichment or kinase-activity inference tools (PhosR, PTM-SEA, RoKAI, etc.), streamlining biological interpretation and network reconstruction.

## 3.4 Differential Expression Analysis using `prolfquapp`

Post-translational modifications (PTMs) play a crucial role in regulating protein function; however, their analysis is complex and challenging. A key challenge is to distinguish changes in PTM abundance that are due to altered protein expression from those that reflect a change in the modification stoichiometry (i.e., the fraction of the protein pool that is modified). This protocol provides a step-by-step computational workflow to analyze PTM data in the context of total protein expression changes. We leverage the R packages `prolfquapp` ([10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911)) for streamlined differential expression analysis and `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) for the integration, analysis, and visualization of PTM and total proteome data. The workflow is divided into two main parts: (1) initial differential expression analysis of the phospho-enriched (site and multisite feature abundances) and total proteome datasets (see **Supplementary Material Section A: Differential expression analysis using `prolfquapp`**), and (2) the integrated analysis of differential PTM usage (see **Supplementary Material Section B: Integration and analysis of PTM features using `prophosqua`**).


### 3.4.1 Data and Software Setup

The first step involves obtaining the FragPipe 22 mass spectrometry output files [10.5281/zenodo.15850770](http://doi.org/10.5281/zenodo.15850770) and setting up the analysis environment in R. The `prolfquapp` ([10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911)) package provides a set of shell scripts that automate the analysis workflow. These scripts are copied into the working directory to be used in subsequent steps. See **Supplementary Material Section A.1-A.2** for complete setup instructions and example dataset download.

### 3.4.2 Sample Annotation

A crucial step is the creation of a detailed sample annotation file. This file maps each raw data file to its experimental conditions . We parse sample names from the FragPipe output files using the `prolfqua_dataset.sh` script (see **Supplementary Material Section A.3.1**). From the sample names we extract experimental factors (e.g., Genotype, Timepoint). If the filenames do not encode the experimental conditions, the explanatory variables would need to be provided by adding columns. We then call the `annotation_add_contrasts` function in the `prolfqua` ([10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441)) to generate the required factorial contrasts automatically.

The `annotation_add_contrasts` function adds two key columns to the sample annotation:
- **`ContrastName`**: A descriptive name for each statistical contrast
- **`Contrast`**: The mathematical formula defining the contrast in terms of group means

For the experimental design in this protocol (2×3 factorial: Genotype [KO/WT] × Time [Uninfected/Early/Late]), the function generates four distinct contrasts:

1. **KO_vs_WT** (Main effect): `( (G_KO_Uninfect + G_KO_Late + G_KO_Early)/3 - (G_WT_Uninfect + G_WT_Late + G_WT_Early)/3 )` - Compares overall genotype effects across all timepoints
2. **KO_vs_WT_at_Uninfect** (Interaction): `G_KO_Uninfect - G_WT_Uninfect` - Compares genotypes specifically in the uninfected condition
3. **KO_vs_WT_at_Late** (Interaction): `G_KO_Late - G_WT_Late` - Compares genotypes specifically in the late infection condition  
4. **KO_vs_WT_at_Early** (Interaction): `G_KO_Early - G_WT_Early` - Compares genotypes specifically in the early infection condition

Note that setting `interaction=FALSE` in the `annotation_add_contrasts` function disables the generation of interaction contrasts between factors. This means the contrasts focus on main effects rather than examining how the effect of one factor depends on the level of another factor.

These contrasts enable comprehensive statistical analysis, distinguishing between general genotype effects and timepoint-specific effects. See **Supplementary Material Section A.3.1-A.3.2** for complete annotation workflow with code examples.


### 3.4.3 Execution of Differential Expression Analysis

Before running the differential expression analysis, we first generate a YAML configuration file with default settings using the `prolfqua_yaml.sh` script. This configuration file `config.yaml` specifies key analysis parameters like normalization method (variance stabilization normalization, vsn, by default), protein-level summarization (Tukey's median polish), and quality-filtering thresholds for peptide or protein observations (e.g., q-value cutoffs). 

With the annotation and configuration file in place, the `prolfqua_dea.sh` script is used to run the differential expression analysis. First, the script reads the sample annotation (sample IDs, experimental groups, blocking factors) alongside the YAML configuration. Intensities are $\log_2$-transformed and variance-stabilized before modeling to ensure homoscedasticity across runs.

Differential abundance is then estimated by fitting linear models via the prolfqua R-package's contrast API, with empirical Bayes variance moderation to improve statistical power. Proteins or PTM features exhibiting excessive missingness are down-weighted rather than imputed, reducing bias in fold-change estimates.

For each data type, `prolfquapp` generates a comprehensive suite of outputs:

- Dynamic HTML reports containing interactive PCA plots, boxplots, volcano plots, and heatmaps for quality control and exploratory analysis
- XLSX tables (DE_<datatype>.xlsx) summarizing fold changes, moderated t-statistics, raw p-values, and Benjamini–Hochberg–adjusted q-values
- GSEA rank files (GSEA_<datatype>.rnk) for downstream gene-set enrichment analysis
- SummarizedExperiment objects (<datatype>.rds) for seamless import into Shiny-based visualization tools

This is done separately for three data types of features derived from the mass spectrometry experiment:

1.  **Total Proteome:** Analysis of protein abundance changes
2.  **Multi-site PTM:** Analysis of multisites features abundances
3.  **Single-site PTM:** Analysis of single site features abundances


Complete execution commands and guidance for interpreting these outputs are detailed in Supplementary Material Sections A.4–A.6. All differential expression analysis results are publicly archived on Zenodo [10.5281/zenodo.15830989](https://doi.org/10.5281/zenodo.15830989).


## 3.5 Integrated PTM Analysis using `prophosqua`

The `prophosqua` package ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) provides tools for integrating and analyzing post-translational modification (PTM) data with total proteome measurements. It enables researchers to distinguish between changes in protein abundance and changes in modification site usage.

This part describes how to load the results of the Differential Expression Analysis, integrate the PTM and total proteome data, and how to determine differential PTM-feature usage. For complete implementation with detailed code examples, see **Supplementary Material Section B** ("Integration and analysis of PTM features using `prophosqua`").

### 3.5.1 Data Loading and Integration

The differential expression results, in Excel format, generated by `prolfquapp`, from the total proteome, and either multi or single-site PTM analyses, are loaded into the R. The two datasets are then integrated by performing a left-join operation on the protein IDs, merging the PTM-level statistics with the corresponding protein-level statistics of the total proteome experiment for each condition. A left join is used because we want to retain all PTM features, even if their parent proteins were not detected in the total proteome analysis - this ensures we don't lose any PTM information while still being able to normalize PTM changes by protein abundance when available. See **Supplementary Material Sections B.1-B.5** for complete data loading and integration workflow.



### 3.5.2 Analysis of Differential PTM-feature Expression (DPE)

**Definition:** DPE tests the raw PTM-feature signal change between conditions, without any correction for its parent protein's expression level. This analysis is used to flag any PTM-feature whose abundance changes, even if the parent protein itself is also up- or down-regulated.

**Method:** DPE results are visualized using N-to-C plots, generated by the `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) `::n_to_c_expression` function. These plots map both the $\log_2$ fold changes of individual modification sites along the primary sequence of the protein and the $\log_2$ fold change of the protein abundance from the total proteome experiment, providing a clear visual summary of both site-specific expression changes and overall protein regulation. See **Supplementary Material Section B.6** for complete DPE analysis implementation and visualization examples.


![N-to-C experssion plot showing differential expression of total protein abundance (light yellow rectangle) and phosphorylation sites (vertical lines, color coded by type of residue) along protein sequences. The vertical lines are either dashed or solid, depending on the site's imputation status. The x-axis represents amino acid position from N to C terminus, while the y-axis shows $\log_2$ fold changes. Significantly regulated sites (FDR < 0.05 and FDR < 0.01) are highlighted with one or two red asterisks, respectively. The plot enables visualization of both overall protein regulation and site-specific phosphorylation changes in a single view.](figures/dpe_n_to_c.png)


### 3.5.3 Analysis of Differential PTM-feature Usage (DPU)

**Definition:** DPU tests the **protein-normalized** changes of PTM-features. This analysis is essential for determining whether the *fraction* of a protein that is modified at a specific site changes between conditions. It isolates changes in modification stoichiometry from changes in overall protein abundance.

**Method:** DPU is calculated by subtracting the protein's $\log_2$ fold change from the PTM-feature's $\log_2$ fold change. The associated p-values are re-calculated using a method that combines the variance from both the PTM and protein-level models, as implemented in `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) `::test_diff` function, or in the MSstatsPTM R-package ([10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)).

The `test_diff` function implements a statistical framework for calculating differential PTM usage (DPU) that accounts for both PTM and protein-level changes. The analysis involves three key mathematical components:

1. **Effect Size Calculation ($\Delta_I$)**
2. **Standard Error Propagation ($SE_I$)**
3. **Degrees of Freedom Calculation ($df_I$)**

The differential PTM usage effect size is calculated as the difference between PTM and protein $\log_2$ fold changes:

$$\Delta_I = \log_2FC_{PTM} - \log_2FC_{Protein}$$

This represents the protein-normalized change in PTM abundance, where:
- **Positive values** indicate increased modification stoichiometry (higher fraction of protein is modified)
- **Negative values** indicate decreased modification stoichiometry (lower fraction of protein is modified)
- **Zero values** indicate no change in modification stoichiometry relative to protein expression

The combined standard error accounts for uncertainty in both PTM and protein measurements:

$$SE_I = \sqrt{SE_{PTM}^2 + SE_{Protein}^2}$$

This error propagation ensures that the statistical significance of usage changes reflects uncertainty from both datasets, providing a more robust assessment of the reliability of the observed changes.

The effective degrees of freedom for the combined test are calculated using Welch's approximation:

$$df_I = \frac{(SE_{PTM}^2 + SE_{Protein}^2)^2}{\frac{SE_{PTM}^4}{DF_{PTM}} + \frac{SE_{Protein}^4}{DF_{Protein}}}$$

This provides appropriate statistical power for the combined analysis while accounting for potentially different sample sizes and variances between PTM and protein measurements. The Welch approximation is particularly important when the two datasets have different experimental designs or sample sizes.

Using these calculated parameters, a t-statistic is computed:

$$t = \frac{\Delta_I}{SE_I}$$

The corresponding p-value is determined using the t-distribution with $df_I$ degrees of freedom, and multiple testing correction (e.g., Benjamini-Hochberg FDR) is applied across all tested sites.


 The PTM-feature usage fold changes, p-values, FDR estimates are then added to the dataframe. The DPU results are then visualized using N-to-C plots generated by `prophosqua::n_to_c_usage`, which display the protein-normalized PTM-feature abundance changes and highlight sites with significant changes in usage. See **Supplementary Material Sections B.7-B.8** for complete DPU analysis implementation and statistical framework details. 

For example, when comparing the N-to-C expression plot shown above, with the N-to-C usage plot, shown below, for Protein SQSTM1 (UniProt id Q64337), we can see that while raw PTM-feature expression is upregulated in the KO vs WT contrast, the actual usage of the phosphorylation sites remains unchanged or even is significantly downregulated, after normalizing for protein abundance (shown by the vertical lines with negative log2 fold changes). This demonstrates how DPU analysis can reveal changes in modification stoichiometry that differ from, or even oppose, changes in overall PTM-feature levels. 

 ![N-to-C usage plot showing differential usage of phosphorylation sites (vertical lines, color coded by type of residue) along the protein sequences. The vertical lines are either dashed or solid, depending on the site's imputation status. The x-axis represents amino acid position from N to C terminus, while the y-axis shows $\log_2$ fold changes. Significantly regulated sites (FDR < 0.05 and FDR < 0.01) are highlighted with one or two red asterisks, respectively. The plot enables visualization of both overall protein regulation and site-specific phosphorylation changes in a single view.](figures/dpu_n_to_c.png)


In the *Atg16l1* macrophage dataset [10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320), SQSTM1 (uniprot id Q64337) protein levels are markedly increased in Atg16l1-knockout (KO) versus wild-type (WT), reflecting loss of autophagy-mediated turnover. Phosphosite analysis reveals that the raw phosphopeptide intensities at $Ser_28$, $Ser_308$, and $Ser_368$ are similarly elevated in cKO compared to WT (DPE), consistent with higher protein abundance. However, when applying DPU normalization, which subtracts the protein $log_2$ fold change from each site's phospho $log_2$ fold change, these apparent phosphorylation increases are substantially reduced, indicating that most of the phosphosite enrichment arises from protein accumulation rather than genuine changes in modification stoichiometry. This example highlights how DPU effectively isolates regulatory events on SQSTM1 from confounding protein-level effects. 


### 3.5.3.1 Integration Report Generation

The integrated analysis results are compiled into a comprehensive HTML report using an R Markdown template `_Overview_PhosphoAndIntegration_site.Rmd`. The report includes:

- Overview statistics on identified phosphorylation sites
- Interactive scatter plots comparing protein vs PTM-feature changes
- Volcano plots highlighting differential PTM-feature usage

The report includes a data tables which allows to search for a specific protein and phosphorylations feature and to highlight it in the scatter and volcano plots. See **Supplementary Material Sections B.9-B.11** for complete report generation workflow and output descriptions.

### 3.5.3.2 Exporting Results to Excel

The integrated analysis results are exported to an Excel file for further analysis and sharing. The Excel file contains multiple worksheets:

- **combinedSiteProteinData**: Contains the merged PTM and protein-level data used for the integrated analysis
- **combinedStats**: Contains the differential PTM-feature usage statistics including fold changes, p-values, and FDR estimates for each phosphorylation site

This Excel format facilitates downstream analysis, data sharing with collaborators, and integration with other bioinformatics tools. The file is automatically saved with a timestamped filename in the results directory.


### 3.5.4 Sequence Motif Analysis

To infer which kinases may be responsible for the observed changes in phosphorylation, a sequence motif analysis is performed on sites with significant DPU.

**Method:** Amino acid sequences surrounding significantly regulated PTM sites (e.g., $FDR < 0.01$) are extracted. These sequences are then grouped by their regulation status (up- or downregulated) and by the experimental contrast. The `ggseqlogo` R package is used to generate sequence logo plots from these groups. These plots visualize conserved amino acid patterns, which can be compared to known kinase recognition motifs to identify potential upstream regulators. See **Supplementary Material Section C** for complete sequence motif analysis implementation and kinase prediction workflow.


# 4 Troubleshooting

## 4.1 Sample Preparation Issues

### 4.1.1 Low Protein Yields
**Problem:** Protein concentration after lysis is unexpectedly low (<0.5 mg/mL)
**Possible causes and solutions:**
- **Insufficient cell lysis:** Increase homogenization time, add additional sonication cycles, or use higher SDS concentration (up to 6%)
- **Sample degradation:** Ensure samples were properly stored at -80°C and minimize freeze-thaw cycles
- **Incomplete tissue disruption:** For tough tissues, increase Tissue Lyzer cycles or use liquid nitrogen grinding before lysis buffer addition
- **Protein precipitation:** Check pH of lysis buffer (should be 7.5±0.2), avoid prolonged heating above 95°C

**Quality control check:** Visual inspection should show complete tissue/cell disruption with no visible particulates after centrifugation

### 4.1.2 Poor Digestion Efficiency
**Problem:** Low peptide yield after trypsin digestion (<50% of expected)
**Possible causes and solutions:**
- **Trypsin activity loss:** Use fresh trypsin, store at recommended conditions, check expiration date
- **Incomplete reduction/alkylation:** Verify TCEP and ClAA concentrations, ensure adequate incubation time
- **Bead-to-protein ratio incorrect:** Maintain 10:1 bead:protein ratio for optimal binding
- **Ethanol concentration wrong:** Verify 60% final ethanol concentration for protein precipitation
- **pH too low for trypsin:** Ensure TEAB buffer pH is 8.0-8.5

**Quality control check:** Digestion efficiency should be >80% as assessed by SDS-PAGE or LC-MS analysis

### 4.1.3 TMT Labeling Problems
**Problem:** Poor labeling efficiency (<90%) or high variability between channels
**Possible causes and solutions:**
- **Residual primary amines:** Ensure complete removal of ammonium-containing buffers during cleanup
- **TMT reagent degradation:** Use fresh reagent, store with desiccant at -20°C, check for precipitates
- **Incorrect pH:** Peptide solution pH should be >7.5 for optimal labeling
- **Water contamination:** Use anhydrous acetonitrile, avoid moisture exposure during labeling
- **Insufficient reagent:** Use 4-fold molar excess of TMT reagent over primary amines

**Quality control check:** Analyze small aliquot by LC-MS/MS to verify >95% labeling before pooling

### 4.1.4 Enrichment Specificity Issues
**Problem:** Low phosphopeptide enrichment specificity (<60% phosphorylated peptides)
**Possible causes and solutions:**
- **Buffer contamination:** Prepare fresh Ti-IMAC buffers daily, use high-purity reagents
- **Incorrect bead amount:** Use 1:4 peptide:bead ratio (w/w) for optimal binding capacity
- **Acidic peptides competing:** Increase glycolic acid concentration to 0.2 M in binding buffer
- **Bead saturation:** Reduce peptide load or increase bead amount
- **Elution too harsh:** Use minimum elution volume and time to avoid non-specific elution

**Quality control check:** Analyze enriched fraction by LC-MS/MS; >70% phosphopeptides indicates good specificity

## 4.2 Instrumental Problems

### 4.2.1 LC Separation Issues
**Problem:** Poor chromatographic separation, peak broadening, or high backpressure
**Possible causes and solutions:**
- **Column degradation:** Replace column if >500 injections or pressure increases >50%
- **Mobile phase contamination:** Prepare fresh buffers, filter through 0.1 μm filters
- **Sample carryover:** Increase wash volume between injections, check autosampler cleaning
- **Blocked frits:** Replace column end fittings, filter samples before injection
- **Temperature fluctuation:** Ensure stable column temperature (±1°C)

**Quality control check:** Monitor backpressure trends and peak symmetry for early problem detection

### 4.2.2 Mass Spectrometer Performance
**Problem:** Low peptide identification rates or poor spectral quality
**Possible causes and solutions:**
- **Source contamination:** Clean ion source components, replace ESI emitter if needed
- **Calibration drift:** Perform mass calibration, check mass accuracy (<5 ppm error)
- **Low spray stability:** Adjust spray voltage, check for air bubbles in lines
- **Inappropriate acquisition parameters:** Optimize collision energy, isolation width, resolution settings
- **Overloading:** Reduce sample injection amount if ion current is too high

**Quality control check:** Monitor peptide identification rates and mass accuracy across runs

### 4.2.3 Data Acquisition Problems
**Problem:** Missing or low-quality MS/MS spectra
**Possible causes and solutions:**
- **Dynamic exclusion too aggressive:** Reduce exclusion time or increase exclusion mass window
- **AGC targets too low:** Increase ion accumulation targets for better spectral quality
- **Cycle time too fast:** Increase maximum injection time for better sensitivity
- **Charge state exclusion:** Ensure 2+ and 3+ charge states are included for TMT
- **Isolation window issues:** Use appropriate isolation width (0.7 Th for TMT)

**Quality control check:** >20,000 peptide identifications expected for typical TMT run

## 4.3 Data Analysis Challenges

### 4.3.1 Software-Specific Issues
**Problem:** FragPipe analysis fails or produces poor results
**Possible causes and solutions:**
- **Memory limitations:** Increase RAM allocation in FragPipe settings (minimum 16GB)
- **Database issues:** Use species-specific database with appropriate contaminants
- **Search parameters:** Verify mass tolerances, modifications, and enzyme specificity
- **File format problems:** Ensure raw files are properly converted or supported
- **Version compatibility:** Use recommended FragPipe version for your data type

**Quality control check:** >1% PSM identification rate indicates successful search

### 4.3.2 Statistical Interpretation Problems
**Problem:** No significant PTM sites detected despite biological expectation
**Possible causes and solutions:**
- **Insufficient statistical power:** Increase sample size or reduce technical variability
- **Inappropriate statistical thresholds:** Consider FDR = 0.05 instead of 0.01 for initial exploration
- **High missing data:** Filter features present in >70% of samples within each group
- **Normalization issues:** Check for batch effects, verify sample loading consistency
- **Experimental design:** Ensure adequate biological replicates (minimum n=3 per group)

**Quality control check:** Examine volcano plots and distribution of p-values for expected patterns

### 4.3.3 Integration Difficulties
**Problem:** Low match rates between PTM and protein datasets (<60%)
**Possible causes and solutions:**
- **Database inconsistencies:** Use identical FASTA databases for both analyses
- **Different filtering criteria:** Apply consistent FDR and quality filters to both datasets
- **Protein inference differences:** Verify identical protein inference settings
- **Sample preparation artifacts:** Check for differential sample loss between protocols
- **Annotation mismatches:** Verify gene/protein ID mapping between datasets

**Quality control check:** >80% match rate indicates good integration potential

## 4.4 Biological Interpretation Issues

### 4.4.1 DPE vs DPU Interpretation
**Problem:** Conflicting results between DPE and DPU analyses
**Biological interpretation:**
- **DPE+ only:** Changes driven by protein abundance alterations
- **DPU+ only:** True signaling changes with stable protein levels
- **DPE+ and DPU+:** Amplified signaling (both protein and modification efficiency increase)
- **DPE- and DPU+:** Compensatory regulation (modification increases despite protein decrease)

**Validation strategies:**
- Cross-reference with known biology and literature
- Examine individual protein examples manually
- Validate key sites by targeted analysis or Western blotting

### 4.4.2 Motif Analysis Limitations
**Problem:** Weak or absent sequence motifs in significantly regulated sites
**Possible causes and solutions:**
- **Too few significant sites:** Lower statistical thresholds for motif analysis only
- **Mixed kinase activities:** Analyze upregulated and downregulated sites separately
- **Non-canonical regulation:** Consider non-kinase mechanisms (phosphatases, binding partners)
- **Sequence window issues:** Try different window sizes (±5, ±7, ±10 amino acids)
- **PTM crosstalk:** Consider other modifications influencing kinase specificity

**Quality control check:** Compare motifs to known databases (PhosphoSitePlus, Kinexus)

## 4.5 Performance Optimization

### 4.5.1 Memory and Processing Issues
**Problem:** Analysis runs slowly or crashes due to resource limitations
**Solutions:**
- **Increase system RAM:** Minimum 32GB recommended for large datasets
- **Use SSD storage:** Faster I/O improves processing speed significantly
- **Process in chunks:** Analyze subsets of fractions separately if needed
- **Parallel processing:** Use multi-core capabilities where available
- **Cloud computing:** Consider high-memory cloud instances for large studies

### 4.5.2 Reproducibility Concerns
**Problem:** Results vary between technical replicates or analysis runs
**Prevention strategies:**
- **Standardize protocols:** Use automated procedures where possible
- **Document parameters:** Keep detailed records of all analysis settings
- **Version control:** Use specific software versions throughout study
- **Batch processing:** Analyze all samples together to minimize technical variation
- **Quality monitoring:** Track quality metrics across all samples

**Quality control check:** Technical replicates should show >0.9 correlation in abundance measurements

# 5. Notes

## 5.1 Critical Parameters

**Note 5.1** (referenced in step 3.1.1.1): **Complete cell disruption** is essential for reproducible protein extraction efficiency. Incomplete lysis leads to variable protein yields and can cause sample-to-sample variability that affects downstream quantitative analysis. Visual inspection should confirm no visible cell debris remains after centrifugation.

**Note 5.2** (referenced in step 3.1.1.13): **Bead-to-protein ratio** of 10:1 (μg beads per μg protein) is critical for efficient protein binding in SP3 protocols. Lower ratios result in protein loss, while higher ratios can cause non-specific binding and increased background. Microsphere-based protein clean-up strategies exploit protein aggregation in organic solvents ([10.1074/mcp.TIR118.001270](https://doi.org/10.1074/mcp.TIR118.001270)).

**Note 5.3** (referenced in step 3.1.1.17): **Digestion temperature and time** are critical for achieving optimal trypsin efficiency. The automated protocol typically yields ~20% digestion efficiency. If efficiency decreases, Lys-C pre-digestion can improve proteolytic efficiency and peptide recovery.

**Note 5.4** (referenced in step 3.1.1.19): **TMT reagent handling** requires careful attention to moisture control. Recent optimizations have reduced required TMT amounts ([10.1074/mcp.TIR119.001385](https://doi.org/10.1074/mcp.TIR119.001385)) by using lower reaction volumes and optimized peptide:TMT ratios. We typically aim for a 4:1 TMT:peptide ratio to ensure high labeling efficiency while minimizing costs.

**Note 5.5** (referenced in step 3.1.2.6): **Antibody bead volume** can be downscaled from manufacturer protocols, but peptide:bead ratios must be empirically tested for each sample type and antibody lot to maintain enrichment efficiency.

**Note 5.6** (referenced in step 3.1.2.7): **Antibody incubation conditions** require careful mixing to keep beads in suspension without splashing to tube walls. Loss of beads during washing significantly reduces enrichment efficiency.

**Note 5.7** (referenced in step 3.1.3.4): **HPLC system equilibration** typically requires 20 minutes for backpressure stabilization. Monitor pressure trends across all injections to ensure optimal column performance and lifetime.

**Note 5.8** (referenced in step 3.1.3.8): **Fraction concatenation strategy** can be programmed to omit fractions containing unreacted TMT reagent, typically observed as increased 260 nm UV absorbance. This optimization improves downstream enrichment efficiency.

**Note 5.9** (referenced in step 3.1.4.1): **Buffer freshness** is crucial for phosphopeptide enrichment efficiency. Ti-IMAC buffer degradation, particularly glycolic acid oxidation, significantly reduces binding specificity and capacity.

**Note 5.10** (referenced in step 3.1.4.3): **Bead preparation** from 20% stock slurry requires thorough mixing before use. Alternative chemistries like Zr-IMAC HP can be substituted and will generate overlapping but complementary peptide pools compared to Ti-IMAC.

## 5.2 Modifications and Variations

**Note 5.11**: **Sample types adaptation** - This protocol has been successfully adapted for various biological systems including yeast, plant tissues, and primary cells. Lysis conditions may require optimization: increase SDS concentration to 6% for tough plant cell walls, or reduce to 2% for fragile primary cells.

**Note 5.12**: **Alternative enrichment strategies** - The antibody-based enrichment step (3.1.2) can be replaced with or supplemented by other PTM enrichment approaches:
- **Tyrosine phosphorylation:** Use pY-specific antibodies (e.g., 4G10, pY100)
- **Motif-specific enrichment:** Employ kinase substrate motif antibodies for pathway-focused analysis
- **Dual enrichment:** Combine antibody and IMAC enrichment for comprehensive coverage

**Note 5.13**: **Fractionation alternatives** - High-pH fractionation can be replaced with:
- **StageTip-based fractionation:** More manual but suitable for smaller sample numbers
- **Strong cation exchange (SCX):** Alternative separation mechanism with different peptide distribution
- **Hydrophilic interaction chromatography (HILIC):** Useful for phosphopeptide-enriched samples

**Note 5.14**: **Instrumentation alternatives** - While optimized for Orbitrap systems, the protocol is compatible with:
- **Q-TOF instruments:** Adjust acquisition parameters for different fragmentation characteristics
- **Tribrid systems:** Take advantage of multiple fragmentation modes for improved site localization
- **Other TMT-compatible platforms:** Adapt reporter ion monitoring to instrument capabilities

## 5.3 Time Considerations

**Note 5.15**: **Protocol timeline** - Complete workflow timeline for 18 samples:
- **Day 1:** Sample lysis and digestion setup (4-6 hours active time)
- **Day 2:** TMT labeling and pooling (3-4 hours active time)
- **Day 3:** Optional antibody enrichment (6-8 hours with overnight incubation)
- **Day 4:** High-pH fractionation (8-10 hours including drying)
- **Day 5:** Phosphopeptide enrichment (6-8 hours)
- **Days 6-8:** LC-MS/MS analysis (36 hours instrument time for 36 fractions)
- **Days 9-10:** Data analysis (24-48 hours depending on dataset size)

**Note 5.16**: **Critical timing steps**:
- **TMT labeling:** 60 minutes at room temperature is optimal; longer incubation does not improve efficiency
- **Phosphopeptide binding:** 5-10 minutes is sufficient; longer incubation may increase non-specific binding
- **Column equilibration:** Allow 20 minutes minimum for stable chromatography

**Note 5.17**: **Storage and stopping points**:
- **After protein extraction:** Samples stable at -80°C for 6 months
- **After digestion:** Peptides stable at -80°C for 1 year
- **After TMT labeling:** Pooled samples stable at -80°C for 1 year
- **After fractionation:** Dried fractions stable at -80°C for 6 months
- **After enrichment:** Loaded Evotips stable at 4°C for 1 week

## 5.4 Safety Considerations

**Note 5.18**: **Chemical hazards**:
- **SDS:** Causes skin and eye irritation; wear appropriate PPE and work in fume hood
- **Chloroacetamide:** Light-sensitive alkylating agent; store in dark, handle with nitrile gloves
- **Trifluoroacetic acid:** Highly corrosive; use in fume hood with appropriate ventilation
- **Acetonitrile:** Flammable organic solvent; store away from ignition sources
- **Formic acid:** Corrosive; causes burns on contact with skin

**Note 5.19**: **Equipment safety**:
- **Tissue Lyzer:** Ensure proper tube sealing to prevent sample spillage at high speeds
- **Sonicator:** Use appropriate personal protective equipment to prevent hearing damage
- **SpeedVac:** Check for proper vacuum seal and ventilation to prevent solvent vapor exposure
- **Mass spectrometer:** Follow manufacturer guidelines for high-voltage safety procedures

**Note 5.20**: **Biological safety**:
- **Sample handling:** Treat all biological samples as potentially hazardous; follow institutional biosafety guidelines
- **Liquid nitrogen:** Use appropriate cryogenic safety equipment when snap-freezing samples
- **Waste disposal:** Follow institutional guidelines for disposal of chemical and biological waste

## 5.6 Quality Control and Validation

**Note 5.21**: **Protein extraction QC**:
- **Expected yields:** 1-5 mg/mL for cultured cells, 0.5-2 mg/mL for tissue samples
- **Purity assessment:** A280/A260 ratio should be >1.5 for protein-enriched extracts
- **Integrity check:** SDS-PAGE should show expected protein size distribution without significant degradation

**Note 5.22**: **TMT labeling QC**:
- **Efficiency verification:** LC-MS/MS analysis should show >95% peptide labeling
- **Channel balance:** Coefficient of variation <20% across TMT channels indicates proper sample loading
- **Completeness check:** <5% free lysine or N-terminal amine signals indicate successful labeling

**Note 5.23**: **Enrichment QC**:
- **Phosphopeptide specificity:** >70% of identified peptides should be phosphorylated
- **Recovery assessment:** Compare total peptide intensity before and after enrichment
- **Localization confidence:** Use PTMProphet scores ≥0.75 for confident site assignment

**Note 5.24**: **Data analysis QC**:
- **Search statistics:** >20,000 peptide identifications expected for typical TMT phospho run
- **Mass accuracy:** <5 ppm for precursor ions, <20 ppm for fragment ions
- **Integration success:** >80% match rate between PTM and protein datasets indicates good data quality

# 6. Alternative Approaches

## 6.1 Comparison with Other PTM Analysis Methods

**Label-free quantification vs. TMT-based approach:**
- **TMT advantages:** Higher throughput, reduced missing values, better quantitative precision
- **Label-free advantages:** No modification artifacts, unlimited sample capacity, lower cost per sample
- **Recommendation:** Use TMT for studies requiring high precision across many conditions; use label-free for discovery-phase studies or when sample numbers exceed multiplex capacity

**MSstatsPTM vs. prophosqua integration:**
- **MSstatsPTM:** More established statistical framework, peptidoform-level analysis
- **prophosqua:** Enhanced visualization tools, protein-centric analysis, N-to-C plots
- **Recommendation:** Consider MSstatsPTM for rigorous statistical validation; use prophosqua for hypothesis generation and biological interpretation

## 6.2 When to Use This vs. Other Approaches

**Use the integrated DPE/DPU approach when:**
- Biological system involves significant protein expression changes
- Need to distinguish signaling from expression effects
- Study design includes protein abundance measurements
- Publication requires rigorous PTM normalization

**Consider alternative approaches when:**
- Protein expression is stable across conditions
- Single PTM class analysis is sufficient
- Limited sample material available
- Budget constraints prohibit TMT labeling

## 6.3 Limitations and Scope

**Technical limitations:**
- Requires both PTM and protein measurements
- Limited to modifications detectable by enrichment
- Cannot analyze modifications with poor ionization
- Requires specialized computational expertise

**Biological scope:**
- Best suited for phosphorylation analysis
- Adaptable to other enrichable modifications (acetylation, ubiquitination)
- Not applicable to non-enrichable modifications
- Limited dynamic range compared to targeted approaches

# 7. Data Interpretation Guidelines

## 7.1 Step-by-Step Result Interpretation

### 7.1.1 Quality Assessment
1. **Check match rates:** >80% between PTM and protein datasets indicates good integration
2. **Examine coverage:** Assess protein and site coverage across experimental conditions
3. **Evaluate reproducibility:** Technical replicates should show >0.9 correlation
4. **Assess missing data:** <30% missing values per group for reliable quantification

### 7.1.2 Statistical Results Review
1. **DPE analysis:** Identifies sites where absolute PTM abundance changes
2. **DPU analysis:** Reveals sites where modification stoichiometry changes
3. **Volcano plot examination:** Look for balanced distribution of up/down regulation
4. **FDR validation:** Confirm appropriate statistical stringency for biological system

### 7.1.3 Biological Interpretation Framework
**DPE+/DPU- sites:** PTM changes driven by protein expression
**DPE-/DPU+ sites:** True signaling changes with stable protein levels
**DPE+/DPU+ sites:** Amplified signaling (protein and modification increase)
**DPE-/DPU- sites:** Coordinated decrease in protein and modification

## 7.2 Common Pitfalls in Biological Interpretation

**Overinterpretation of statistical significance:**
- Statistical significance does not guarantee biological relevance
- Consider effect sizes alongside p-values
- Validate key findings with orthogonal methods

**Ignoring protein context:**
- Always examine protein-level changes before interpreting PTM data
- Consider protein function and localization
- Check for known regulatory sites in literature

**Motif analysis overconfidence:**
- Sequence motifs suggest but do not prove kinase involvement
- Multiple kinases may share similar motifs
- Consider cellular context and kinase expression

## 7.3 Validation Recommendations

**High-priority sites for validation:**
- Sites with large effect sizes (|log2FC| > 1.5)
- Sites in proteins relevant to study hypothesis
- Novel regulatory sites not in databases
- Sites with contradictory DPE/DPU patterns

**Validation approaches:**
- Targeted mass spectrometry (PRM/SRM)
- Phospho-specific antibody validation
- Site-directed mutagenesis studies
- Kinase inhibitor experiments

# 8. Expected Results and Outcomes

## 8.1 Representative Dataset Characteristics

**Typical identification numbers (TMT 11-plex, 36 fractions):**
- Total proteome: 8,000-12,000 proteins
- Phosphoproteome: 15,000-25,000 phosphosites
- Integration success: 80-90% of phosphosites matched to proteins
- Significant changes: 5-15% of quantified sites (FDR < 0.05)

**Quality metrics:**
- Phosphopeptide specificity: >70%
- TMT labeling efficiency: >95%
- Missing data rate: <20% per experimental group
- Mass accuracy: <5 ppm precursor, <20 ppm fragment

## 8.2 Publication-Quality Outputs

**Generated files for publication:**
- Interactive HTML reports with quality control metrics
- Excel spreadsheets with complete statistical results
- Publication-ready figures (volcano plots, N-to-C plots, motif logos)
- Supplementary data tables formatted for journal submission

**Recommended data sharing:**
- Raw mass spectrometry files (ProteomeXchange)
- Search results and quantification matrices
- Complete R analysis scripts and parameters
- Supplementary material with detailed methods

## **Supplementary Materials**

The following materials are provided to support protocol implementation:

- **Supplementary Material Section A:** Complete differential expression analysis workflow using `prolfquapp` with example code and outputs
- **Supplementary Material Section B:** Integrated PTM analysis using `prophosqua` with detailed visualization examples  
- **Supplementary Material Section C:** Sequence motif analysis and kinase prediction workflows
- **Configuration files:** FragPipe search parameters and R package settings optimized for TMT phosphoproteomics
- **Example datasets:** Processed results from Atg16l1 macrophage study for protocol validation
- **Quality control reports:** Template reports for assessing data quality and integration success

All supplementary materials are available through the following repositories:
- **Protocol files:** Zenodo [10.5281/zenodo.15850770](http://doi.org/10.5281/zenodo.15850770)
- **Example datasets:** GitLab [gitlab.bfabric.org/wolski/PTM_example](https://gitlab.bfabric.org/wolski/PTM_example)
- **Analysis results:** Zenodo [10.5281/zenodo.15830989](https://doi.org/10.5281/zenodo.15830989)

**References**

