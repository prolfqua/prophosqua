---
title: "Integrated Analysis of Post-Translational Modifications and Total Proteome: Methods for Distinguishing Expression from Usage Changes"
author:
  - Witold Wolski
  - Antje Dittmann
  - (Christian Panse)
  - (Laura Kunz)
  - Jonas Grossmann
date: 2025-07-15
running_head: "Integrated PTM and Proteome Analysis"
output:
  pdf_document:
    toc: true
    toc_depth: 2
    number_sections: true
    toc_title: "Table of Contents"
  word_document:
    toc: true
    toc_depth: 2
    number_sections: true
    toc_title: "Table of Contents"
  html_document:
bibliography: doi_references_key.bib
---


| Author             | Affiliation(s)    | Corresponding Author Email                             |
|--------------------|-------------------|--------------------------------------------------------|
| Jonas Grossmann    | FGCZ, SIB         | [jonas.grossmann@fgcz.uzh](mailto:jonas.grossmann@fgcz.uzh)      |
| Antje Dittmann     | FGCZ              | [antje.dittmann@fgcz.uzh.ch](mailto:antje.dittmann@fgcz.uzh.ch)  |
| Laura Kunz         | FGCZ              |                                                        |
| Christian Panse    | FGCZ, SIB         |                                                        |
| Witold Wolski      | FGCZ, SIB         | [witold.wolski@fgcz.uzh.ch](mailto:witold.wolski@fgcz.uzh.ch)   |

- FGCZ - Functional Genomics Center Zurich, Winterthurerstrasse 190, CH-8057 Zurich, 
- SIB - Swiss Institute of Bioinformatics - Amphipole, Quartier UNIL-Sorge, CH-1015 Lausanne


## **Abstract**

**Background:** Post-translational modifications (PTMs), particularly phosphorylation, regulate protein function at substoichiometric levels, making their quantitative analysis technically challenging. A critical limitation in PTM studies is distinguishing between changes in modification abundance due to altered protein expression versus genuine changes in modification stoichiometry (the fraction of protein molecules that are modified).

**Methods:** This chapter presents a comprehensive workflow integrating TMT-based quantitative proteomics with specialized computational analysis to distinguish differential PTM expression (DPE) from differential PTM usage (DPU). The protocol encompasses automated sample preparation using Single-Pot Solid-Phase-Enhanced Sample Preparation (SP3)-based digestion, phosphopeptide enrichment with Ti-IMAC, high-pH offline fractionation, LC-MS/MS analysis, and integrated statistical analysis using the `prolfqua`, `prolfquapp` and `prophosqua` R packages.

**Key Features:** The workflow provides automated, scalable protocols for large sample cohorts, robust statistical frameworks for integrated PTM-protein analysis, comprehensive visualization tools including N-to-C plots for protein-centric PTM mapping, and sequence motif analysis for kinase prediction. All protocols are optimized for reproducibility and include extensive quality control measures.

**Applications:** This integrated approach enables researchers to identify true signaling changes in PTM studies, prioritize functionally relevant modification sites, understand pathway-level regulation in disease contexts, and generate robust datasets suitable for publication in high-impact journals.

## Key Words

Post-translational modifications, phosphoproteomics, TMT labeling, differential expression analysis, protein stoichiometry, mass spectrometry, prolfquapp, prophosqua, data integration, kinase activity

# 1. Introduction

## 1.1 Biological Significance of Post-Translational Modifications

Protein phosphorylation is a reversible post-translational modification (PTM) that modulates a protein’s conformation, enzymatic activity, subcellular localization, and binding interactions. Because this modification can be reversed rapidly, it acts as a molecular switch in cellular signaling networks. As a result, cells respond to environmental cues without requiring new protein synthesis. Dysregulation of phosphorylation is implicated in cancer, neurodegeneration, and metabolic disorders [10.1080/14789450.2021.1976152](https://doi.org/10.1080/14789450.2021.1976152). Large-scale phosphoproteomic profiling by high-throughput mass spectrometry (e.g., TMT-based workflows) is therefore essential to dissect disease-associated signaling alterations and identify candidate therapeutic targets.

## 1.2 Technical Challenges in PTM Analysis

The analysis of post-translational modifications (PTMs), such as phosphorylation, presents significant analytical challenges due to the low abundance, or substoichiometry, of modified peptides. Typically, phosphorylated peptides represent less than 1–5% of the total peptide population ([10.1186/1477-5956-4-15](https://doi.org/10.1186/1477-5956-4-15)). Additionally, phosphopeptides generally exhibit poorer ionization efficiency in mass spectrometry-based proteomics compared to non-modified peptides, further complicating their detection and quantification ([10.1021/acs.jproteome.0c00613](https://doi.org/10.1021/acs.jproteome.0c00613)). This low abundance and challenging ionization behavior necessitate specialized experimental workflows and computational strategies to achieve accurate quantification ([10.1016/j.mcpro.2024.100754](https://doi.org/10.1016/j.mcpro.2024.100754)).

A key analytical challenge is distinguishing genuine changes in PTM levels, due to altered kinase activity or signaling states, from changes caused by altered protein abundance ([10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477); [10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708)). Without proper correction, this confounding factor can lead researchers to misinterpret PTM data, incorrectly attributing higher phosphorylation signals solely to increased modification efficiency or pathway activation, rather than to elevated protein expression ([10.1039/c5mb00024f](https://doi.org/10.1039/c5mb00024f)). Controlling this confounder is critical to discern rapid, signaling-driven phosphorylation events from slower, transcriptionally or translationally regulated changes in protein abundance.


Because phosphorylation is substoichiometric, only a small fraction of protein molecules carry the modification, enrichment to isolate phosphopeptides is required. Enrichment steps (such as IMAC or TiO₂) introduce additional sample handling and technical variability (10.1039/c5mb00024f; 10.1016/j.aca.2021.338716). Therefore, each stage—from digestion to cleanup to MS acquisition must be carefully optimized to improve reproducibility (10.1038/s41467-018-03309-6). Effective experimental design combined with tailored computational approaches (e.g. normalization strategies) are therefore essential to recover reliable, biologically meaningful PTM signals.

In TMT-based phosphoproteomics, sample multiplexing improves reproducibility and reduces missing values across experiments. Adding offline fractionation before LC–MS boosts depth—revealing thousands more phosphosites—while also distributing peptide load and minimizing matrix effects in each run. Together, these strategies control technical variance and ensure consistent quantification across fractions and conditions [10.1002/pmic.202100245](https://doi.org/10.1002/pmic.202100245).

In TMT-based proteomics, a “plex” refers to the complete set of isobarically tagged samples analyzed together in a single LC–MS/MS run. When the total number of samples exceeds the available TMT channels, each plex should include at least one biological replicate from every experimental group and a common reference or pooled channel for normalization [10.1021/acs.jproteome.0c00536](https://doi.org/10.1021/acs.jproteome.0c00536). Replicates of the same group are then distributed across different plexes rather than clustered within a single run, preventing any group’s replicates from being confined to one batch. If all samples fit within the channel capacity of a single plex, these balancing measures are inherently satisfied in that run.

## 1.3 Current Methodological Approaches

### 1.3.1 Mass Spectrometry-Based Quantification Strategies

Isobaric tagging with TMT reagents has emerged as the gold standard for multiplexed PTM analysis, enabling simultaneous analysis of up to 35 samples in a single combined run and significantly reducing inter-sample technical variability [10.1021/acs.jproteome.4c00668](https://doi.org/10.1021/acs.jproteome.4c00668). This multiplexing approach provides enhanced statistical power while minimizing batch effects that can confound biological interpretation, creating closed systems with unique statistical properties that allow for greater confidence in comparative results [10.1021/ac0262560](https://doi.org/10.1021/ac0262560). Enrichment of phosphorylated peptides is typically achieved through metal affinity chromatography (IMAC) or metal oxide affinity chromatography (MOAC), with Fe³⁺, Zr⁴⁺, Ti⁴⁺ ions, or TiO₂ providing complementary selectivity profiles that enable coverage of the phosphoproteome ([10.1074/mcp.m114.045609](https://doi.org/10.1074/mcp.m114.045609); [10.1101/2020.04.13.038810](https://doi.org/10.1101/2020.04.13.038810)). The practical implementation of these strategies, including automated sample preparation protocols using SP3-based digestion and Fe-NTA enrichment, is detailed in the experimental methods section of this chapter.

### 1.3.2 Computational Tools and Workflows

The computational landscape for quantitative protein post-translational modification (PTM) analysis has evolved significantly, with numerous specialized software platforms supporting bottom-up proteomics workflows. Following mass spectrometry data acquisition, spectra undergo database searching to identify peptide sequences, which is subsequently followed by protein inference to determine protein identities based on these identified peptides. The ecosystem of PTM analysis includes several well-established DDA-TMT compatible software suites, both free and commercial, such as Andromeda integrated within MaxQuant [10.1021/acs.jproteome.4c00869](https://doi.org/10.1021/acs.jproteome.4c00869), Proteome Discoverer (Thermo Fisher Scientific), FragPipe [10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7), and PeptideShaker [10.1021/acs.jproteome.1c00678](https://doi.org/10.1021/acs.jproteome.1c00678).

Among the analytical steps following peptide identification, PTM site localization scoring represents a particularly crucial computational challenge. Site localization determines the confidence with which a modification can be assigned to a specific amino acid residue within a peptide sequence. This process is especially challenging when peptides contain multiple potential modification sites, necessitating algorithms to evaluate the quality and specificity of site-determining fragment ions. Specialized tools have been developed explicitly for robust PTM site localization, including PTMProphet [10.1038/s41467-020-17914-x](https://doi.org/10.1038/s41467-020-17914-x), PhosphoRS [10.1021/pr200611n](https://doi.org/10.1021/pr200611n), and Ascore [10.1038/nbt1240](https://doi.org/10.1038/nbt1240). Accurate localization is essential, as incorrect site assignments can lead to biological misinterpretations, particularly in signaling pathway analyses and kinase-substrate predictions.

In data-dependent acquisition (DDA) workflows, FragPipe [10.1101/2025.05.27.656447](https://doi.org/10.1101/2025.05.27.656447) represents a widely adopted and integrated platform that comprehensively addresses PTM analysis, including site localization challenges. FragPipe integrates the MSFragger search engine [10.1038/nmeth.4256](https://doi.org/10.1038/nmeth.4256), enabling highly sensitive peptide identification, and employs PTMProphet functionality for robust site localization scoring. It also generates detailed quantitative outputs, including site-level quantification reports and multisite feature reports grouping peptides that share identical modification patterns. A detailed description of the FragPipe workflow is provided in the experimental methods section of this chapter.

While the above discussion has focused primarily on DDA workflows, data-independent acquisition (DIA)-based approaches represent alternative strategies with distinct computational considerations for PTM analysis. Specialized DIA platforms such as Spectronaut (Biognosys) and DIA-NN [10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x) offer PTM support through dedicated site-localization scoring algorithms and site-level quantification capabilities. Unlike DDA, DIA workflows require specialized computational approaches to extract and quantify modification-specific information from highly multiplexed fragmentation spectra. Both Spectronaut and DIA-NN incorporate advanced statistical frameworks specifically designed for robust PTM identification and quantification in DIA workflows.

 Data Analysis Frameworks

Proteomic data analysis can be conducted at multiple levels of granularity, each providing distinct analytical perspectives. A peptidoform represents a specific peptide sequence with a particular set of modifications, while site-level analysis aggregates signals for each residue position across multiple peptides. Differential analysis can therefore target either peptidoforms or individual modification sites, with each approach offering unique advantages for different research questions.

FragPipe introduces an additional analytical concept through multisite features, which refers to the set of all identified peptide-forms sharing the same set of modification sites under investigation, regardless of sequence derivatives, cleavage state, or whether the modifications are unambiguously localized. This approach groups peptidoforms with different sequences, which can arise from missed cleavages or alternative cleavage patterns, but identical possible modification sites into the same multisite feature, providing a balanced approach between peptidoform specificity and site-level aggregation. Discarding all peptidoforms with ambiguous localization would result in the loss of a large fraction of PTM information. Furthermore, multisite features preserve information about phosphorylation sites that occur together on the same peptide molecule, revealing co-modification patterns that indicate coordinated regulation by kinases or functional relationships between sites. Traditional site-level analysis treats each phosphorylation site independently, losing this valuable information about which specific combinations of sites are modified together in biological samples.

Site-level reports collapse signals from all peptidoforms mapping to the same PTM at an amino acid residue into a single quantitative value, ensuring each site has exactly one intensity measurement per sample. In downstream analysis, site-level intensities integrate seamlessly with site-centric enrichment tools and kinase-activity inference platforms (PhosR [10.1016/j.cels.2021.04.007](https://doi.org/10.1016/j.cels.2021.04.007), PTM-SEA [10.1016/j.molcel.2019.05.030](https://doi.org/10.1016/j.molcel.2019.05.030), RoKAI [10.1038/s41467-021-21211-6](https://doi.org/10.1038/s41467-021-21211-6), etc.), streamlining biological interpretation and network reconstruction analyses [10.1186/s12014-020-09290-x](https://doi.org/10.1186/s12014-020-09290-x).

### 1.3.4 Statistical Analysis of PTM Data

Statistical analysis of PTM data requires specialized methods due to confounding effects between PTM levels and protein abundances. Several statistical packages have emerged to address the critical issue of integrating PTM-feature quantifications with protein abundance data. Notable examples include `MSstatsPTM` [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477) and `msqrob2PTM` [10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708), which explicitly model both modified peptide and protein-level changes to accurately distinguish genuine PTM regulation from protein abundance effects.

The `msqrob2PTM` framework exemplifies the statistical approaches now available for PTM analysis, defining two complementary analytical strategies: Differential Peptidoform Abundance (DPA) and Differential Peptidoform Usage (DPU) [10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708). DPA directly models the log₂ intensities of each modified peptide (peptidoform), detecting absolute changes in PTM levels between experimental conditions. In contrast, DPU adjusts the PTM intensities by the corresponding protein-level changes, effectively testing for changes in the relative usage or stoichiometry of a modification site. This dual approach enables researchers to distinguish between PTM changes driven by protein abundance differences (DPA) versus those reflecting genuine changes in modification efficiency or regulatory activity (DPU).

These two frameworks differ in their analytical approaches for addressing confounding between PTM and protein abundance changes. `MSstatsPTM` employs a "model then correct" strategy, separately fitting linear models to modified and unmodified peptide data, and then combining statistical inferences post-modeling. In contrast, `msqrob2PTM` follows a "correct then model" approach, first adjusting PTM intensities by their corresponding protein abundances at the data level, and subsequently performing statistical testing on these normalized values. While both methodologies ultimately estimate PTM abundance changes corrected for protein-level effects, differences in workflow design can influence practical considerations such as computational complexity, ease of interpretation, and user preference.

To further extend these capabilities and offer researchers additional flexibility, we developed `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272), an R package designed to streamline statistical analysis of phosphoproteomics and other PTM datasets. Similar to `MSstatsPTM` and `msqrob2PTM`, `prophosqua` explicitly integrates protein-level abundance data into PTM analyses to accurately distinguish genuine PTM regulation from protein-level changes. `prophosqua` implements a robust statistical framework, based on the `prolfqua` package [10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441), that leverages linear modeling of protein-corrected PTM data, generates adjusted significance estimates, and provides intuitive visualizations for comprehensive exploratory analysis.

### 1.3.5 The Protein Assignment Problem in Differential PTM-feature Usage 

Identification and quantification tools need to assign peptides to proteins, i.e., find the leading or representative protein ID. This assignment may be dataset-specific and may differ between the PTM and the total proteome dataset. Therefore, using protein IDs to match PTM features with proteins from the total run might lead to a loss of information because of differing representative IDs [10.1093/bib/bbr065](https://doi.org/10.1093/bib/bbr065). For example, a phosphopeptide might be assigned to protein isoform A in the PTM dataset while the corresponding total protein abundance is quantified under isoform B in the proteome dataset, preventing proper integration and DPU analysis.

A robust solution involves matching stripped peptide sequences from the PTM dataset to peptide sequences of proteins quantified in the total proteome dataset, bypassing protein ID dependency entirely. When peptides are shared among multiple protein isoforms, the averaged protein abundance across all matching isoforms provides a more stable reference for normalization. While this sequence-based matching approach requires additional computational implementation, it significantly improves integration success rates and reduces the loss of valuable PTM information due to arbitrary protein assignment differences between datasets.


## 1.5 Chapter Overview and Learning Objectives

This chapter introduces an analytical workflow that addresses the biological challenge of accurately distinguishing Differential PTM Expression (DPE) from Differential PTM Usage (DPU). Distinguishing between protein-level effects and genuine modification-specific changes is critical for meaningful interpretation of PTM data.

**Differential PTM Expression (DPE)** tests raw PTM signal changes between experimental conditions, identifying any modification abundance changes regardless of underlying protein-level effects. This analysis captures the total effect of experimental perturbations on PTM levels, including both direct modification effects and indirect effects mediated through protein abundance changes.

**Differential PTM Usage (DPU)** evaluates protein-normalized PTM changes, explicitly identifying modification sites where stoichiometry is genuinely altered independent of protein abundance. This analysis pinpoints sites where experimental conditions directly impact modification efficiency, kinase activity, or phosphatase activity, thereby highlighting specific regulatory events.

The integrated analytical approach combines optimized sample preparation protocols, mass spectrometry acquisition strategies, and computational tools (including the `prolfquapp` and `prophosqua` R packages) to provide a complete solution spanning from sample processing to biological interpretation. This framework enables researchers to extract maximum biological insight from PTM datasets while avoiding common analytical pitfalls that have historically complicated PTM data interpretation.

This protocol presents a complete experimental and analytical workflow for integrated PTM studies. The wet lab component describes the optimized protocols used at the Functional Genomics Center Zurich for PTM services, covering automated sample preparation, enrichment strategies, and LC-MS/MS acquisition. The computational analysis is demonstrated using the Atg16l1 macrophage dataset (Maculins *et al.*, eLife 2021, [10.7554/elife.62320](https://doi.org/10.7554/elife.62320)), which comprises TMT-11-plex measurements from six conditions (WT/KO × uninfected/early/late infection) with both total proteome and phospho-enriched samples.

**Learning objectives include:**
1. Implementing automated, scalable sample preparation protocols using SP3-based digestion and KingFisher Flex automation
2. Performing phosphopeptide enrichment with Fe-NTA and optional antibody-based enrichment strategies
3. Executing high-pH offline fractionation and optimized LC-MS/MS acquisition on Evosep-Orbitrap systems
4. Conducting integrated statistical analysis of PTM and protein data using `prolfquapp` and `prophosqua` R packages
5. Interpreting DPE vs. DPU results for biological insights and generating publication-quality visualizations
6. Performing sequence motif analysis for kinase prediction and pathway reconstruction

The workflow emphasizes reproducibility, automation, and biological interpretation, spanning from sample lysis through data analysis and making it suitable for both experienced researchers and newcomers to integrated PTM analysis. All protocols include extensive quality control measures and troubleshooting guidelines for large-scale studies.

# 2. Materials

## 2.1 Biological Samples

### 2.1.1 Demonstration Dataset

This protocol uses the Atg16l1 macrophage dataset from Maculins et al. [10.7554/elife.62320](https://doi.org/10.7554/elife.62320) to demonstrate the bioinformatics workflow and integrated statistical analysis. The dataset comprises TMT-11-plex measurements from six conditions (WT/KO × uninfected/early/late infection) with both phospho-enriched and total proteome samples, making it ideal for illustrating the principles of integrated PTM analysis. In the original publication, the authors also performed KGG-enrichment (ubiquitin remnant), though we focus on the phosphoproteome and total proteome datasets for demonstration purposes. This same dataset was previously used in the `MSstatsPTM` publication (Kohler et al., MCP 2023, [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)), providing additional validation of the analytical approaches.

### 2.1.2 Sample Types Compatible with FGCZ Protocol

The wet lab protocols described in this chapter are optimized for high-throughput processing and are compatible with a wide range of biological sample types:

- **Fresh-frozen tissue samples** (≤5 mg per sample): Optimal for preserving PTM states and minimizing degradation
- **Cultured cell pellets** (≥1×10⁶ cells per sample): Standard format for most cell culture experiments  
- **Primary cells isolated from tissue**: Requires careful handling to maintain viability during isolation
- **Organoids and spheroid cultures**: Emerging model systems requiring specialized collection protocols

### 2.1.3 Sample Storage and Handling Requirements

Proper sample handling is critical for PTM preservation and protocol success:

**Storage requirements:**
- Store samples at -80°C immediately after collection to prevent PTM degradation
- Avoid freeze-thaw cycles, which can lead to protein degradation and PTM loss
- Process samples within 6 months of collection for optimal results
- Use appropriate sample labeling and tracking systems for large studies

**Pre-processing considerations:**
- Flash-freeze samples in liquid nitrogen when possible
- Minimize time between sample collection and freezing (target <5 minutes)
- Consider phosphatase and kinase inhibitor cocktails for certain sample types
- Document collection conditions and timing for batch effect analysis

### 2.1.4 Protocol Compatibility and Adaptations

The experimental methods presented here represent optimized protocols developed at the Functional Genomics Center Zurich. While the demonstration uses the Maculins et al. dataset, the analytical approaches are broadly applicable to TMT-based phosphoproteomics studies regardless of the specific sample preparation method used, provided that:

- Chemical labeling (TMT or similar isobaric tags) was employed
- Both PTM-enriched and total proteome datasets are available
- Consistent sample preparation was applied across experimental conditions
- Appropriate biological and technical replicates are included in the experimental design

The protocols emphasize automation, scalability, and reproducibility, making them suitable for large cohort studies while maintaining the analytical rigor necessary for high-impact publications.


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
   - `prolfquapp` [10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911): `install.packages("prolfquapp")`
   - `prolfqua` [10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441): dependency of `prolfquapp`
   - `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272): `devtools::install_github("fgcz/prophosqua")`
   - Additional dependencies: `tidyverse`, `ggseqlogo`, `writexl`, `rmarkdown`

**System requirements:**
- **Minimum:** 16GB RAM, 100GB free storage, 4-core processor
- **Recommended:** 32GB+ RAM, 500GB+ SSD storage, 8+ core processor
- **Operating system:** Windows 10+, macOS 10.15+, or Linux Ubuntu 18.04+

# 3. Methods

## 3.1 Common Sample Preparation Protocol

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
22. Pool a small fraction of each sample (e.g.,1 μL, ensure precise pipetting) and pool at equal ratio to check for sample loading variation and labeling efficiency using LC-MS/MS in QC step (**Critical QC step:** Assess before proceeding; see Note 4.11)
23. If needed, adjust the pooling ratio based on QC LC-MS/MS, otherwise, proceed with pooling all samples at an equal ratio (**Expected CV:** <20% across channels)
24. Optional: if antibody-based enrichment is required, split the sample in half
25. Dry down in SpeedVac (**Stopping point:** can be stored frozen at -80°C for up to 1 year)
26. Optional: For final clean-up, subject sample to C18 solid-phase extraction (SPE) (e.g., Sep-Pak, Waters) (**Recovery:** >85% of peptides)

**Expected outcomes:**
- **Total peptide yield:** 20-40 μg per 50 μg input protein
- **TMT labeling efficiency:** >95% (verify by LC-MS/MS)
- **Sample purity:** Suitable for downstream enrichment and fractionation
- **Storage stability:** 1 year at -80°C without significant degradation

### 3.1.2 Sample splitting and workflow assignment

After TMT labeling and pooling, split the sample into two workflows:

- **Aliquot A (5-10% of total):** Total proteome analysis (proceed to section 3.2)
- **Aliquot B (90-95% of total):** Phospho-enrichment analysis (proceed to section 3.3)


## 3.2 Total Proteome Workflow

### 3.2.1 Sample cleanup and fractionation for total proteome

1. Take 5-10% of the pooled TMT-labeled peptides for total proteome analysis
2. Subject sample to C18 solid-phase extraction (SPE) for cleanup if not already performed
3. Proceed with high-pH offline fractionation following the protocol described in section 3.2.2
4. Collect fewer fractions (12-24) compared to the phospho-enrichment workflow to account for lower sample complexity

### 3.2.2 High pH offline fractionation for total proteome

Follow the general fractionation protocol with modifications for the total proteome:

1. For total proteome analysis, aim for 2-5 µg of peptide per fraction
2. Generate 12-24 fractions (adjust FC.FractionRange accordingly)
3. Follow steps 2-11 from section 3.3.1 with appropriate scaling

### 3.2.3 LC-MS/MS analysis - Total proteome

1. Peptides are separated on an Aurora Elite XT column using the 40SPD Evosep method, interfaced with the mass spectrometer through an EASY-spray source
2. On MS1 level, peptide masses were detected at a resolution of 120,000, with an ion target of 3×10⁶, maximal injection time of 45 ms, and RF lens set to 40%
3. MS2 spectra were recorded for the top 12 precursors with an intensity threshold of 10³, which were isolated at 0.7 Th and subjected to dynamic exclusion for 16 s. Normalized collision energy was set to 32% and spectra were recorded with a resolution of **30,000 with TMTpro on**, ion target of **1×10⁵** and a maximal injection time of **Auto**

## 3.3 Phospho-enrichment Workflow

### 3.3.1 High pH offline fractionation for phospho-enrichment

This protocol describes the general flow of events for high-pH offline fractionation of peptides on the Thermo Vanquish Flex System, featuring automatic fraction collection and concatenation for subsequent phospho-peptide enrichment on the KingFisher Flex.

1. For robust and sensitive performance of phosphoproteomics experiments, we typically aim for 10-25 µg of peptide per fraction. These amounts have been demonstrated to be effective for automated protocols on the KingFisher Flex ([10.1016/j.mcpro.2024.100754](https://doi.org/10.1016/j.mcpro.2024.100754), [10.1002/pmic.202100245](https://doi.org/10.1002/pmic.202100245)). We generate 36 fractions, corresponding to a total of approximately 180-540 µg of pooled peptide before fractionation
2. Resuspend the peptides in 100 µL of high-pH buffer A, and sonicate for 5 min to ensure complete resuspension. Pellet insoluble material and only transfer supernatant into the HPLC vial
3. Purge and equilibrate HPLC with high pH buffer A and 2.1% high pH buffer B
4. Connect the XBridge peptide BEH column and adjust the flow rate to 0.75 mL/min, column compartment heater to 40 °C, and equilibrate the system until pressure and UV readings stabilize (see Note 4.7)
5. Place the KingFisher deep-well plate in the fraction collector (FC) and check the 'Use Safe Needle Height' setting in Chromeleon FC General Settings
6. Set up a sample queue in Chromeleon and add the custom variable specifying the fractionation starting position (needs to be specified in the Script editor)
7. Always run QC without fractionation using BSA peptide standard and check the elution profile for proper peak separation before an actual sample
8. For sample fractionation, adjust the FC.FractionRange in the Script editor to 36 and the MaxTubesPerFraction to 72, collecting each fraction for 60 seconds. The needle will start again at the first position after 36 fractions, i.e., concatenating fraction 37 into fraction 1, fraction 38 into fraction 2, etc. (see Note 4.8)
9. Set up the flow gradient to run from 2.1% B to 42.1% B in 54 minutes, followed by 5 minutes at 95% B and another 10 minutes equilibrating back to 2.1% B
10. After the run, equilibrate the column with 100% acetonitrile and then return the system to 20% methanol
11. Remove the deep-well plate from the fraction collector, snap-freeze, and dry down the samples in a SpeedVac. Peptides can be stored at -80°C

### 3.3.2 Phosphopeptide enrichment and Evotip loading

1. Resuspend each fraction in 150 µL Ti-IMAC buffer I (see Note 4.9), sonicate for 10 min
2. If total proteome measurements are required as well, take a small aliquot of either the input before phospho-peptide binding or the non-bound fraction after phospho-peptide binding
3. Prepare an appropriate volume of Ti-IMAC HP beads (see Note 4.10) from 20% bead slurry, aiming for a peptide-to-bead ratio of 1:4
4. Wash the beads with 200 µL Ti-IMAC buffer I, collecting the beads on a magnetic rack each time and discarding the wash buffer. Repeat for a total of three washes, then add beads to the microwell plate in a final volume of 150 µL per well
5. Prepare three additional microwell plates, each with 150 µL of Ti-IMAC buffer I, II, and III
6. Start enrichment protocol from the BindIT control software as described in Leutert et al. ([10.15252/msb.20199021](https://doi.org/10.15252/msb.20199021))
7. Add 80 µL Ti-IMAC elution buffer to each well of a microwell plate and add the plate to the indicated KingFisher slot during the pause after the last wash step, and continue with the peptide elution
8. Immediately after elution, neutralize peptides by adding 10 μL of neutralization buffer
9. Follow the Evotip pure loading protocol for peptide binding to Evotips

### 3.3.3 Optional antibody-based enrichment

This protocol describes the targeted enrichment of peptides harboring specific kinase substrate motifs or phosphorylated tyrosines using PTMScan HS antibody beads from Cell Signaling Technology. Binding and wash buffers have been modified, and an additional second enrichment step using Fe-NTA spin columns is used to further increase phospho-peptide enrichment specificity (see Note 4.5).

1. Resuspend one half of the pooled labeled peptides in 400 µL cold IP binding buffer, vortex, and incubate rotating at 4°C for 5 min, check pH (should be neutral or slightly basic, not below pH 7)
2. Spin down solution 5 min at 20,000×g to clear the solution and pellet insoluble material
3. Gently mix the antibody bead slurry to obtain a uniform bead suspension and transfer 10 µL of the antibody slurry (see Note 4.6) into an Eppendorf tube
4. Wash the beads with 500 µL of cold IP binding buffer by inverting the tube to resuspend the beads in the buffer, and then remove the buffer by placing the tube on a magnetic rack. Repeat for a total of four washes
5. Transfer the soluble peptide solution to the washed beads
6. Incubate on an end-over-end shaker for 2 hours at 4°C (see Note 4.7)
7. Collect buffer solution by spinning at no more than 1000×g for 4-5 seconds
8. Collect beads on magnetic rack and transfer supernatant into a new tube (can be saved and stored at -80°C)
9. Wash the beads with 400 µL of cold IP binding buffer by inverting the tube to resuspend the beads. Collect the beads on a magnetic rack and discard the buffer. Repeat for a total of two washes
10. Wash the beads with 400 µL of cold IP wash buffer by inverting the tube to resuspend the beads. Collect the beads on a magnetic rack and discard the buffer. Repeat for a total of two washes
11. Wash the beads with 400 µL of cold water by inverting the tube to resuspend the beads. Collect the beads on a magnetic rack and discard the water. Repeat for a total of two washes
12. Add 50 µL High Select Fe-NTA binding buffer and incubate beads for 10 min at room temperature on a shaker or Thermomixer. Ensure beads stay in suspension but are not splashed to the sides of the tube
13. Collect beads and transfer elution to High Select Fe-NTA spin column conditioned as per the manufacturer's instructions
14. Repeat the elution step and combine with the first elution in the same spin column
15. Follow the manufacturer's instructions for binding and washing routines
16. Add 30 µL elution buffer and collect elution into a fresh Eppendorf tube. Repeat elution once into the same collection tube
17. Dry peptides down to near-completeness, leaving ~5-10 µL of liquid in the tube if continuing with data acquisition immediately. Otherwise, dry and store peptides at -80°C
18. Acidify with 40 µL buffer A, check pH (should be at or below pH 4, see section 3.4.4)

### 3.3.4 LC-MS/MS analysis - Phosphoproteome

1. Peptides are separated on an Aurora Elite XT column using the 40SPD Evosep method, interfaced with the mass spectrometer through an EASY-spray source
2. On MS1 level, peptide masses were detected at a resolution of 120,000, with an ion target of $3\times10^6$, maximal injection time of 45 ms, and RF lens set to 40%
3. MS2 spectra were recorded for the top 12 precursors with an intensity threshold of $10^3$, which were isolated at 0.7 Th and subjected to dynamic exclusion for 16 s. Normalized collision energy was set to 32 $\%$ and spectra were recorded with a resolution of **45,000**, ion target of **$1\times10^5$** and a maximal injection time of **250 ms**

## 3.4 Mass Spectrometry Data Processing

Raw mass spectrometry (MS) data acquired using data-dependent acquisition (DDA) combined with tandem mass tag (TMT) labeling are processed using specialized software to:
1. Database searching: Spectra are searched against a protein sequence database including species-specific sequences, common contaminants, and decoy sequences
2. Reporter ion quantification: Reporter ion intensities are extracted
3. Site localization scoring: Post-translational modification (PTM) site localization probabilities are calculated


### 3.3.1 FragPipe Method

The TMT 16-plex phospho workflow in FragPipe [10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7) 22.0 with MSFragger [10.1038/nmeth.4256](https://doi.org/10.1038/nmeth.4256) (version 4.1) was utilized.
Database searching was performed against a species-specific protein sequence database supplemented with common contaminants and reversed decoy sequences. The search parameters included:

- Fixed modifications:
  - Carbamidomethylation of cysteine (+57.0215 Da)
  - TMT labeling of lysine residues and peptide N-termini (+304.2071 Da)
- Variable modifications:
  - Phosphorylation on serine, threonine, and tyrosine residues (+79.9663 Da)
  - Oxidation of methionine (+15.9949 Da)
  - Acetylation at protein N-termini (+42.0106 Da)

Reporter ions are quantified by IonQuant (version 1.10.27), with downstream processing in TMTIntegrator (version 1.10.27). TMTIntegrator [10.1101/2025.05.27.656447](https://doi.org/10.1101/2025.05.27.656447) performs quantification and normalization specifically tailored for multiplexed TMT experiments by:

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

For phospho-enriched samples, PTM sites with PTMProphet [10.1021/acs.jproteome.9b00205](https://doi.org/10.1021/acs.jproteome.9b00205) localization scores ≥ 0.75 are retained.

Note: FragPipe workflow parameters are stored in configuration files, which can be customized for different experimental designs. The exact configuration file used to process the data in this protocol is available for download from Zenodo [10.5281/zenodo.15850770](https://doi.org/10.5281/zenodo.15850770) or [gitlab.bfabric.org/wolski/PTM_example](https://gitlab.bfabric.org/wolski/PTM_example).

FragPipe provides two types of PTM-features reports for the phospho-enriched samples:

- Multisite features - refers to the set of all identified peptide-forms sharing the same set of modification sites under investigation, regardless of peptide length, its derivatives, cleavage state, or whether the modifications are unambiguously localized. Hence, also peptides with different sequence, ambiguous localization but identical possible modification sites are grouped into the same **multisite** feature.

- Site-level reports collapse signals from all peptidoforms mapping to the same amino acid position into a single quantitative value per PTM site, ensuring each site has exactly one intensity per sample. In downstream analysis, site-level intensities plug into site-centric enrichment or kinase-activity inference tools (PhosR, PTM-SEA, RoKAI, etc.), streamlining biological interpretation and network reconstruction.

## 3.5 Differential Expression Analysis using `prolfquapp`

Post-translational modifications (PTMs) play a crucial role in regulating protein function; however, their analysis is complex and challenging. A key challenge is to distinguish changes in PTM abundance that are due to altered protein expression from those that reflect a change in the modification stoichiometry (i.e., the fraction of the protein pool that is modified). This protocol provides a step-by-step computational workflow to analyze PTM data in the context of total protein expression changes. We leverage the R packages `prolfquapp` [10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911) for streamlined differential expression analysis and `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272) for the integration, analysis, and visualization of PTM and total proteome data. The workflow is divided into two main parts: (1) initial differential expression analysis of the phospho-enriched (site and multisite feature abundances) and total proteome datasets (see **Supplementary Material Section A: Differential expression analysis using `prolfquapp`**), and (2) the integrated analysis of differential PTM usage (see **Supplementary Material Section B: Integration and analysis of PTM features using `prophosqua`**).


### 3.4.1 Data and Software Setup

The first step involves obtaining the FragPipe 22 mass spectrometry output files [10.5281/zenodo.15850770](https://doi.org/10.5281/zenodo.15850770) and setting up the analysis environment in R.

- For total proteome samples `psm.tsv` files containing peptide-spectrum match data
- For phospho-enriched samples:
  - abundance_multi-site.tsv files for multisite feature analysis
- abundance_single-site.tsv files for individual modification site analysis


 The `prolfquapp` [10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911) package provides a set of shell scripts that automate the analysis workflow. These scripts are copied into the working directory to be used in subsequent steps. See **Supplementary Material Section A.1-A.2** for complete setup instructions and example dataset download.

### 3.4.2 Sample Annotation

A crucial step is the creation of a detailed sample annotation file. This file maps each raw data file to its experimental conditions. We parse sample names from the FragPipe output files using the `prolfqua_dataset.sh` script (see **Supplementary Material Section A.3.1**). From the sample names, we extract experimental factors (e.g., Genotype, Timepoint). If the filenames do not encode the experimental conditions, the explanatory variables would need to be provided by adding columns. We then call the `annotation_add_contrasts` function in the `prolfqua` [10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441) to generate the required factorial contrasts automatically.

The `annotation_add_contrasts` function adds two key columns to the sample annotation:
- **`ContrastName`**: A descriptive name for each statistical contrast
- **`Contrast`**: The mathematical formula defining the contrast in terms of group means

For the experimental design in this protocol (2×3 factorial: Genotype [KO/WT] × Time [Uninfected/Early/Late]), the function generates four distinct contrasts:

1. **KO_vs_WT** (Main effect): `( (G_KO_Uninfect + G_KO_Late + G_KO_Early)/3 - (G_WT_Uninfect + G_WT_Late + G_WT_Early)/3 )` - Compares overall genotype effects across all timepoints
2. **KO_vs_WT_at_Uninfect** (Interaction): `G_KO_Uninfect - G_WT_Uninfect` - Compares genotypes specifically in the uninfected condition
3. **KO_vs_WT_at_Late** (Interaction): `G_KO_Late - G_WT_Late` - Compares genotypes specifically in the late infection condition  
4. **KO_vs_WT_at_Early** (Interaction): `G_KO_Early - G_WT_Early` - Compares genotypes specifically in the early infection condition

Note that setting `interaction=FALSE` in the `annotation_add_contrasts` function disables the generation of interaction contrasts between factors. Interaction contrasts assess whether the effect of one factor depends on the level of another, thereby capturing non‑additive interplay between factors.

These contrasts enable comprehensive statistical analysis, distinguishing between general genotype effects and timepoint-specific effects. See **Supplementary Material Section A.3.1-A.3.2** for complete annotation workflow with code examples.


### 3.4.3 Execution of Differential Expression Analysis

Before running the differential expression analysis, we first generate a YAML configuration file with default settings using the `prolfqua_yaml.sh` script. This configuration file `config.yaml` specifies key analysis parameters like normalization method (variance stabilization normalization, vsn, by default), protein-level summarization (Tukey's median polish), and quality-filtering thresholds for peptide or protein observations (e.g., q-value cutoffs). 

With the annotation and configuration file in place, the `prolfqua_dea.sh` script is used to run the differential expression analysis. First, the script reads the sample annotation (sample IDs, experimental groups, blocking factors) alongside the YAML configuration. Intensities are $\log_2$-transformed and variance-stabilized before modeling to ensure homoscedasticity across runs.

Differential abundance is then estimated by fitting linear models via the `prolfqua` R-package's contrast API, with empirical Bayes variance moderation to improve statistical power. Proteins or PTM features exhibiting excessive missingness are down-weighted rather than imputed, reducing bias in fold-change estimates.

For each data type, `prolfquapp` generates a comprehensive suite of outputs:

- Dynamic HTML reports containing interactive PCA plots, boxplots, volcano plots, and heatmaps for quality control and exploratory analysis
- XLSX tables (DE_<datatype>.xlsx) summarizing fold changes, moderated t-statistics, raw p-values, and Benjamini–Hochberg–adjusted q-values
- GSEA rank files (GSEA_<datatype>.rnk) for downstream gene-set enrichment analysis
- SummarizedExperiment objects (<datatype>.rds) for seamless import into Shiny-based visualization tools

This is done separately for three data types of features derived from the mass spectrometry experiment:

1.  **Total Proteome:** Analysis of protein abundance changes
2.  **Multi-site PTM:** Analysis of multisites features abundances
3.  **Single-site PTM:** Analysis of single site features abundances


Complete execution commands and guidance for interpreting these outputs are detailed in Supplementary Material Sections A.4–A.6. All differential expression analysis results are publicly archived on Zenodo [10.5281/zenodo.15830988](https://doi.org/10.5281/zenodo.15830988).


## 3.6 Integrated PTM Analysis using `prophosqua`

The `prophosqua` package [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272) provides tools for integrating and analyzing post-translational modification (PTM) data with total proteome measurements. It enables researchers to distinguish between changes in protein abundance and changes in modification site usage.

This part describes how to load the results of the Differential Expression Analysis, integrate the PTM and total proteome data, and how to determine differential PTM-feature usage. For complete implementation with detailed code examples, see **Supplementary Material Section B** ("Integration and analysis of PTM features using `prophosqua`").

### 3.6.1 Data Loading and Integration

The differential expression results, in Excel format, generated by `prolfquapp`, from the total proteome, and either multi or single-site PTM analyses, are loaded into the R. The two datasets are then integrated by performing a left-join operation on the protein IDs, merging the PTM-level statistics with the corresponding protein-level statistics of the total proteome experiment for each condition. A left join is used because we want to retain all PTM features, even if their parent proteins were not detected in the total proteome analysis - this ensures we don't lose any PTM information while still being able to normalize PTM changes by protein abundance when available. See **Supplementary Material Sections B.1-B.5** for complete data loading and integration workflow.

The current prolfquapp implementation uses protein ID matching. This is a limitation, as it may lead to a loss of information because of differing representative IDs in the PTM-enriched and total proteome dataset (Note 4.12).



### 3.6.2 Analysis of Differential PTM-feature Expression (DPE)

**Definition:** DPE tests the raw PTM-feature signal change between conditions, without any correction for its parent protein's expression level. This analysis is used to flag any PTM-feature whose abundance changes, even if the parent protein itself is also up- or down-regulated.

**Method:** DPE results are visualized using N-to-C plots, generated by the `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272) `::n_to_c_expression` function. These plots map both the $\log_2$ fold changes of individual modification sites along the primary sequence of the protein and the $\log_2$ fold change of the protein abundance from the total proteome experiment, providing a clear visual summary of both site-specific expression changes and overall protein regulation. See **Supplementary Material Section B.6** for complete DPE analysis implementation and visualization examples.


![N-to-C expression plot showing differential expression of total protein abundance (light yellow rectangle) and phosphorylation sites (vertical lines, color coded by type of residue) along protein sequences. The vertical lines are either dashed or solid, depending on the site's imputation status. The x-axis represents amino acid position from N to C terminus, while the y-axis shows $\log_2$ fold changes. Significantly regulated sites (FDR < 0.05 and FDR < 0.01) are highlighted with one or two red asterisks, respectively. The plot enables visualization of both overall protein regulation and site-specific phosphorylation changes in a single view.](figures/dpe_n_to_c.png)


### 3.6.3 Analysis of Differential PTM-feature Usage (DPU)

**Definition:** DPU tests the **protein-normalized** changes of PTM-features. This analysis is essential for determining whether the *fraction* of a protein that is modified at a specific site changes between conditions. It isolates changes in modification stoichiometry from changes in overall protein abundance.

**Method:** DPU is calculated by subtracting the protein's $\log_2$ fold change from the PTM-feature's $\log_2$ fold change. The associated p-values are recalculated using a method that combines the variance from both the PTM and protein-level models, as implemented in `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272) `::test_diff` function, or in the `MSstatsPTM` R-package [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477).

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


The PTM-feature usage fold changes, p-values, and FDR estimates are then added to the dataframe. The DPU results are then visualized using N-to-C plots generated by `prophosqua::n_to_c_usage`, which display the protein-normalized PTM-feature abundance changes and highlight sites with significant changes in usage. See **Supplementary Material Sections B.7-B.8** for complete DPU analysis implementation and statistical framework details.


For example, when comparing the N-to-C expression plot shown above, with the N-to-C usage plot, shown below, for Protein SQSTM1 (UniProt id Q64337), we can see that while raw PTM-feature expression is upregulated in the KO vs WT contrast, the actual usage of the phosphorylation sites remains unchanged or even is significantly downregulated, after normalizing for protein abundance (shown by the vertical lines with negative log2 fold changes). This demonstrates how DPU analysis can reveal changes in modification stoichiometry that differ from, or even oppose, changes in overall PTM-feature levels. 


![N-to-C usage plot showing differential usage of phosphorylation sites (vertical lines, color coded by type of residue) along the protein sequences. The vertical lines are either dashed or solid, depending on the site's imputation status. The x-axis represents amino acid position from N to C terminus, while the y-axis shows $\log_2$ fold changes. Significantly regulated sites (FDR < 0.05 and FDR < 0.01) are highlighted with one or two red asterisks, respectively. The plot enables visualization of both overall protein regulation and site-specific phosphorylation changes in a single view.](figures/dpu_n_to_c.png)


Because generating N-to-C plots is computationally intensive, we generate plots only for proteins with at least one significant phosphorylation site or multi-site feature. The FDR threshold for determining significance can be adjusted in the `n_to_c_expression` function, allowing researchers to balance computational efficiency with comprehensive visualization based on their specific statistical considerations and interpretation guidelines. When no significant PTM sites are detected despite biological expectations, several factors should be considered, which include insufficient statistical power, which can be addressed by increasing sample size (see Note 4.13). 

In the *Atg16l1* macrophage dataset [10.7554/elife.62320](https://doi.org/10.7554/elife.62320), SQSTM1 (uniprot id Q64337) protein levels are markedly increased in Atg16l1-knockout (KO) versus wild-type (WT), reflecting loss of autophagy-mediated turnover. Phosphosite analysis reveals that the raw phosphopeptide intensities at $Ser_28$, $Ser_308$, and $Ser_368$ are similarly elevated in cKO compared to WT (DPE), consistent with higher protein abundance. However, when applying DPU normalization, which subtracts the protein $log_2$ fold change from each site's phospho $log_2$ fold change, these apparent phosphorylation increases are substantially reduced, indicating that most of the phosphosite enrichment arises from protein accumulation rather than genuine changes in modification stoichiometry. This example highlights how DPU effectively isolates regulatory events on SQSTM1 from confounding protein-level effects. 

When evaluating statistical significance in DPU analysis, it is important to avoid overinterpretation of results (see Note 4.13). Statistical significance does not guarantee biological relevance, and effect sizes should be considered alongside p-values when prioritizing sites for further investigation. Key findings should be validated through orthogonal methods such as targeted mass spectrometry or phospho-specific antibodies to confirm biological relevance.


### 3.6.4 Integration Report Generation

The integrated analysis results are compiled into a comprehensive HTML report using an R Markdown template `_Overview_PhosphoAndIntegration_site.Rmd`. The report includes:

- Overview statistics on identified phosphorylation sites
- Interactive scatter plots comparing protein vs PTM-feature changes
- Volcano plots highlighting differential PTM-feature usage

The report includes data tables that allow searching for a specific protein and phosphorylation feature and highlighting it in the scatter and volcano plots. See **Supplementary Material Sections B.9-B.11** for complete report generation workflow and output descriptions.

When analyzing the interactive reports, different combinations of DPE and DPU results provide distinct biological insights that guide interpretation. For instance, **DPE+ only** indicates changes driven primarily by protein abundance alterations, where PTM increases reflect higher protein levels rather than enhanced modification efficiency (see Note 4.15).

### 3.6.5 Exporting Results to Excel

The integrated analysis results are exported to an Excel file for further analysis and sharing. The Excel file contains multiple worksheets:

- **combinedSiteProteinData**: Contains the merged PTM and protein-level data used for the integrated analysis
- **combinedStats**: Contains the differential PTM-feature usage statistics, including fold changes, p-values, and FDR estimates for each phosphorylation site

This Excel format facilitates downstream analysis, data sharing with collaborators, and integration with other bioinformatics tools. The file is automatically saved with a timestamped filename in the results directory.


### 3.6.6 Sequence Motif Analysis

To infer which kinases may be responsible for the observed changes in phosphorylation, a sequence motif analysis is performed on sites with significant DPU.

Amino acid sequences surrounding significantly regulated PTM sites (e.g., $FDR < 0.01$) are extracted. These sequences are then grouped by their regulation status (up- or downregulated) and by the experimental contrast. The `ggseqlogo` R package is used to generate sequence logo plots from these groups. These plots visualize conserved amino acid patterns, which can be compared to known kinase recognition motifs to identify potential upstream regulators. When motif analysis yields weak or absent sequence patterns, several factors may require optimization including the number of significant sites, mixed kinase activities, or sequence window parameters (see Note 4.16). See **Supplementary Material Section C** for complete sequence motif analysis implementation and kinase prediction workflow.



![Sequence logo plots showing amino acid motifs surrounding significantly regulated phosphorylation at Serine. The height of each letter represents the frequency of that amino acid at each position relative to the phosphorylation site (position 8). Separate logos are generated for down- and up-regulated sites in each experimental contrast. Common kinase recognition motifs can be identified from these patterns to suggest upstream regulators.](figures/SequenceLogoPlot.png)


# 4. Notes


## 4.1 Complete cell disruption

Complete cell disruption is essential for reproducible protein extraction efficiency. Incomplete lysis leads to variable protein yields and can cause sample-to-sample variability that affects downstream quantitative analysis. Visual inspection should confirm no visible cell debris remains after centrifugation.

## 4.2 Bead-to-protein ratio

Bead-to-protein ratio of 10:1 (μg beads per μg protein) is critical for efficient protein binding in SP3 protocols. Lower ratios result in protein loss, while higher ratios can cause non-specific binding and increased background. Microsphere-based protein clean-up strategies exploit protein aggregation in organic solvents [10.1074/mcp.tir118.001270](https://doi.org/10.1074/mcp.tir118.001270).

## 4.3 Digestion temperature and time

Digestion temperature and time are critical for achieving optimal trypsin efficiency. The automated protocol typically yields ~20% digestion efficiency. If efficiency decreases, Lys-C pre-digestion can improve proteolytic efficiency and peptide recovery.

## 4.4 TMT reagent handling

TMT reagent handling requires careful attention to moisture control. Recent optimizations have reduced required TMT amounts [10.1074/mcp.tir119.001385](https://doi.org/10.1074/mcp.tir119.001385) by using lower reaction volumes and optimized peptide:TMT ratios. We typically aim for a 4:1 TMT:peptide ratio to ensure high labeling efficiency while minimizing costs.

## 4.5 Antibody bead volume

Antibody bead volume can be downscaled from manufacturer protocols, but peptide:bead ratios must be empirically tested for each sample type and antibody lot to maintain enrichment efficiency.

## 4.6 Antibody incubation conditions

Antibody incubation conditions require careful mixing to keep beads in suspension without splashing to tube walls. Loss of beads during washing significantly reduces enrichment efficiency.

## 4.7 HPLC system equilibration

HPLC system equilibration typically requires 20 minutes for backpressure stabilization. Monitor pressure trends across all injections to ensure optimal column performance and lifetime.

## 4.8 Fraction concatenation strategy

Fraction concatenation strategy can be programmed to omit fractions containing unreacted TMT reagent, typically observed as increased 260 nm UV absorbance. This optimization improves downstream enrichment efficiency.

## 4.9 Buffer freshness

Buffer freshness is crucial for phosphopeptide enrichment efficiency. Ti-IMAC buffer degradation, particularly glycolic acid oxidation, significantly reduces binding specificity and capacity.

## 4.10 Bead preparation

Bead preparation from 20% stock slurry requires thorough mixing before use. Alternative chemistries like Zr-IMAC HP can be substituted and will generate overlapping but complementary peptide pools compared to Ti-IMAC.


## 4.11 Automated QC reporting with FragPipe TMT QC script

An automated quality control script generates comprehensive TMT labeling efficiency reports directly from FragPipe PSM output files. The script evaluates labeling completeness at both peptide and PSM levels for N-terminal and lysine modifications, calculates missed cleavage rates for tryptic digestion efficiency, and provides quantitative channel balance assessment across all TMT channels. Key metrics include percentage of modified N-termini and lysine residues (expected >95%), missed cleavage rates for lysine and arginine residues (typically 5-15%), and relative abundance distributions across channels to detect loading imbalances. The script outputs interactive HTML reports with visualization of identification numbers per channel, modification frequencies, total abundance distributions, and density plots for rapid assessment of data quality before proceeding with downstream analysis. The QC script is part of the `prophosqua` R-package vignettes [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272).


## 4.12 Match Rates

Check match rates:** >80% between PTM and protein datasets indicates good integration. Lower rates may be due to protein assignment differences in the PTM and total proteome datasets, and sequence-based matching could improve the match rates.

## 4.13 No significant PTM sites detected despite biological expectation

No significant PTM sites detected despite biological expectation can be due to:

- **Insufficient statistical power:** Increase sample size or reduce technical variability
- **Inappropriate statistical thresholds:** Consider FDR = 0.25 instead of $0.05$ for initial exploration
- **Normalization issues:** Check for batch effects, verify sample loading consistency
- **Outlier detection:** Check for outliers in the data, e.g., using density or scatter plots
- **Experimental design:** Ensure adequate biological replicates (minimum n=3 per group)
- **Quality control check:** Examine volcano plots and distribution of p-values for expected patterns

## 4.14 Overinterpretation of statistical significance

Overinterpretation of statistical significance can be due to:

- Statistical significance does not guarantee biological relevance
- Consider effect sizes alongside p-values
- Validate key findings with orthogonal methods


## 4.15 DPE vs DPU Interpretation

Conflicting results between DPE and DPU analyses can be due to:

- **DPE+ only:** Changes driven by protein abundance alterations
- **DPU+ only:** True signaling changes with stable protein levels
- **DPE+ and DPU+:** Amplified signaling (both protein and modification efficiency increase)
- **DPE- and DPU+:** Compensatory regulation (modification increases despite protein decrease)


## 4.16 Motif Analysis Limitations

Weak or absent sequence motifs in significantly regulated sites can be due to:

- **Too few significant sites:** Lower statistical thresholds for motif analysis only
- **Mixed kinase activities:** Analyze upregulated and downregulated sites separately
- **Non-canonical regulation:** Consider non-kinase mechanisms (phosphatases, binding partners)
- **Sequence window issues:** Try different window sizes (±5, ±7, ±10 amino acids)
- **PTM crosstalk:** Consider other modifications influencing kinase specificity

## 4.17 Label-free quantification vs. TMT-based approach

- **TMT advantages:** Higher throughput, reduced missing values, better quantitative precision; enables efficient use of fractionation with all samples in each fraction
- **Label-free advantages:** No modification artifacts, unlimited sample capacity, lower cost per sample
- **Fractionation synergy:** TMT allows all samples to be compared within each fraction, with cross-fraction integration; DIA preferred when avoiding fractionation entirely
- **Recommendation:** Use TMT for studies requiring high precision across many conditions, especially when fractionation can improve depth; use label-free/DIA for discovery-phase studies or when sample numbers exceed multiplex capacity

## 4.18 MSstatsPTM vs. prophosqua vs. msqrob2PTM

- **MSstatsPTM:** Model-first approach (separate modeling then adjustment), established workflow
- **prophosqua:** Model-first approach (with easy correct-first option), enhanced visualization tools, flexible input handling (multi-site, single-site, FragPipe, DIA-NN outputs), N-to-C plots
- **msqrob2PTM:** Correct-first approach (normalize then model), built on QFeatures infrastructure, handles both peptidoform and PTM-level analysis
- **Statistical framework:** All use equivalent statistical rigor and correct for protein abundance; key difference is modeling strategy (separate modeling vs. direct normalization)
- **Recommendation:** Consider `MSstatsPTM` for the model-first approach, `msqrob2PTM` for the correct-first approach, and `prophosqua` for either approach plus flexible input handling and enhanced visualization


## 4.19 Use the integrated DPE/DPU approach when

Use the integrated DPE/DPU approach when:

- Biological system involves significant protein expression changes
- Need to distinguish signaling from expression effects
- Study design includes protein abundance measurements
- Publication requires rigorous PTM normalization

Consider enriched phosphoproteome measurement only when:

- Protein expression is stable across conditions, e.g., timescale is too short for protein expression changes to occur
- Budget constraints prohibit comprehensive analysis

# 5. Supplementary Materials

| Header | Description | File(s) | Web Archive |
|--------|-------------|-------|----------------|
| Sup. Material Sect. A: | Differential expression analysis workflow using `prolfquapp` and `prolfqua` | Supplementary_Material_v2.pdf | [10.5281/zenodo.15830988](https://doi.org/10.5281/zenodo.15830988) |
| Sup. Material Sect. B: | Integrated PTM analysis using `prophosqua` | Supplementary_Material_v2.pdf | [10.5281/zenodo.15830988](https://doi.org/10.5281/zenodo.15830988) |
| Sup. Material Sect. C: | Sequence motif analysis for kinase prediction workflows | Supplementary_Material_v2.pdf | [10.5281/zenodo.15830988](https://doi.org/10.5281/zenodo.15830988) |
| Quality control for TMT labelling | Template reports for assessing data quality and integration success | QCReport.html | [10.5281/zenodo.15830988](https://doi.org/10.5281/zenodo.15830988) |
| FragPipe Configuration files | FragPipe search parameters optimized for TMT phosphoproteomics | fragpipe.workflow | [10.5281/zenodo.15850770](http://doi.org/10.5281/zenodo.15850770) |
| FragPipe results total proteome | Processed results from Atg16l1 macrophage study - total proteome | p1/psm.tsv; p2/psm.tsv | [10.5281/zenodo.15850770](http://doi.org/10.5281/zenodo.15850770) |
| FragPipe results phospho enriched | Processed results from Atg16l1 macrophage study - phospho enriched samples | abundance_multi-site_None.tsv; abundance_single-site_None.tsv | [10.5281/zenodo.15850770](http://doi.org/10.5281/zenodo.15850770) |



**References**

