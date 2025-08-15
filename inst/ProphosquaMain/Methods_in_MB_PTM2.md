---
title: "Integrated Analysis of Post-Translational Modifications and Total Proteome: Methods for Distinguishing Expression from Usage Changes"
author:
  - Witold Wolski
  - Antje Dittmann
  - LK
  - CP
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
| Witold Wolski      | FGCZ, SIB         | [witold.wolski@fgcz.uzh.ch](mailto:witold.wolski@fgcz.uzh.ch)   |
| Antje Dittmann     | FGCZ              | [antje.dittmann@fgcz.uzh.ch](mailto:antje.dittmann@fgcz.uzh.ch)  |
| Laura Kunz         | FGCZ              |                                                        |
| Christian Panse    | FGCZ, SIB         |                                                        |
| Jonas Grossmann    | FGCZ, SIB         | [jonas.grossmann@fgcz.uzh](mailto:jonas.grossmann@fgcz.uzh)      |


- FGCZ - Functional Genomics Center Zurich, Winterthurerstrasse 190, CH-8057 Zurich, 
- SIB - Swiss Institute of Bioinformatics - Amphipole, Quartier UNIL-Sorge, CH-1015 Lausanne


## **Abstract**

**Background:** Post-translational modifications (PTMs), particularly phosphorylation, regulate protein function at substoichiometric levels, making their quantitative analysis technically challenging. A critical limitation in PTM studies is distinguishing between changes in modification abundance due to altered protein expression versus genuine changes in modification stoichiometry (the fraction of protein molecules that are modified).

**Methods:** This chapter presents a comprehensive workflow integrating TMT-based quantitative proteomics with specialized computational analysis to distinguish differential PTM expression (DPE) from differential PTM usage (DPU). The protocol encompasses automated sample preparation using Single-Pot Solid-Phase-Enhanced Sample Preparation (SP3)-based digestion, phosphopeptide enrichment with Ti-IMAC, high-pH offline fractionation, LC-MS/MS analysis, and integrated statistical analysis using the `prolfqua`, `prolfquapp` and `prophosqua` `R` packages.

**Key Features:** The workflow provides automated, scalable protocols for large sample cohorts, robust statistical frameworks for integrated PTM-protein analysis, comprehensive visualization tools including N-to-C plots for protein-centric PTM mapping, and sequence motif analysis for kinase prediction. All protocols are optimized for reproducibility and include extensive quality control measures.

**Applications:** This integrated approach enables researchers to identify true signaling changes in PTM studies, prioritize functionally relevant modification sites, understand pathway-level regulation in disease contexts, and generate robust datasets suitable for publication in high-impact journals.

## Key Words

Post-translational modifications, phosphoproteomics, TMT labeling, differential expression analysis, protein stoichiometry, mass spectrometry, prolfquapp, prophosqua, data integration, kinase activity

# 1. Introduction

## 1.1 Biological Significance of Post-Translational Modifications

Protein phosphorylation is a reversible post-translational modification (PTM) that modulates a protein’s conformation, enzymatic activity, subcellular localization, and binding interactions. Because this modification can be reversed rapidly, it acts as a molecular switch in cellular signaling networks. As a result, cells respond to environmental cues without requiring new protein synthesis. Dysregulation of phosphorylation is implicated in cancer, neurodegeneration, and metabolic disorders [10.1080/14789450.2021.1976152](https://doi.org/10.1080/14789450.2021.1976152). Large-scale phosphoproteomic profiling by high-throughput mass spectrometry (e.g., TMT-based workflows) is therefore essential to dissect disease-associated signaling alterations and identify candidate therapeutic targets.

## 1.2 Technical Challenges in PTM Analysis

The analysis of post-translational modifications (PTMs), such as phosphorylation, presents significant analytical challenges due to the low abundance, or substoichiometry, of modified peptides. Typically, phosphorylated peptides represent less than 1–5% of the total peptide population ([10.1186/1477-5956-4-15](https://doi.org/10.1186/1477-5956-4-15)). Additionally, phosphopeptides generally exhibit poorer ionization efficiency in mass spectrometry-based proteomics compared to non-modified peptides, further complicating their detection and quantification. This low abundance and challenging ionization behavior necessitate specialized experimental workflows and computational strategies to achieve accurate quantification ([10.1016/j.mcpro.2024.100754](https://doi.org/10.1016/j.mcpro.2024.100754)).

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

### 1.3.3 Data Analysis Frameworks

Proteomic data analysis can be conducted at multiple levels of granularity, each providing distinct analytical perspectives. A peptidoform represents a specific peptide sequence with a particular set of modifications, while site-level analysis aggregates signals for each residue position across multiple peptides. Differential analysis can therefore target either peptidoforms or individual modification sites, with each approach offering unique advantages for different research questions.

FragPipe introduces an additional analytical concept through multisite features, which refers to the set of all identified peptide-forms sharing the same set of modification sites under investigation, regardless of sequence derivatives, cleavage state, or whether the modifications are unambiguously localized. This approach groups peptidoforms with different sequences, which can arise from missed cleavages or alternative cleavage patterns, but identical possible modification sites into the same multisite feature, providing a balanced approach between peptidoform specificity and site-level aggregation. Discarding all peptidoforms with ambiguous localization would result in the loss of a large fraction of PTM information. Furthermore, multisite features preserve information about phosphorylation sites that occur together on the same peptide molecule, revealing co-modification patterns that indicate coordinated regulation by kinases or functional relationships between sites. Traditional site-level analysis treats each phosphorylation site independently, losing this valuable information about which specific combinations of sites are modified together in biological samples.

Site-level reports collapse signals from all peptidoforms mapping to the same PTM at an amino acid residue into a single quantitative value, ensuring each site has exactly one intensity measurement per sample. In downstream analysis, site-level intensities integrate seamlessly with site-centric enrichment tools and kinase-activity inference platforms (PhosR [10.1016/j.cels.2021.04.007](https://doi.org/10.1016/j.cels.2021.04.007), PTM-SEA [10.1016/j.molcel.2019.05.030](https://doi.org/10.1016/j.molcel.2019.05.030), RoKAI [10.1038/s41467-021-21211-6](https://doi.org/10.1038/s41467-021-21211-6), etc.), streamlining biological interpretation and network reconstruction analyses [10.1186/s12014-020-09290-x](https://doi.org/10.1186/s12014-020-09290-x).

### 1.3.4 Statistical Analysis of PTM Data

Statistical analysis of PTM data requires specialized methods due to confounding effects between PTM levels and protein abundances. Several statistical packages have emerged to address the critical issue of integrating PTM-feature quantifications with protein abundance data. Notable examples include `MSstatsPTM` [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477) and `msqrob2PTM` [10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708), which explicitly model both modified peptide and protein-level changes to distinguish genuine PTM regulation from protein abundance effects accurately.

The `msqrob2PTM` framework exemplifies the statistical approaches now available for PTM analysis, defining two complementary analytical strategies: Differential Peptidoform Abundance (DPA) and Differential Peptidoform Usage (DPU) [10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708). DPA directly models the $log^2$ intensities of each modified peptide (peptidoform), detecting absolute changes in PTM levels between experimental conditions. In contrast, DPU adjusts the PTM intensities by the corresponding protein-level changes, effectively testing for changes in the relative usage or stoichiometry of a modification site. This dual approach enables researchers to distinguish between PTM changes driven by protein abundance differences (DPA) versus those reflecting genuine changes in modification efficiency or regulatory activity (DPU).

These two frameworks differ in their analytical approaches for addressing confounding between PTM and protein abundance changes. `MSstatsPTM` employs a "model then correct" strategy, separately fitting linear models to modified and unmodified peptide data, and then combining statistical inferences post-modeling. In contrast, `msqrob2PTM` follows a "correct then model" approach, first adjusting PTM intensities by their corresponding protein abundances at the data level, and subsequently performing statistical testing on these normalized values. 

A further option, discussed in [10.1007/978-1-0716-1967-4_12](https://doi.org/10.1007/978-1-0716-1967-4_12) is to model the PTM and total proteome data as a factorial design, that jointly includes time, or another biological factor, and a "component" indicator, i.e. modified peptide and its parent protein, fits a linear model with their interaction, and then tests interaction term, to flag modified peptides whose dynamics deviate from the corresponding protein. Here, the package `limma` is used to fit the model and perform the testing.

While all three methodologies ultimately estimate PTM abundance changes corrected for protein-level effects, each of the methods can lead to different test results. 

There are also different assumptions about the relationship between the PTM and the protein abundances within a sample and differences in how variance and degrees of freedom are estimated.

- `MSstatsPTM` fits the PTM and the protein in separate models, then computes an adjusted effect as a difference of estimated contrasts: $\Delta_{adj} = 
\Delta_{PTM} − \Delta_{protein}$ with standard error combined in quadrature and Satterthwaite $df$ for the test statistic, which yields p-values for the adjusted fold change directly on the contrast scale. `MSstatsPTM` assumes that there is no correlation between PTM and the protein abundances within the samples.
- The single-model approach [10.1007/978-1-0716-1967-4_12](https://doi.org/10.1007/978-1-0716-1967-4_12) uses `limma` moderated (empirical-Bayes) variances for the interaction contrasts within one fit, which shrinks feature-wise variances toward a common prior; $df$ are moderated accordingly. The single-model formulation can, in principle, capture shared covariance (e.g., by adding a sample blocking/random effect), but as illustrated in the chapter [10.1007/978-1-0716-1967-4_12](https://doi.org/10.1007/978-1-0716-1967-4_12), the baseline model omits a sample term, effectively treating PTM and protein measurements as independent.
- `msqrob2PTM` assumes that the PTM and the protein abundances within a sample are correlated and therefore corrects PTM abundances with the matched protein abundances. Then, for modelling, only the corrected PTM intensities are used.

The three pipelines target the PTM-vs-protein effect through different modeling paths and make different assumptions about within-sample covariance, and differ in variance moderation. Therefore, even with the same raw data and contrasts, we can expect differences in estimated effects, standard errors, and p-values.

To offer researchers additional flexibility, we developed `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272), an `R` package designed to streamline statistical analysis of phosphoproteomics and other PTM datasets. Using `prophosqua' and the `prolfqua` package [10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441), the `MStatsPTM` "model first" and `msqrob2PTM` "correct first" approach can be executed to obtain protein-corrected PTM data and provide intuitive visualizations for comprehensive exploratory analysis.
 
 
### 1.3.5 The Protein Assignment Problem in Differential PTM-feature Usage 

Identification and quantification tools need to assign peptides to proteins, i.e., find the leading or representative protein ID. This assignment may be dataset-specific and may differ between the PTM and the total proteome dataset. Therefore, using protein IDs to match PTM features with proteins from the total run might lead to a loss of information because of differing representative IDs [10.1093/bib/bbr065](https://doi.org/10.1093/bib/bbr065). For example, a phosphopeptide might be assigned to protein isoform A in the PTM dataset. In contrast, in the proteome dataset, the corresponding total protein abundance is quantified under isoform B, preventing proper integration and DPU analysis.

A robust solution involves matching stripped peptide sequences from the PTM dataset to peptide sequences of proteins quantified in the total proteome dataset, bypassing protein ID dependency entirely. When peptides are shared among multiple protein isoforms, the averaged protein abundance across all matching isoforms provides a more stable reference for normalization. While this sequence-based matching approach requires additional computational implementation, it improves integration success rates and reduces the loss of valuable PTM information due to protein assignment differences between datasets.


## 1.5 Chapter Overview and Learning Objectives

This chapter introduces an analytical workflow that addresses the biological challenge of accurately distinguishing Differential PTM Expression (DPE) from Differential PTM Usage (DPU). Distinguishing between protein-level effects and genuine modification-specific changes is critical for meaningful interpretation of PTM data.

**Differential PTM Expression (DPE)** tests raw PTM signal changes between experimental conditions, identifying any modification abundance changes regardless of underlying protein-level effects. This analysis captures the total effect of experimental perturbations on PTM levels, including both direct modification effects and indirect effects mediated through protein abundance changes.

**Differential PTM Usage (DPU)** evaluates protein-normalized PTM changes, explicitly identifying modification sites where stoichiometry is genuinely altered independent of protein abundance. This analysis pinpoints sites where experimental conditions directly impact modification efficiency, kinase activity, or phosphatase activity, thereby highlighting specific regulatory events.

The integrated analytical approach combines optimized sample preparation protocols, mass spectrometry acquisition strategies, and computational tools (including the `prolfquapp` and `prophosqua` `R` packages) to provide a complete solution spanning from sample processing to biological interpretation. This framework enables researchers to extract maximum biological insight from PTM datasets while avoiding common analytical pitfalls that have historically complicated PTM data interpretation.

This protocol presents a complete experimental and analytical workflow for integrated PTM studies. The wet lab component describes the optimized protocols used at the Functional Genomics Center Zurich for PTM services, covering automated sample preparation, enrichment strategies, and LC-MS/MS acquisition. The computational analysis is demonstrated using the Atg16l1 macrophage dataset (Maculins *et al.*, eLife 2021, [10.7554/elife.62320](https://doi.org/10.7554/elife.62320)), which comprises two TMT-11-plex measurements from six conditions (WT/KO × uninfected/early/late infection) with both total proteome and phospho-enriched samples.

**Learning objectives include:**

1. Implementing automated, scalable sample preparation protocols using SP3-based digestion and KingFisher Flex automation
2. Performing phosphopeptide enrichment with Fe-NTA and optional antibody-based enrichment strategies
3. Executing high-pH offline fractionation and optimized LC-MS/MS acquisition on Evosep-Orbitrap systems
4. Conducting integrated statistical analysis of PTM and protein data using `prolfquapp` and `prophosqua` `R` packages
5. Interpreting DPE vs. DPU results for biological insights and generating publication-quality visualizations
6. Performing sequence motif analysis for kinase prediction and pathway reconstruction

The workflow emphasizes reproducibility, automation, and biological interpretation, spanning from sample lysis through data analysis and making it suitable for both experienced researchers and newcomers to integrated PTM analysis. All protocols include extensive quality control measures and troubleshooting guidelines for large-scale studies.

# 2. Materials

## 2.1 Biological Samples

### 2.1.1 Demonstration Dataset

This protocol uses the Atg16l1 macrophage dataset from Maculins et al. [10.7554/elife.62320](https://doi.org/10.7554/elife.62320) to demonstrate the bioinformatics workflow and integrated statistical analysis. The dataset comprises TMT-11-plex measurements from six conditions (WT/KO × uninfected/early/late infection) with both phospho-enriched and total proteome samples, making it ideal for illustrating the principles of integrated PTM analysis. In the original publication, the authors also performed ubiquitin remnant-enrichment, though we focus on the phosphoproteome and total proteome datasets for demonstration purposes. This same dataset was previously used in the `MSstatsPTM` publication (Kohler et al., MCP 2023, [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)), providing additional validation of the analytical approaches.

### 2.1.2 Sample Preparation

Overview: This protocol utilizes optimized sample processing methods for high-throughput phosphoproteomics analyses, featuring automated sample clean-up, digestion and enrichment steps, TMT labeling, and dual-workflow processing of both input and phosphopeptide-enriched fractions. The approach is applicable to diverse sample types and compatible with commonly used laboratory reagents (see **Supplementary Material Wetlab Section 2.1**). Key workflow components include: 

1. Sample preparation using denaturing lysis, automated protein cleanup/digestion, TMT labeling and offline-fractionation (see **Supplementary Material Wetlab Section 3.1**); 
2. Dual-workflow splitting into total proteome analysis (see Supplementary Material Wetlab Section 3.2) and phosphopeptide-enrichment (see **Supplementary Material Wetlab Section 3.3**)
2. LC-MS/MS analysis with optimized parameters for each workflow (see **Supplementary Material Wetlab Section 3.3.5  and Section 3.3.6** )

The protocols emphasize automation, scalability, and reproducibility, making them suitable for large cohort studies while maintaining analytical rigor. Expected outcomes: 80-95% protein recovery, >95% TMT labeling efficiency, 20-40 μg peptide yield per 50 μg input protein. Complete wet lab protocols, including detailed materials lists and step-by-step methods, are provided in **Supplementary Material Wetlab: Automated Sample Preparation and Wet Lab Protocols**.

## 2.2 Data analysis

1. **Search software:** FragPipe 22.0 (free download for academic usefrom [fragpipe.nesvilab.org](https://fragpipe.nesvilab.org); system requirements: 32GB+ RAM, 50GB+ storage)
2. **Statistical platform:** `R` version ≥4.0.0 ([https://www.r-project.org/](https://www.r-project.org/)) and RStudio ([https://posit.co/](https://posit.co/))
3. **R packages:**
   - `prolfqua` [10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441)
   - `prolfquapp` [10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911)
   - `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)
   - Additional dependencies: `tidyverse`, `ggseqlogo`, `writexl`, `rmarkdown`

For detailed installation instructions, see the package documentation on [github.com/fgcz/prolfqua](https://github.com/fgcz/prolfqua), [github.com/prolfqua/prolfquapp](https://github.com/prolfqua/prolfquapp), [github.com/prolfqua/prophosqua](https://github.com/prolfqua/prophosqua).

# 3. Methods

## 3.1 Mass Spectrometry Data Processing

Raw mass spectrometry (MS) data acquired using data-dependent acquisition (DDA) combined with tandem mass tag (TMT) labeling are processed using specialized software to:

1. Database searching: Spectra are searched against a protein sequence database including species-specific sequences, common contaminants, and decoy sequences (see 3.1.1)
2. Reporter ion quantification: Reporter ion intensities are extracted
3. Site localization scoring: Post-translational modification (PTM) site localization probabilities are calculated


### 3.1.1 FragPipe Method

The TMT 11-plex phospho workflow in FragPipe [10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7) 22.0 with MSFragger [10.1038/nmeth.4256](https://doi.org/10.1038/nmeth.4256) (version 4.1) was utilized.
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

For phospho-enriched samples, PTM sites with PTMProphet [10.1021/acs.jproteome.9b00205](https://doi.org/10.1021/acs.jproteome.9b00205) localization scores ≥ 0.75 are retained. FragPipe workflow parameters are stored in configuration files (`fragpipe.workflow`) (see Note 4.1).

FragPipe provides two types of PTM-features reports for the phospho-enriched samples:

- Multisite features - refers to the set of all identified peptide-forms sharing the same set of modification sites under investigation, regardless of peptide length, its derivatives, cleavage state, or whether the modifications are unambiguously localized. Hence, also peptides with different sequence, ambiguous localization but identical possible modification sites are grouped into the same **multisite** feature.

- Site-level reports collapse signals from all peptidoforms mapping to the same amino acid position into a single quantitative value per PTM site, ensuring each site has exactly one intensity per sample. In downstream analysis, site-level intensities plug into site-centric enrichment or kinase-activity inference tools (PhosR, PTM-SEA, RoKAI, etc.), streamlining biological interpretation and network reconstruction.

## 3.2 Differential Expression Analysis using `prolfquapp`

A key challenge, when analysing post-translational modifications (PTMs), is to distinguish changes in PTM abundance that are due to altered protein expression from those that reflect a change in the modification stoichiometry (i.e., the fraction of the protein pool that is modified). This protocol provides a step-by-step computational workflow to analyze PTM data in the context of total protein expression changes. We leverage the `R` packages `prolfquapp` [10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911) for streamlined differential expression analysis [10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441) and `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272) for the integration, analysis, and visualization of PTM and total proteome data. The workflow is divided into two main parts: 

1. initial differential expression analysis of the phospho-enriched (site and multisite feature abundances) and total proteome datasets (see **Supplementary Material Section A: Differential expression analysis using `prolfquapp`**)
2. the integrated analysis of differential PTM usage (see **Supplementary Material Section B: Integration and analysis of PTM features using `prophosqua`**).


### 3.2.1 Data and Software Setup

The first step involves obtaining the FragPipe 22 mass spectrometry output files [10.5281/zenodo.15850770](https://doi.org/10.5281/zenodo.15850770) and setting up the analysis environment in `R`. Depending on the sample type, different output files are required:

- **total proteome** samples `psm.tsv` files, containing peptide-spectrum match (PSM) information
- **phospho-enriched** samples:
  - `abundance_multi-site.tsv` files for multisite feature analysis
  - `abundance_single-site.tsv` files for individual modification site analysis
  

 The `prolfquapp` [10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911) package provides a set of shell scripts that automate the analysis workflow. These scripts are copied into the working directory to be used in subsequent steps. See **Supplementary Material Section A.1-A.2** for complete setup instructions and example dataset download.

### 3.2.2 Sample Annotation

A crucial step is the creation of a detailed sample annotation file. This file maps each raw data file to its experimental conditions. We parse sample names from the FragPipe output files using the `prolfqua_dataset.sh` script (see **Supplementary Material Section A.3.1**). From the sample names, we use `R` to extract experimental factors (e.g., Genotype, Timepoint). If the filenames do not encode the experimental conditions, the explanatory variables would need to be provided by adding columns to the annotation file. We then call the `annotation_add_contrasts` function in the `prolfqua` [10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441) to generate the required factorial contrasts automatically.

The `annotation_add_contrasts` function adds two key columns to the sample annotation:

- **`ContrastName`**: A descriptive name for each statistical contrast
- **`Contrast`**: The mathematical formula defining the contrast in terms of group means

For the experimental design in this protocol (2×3 factorial: Genotype [KO/WT] × Time [Uninfected/Early/Late]), the function generates four distinct contrasts:

1. **KO_vs_WT** (Main effect): `( (G_KO_Uninfect + G_KO_Late + G_KO_Early)/3 - (G_WT_Uninfect + G_WT_Late + G_WT_Early)/3 )` - Compares overall genotype effects across all timepoints
2. **KO_vs_WT_at_Uninfect** (Interaction): `G_KO_Uninfect - G_WT_Uninfect` - Compares genotypes specifically in the uninfected condition
3. **KO_vs_WT_at_Late** (Interaction): `G_KO_Late - G_WT_Late` - Compares genotypes specifically in the late infection condition  
4. **KO_vs_WT_at_Early** (Interaction): `G_KO_Early - G_WT_Early` - Compares genotypes specifically in the early infection condition

Note that setting the argument`interaction=FALSE` in the `annotation_add_contrasts` function disables the generation of interaction contrasts for factor level combinations. Interaction contrasts assess whether the effect of one factor depends on the level of another, thereby capturing non‑additive interplay between factors.

These contrasts enable comprehensive statistical analysis, distinguishing between general genotype effects and timepoint-specific effects. See **Supplementary Material Section A.3.1-A.3.2** for complete annotation workflow with code examples.


### 3.2.3 Execution of Differential Expression Analysis

Before performing differential expression analysis, a YAML configuration file is generated using the `prolfqua_yaml.sh` script. This file (`config.yaml`) defines key analysis parameters, including the normalization method (variance stabilization normalization, VSN, by default), protein summarization strategy (Tukey's median polish), and quality-filtering thresholds such as q-value cutoffs for peptide or protein-level evidence.

Once the sample annotation and configuration files are prepared, differential expression analysis is executed via the `prolfqua_dea.sh` script. This script reads the experimental design (sample IDs, group labels, and blocking factors) alongside the configuration file. Peptide or protein intensities are $\log_2$-transformed and variance-stabilized prior to modeling, improving homoscedasticity across conditions.

Linear models are fitted and contrasts are tested using the `prolfqua` `R`-package’s contrast API. Rather than imputing missing values, the model adjusts the degrees of freedom based on the number of actual observations per feature, which naturally down-weights features with high missingness and improves the reliability of the fold-change and p-value estimates. Empirical Bayes variance moderation is applied to stabilize variance estimates, particularly when sample size is limited. 

For each data type, `prolfquapp` generates a comprehensive suite of outputs:

- Dynamic HTML reports containing interactive PCA plots, boxplots, volcano plots, and heatmaps for quality control and exploratory analysis
- XLSX tables (DE_<datatype>.xlsx) summarizing fold changes, moderated t-statistics, raw p-values, and Benjamini–Hochberg–adjusted q-values
- GSEA rank files (GSEA_<datatype>.rnk) for downstream gene-set enrichment analysis
- SummarizedExperiment objects (<datatype>.rds) for seamless import into Shiny-based visualization tools

This is done separately for three data types of features derived from the mass spectrometry experiment:

1.  **Total Proteome:** Analysis of protein abundance changes
2.  **Multi-site PTM:** Analysis of multisites features abundances
3.  **Single-site PTM:** Analysis of single site features abundances


Complete execution commands and guidance for interpreting these outputs are detailed in **Supplementary Material Sections A.4–A.6**. All differential expression analysis results are publicly archived on Zenodo [10.5281/zenodo.15830988](https://doi.org/10.5281/zenodo.15830988).


## 3.3 Integrated PTM Analysis using `prophosqua`

The `prophosqua` package [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272) provides tools for integrating and analyzing post-translational modification (PTM) data with total proteome measurements. It enables researchers to distinguish between changes in protein abundance and changes in modification site usage.

This section explains how to load the results of the Differential Expression Analysis, integrate the PTM and total proteome data, and how to determine differential PTM-feature usage. For complete implementation with detailed code examples, see **Supplementary Material Section B** ("Integration and analysis of PTM features using `prophosqua`").

### 3.3.1 Data Loading and Integration

Differential expression results from `prolfquapp`—provided in Excel format—for the total proteome and the multi- or single-site PTM analyses are imported into `R`. These datasets are then integrated by performing a left join on protein identifiers, merging PTM-level statistics with the corresponding protein-level values for each condition. A left join is specifically chosen to retain all PTM features, even if their parent proteins were not quantified in the total proteome dataset. This ensures that no PTM information is lost, while still allowing normalization by protein abundance when available. See **Supplementary Material Sections B.1–B.5** for the full data loading and integration workflow.

Currently, `prolfquapp` matches PTM and protein features using exact protein IDs. This introduces a limitation: discrepancies in representative protein identifiers between the PTM-enriched and total proteome datasets may lead to loss of mapping information (see Note 4.12).

### 3.3.2 Analysis of Differential PTM-feature Expression (DPE)

**Definition:** DPE tests the raw PTM-feature signal change between conditions, without any correction for its parent protein's expression level. This analysis is used to flag any PTM-feature whose abundance changes, even if the parent protein itself is also up- or down-regulated.

**Method:** DPE results are visualized using N-to-C plots, generated by the `prophosqua` [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272) `n_to_c_expression` function. These plots map both the $\log_2$ fold changes of individual modification sites along the primary sequence of the protein and the $\log_2$ fold change of the protein abundance from the total proteome experiment. This provides a clear visual summary of both site-specific and overall protein expression changes. See **Supplementary Material Section B.6** for complete DPE analysis implementation and visualization examples. The N-to-C expression plot shown below, for example, shows that the protein SQSTM1 (UniProt id Q64337) and the phosphorylation sites are upregulated in the KO vs WT comparison.


![N-to-C expression plot showing differential expression of total protein abundance (light yellow rectangle) and phosphorylation sites (vertical lines, color coded by type of residue) along protein sequences. The vertical lines are either dashed or solid, depending on the site's imputation status. The x-axis represents amino acid position from N to C terminus, while the y-axis shows $\log_2$ fold changes. Significantly regulated sites (FDR < 0.05 and FDR < 0.01) are highlighted with one or two red asterisks, respectively. The plot enables visualization of both overall protein regulation and site-specific phosphorylation changes in a single view.](figures/dpe_n_to_c.png)


### 3.3.3 Analysis of Differential PTM-feature Usage (DPU)

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


The PTM-feature usage fold changes, p-values, and FDR estimates are then added to the dataframe. The DPU results are visualized using N-to-C plots generated by `prophosqua::n_to_c_usage`, which display the protein-normalized PTM-feature abundance changes and highlight sites with significant changes in usage. See **Supplementary Material Sections B.7-B.8** for complete DPU analysis implementation and statistical framework details.


For example, when comparing the N-to-C expression plot shown above, with the N-to-C usage plot, shown below, for Protein SQSTM1 (UniProt id Q64337), we can see that while raw PTM-feature expression is upregulated in the KO vs WT contrast, the actual usage of the phosphorylation sites remains unchanged or even is significantly downregulated, after normalizing for protein abundance (shown by the vertical lines with negative log2 fold changes). This demonstrates how DPU analysis can reveal changes in modification stoichiometry that differ from, or even oppose, changes in overall PTM-feature levels. 


![N-to-C usage plot showing differential usage of phosphorylation sites (vertical lines, color coded by type of residue) along the protein sequences. The vertical lines are either dashed or solid, depending on the site's imputation status. The x-axis represents amino acid position from N to C terminus, while the y-axis shows $\log_2$ fold changes. Significantly regulated sites (FDR < 0.05 and FDR < 0.01) are highlighted with one or two red asterisks, respectively. The plot enables visualization of both overall protein regulation and site-specific phosphorylation changes in a single view.](figures/dpu_n_to_c.png)


Because generating N-to-C plots is computationally intensive, we generate plots only for proteins with at least one significant phosphorylation site or multi-site feature. The FDR threshold for determining significance can be adjusted in the `n_to_c_expression` or `n_to_c_usage` functions, allowing researchers to balance computational efficiency with comprehensive visualization based on their specific statistical considerations and interpretation guidelines. When no significant PTM sites are detected despite biological expectations, one possible reason is insufficient statistical power, which can be mitigated by increasing the sample size (see Note 4.13).

In the *Atg16l1* macrophage dataset [10.7554/elife.62320](https://doi.org/10.7554/elife.62320), SQSTM1 (uniprot id Q64337) protein levels are markedly increased in Atg16l1-knockout (KO) versus wild-type (WT), reflecting loss of autophagy-mediated turnover. Phosphosite analysis reveals that the raw phosphopeptide intensities at $Ser_28$, $Ser_308$, and $Ser_368$ are similarly elevated in cKO compared to WT (DPE), consistent with higher protein abundance. However, when applying DPU normalization, which subtracts the protein $log_2$ fold change from each site's phospho $log_2$ fold change, these apparent phosphorylation increases are substantially reduced, indicating that most of the phosphosite enrichment arises from protein accumulation rather than genuine changes in modification stoichiometry. This example highlights how DPU effectively isolates regulatory events on SQSTM1 from confounding protein-level effects. 

When evaluating statistical significance in DPU analysis, it is important to avoid overinterpretation of results (see Note 4.13). Statistical significance does not guarantee biological relevance, and effect sizes should be considered alongside FDR values, when prioritizing sites for further investigation. Key findings should be validated through orthogonal methods, to confirm biological relevance.


### 3.3.4 Integration Report Generation

The integrated analysis results are compiled into a comprehensive HTML report using an R Markdown template `_Overview_PhosphoAndIntegration_site.Rmd`. The report includes:

- Overview statistics on identified phosphorylation sites
- Interactive scatter plots comparing protein vs PTM-feature changes
- Volcano plots highlighting differential PTM-feature usage

The report includes data tables that allow searching for a specific protein and phosphorylation feature and highlighting it in the scatter and volcano plots. See **Supplementary Material Sections B.9-B.11** for complete report generation workflow and output descriptions.

When analyzing the interactive reports, different combinations of DPE and DPU results provide distinct biological insights that guide interpretation. For instance, **DPE+ only** indicates changes driven primarily by protein abundance alterations, where PTM increases reflect higher protein levels rather than enhanced modification efficiency (see Note 4.15).

### 3.3.5 Exporting Results to Excel

The integrated analysis results are exported to an Excel file for further analysis and sharing. The Excel file contains multiple worksheets:

- **combinedSiteProteinData**: Contains the merged PTM and protein-level data used for the integrated analysis
- **combinedStats**: Contains the differential PTM-feature usage statistics, including fold changes, p-values, and FDR estimates for each phosphorylation site

This Excel format facilitates downstream analysis, data sharing with collaborators, and integration with other bioinformatics tools. The file is automatically saved with a timestamped filename in the results directory.


### 3.3.6 Sequence Motif Analysis

To identify potential kinases responsible for the observed phosphorylation changes, sequence motif analysis is performed on significantly regulated PTM sites identified in the DPU analysis.

For each contrast, amino acid sequences flanking the significantly regulated sites (e.g., $FDR < 0.01$) are extracted and grouped by regulation status (upregulated or downregulated). The `ggseqlogo` `R` package is then used to generate sequence logo plots for each group. These plots visualize conserved amino acid patterns around the modification site, facilitating comparison to known kinase recognition motifs and enabling inference of upstream regulatory kinases.

When motif analysis yields weak or absent sequence patterns, several factors may require optimization including the number of significant sites, or the sequence window size used for motif extraction (see Note 4.16). See **Supplementary Material Section C** for complete sequence motif analysis implementation and kinase prediction workflow.

![Sequence logo plots showing amino acid motifs surrounding significantly regulated phosphorylation at serine residues. The height of each letter indicates the relative frequency of each amino acid at each position surrounding the phosphorylation site (position 8). Separate logos are generated for up- and downregulated sites within each experimental contrast. Common kinase recognition motifs can be identified from these patterns to suggest upstream regulators.](figures/sequence_logo_plot.png)


# 4. Notes

## 4.1 FragPipe workflow parameters

FragPipe workflow parameters are stored in configuration files (`fragpipe.workflow`), which can be customized for different experimental designs. The exact configuration file used to process the data in this protocol is available for download from Zenodo [10.5281/zenodo.15850770](https://doi.org/10.5281/zenodo.15850770) or [gitlab.bfabric.org/wolski/PTM_example](https://gitlab.bfabric.org/wolski/PTM_example).

## 4.11 Automated QC reporting with FragPipe TMT QC script

An automated quality control script generates comprehensive TMT labeling efficiency reports directly from FragPipe PSM output files. The script evaluates labeling completeness at both peptide and PSM levels for N-terminal and lysine modifications, calculates missed cleavage rates for tryptic digestion efficiency, and provides quantitative channel balance assessment across all TMT channels. Key metrics include percentage of modified N-termini and lysine residues (expected >95%), missed cleavage rates for lysine and arginine residues (typically 5-15%), and relative abundance distributions across channels to detect loading imbalances. The script outputs interactive HTML reports with visualization of identification numbers per channel, modification frequencies, total abundance distributions, and density plots for rapid assessment of data quality before proceeding with downstream analysis. The QC script is part of the `prophosqua` `R`-package vignettes [10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272).


## 4.12 Match Rates

**Check match rates:** $>80\%$ between PTM and protein datasets indicates good integration. Lower rates may be due to protein assignment differences in the PTM and total proteome datasets, and sequence-based matching could improve the match rates.

## 4.13 No significant PTM sites detected despite biological expectation

No significant PTM sites detected despite biological expectation can be due to:

- **Insufficient statistical power:** Increase sample size or reduce technical variability
- **Inappropriate statistical thresholds:** Consider $FDR = 0.25$ instead of $0.05$ for initial exploration
- **Normalization issues:** Check for batch effects, verify sample loading consistency
- **Outlier detection:** Check for outliers in the data, e.g., using density or scatter plots in `prolfquapp` generated reports.
- **Experimental design:** Ensure adequate biological replicates (minimum $n=3$ per group)
- **Quality control check:** Examine volcano plots and distribution of p-values for expected patterns in `prolfquapp` qc reports.

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



# 6. Abbreviations

| Abbreviation | Definition |
|--------------|------------|
| BSA | Bovine Serum Albumin |
| ClAA | Chloroacetamide |
| CV | Coefficient of Variation |
| DDA | Data-Dependent Acquisition |
| DEA | Differential Expression Analysis |
| DIA | Data-Independent Acquisition |
| DPA | Differential Peptidoform Abundance |
| DPE | Differential Post-Translational Modification Expression |
| DPU | Differential Post-Translational Modification Usage |
| EtOH | Ethanol |
| FC | Fraction Collector |
| FDR | False Discovery Rate |
| FGCZ | Functional Genomics Center Zurich |
| HPLC | High-Performance Liquid Chromatography |
| ID | Identification |
| IMAC | Immobilized Metal Affinity Chromatography |
| IP | Immunoprecipitation |
| LC | Liquid Chromatography |
| MOAC | Metal Oxide Affinity Chromatography |
| MS | Mass Spectrometry |
| N-to-C | N-terminus to C-terminus |
| NTA | Nitrilotriacetic Acid |
| PSM | Peptide-Spectrum Match |
| PTM | Post-Translational Modification |
| QC | Quality Control |
| RT | Room Temperature / Retention Time |
| SDS | Sodium Dodecyl Sulfate |
| SIB | Swiss Institute of Bioinformatics |
| SP3 | Single-Pot Solid-Phase-Enhanced Sample Preparation |
| SPE | Solid-Phase Extraction |
| TCEP | Tris(2-carboxyethyl)phosphine |
| TEAB | Tetraethylammonium Bicarbonate |
| TMT | Tandem Mass Tag |
| UHPLC | Ultra-High-Performance Liquid Chromatography |
| UV | Ultraviolet |

**References**

