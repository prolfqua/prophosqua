# Outline for Methods in Molecular Biology Article

## **Abstract**

Phosphoproteomics studies are crucial for our understanding of cellular signaling pathways and the molecular mechanisms underlying biological processes and their deregulation in disease. Utilizing Tandem Mass Tag (TMT) labeling significantly enhances quantification sensitivity by allowing the simultaneous analysis of multiple samples, improving compatibility with phosphopeptide enrichment and fractionation strategies, and minimizing issues with missing data. Offline fractionation before LC-MS/MS is crucial for achieving high sensitivity in phosphopeptide analysis. Our workflow encompasses both enriched phosphopeptide samples and whole-proteome analysis to accurately normalize phosphorylation data, correcting fold changes via Differential Protein Usage (DPU) to distinguish between PTM-specific changes and overall variations in protein abundance. We leverage FragPipe software tools for robust peptide identification, quantification, and bioinformatics analysis. Our approach, emphasizing flexibility, reproducibility, and automation, streamlines phosphoproteomics studies, facilitating rapid, accurate, and scalable PTM quantification.

# 1. Introduction

Protein phosphorylation is a reversible post‐translational modification (PTM) that modulates protein conformation, activity, localization, and interactions, acting as a key switch in cellular signaling networks ([10.1080/14789450.2021.1976152](https://doi.org/10.1080/14789450.2021.1976152)). Dysregulated phosphorylation underlies many diseases; therefore, large-scale phosphoproteomic profiling by mass spectrometry provides critical insight into dynamic signaling pathways. However, the low stoichiometry of phosphopeptides necessitates the use of specific workflows for robust quantification. Isobaric tagging (TMT) has become a standard approach for multiplexed PTM analysis. TMT reagents label peptides from up to 35 samples in a single combined run, thereby significantly reducing inter-sample technical variability ([10.1021/acs.jproteome.4c00668](https://doi.org/10.1021/acs.jproteome.4c00668)). In practice, each sample is labeled with a distinct TMT reagent, and the labeled digests are then pooled and co-analyzed by LC-MS/MS. In this way, one can achieve highly reproducible relative quantification across many conditions. For example, Maculins *et al.* (eLife 2021\) performed two 11‐plex TMT experiments on wild‐type and *Atg16l1*‐knockout murine macrophages (6 conditions: uninfected or early/late *Shigella* infection, 3–4 replicates each) ([10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320)). This dataset includes both total proteome and phosphoproteome measurements, enabling integrated PTM analyses. Enrichment of phosphorylated peptides is typically achieved through metal affinity chromatography, which exploits the high affinity of phosphate groups for positively charged metal cations or metal oxides. Standard protocols include immobilized metal affinity chromatography (IMAC) or metal oxide affinity chromatography (MOAC), mainly employing Fe³⁺, Zr4+, or Ti4+ ions (IMAC) or TiO2 (MOAC). The different chemistries generally yield similar depths of phosphopeptide coverage but also show some degree of complementarity. Therefore, the choice between them often depends on project-specific or practical considerations, with either approach being compatible with large-scale phosphoproteomic studies. 

Computationally, numerous software tools support the analysis of quantitative protein post-translational modifications (PTMs). For DDA (data‐dependent acquisition) data, platforms like FragPipe ([10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7)) (Yu *et al.*, *Nat Commun* 2023\) integrate the MSFragger search engine ([10.1038/nmeth.4256](https://doi.org/10.1038/nmeth.4256)) and are widely used for TMT and label‐free workflows. FragPipe can process TMT‐labeled data and output site‐localized quantifications. For DIA (data-independent acquisition) data, tools such as Spectronaut (Biognosys) and DIA-NN (Demichev *et* *al.*, *Nat Methods* 2020, [10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x)) provide PTM support. Both of these tools implement site-localization scoring and generate site-level quantification reports ([10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x)). 

After peptide identification and quantification, statistical packages are used to detect regulated modifications. MSstatsPTM (Kohler *et al.*, *Mol. Cell. Proteomics 2023, [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)) models both modified and unmodified peptides to account for protein-level* changes. Likewise, the R-based msqrob2PTM workflow (Demeulemeester *et al.*, *Mol. Cell. Proteomics 2024, [10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708)) explicitly performs differential analysis at the peptidoform level. These approaches correct for confounding between PTM and protein abundance, providing statistical tests for modification changes. Proteomic data can be analyzed at different levels: a *peptidoform* is a specific peptide sequence with a particular set of modifications, whereas a *site*\-level analysis aggregates signals for each residue position across peptides. Differential analysis can thus target peptidoforms or sites. For example, msqrob2PTM defines two complementary tests: Differential Peptidoform Abundance (DPA) and Differential Peptidoform Usage (DPU) ([10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708)). DPA directly models the log₂ intensities of each modified peptide (peptidoform), detecting absolute changes in PTM levels. In contrast, DPU adjusts the PTM intensities by the corresponding protein-level changes, effectively testing for changes in the relative usage of a site. In short, DPA flags peptidoforms whose abundance changes, while DPU flags modifications whose change exceeds that expected from the parent protein’s change.

* Differential Peptidoform Abundance (DPA): tests raw PTM peptide intensity differences  
* Differential Peptidoform Usage (DPU): corrects PTM intensities for protein-level variation to highlight proper site-specific regulation

In this chapter, we apply these methods to the *Atg16l1* macrophage dataset. This MassIVE MSV000085565 dataset (Maculins *et al.*, eLife 2021, [10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320)) comprises two TMT-11-plex runs of six conditions (WT/KO × \[uninfected, early, late\] infection, with 3–4 replicates per condition). It includes both global proteome and phosphopeptide‐enriched runs. This dataset was also analyzed in the development of MSstatsPTM. Notably, the *Atg16l1* knockout results in a near-complete loss of the Atg16l1 protein itself, illustrating that extreme protein-level changes (e.g., a null allele) must be considered when interpreting PTM results. By integrating TMT‐based quantification, phospho‐enrichment, and advanced statistics (e.g,. DPA/DPU tests), we aim to robustly identify actual phosphorylation changes across conditions without bias from protein abundance shifts. 

References: Key methods and tools are detailed in \[Maculins *et al.* eLife 2021 ([10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320))\], \[Kohler *et al.* MCP 2023 ([10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477))\], \[Demeulemeester *et al.* MCP 2024 ([10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708))\], and software documentation for FragPipe and DIA-NN

# 2. Materials

## 2.1 Biological Samples

As described above, we are focusing on the dataset from Maculins et al. ([10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320)) in this section. We are using the phospho-enriched and the total proteome datasets for demonstration purposes. In the original publication, the authors also performed a KGG-enrichment (ubiquitin remnant). Furthermore, the same dataset was also used in the MSstatsPTM publication (MCP, Kohler et. al, 2023, 10.1016/j.mcpro.2022.100477). The material and methods that were used in this publication are not identical to what we suggest here. Nevertheless, the most important parts are that chemical labeling was used and there is a PTM-enriched dataset with potentially multiple fractions (phospho-enriched here) as well as a total proteome part that can be integrated. The main differences in sample preparation between the proposed method and the method used by Maculins et al. are summarized in the table. 

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

1. Lysis buffer: 50 mM Tris-HCl, pH 7.5, 4 % (w/v) SDS, optional: PhosSTOP (Merck, Germany)  
2. Optional: 5 U/µL Benzonase Nuclease (Merck, Germany)  
3. Tissue Lyzer II (Qiagen, Germany)  
4. Probe-type sonicator UP200St (Hielscher, Germany)  
5. Orbital benchtop mixer for microtubes, e.g., Eppendorf Thermomixer C  
6. UV/Vis spectrophotometer, e.g,. Lunatic (Unchained Labs, USA)  
7. Speedmag carboxylated magnetic beads (hydrophilic and hydrophobic, Merck, Germany)  
8. Multichannel pipette  
9. Microtube centrifuge, e.g., Eppendorf 5424  
10. KingFisher Flex with deep-well head (Thermo, USA)  
11. KingFisher plasticware (deep-well plate, microwell plate, tip comb)  
12. BindIT 4.0 Software  
13. 250 mM tris(2-carboxyethyl)phosphine (TCEP)  
14. 500 mM chloroacetamide (ClAA)  
15. Trypsin, mass spectrometry grade  
16. 80 % and 100 % Ethanol (EtOH)  
17. LC-grade water  
18. Anhydrous acetonitrile  
19. 50 mM Tetraethylammonium bicarbonate (TEAB), pH 8.5  
20. Stop solution: 5% hydroxylamine  
21. TMTpro reagent, e.g., 18plex (Thermo, USA)  
22. SpeedVac

    

### 2.2.2 Optional antibody-based enrichment

1. PTMscan antibody targeting site and/or kinase substrate motif of interest (e.g., PTMScan HS Phospho-Tyrosine, Cell Signaling, USA)  
2. High Select Fe-NTA Phosphopeptide Enrichment Kit (Thermo, USA)  
3. IP binding buffer: 50 mM Tris-HCl, pH 7.5, 1% NP40 (IGEPAL CA-630)  
4. IP wash buffer: 50 mM Tris-HCl  
5. End-over-end shaker

   

### 2.2.3 High pH offline fractionation

1. Vanquish Flex UHPLC with fraction collector and 100 µL sample loop  
2. Chromeleon 7.3.1 instrument control software  
3. BSA peptide standard (500 μmol dry, diluted in 600 μl, 83.3 μmol on column, in 100 μl)  
4. High pH buffer A: 2.5mM ammonium hydroxide in LC-grade water  
5. High pH buffer B: 95% acetonitrile/10mM ammonium hydroxide in LC-grade water  
6. XBridge peptide BEH C18 130Å, 3.5 µm, 4.6 mm X 250 mm (Waters)  
7. Storage solution: 20% Methanol  
8. SpeedVac compatible with deep-well plates (e.g., Savant SpeedVac Vacuum Concentrator SPD210, Thermo)


### 2.2.4 Phosphopeptide enrichment and Evotip loading

1. KingFisher Flex as specified in 2.2.1  
2. Ti-IMAC buffer I: 80% acetonitrile, 5% trifluoroacetic acid, 0.1 M Glycolic acid  
3. Ti-IMAC buffer II: 80% acetonitrile, 1% trifluoroacetic acid  
4. Ti-IMAC buffer III: 10% acetonitrile, 0.2% trifluoroacetic acid  
5. Ti-IMAC elution buffer: 1% ammonium hydroxide  
6. Ti-IMAC neutralization buffer: 20% formic acid  
7. ReSyn Biosciences Ti-IMAC HP beads  
8. Evotip pure (Evosep, Denmark)

   

### 2.2.5 LC-MS/MS analysis

1. Evosep One (Evosep, Denmark)  
2. Buffer A: 0.1% formic acid in LC-grade water  
3. Buffer B: 99.9% acetonitrile/0.1% formic acid in LC-grade water  
4. Separation column, e.g., Aurora Elite XT 15x75 C18 UHPLC column (IonOpticks, Australia)  
5. Mass spectrometer equipped with nano-spray ion source, e.g. EASY-Spray Source, Orbitrap Exploris 480 (Thermo, USA)  
6. Column heater and heater controller (IonOpticks, Australia)

   

### 2.2.6 Data analysis

1. FragPipe  22.0  [10.1101/2025.05.27.656447](https://doi.org/10.1101/2025.05.27.656447)
2. R [10.1080/10618600.1996.10474713](https://doi.org/10.1080/10618600.1996.10474713) and RStudio  [https://posit.co/](https://posit.co/)
3. R packages: `prolfquapp` ([10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911)), `prolfqua` ([10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441)), `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) (and its dependencies)

# 3. Methods

## 3.1 Sample Preparation

### 3.1.1 Cell lysis, digestion, and labeling of tryptic digests with TMT reagents

This protocol describes a generic lysis approach for fresh-frozen tissue or frozen cell pellets, utilizing denaturing SDS for enzymatic deactivation and efficient membrane disruption. It has been successfully used for most eukaryotic systems, including yeast and plants (see Note 4.1). Protein clean-up and digestion are based on the automated SP3 protocol, as described in Leutert et al. ([10.15252/msb.20199021](https://doi.org/10.15252/msb.20199021)). Still, they can be performed with any protocol that efficiently removes all components incompatible with downstream processing steps, including buffers containing primary amines.

      

1. Add up to 50 μL lysis buffer for small pellets (≤10 μL pellet volume) and more for larger pellets, using 3x – 5x the pellet volume and completely resuspend and disrupt cells in buffer by pipetting up and down; apply 1 min Tissue Lyzer (30 Hz) cycles with glass beads if cell pellet is not completely disrupted and lysed  
2. For tissue samples, add 300 μL lysis buffer to up to 5 mg of tissue (optimal volume depends on tissue origin and protein content) and disrupt the sample using the Tissue Lyzer II with 2x 1 min cycles (30 Hz)  
3. Apply a 1-minute sonication pulse to the sample (probe-type sonicator)  
4. Incubate 5min at 95 °C  
5. Apply a 1-minute sonication pulse to the sample (probe-type sonicator)  
6. If viscosity is still high due to a high DNA content, dilute the sample to 1 % SDS and add 5 – 10 U of Benzonase and incubate for 15 min at room temperature  
7. Centrifuge 5 min at 20’000 x g and proceed with supernatant  
8. Take a small aliquot for protein concentration determination and dilute at least 1:20 for 280 nm reading on the Lunatic  
9. Snap-freeze samples in liquid nitrogen and store at \-20°C or proceed with the protocol for 30 – 50 µg of total protein with 0.5 – 2 µg/µL   
10. Add TCEP to a final concentration of 5 mM  
11. Add ClAA to a final concentration of 15 mM   
12. Incubate 30 min at 30 °C in the dark, shaking at 500 rpm  
13. Prepare carboxylated magnetic beads (see Note 4.2)  
    1. Use 10 ug of beads for each ug of protein  
    2. Beads are available at a stock concentration of 50 μg/μL. Take the required volume of beads and mix hydrophilic and hydrophobic beads at a 1:1 ratio  
    3. Wash beads three times with water and add the respective amount in 250 µL to each well  
14. Add 100% EtOH to the reduced and alkylated sample to reach 60% EtOH (v/v) (1.5x of reduced and alkylated sample volume)  
15. Fill three deep-well plates with 500 µL wash buffer (80 % EtOH)  
16. Add 100 µL of trypsin in 20 mM TEAB at an enzyme:protein ratio of 1:50 to a microwell plate  
17. Process samples with KingFisher protocol as described in Leutert et al. (https://doi.org/10.15252/msb.20199021) with minor changes  
    1. After transferring beads into the trypsin digestion buffer, pause the KingFisher protocol, transfer the plate to the Thermomixer heated to 37°C, and incubate overnight while shaking at 500 rpm (see Note 4.3)  
    2. The next day, return the digestion plate to its original KingFisher slot, add 100 µL water to a second microwell plate, and place it in the empty KingFisher slot   
    3. Continue with the  KingFisher protocol for bead collection and transfer into the water elution plate  
    4. Pool digest, water elution, and dry samples in SpeedVac (peptides can be stored frozen at \-80 °C)  
18. Resuspend peptides in 45 µL 50 mM TEAB and sonicate in a water bath for 10 min   
19. Resuspend 100-200 µg TMTpro reagent in 15 µL anhydrous acetonitrile, vortex, and spin down (see Note 4.4)  
20. Add TMTpro reagents to the respective sample, mix, quickly spin down, and incubate for  60 min at RT in the dark  
21. Add 3.5 µL of 5 % hydroxylamine to quench the reaction  
22. Pool a small fraction of each sample (e.g. 1uL, ensure precise pipetting) and pool at equal ratio to check for sample loading variation and labeling efficiency using LC-MS/MS in QC step  
23. If needed, adjust the pooling ratio based on QC LC-MS/MS, otherwise, proceed with pooling all samples at an equal ratio  
24. Optional: if antibody-based enrichment is required, split the sample in half  
25. Dry down in SpeedVac (can be stored frozen at \-80°C)  
26. Optional: For final clean-up, subject sample to C18 solid-phase extraction (SPE) (e.g., Sep-Pak, Waters) 

### 3.1.2 Optional antibody-based enrichment

This protocol describes the targeted enrichment of peptides harboring specific kinase substrate motifs or phosphorylated tyrosines using PTMScan HS antibody beads from Cell Signaling Technology. Binding and wash buffers have been modified, and an additional second enrichment step using Fe-NTA spin columns is used to further increase phospho-peptide enrichment specificity (see Note 4.5).

1. Resuspend one half of the pooled labeled peptides in 400 µL cold IP binding buffer, vortex, and incubate rotating at 4 °C for 5 min, check pH (should be neutral or slightly basic, not below pH 7\)  
2. Spin down solution 5 min at 20’000 x g to clear the solution and pellet insoluble material  
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
13. Collect beads and transfer elution to High Select Fe-NTA spin column conditioned as per the manufacturer’s instructions  
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
11. Remove the deep-well plate from the fraction collector, snap-freeze, and dry down the samples in a SpeedVac. Peptides can be stored at \-80 °.C

    

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
2. On MS1 level, peptide masses were detected at a resolution of 120’000, with an ion target of 3E6, maximal injection time of 45 ms, and RF lens set to 40 %  
3. MS2 spectra were recorded for the top 12 precursors with an intensity threshold of 10^3, which were isolated at 0.7 Th and subjected to dynamic exclusion for 16 s. Normalized collision energy was set to 32 % and spectra were recorded with a resolution of 45’000 (phosphopeptides) or 30’000 and TMTpro on (total proteome), ion target of 1E5 and a maximal injection time of 250 ms (phosphopeptides) or Auto (total proteome)


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
- Normalizing intensities by retention time bins and applying log₂ transformations.
- Removing outliers based on interquartile range filtering.
- Aggregating quantified values to protein or PTM-site levels.

Selected TMTIntegrator parameters are:

- Best PSM selection
- Outlier removal enabled
- Allow overlabel/underlabel: OFF
- Use MS1 intensity: OFF
- Aggregation method: Median

For phospho-enriched samples, PTM sites with PTMProphet ([10.1038/s41467-020-17914-x](https://doi.org/10.1038/s41467-020-17914-x)) localization scores ≥ 0.75 are retained.

Note: FragPipe workflow parameters are stored in configuration files, which can be customized for different experimental designs. The exact configuration file (fragpipe_TMT16_phospho.config) is available as supplementary material in PTM_experiment.zip.

## 3.4 Differential Expression Analysis using `prolfquapp`

Post-translational modifications (PTMs) play a crucial role in regulating protein function; however, their analysis is complex and challenging. A key challenge is to distinguish changes in PTM abundance that are due to altered protein expression from those that reflect a change in the modification stoichiometry (i.e., the fraction of the protein pool that is modified). This protocol provides a step-by-step computational workflow to analyze PTM data in the context of total protein expression changes. We leverage the R packages `prolfquapp` ([10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911)) for streamlined differential expression analysis and `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) for the integration, analysis, and visualization of PTM and total proteome data. The workflow is divided into two main parts: (1) initial differential expression analysis of separate datasets (see **Supplementary Material Section A: Differential expression analysis using `prolfquapp`**), and (2) the integrated analysis of differential PTM usage (see **Supplementary Material Section B: Integration and analysis of PTM features using `prophosqua`**).


### 3.4.1 Data and Software Setup

The first step involves obtaining the mass spectrometry output files and setting up the analysis environment. The `prolfquapp` ([10.1021/acs.jproteome.4c00911](https://doi.org/10.1021/acs.jproteome.4c00911)) package provides a set of shell scripts that automate the analysis workflow. These scripts are copied into the working directory to be used in subsequent steps. See **Supplementary Material Section A.1-A.2** for complete setup instructions and example dataset download.

### 3.4.2 Sample Annotation

A crucial step is the creation of a detailed sample annotation file. This file maps each raw data file to its experimental conditions and defines the statistical comparisons to be made. For this protocol, we parse sample names to extract experimental factors (e.g., Genotype, Timepoint). The `prolfqua` ([10.1021/acs.jproteome.2c00441](https://doi.org/10.1021/acs.jproteome.2c00441)) `::annotation_add_contrasts` function is then used to automatically generate the statistical contrasts required for a factorial analysis design. See **Supplementary Material Section A.3** for complete annotation workflow with code examples.

### 3.4.3 Execution of Differential Expression Analysis

With the annotation file in place, we execute the `prolfqua_dea.sh` command-line script to perform the differential expression analysis. This is done separately for three data types derived from the mass spectrometry experiment:

1.  **Total Proteome:** Analysis of protein abundance changes
2.  **Multi-site PTM:** Analysis of peptides containing multiple modification sites
3.  **Single-site PTM:** Analysis of peptides containing single, localized modification sites

This step generates a set of detailed reports, including statistical tables (p-values, fold changes) and quality control plots for each data type. See **Supplementary Material Section A.4-A.6** for complete execution commands and output interpretation.

## 3.5 Integrated PTM Analysis using `prophosqua`

This part describes how to load the results from Part 1, integrate the PTM and total proteome data, and perform an in-depth analysis of PTM dynamics. For complete implementation with detailed code examples, see **Supplementary Material Section B** ("Integration and analysis of PTM features using `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272))").

### 3.5.1 Data Loading and Integration

The differential expression results, typically in Excel format, from the total proteome and single-site PTM analyses are loaded into the R environment. The two datasets are then integrated by performing a left-join operation, merging the PTM-level statistics with the corresponding protein-level statistics for each condition. See **Supplementary Material Sections B.1-B.5** for complete data loading and integration workflow.

### 3.5.2 Analysis of Differential PTM-feature Expression (DPE)

**Definition:** DPE tests the raw PTM-feature signal change between conditions, without any correction for its parent protein's expression level. This analysis is used to flag any PTM-feature whose abundance changes, even if the parent protein itself is also up- or down-regulated.

**Method:** DPE results are visualized using N-to-C plots, generated by the `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) `::n_to_c_expression` function. These plots map the log2 fold changes of individual modification sites along the primary sequence of the protein, providing a clear visual summary of site-specific expression changes. See **Supplementary Material Section B.6** for complete DPE analysis implementation and visualization examples.

### 3.5.3 Analysis of Differential PTM-feature Usage (DPU)

**Definition:** DPU tests the **protein-normalized** changes of PTM-features. This analysis is essential for determining whether the *fraction* of a protein that is modified at a specific site changes between conditions. It isolates changes in modification stoichiometry from changes in overall protein abundance.

**Method:** DPU is calculated by subtracting the protein's log2 fold change from the PTM-feature's log2 fold change. The associated p-values are re-calculated using a method that combines the variance from both the PTM and protein-level models, as implemented in `prophosqua` ([10.5281/zenodo.15845272](https://doi.org/10.5281/zenodo.15845272)) `::test_diff`. The DPU results are then visualized using N-to-C plots generated by `prophosqua::n_to_c_usage`, which display the protein-normalized changes and highlight sites with significant changes in usage. See **Supplementary Material Sections B.7-B.8** for complete DPU analysis implementation and statistical framework details.

### 3.5.3.1 Integration Report Generation

The integrated analysis results are compiled into a comprehensive HTML report using an R Markdown template `_Overview_PhosphoAndIntegration_site.Rmd`. The report includes:

- Overview statistics on identified phosphorylation sites
- Interactive scatter plots comparing protein vs PTM-level changes
- Volcano plots highlighting significant DPU sites

The report is generated using the `_Overview_PhosphoAndIntegration_site.Rmd` template, which creates publication-ready figures and tables while maintaining reproducibility. The report includes both raw data tables and interactive visualizations to facilitate data exploration and biological interpretation. See **Supplementary Material Sections B.9-B.11** for complete report generation workflow and output descriptions.

### 3.5.4 Sequence Motif Analysis

To infer which kinases may be responsible for the observed changes in phosphorylation, a sequence motif analysis is performed on sites with significant DPU.

**Method:** Amino acid sequences surrounding significantly regulated PTM sites (e.g., FDR < 0.01) are extracted. These sequences are then grouped by their regulation status (up- or downregulated) and by the experimental contrast. The `ggseqlogo` R package is used to generate sequence logo plots from these groups. These plots visualize conserved amino acid patterns, which can be compared to known kinase recognition motifs to identify potential upstream regulators. See **Supplementary Material Section C** for complete sequence motif analysis implementation and kinase prediction workflow.

## Results and Reporting

The entire workflow is captured in an R Markdown document, which generates a comprehensive HTML report. This report includes all statistical tables, quality control plots, DPE and DPU visualizations, and sequence logo analyses, providing a complete and reproducible record of the integrated PTM analysis. See **Supplementary Material Sections B.9-B.11** for complete report generation workflow and output file descriptions. 
## 3.5 Validation and Quality Control

* QC checks for TMT labeling efficiency
* Reproducibility assessments (technical replicates)
* Data quality metrics and visualization
* Cross-validation with orthogonal methods

# 4 Troubleshooting

* Common issues in phosphopeptide enrichment
* Tips for improving MS sensitivity and quantification accuracy
* Optimization of enrichment protocols
* Data quality improvement strategies

# 4. Notes

   1. It is critical to minimize the time of sample exposure to non-denaturing conditions after sample harvesting, as enzymes such as phosphatases and kinases can significantly confound phosphoproteomic studies, leading to artifacts and biases ([10.1158/0008-5472.CAN-14-2309](https://doi.org/10.1158/0008-5472.CAN-14-2309)) by affecting the stability of phosphorylation marks. Studies that include sampling of early signaling dynamics should include fast lysing and denaturing strategies ([10.1073/pnas.1521288113](https://doi.org/10.1073/pnas.1521288113)).  
   2. Microsphere-based protein clean-up and digestion strategies can be based on a variety of bead surface chemistries that exploit protein aggregation ([10.1074/mcp.TIR118.001270](https://doi.org/10.1074/mcp.TIR118.001270)).   
   3. The conditions described result in typical digestion efficiencies of around 20%. If there is a decrease in efficiency, Lys C can be used to improve proteolytic efficiency.    
   4. Recent optimizations of the TMT labeling protocol have resulted in a reduction of the required TMT reagent amounts and, consequently, the associated costs ([10.1074/mcp.TIR119.001385](https://doi.org/10.1074/mcp.TIR119.001385)). This is achieved by using lower reaction volumes, narrow peptide and TMT concentration ranges, and careful fine-tuning of peptide: TMT ratios. We typically aim for a TMT:peptide ratio of 4:1 to ensure a high degree of labeling efficiency.  
   5. Antibody-based enrichment specificity can be increased by subjecting the acid-eluted peptides to a second enrichment using IMAC beads. We describe the enrichment using commercial ready-to-go spin columns based on nitrilotriacetic acid (NTA), a chelating agent. However, products based on alternative chelating chemistries can also be used.    
   6. It is possible to modify the manufacturer's protocol and downscale the bead volume and peptide input, but specific ratios must be tested for the given sample.  
   7. Always ensure mixing of beads and peptide solution during the incubation while avoiding beads sticking to the tube walls and caps. Be careful not to remove any beads during the washing steps.  
   8. Typically, system backpressure stabilizes within 20 min and should be monitored over all injections to ensure optimal performance and column lifetime.   
   9. Different concatenation methods can be programmed and may include omitting fractions that contain the bulk of unreacted TMT reagent, as observed by an increase in the 260 nm UV reading  
   10. It is crucial to work with fresh binding and wash buffers, as the enrichment efficiency can be negatively affected.  
   11. Alternatively, Zr-IMAC HP can be use,d which will generate overlapping as well as unique/complementary peptide pools. 

	

* Critical points and considerations for the experimental setup

* Tips for data interpretation and limitations

## **References**

* Relevant and recent methodological publications, including FragPipe ([10.1038/s41467-023-39891-7](https://doi.org/10.1038/s41467-023-39891-7)), DIA-NN ([10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x)), msqrob2PTM ([10.1016/j.mcpro.2023.100708](https://doi.org/10.1016/j.mcpro.2023.100708)), and MSstatsPTM ([10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477))

* MSstatsPTM: MCP, Kohler et. al, 2023, [10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)
* TMT Mouse2plex Dataset: Maculins et. al 2021, eLife, [10.7554/eLife.62320](https://doi.org/10.7554/eLife.62320)  
* Zuniga et al, [10.1021/acs.jproteome.4c00668](https://doi.org/10.1021/acs.jproteome.4c00668)  
* Leutert et al, [10.15252/msb.20199021](https://doi.org/10.15252/msb.20199021)  
* Bortel et al, [10.1016/j.mcpro.2024.100754](https://doi.org/10.1016/j.mcpro.2024.100754)  
* Koenig et al, [10.1002/pmic.202100245](https://doi.org/10.1002/pmic.202100245)  
* Gajadhar et al, [10.1158/0008-5472.CAN-14-2309](https://doi.org/10.1158/0008-5472.CAN-14-2309)  
* Reddy et al, [10.1073/pnas.1521288113](https://doi.org/10.1073/pnas.1521288113)  
* Batth et al, [10.1074/mcp.TIR118.001270](https://doi.org/10.1074/mcp.TIR118.001270)  
* Zecha et al, [10.1074/mcp.TIR119.001385](https://doi.org/10.1074/mcp.TIR119.001385)


## **Supplementary Materials**

* Additional data tables
* Detailed protocols and software parameter settings
* Complete workflow scripts and parameter files
* Quality control reports and validation data

