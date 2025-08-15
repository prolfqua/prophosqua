---
title: "Supplementary Methods: Automated Sample Preparation and Wet Lab Protocols for Integrated PTM Analysis"

author:
  - Witold Wolski
  - Antje Dittmann
  - Laura Kunz
  - Christian Panse
  - Jonas Grossmann
date: 2025-08-13
running_head: "Supplement 1: Wet Lab Protocols"
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
| Witold Wolski      | FGCZ, SIB         | [witold.wolski@fgcz.uzh.ch](mailto:witold.wolski@fgcz.uzh.ch)   |
| Antje Dittmann     | FGCZ              | [antje.dittmann@fgcz.uzh.ch](mailto:antje.dittmann@fgcz.uzh.ch)  |
| Laura Kunz         | FGCZ              |                                                        |
| Christian Panse    | FGCZ, SIB         |                                                        |
| Jonas Grossmann    | FGCZ, SIB         | [jonas.grossmann@fgcz.uzh](mailto:jonas.grossmann@fgcz.uzh)      |

- FGCZ - Functional Genomics Center Zurich, Winterthurerstrasse 190, CH-8057 Zurich, 
- SIB - Swiss Institute of Bioinformatics - Amphipole, Quartier UNIL-Sorge, CH-1015 Lausanne

# 1. Introduction

Identification and quantification of phosphorylated peptides using quantitative proteomics enables a largely unbiased profiling of functionally relevant protein states. Among others, developments in sample preparation methods aimed at improving robustness, reproducibility, and throughput in bottom-up phosphoproteomics have enabled researchers to generate data sets with greater analytical depth and consistent quantification within larger study cohorts than ever before. Here, we describe the implementation and integration of largely automated processes, including protein digestion, high pH offline fractionation, phosphopeptide enrichment, and LC-MS/MS analyses of both input and phospho-enriched peptide fractions. The workflow leverages consistent sample processing and use of smaller input amounts (20-50 µg per sample) afforded by TMTpro-labeling and multiplexing. It offers various entry points to modify specific steps flexibly. These include but are not limited to the number of fractions, which can range from 12 to 36, the type of phosphopeptide enrichment matrix (e.g, Ti-IMAC, Zr-IMAC), and the potential use of complementary antibody-based affinity enrichments. 


# 2. Materials

## 2.1 Biological Samples

### 2.1.1 Sample Types Compatible with FGCZ Protocol

The wet lab protocols described in this chapter are optimized for high-throughput processing and are compatible with a wide range of biological sample types:

- **Fresh-frozen tissue samples** ($\le 5$ mg per sample): Optimal for preserving PTM states and minimizing degradation
- **Cultured cell pellets** ($\ge 1×10^{6}$ cells per sample): Standard format for most cell culture experiments  
- **Primary cells isolated from tissue**: Requires careful handling to maintain viability during isolation
- **Organoids and spheroid cultures**: Emerging model systems requiring specialized collection protocols

### 2.1.2 Sample Storage and Handling Requirements

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

### 2.1.3 Protocol Compatibility and Adaptations

The experimental methods presented here represent optimized protocols implemented at the Functional Genomics Center Zurich. While the demonstration uses the Maculins et al. dataset, the analytical approaches are broadly applicable to TMT-based phosphoproteomics studies regardless of the specific sample preparation method used, provided that:

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
3. MS2 spectra were recorded for the top 12 precursors with an intensity threshold of 10³, which were isolated at 0.7 Th and subjected to dynamic exclusion for 16 s. Normalized collision energy was set to 32% and spectra were recorded with a resolution of **30,000 with turboTM set to TMTpro**, ion target of **1×10⁵** and a maximal injection time of **Auto**

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
