# Untergasser Data
This file and the mentioned data are part of the test data set of RDML tools ([Untergasser Data](https://github.com/RDML-consortium/rdmlpython/tree/main/experiments/untergasser)).

The Untergasser data were created by Andreas Untergasser in collaboration with Maurice van der Hoff (AMC, Amsterdam, Netherlands) and Vladimir Benes (EMBL, Heidelberg, Germany). These experiments should provide the data to optimize the algorithms of the RDML-Tools and should serve as a large technical data collection to compare different machines, reaction mixes, amplicon sizes, primer concentrations, detection methods (SYBR-Green, LC-Green, Eva-Green and hydrolysis probes) and DNA standards. Andreas Untergasser did all pipetting.

This folder contains the pipetting schemes for each run and the corresponding RDML files with the raw fluorescence data. The RDML files can be viewed with [RDML-Edit](https://www.rdml-tools.com/edit.html). The raw data can be exported as spreadsheet table in RDES format on the experiment tab. [RDML-RunAnalysis](https://www.rdml-tools.com/runanalysis.html) can analyze each run using the AmplificationCurveAnalysis tab. For convenience, an edited version combining the results of each experiment is provided as results.tsv. Additionally, the data can be analyzed by running the [python test scripts](../../test/README.md), creating analysis results files starting with temp_ and providing an overview of the results as test output. As this output is compared with the expected results, it can be used to confirm the reliability of the calculations.

## The reaction mixes used in this study
|Mix Name   |Type        | Order Information                                                            |
|-----------|------------|------------------------------------------------------------------------------|
|Roche SYBR |SYBR Green I|LightCycler 480 SYBR Green I Master from Roche (04887352001)                  |
|SensiFast  |SYBR Green I|SensiFast SYBR No-ROX Mix from bioline (CSA-01190 or BIO-98005)               |
|LCGreen    |SYBR Green I|5x LCG PCR Master from the Wittwer Lab                                        |
|Roche Probe|Probe       |LightCycler 480 Probes Master from Roche (04707494001)                        |
|IDT Probe  |Probe       |PrimeTime Gene Expression 2X Master from Integrated DNA Technologies (1055770)|
|Roche PCR  |PCR         |                          |

## The qPCR machines used in this study
|Machine Name      |Wells|Location| Machine Details                          |
|------------------|-----|--------|------------------------------------------|
|Lightcycler480    |  384|  AMC NL|Lightcycler 480 II from Roche             |
|QuantStudio6Flex  |  384| EMBL DE|QuantStudio 6 Flex from applied biosystems|
|StepOne           |   96| EMBL DE|StepOne from applied biosystems           |
|Lightcycler480DKFZ|  384| DKFZ DE|Lightcycler 480 II from Roche             |

## The human genomic DNA used in this study
The human genomic DNA was ordered as Human Genomic DNA (11691112001, Roche Diagnostics GmbH, Mannheim, Germany). 

### A human genome has the size of:
|        |  GB  |  pg |
|:-------|-----:| ---:|
|male    | 6.27 | 6.41|
|female  | 6.37 | 6.51|
|both    | 6.32 | 6.46|
|haplo   | 3.16 | 3.23|

1pg == 0.978 * 10^9 bp

<pre>
 1.938ng ==  600 copies  
 4.845ng == 1500 copies  
 9.690ng == 3000 copies  
19.380ng == 6000 copies  
  
     1ng ==  309 copies  
     2ng ==  619 copies  
    20ng == 6192 copies  
</pre>
Values originate from [On the length, weight and GC content of the human genome. BMC Res Notes. 2019 Feb 27;12(1):106. doi: 10.1186/s13104-019-4137-z](https://doi.org/10.1186/s13104-019-4137-z)

<a id="experiment_volume"></a>

## Experiment Volumes and Machines (14.02.2024, 06.03.2024)
This experiment investigates the impact of different qPCR machines and different volumes on the measured Cq, TD0 and final Ncopy numbers. Care was taken that one large volume of qPCR mix was created using primer pair FSTL-1-*-259 and spread out in different volumes on each plate. The DNA dilutions and primer dilutions were created in AMC and were frozen and transported to EMBL after the experiment "AMC Amsterdam NL - 384 Wells" on the Lightcycler480 (14.02.2024). Then identical reactions were prepared in larger volumes and pipetted from the same tube to the plates for QuantStudio6Flex and StepOne in EMBL (06.03.2024), creating the experiments "EMBL Heidelberg DE - 384 Wells" and "EMBL Heidelberg DE - 96 Wells". The pipetting schemes are available as [volume_pipett_AMC.tsv](volume_pipett_AMC.tsv) and [volume_pipett_EMBL.tsv](volume_pipett_EMBL.tsv). The raw data are stored in [volume_machine.rdml](volume_machine.rdml), the analyzed results in [volume_results.tsv](volume_results.tsv) and the Cq values called by the machines in [volume_machine_cq.tsv](volume_machine_cq.tsv).

<a id="primer_test"></a>

## Experiment Primer Test (20.11.2024)
As genomic DNA should be used as a standard, the single copy gene [FSTL1](primers_FSTL1.fa) was selected as target. Then a big set of primers with probes was calculated with amplicon sizes of 50bp, 100bp, 200bp, 400bp and 800bp using [NCBI PrimerBLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/), see file [primers_evaluated.tsv](primers_evaluated.tsv). For each amplicon size a set of 3-4 primer pairs were selected. Care was taken that one hydrolysis probe would bind to several primer pairs [primers_used.tsv](primers_used.tsv). The default primer pair from AMC was included as FSTL-1-*-259. The primers were tested on QuantStudio6Flex in EMBL, using the LightCycler 480 SYBR Green I Master Mix from Roche (04887352001) and 4ng/1200 copies human genomic DNA (11691112001, Roche Diagnostics GmbH, Mannheim, Germany). The pipetting scheme is available as [amplicon_primer_mix_pipett_primer_test.tsv](amplicon_primer_mix_pipett_primer_test.tsv). The raw data are stored as experiment "Primer Test - EMBL" in [amplicon_primer_mix.rdml](amplicon_primer_mix.rdml) and the analyzed results in [amplicon_primer_mix_results.tsv](amplicon_primer_mix_results.tsv).

<a id="amc_experiments"></a>

## The AMC Experiments (09.12.2024 - 11.12.2024)
The following experiments were performed within three days in AMC using the Lightcycler480. Although the DNA dilutions are given in each pipetting scheme, only one DNA dilution series of human DNA was created once and was used for all the remaining experiments. During the three days in AMC, these DNA dilutions and the primers were stored at 4°C and not frozen. For the later experiments, DNAs and primers were stored at -20°C.

<a id="primer_concentation"></a>

## Experiment Primer Concentration (09.12.2024, 10.12.2024)
This experiment evaluates the influence of the primer concentration and the reaction mix on amplification of 4ng human DNA using 100nM, 250nM and 750nM primers in combination with Roche SYBR, SensiFast and LCGreen. The pipetting schemes are available as [amplicon_primer_mix_pipett_p1_primer_conc.tsv](amplicon_primer_mix_pipett_p1_primer_conc.tsv) for Roche SYBR (09.12.2024), [amplicon_primer_mix_pipett_p2_primer_conc.tsv](amplicon_primer_mix_pipett_p2_primer_conc.tsv) for LCGreen (09.12.2024) and [amplicon_primer_mix_pipett_p5_primer_conc.tsv](amplicon_primer_mix_pipett_p5_primer_conc.tsv) for SensiFast (10.12.2024). The raw data are stored as experiment "Primer Conc - AMC" in [amplicon_primer_mix.rdml](amplicon_primer_mix.rdml) and the analyzed results in [amplicon_primer_mix_results.tsv](amplicon_primer_mix_results.tsv).

<a id="dna_concentation"></a>

## Experiment DNA Concentration (09.12.2024 - 11.12.2024)
This experiment uses the DNA dilutions to evaluate PCR efficiency using 250nM primers in combination with Roche SYBR, SensiFast and LCGreen. The pipetting schemes are available as [amplicon_primer_mix_pipett_p3_dna_conc.tsv](amplicon_primer_mix_pipett_p3_dna_conc.tsv) for Roche SYBR (09.12.2024), [amplicon_primer_mix_pipett_p6_dna_conc.tsv](amplicon_primer_mix_pipett_p6_dna_conc.tsv) for SensiFast (10.12.2024) and [amplicon_primer_mix_pipett_p9_dna_conc.tsv](amplicon_primer_mix_pipett_p9_dna_conc.tsv) for LCGreen (11.12.2024). The raw data are stored as experiment "DNA Dilution - AMC" in [amplicon_primer_mix.rdml](amplicon_primer_mix.rdml) and the analyzed results in [amplicon_primer_mix_results.tsv](amplicon_primer_mix_results.tsv).

<a id="hydrolysis_probe"></a>

## Experiment Hydrolysis Probe (10.12.2024, 11.12.2024)
The first plate uses a hydrolysis probe to quantify DNA dilutions with 250nM primers and 150nM probe in combination with Roche Probe and IDT Probe master mixes (10.12.2024). The pipetting scheme is available as [probes_pipett_p4_probe_dna_conc.tsv](probes_pipett_p4_probe_dna_conc.tsv). The raw data are stored as run "P4 - Probe Mixes - DNA Dilution" in [probes.rdml](probes.rdml) and the analyzed results in [probes_results.tsv](probes_results.tsv). The second plate uses different primer and probe concentrations with 4ng human DNA and Roche Probe mix (11.12.2024). The pipetting scheme is available as [probes_pipett_p8_primer_probe_conc.tsv](probes_pipett_p8_primer_probe_conc.tsv). The raw data are stored as run "P8 - Primer conc - Probe conc" in [probes.rdml](probes.rdml) and the analyzed results in [probes_results.tsv](probes_results.tsv).

<a id="sybr_concentation"></a>

## Experiment SYBR Green I Concentration (09.12.2024 - 11.12.2024)
Initially, SYBR Green I was added to regular reactions with the Roche SYBR and SensiFast mix. The additional reactions were added to the plates of the experiment for DNA concentration. The pipetting schemes are available as [amplicon_primer_mix_pipett_p3_dna_conc.tsv](amplicon_primer_mix_pipett_p3_dna_conc.tsv) for Roche SYBR (09.12.2024) and [amplicon_primer_mix_pipett_p6_dna_conc.tsv](amplicon_primer_mix_pipett_p6_dna_conc.tsv) for SensiFast (10.12.2024). The raw data are stored as experiment "Add SYBR - AMC" in [sybr_conc.rdml](sybr_conc.rdml) and the analyzed results in [sybr_conc_results.tsv](sybr_conc_results.tsv). Then dilutions of SYBR Green I were added to Roche Probe mix using 100nM, 250nM and 750nM primers (11.12.2024). The pipetting scheme is available as [sybr_conc_pipett_p7_probe_mix.tsv](sybr_conc_pipett_p7_probe_mix.tsv). The raw data are stored as experiment "Probe Mix SYBR" in [sybr_conc.rdml](sybr_conc.rdml) and the analyzed results in [sybr_conc_results.tsv](sybr_conc_results.tsv).

<a id="embl_experiments"></a>

## The AMC Experiments (09.12.2024 - 11.12.2024)
The remaining experiments were performed in EMBL using the QuantStudio6Flex. The DNA dilution series of human DNA prepared in AMC was used for all the remaining experiments and was stored at -20°C between experiments.









