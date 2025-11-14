# Untergasser Data
This file and the mentioned data are publically avaiable and part of the test data set of RDML tools ([Untergasser Data](https://github.com/RDML-consortium/rdmlpython/tree/main/experiments/untergasser)).

The Untergasser data were created by Andreas Untergasser in collaboration with Maurice van der Hoff (AMC, Amsterdam, Netherlands) and Vladimir Benes (EMBL, Heidelberg, Germany). These experiments should provide the data to optimize the algorithms of the RDML-Tools and should serve as a large technical data collection to compare different machines, reaction mixes, amplicon sizes, primer concentrations, detection methods (SYBR-Green, LC-Green, Eva-Green and hydrolysis probes) and DNA standards. All pipetting was done by Andreas Untergasser.

The data can be analyzed by running the [python test scripts](../../test/README.md).

### The reaction mixes used in this study
|Mix Name   |Type        | Order Information                                              |
|-----------|------------|----------------------------------------------------------------|
|Roche SYBR |SYBR Green I|LightCycler 480 SYBR Green I Master Mix from Roche (04887352001)|
|SensiFast  |SYBR Green I|SensiFast SYBR No-ROX Mix from bioline (CSA-01190)              |
|LCGreen    |SYBR Green I|5x LCG PCR Master from the Wittwer Lab                          |
|Roche Probe|Probe       |LightCycler 480 Probes Master from Roche (04707494001)          |
|IDT Probe  |Probe       |                          |
|Roche PCR  |PCR         |                          |


### The qPCR machines used in this study
|Machine Name      |Wells|Location| Machine Details                          |
|------------------|-----|--------|------------------------------------------|
|Lightcycler480    |  384|  AMC NL|Lightcycler 480 II from Roche             |
|QuantStudio6Flex  |  384| EMBL DE|QuantStudio 6 Flex from applied biosystems|
|StepOne           |   96| EMBL DE|StepOne from applied biosystems           |
|Lightcycler480DKFZ|  384| DKFZ DE|Lightcycler 480 II from Roche             |

### The human genomic DNA used in this study
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
This experiment investigates the impact of different qPCR machines and different volumes on the measured Cq, TD0 and final Ncopy numbers. Care was taken that one large volume of qPCR mix was created using primer pair FSTL-1-*-259 and spread out in different volumes on each plate. The DNA dilutions and primer dilutions were created in AMC and were frozen and transported to EMBL after the experiment "AMC Amsterdam NL - 384 Wells" on the Lightcycler480 (14.02.2024). Then identical reactions were prepared in larger volumes and pipetted from the same tube to the plates for QuantStudio6Flex and StepOne in EMBL (06.03.2024), creating the experiments "EMBL Heidelberg DE - 384 Wells" and "EMBL Heidelberg DE - 96 Wells". The pipetting schemes are availabe as [volume_pipett_AMC.tsv](volume_pipett_AMC.tsv) and [volume_pipett_EMBL.tsv](volume_pipett_EMBL.tsv). The raw data are available in [volume_machine.rdml](volume_machine.rdml), the analyzed results in [volume_results.tsv](volume_results.tsv) and the Cq values called by the machines in [volume_machine_cq.tsv](volume_machine_cq.tsv).

<a id="primer_test"></a>

## Experiment Primer Test (20.11.2024)
As genomic DNA should be used as a standard, the single copy gene [FSTL1](primers_FSTL1.fa) was selected as target. Then a big set of primers with probes was calculated with amplicon sizes of 50bp, 100bp, 200bp, 400bp and 800bp using [NCBI PrimerBLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/), see file [primers_evaluated.tsv](primers_evaluated.tsv). For each amplicon size a set of 3-4 primer pairs were selected. Care was taken that one hydrolysis probe would bind to several primer pairs [primers_used.tsv](primers_used.tsv). The default from AMC primer was included as FSTL-1-*-259. The primers were tested on QuantStudio6Flex in EMBL, using the LightCycler 480 SYBR Green I Master Mix from Roche (04887352001) and 4ng/1200 copies human genomic DNA (11691112001, Roche Diagnostics GmbH, Mannheim, Germany). The pipetting scheme is availabe as [amplicon_primer_mix_pipett_primer_test.csv](amplicon_primer_mix_pipett_primer_test.csv). The raw data are included in [amplicon_primer_mix.rdml](amplicon_primer_mix.rdml) and the analyzed results in [amplicon_primer_mix_results.tsv](amplicon_primer_mix_results.tsv).







In the file primers.csv all designed primers and probe combinations are given. In top you find the selected set.


