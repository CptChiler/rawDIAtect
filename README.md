# [rawDIAtect](https://pubs.acs.org/doi/10.1021/acs.analchem.1c00594)
This script process the [DIA-NN](https://github.com/vdemichev/DiaNN) output of any DIA LC-MS run analyzed with the [CARD](https://card.mcmaster.ca/) AMR databases and a species background from [UniProt](https://www.uniprot.org/) (modified databases can be found in ``/database``). It filters for unique AMR related peptides and only takes a selected amount of them. Furthermore it calculates the TOP-3 precursor quantity and potiental Drug Class spectra. The results are then collected into ``pdf reports``. It can analyze multiple DIA-NN outputs at once and reports every sample individual.

Paper : https://pubs.acs.org/doi/10.1021/acs.analchem.1c00594

### Dependecies
 - readr
 - tidyverse
 - profvis
 - tidyr
 - ggplot2
 - dplyr
 - cowplot
 - plyr
 - devtools
 - gameofthrones
 - tm
 - crayon
 - progress
 - Hmisc
 
## Installation

**[Tested on Windwos 10, RStudio v1.2.1335, R v4.0.3 & v3.6.1(2020-11-26)**

Run with RStudio
```
install.packages("devtools")
library("devtools")
devtools::install_github("https://github.com/CptChiler/rawDIAtect")
```
or from release https://github.com/CptChiler/rawDIAtect/releases

You may need to restart RStudio to see the new package.

## Databases for DIA-NN and settings
In the folder ``/database`` are two zip files, those are the two protein CARD fastas and some proteomic background fastas. If your want to run DIA-NN with your MS data you need to unzip them to use them.

#### AMR Databases
>1. “Full CARD” (196.072 entries) -> Modified [CARD](https://card.mcmaster.ca/)
>2. “Whitelisted CARD” (22 entries) -> PCR (CMY-2, CTX-M-9/15, KPC-2/3, NDM-1, OXA-1/48, SHV-1/12, TEM-1/52, VIM-1, AAC(6)-lb-cr, VanA-R/S and VanB-R/S)

#### DIA-NN setup
![DIA-NN setup example](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/Folie1.png )

**1**=Add your raw DIA file, **2**=Add a species background and a CARD target fasta, **3**= Check the boxes for Fasta digestion and Deep learning spectra prediction. **4**= Uncheck batch process and check Unrelated runs. Maunual set Mass accuracy and MS1 accuracy to 10 and 20 ppm.

## Usage example

Input needs to be a folder with "report".tsv files from DIA-NN and runs should have unique names. You can put multiple reports in one folder.

<ins>Minimal input:<ins>

```
rawDIAtect(path_in = "/Quants")
```

<ins>Options:<ins>

iso_diff = 2 (0 = Off), _Allows a max difference to other protein isoforms_ <br/>
pep_filter = 3 (<1 = Off), _Sets threshold for unique peptide per hit_ <br/>
path_out = (Default = path_in), _Path were report and split data will be stored_ <br/>
Exp_name = "Induced", _Experiment name (all files will get this name)_ <br/>

<ins>Expert input:<ins>

```
rawDIAtect(path_in = "/Quants", pep_filter = 3 ,iso_diff = 2, path_out = "/AMRs", Exp_name = "Induced")
```

## Output
![Main-Page](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/86-09_main_page.png )
**A**=AMR gene families, **B**=expected Drug classs spectra, **C**= On top possibly protein isoform and bottom the count of unique peptides.

![Details-1](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/86-09_overview_page.png)
**A**=AMR gene families, **B**= All possible protein isoforms in relation.

![Details-2](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/86-09_PQ_top3.png)
**A**=TOP-3 preqursor quantity per AMR gene family, **B**= Drug spectra proportion of all peptides found.

## Release History
* v0.4 (20.01.2021)

## Meta
Christian Blumenscheit – [@ChrisMiBiFlower](https://twitter.com/chrismibiflower)
Distributed under the GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007. See ``LICENSE`` for more information.
[https://github.com/CptChiler](https://github.com/CptChiler)

