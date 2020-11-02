# rawDIAtect
This script process the [DIA-NN](https://github.com/vdemichev/DiaNN) output of any DIA LC-MS run analyzed with the [CARD](https://card.mcmaster.ca/) AMR databases and a species background from [UniProt](https://www.uniprot.org/) (modified databases can be found in ``/database``). It filters for unique AMR related peptides and only takes a selected amount of them. Furthermore it calculates the TOP-3 precursor quantity and potiental Drug Class spectra. The results are then collected into ``pdf reports``. It can analyze multiple DIA-NN outputs at once and reports every sample individual.  

## Installation

**[Only tested on Windwos 10, RStudio v1.2.1335, R v3.6.1 (2019-07-05)]**

Run script with RStudio
```
devtools::install_github("https://github.com/CptChiler/rawDIAtect")
```
You may need to restart RStudio to see the new package

### Databases for DIA-NN
In the folder ``/database`` are two zip files, those are the two protein CARD fastas and some proteomic background fastas. If your want to run DIA-NN with your MS data you need to unzip them to use them.

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

## Usage example

Input needs to be a folder with .csv files from DIA-NN!

Minimal input:

```
rawDIAtect(path_in = "/Quants")
```

Options:

iso_diff = 2 (0 = Off) 
pep_filter = 3 (<1 = Off)
path_out = (Default = path_in)

Expert input:

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
* v0.3 (11.02.2020)

## Meta
Christian Blumenscheit – [@ChrisMiBiFlower](https://twitter.com/chrismibiflower)
Distributed under the GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007. See ``LICENSE`` for more information.
[https://github.com/CptChiler](https://github.com/CptChiler)

