# rawDIAtect
This script process the [DIA-NN](https://github.com/vdemichev/DiaNN) output of any DIA LC-MS run analyzed with the [CARD](https://card.mcmaster.ca/) AMR databases and a species background from [UniProt](https://www.uniprot.org/) (modified databases can be found in ``/Database``). It filters for unique AMR related peptides and only takes a selected amount of them. Furthermore it calculates the TOP-3 precursor quantity and potiental Drug Class spectra. The results are then collected into ``pdf reports``. It can analyze multiple DIA-NN outputs at once and reports every sample individual.  

## Installation

**[Only tested on Windwos 10, RStudio v1.2.1335, R v3.6.1 (2019-07-05)]**

Copy the repo folder ``/rawDIAtect`` including the script ``rawDIAtect.R`` (*you can make a shortcut of the script on the desktop*) anywhere on your computer an run it with RStudio. Best read setup before first run.

Run script with RStudio
> Run Source or Ctrl+Shift+S
 
 or the inside the bash by 
```
sudo Rscript path/to/your/file/rawDIAtect.R
```
### Databases for DIA-NN
In the folder ``/Database`` are two zip files, those are the two protein CARD fastas and some proteomic background fastas. If your want to run DIA-NN with your MS data you need to unzip them to use them.

### Test run
The script will install all dependecies upon the first run.

In ``/test_files`` is a DIA-NN test output. Run correctly it should result in 3 ``pdf reports``. One is ``empty`` the other two contain either ``VIM-1`` or ``CMY-2``.

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
![options.PNG](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/options.PNG)

### Setup script and options
Before the first run change the path to the ``aro_index_DIA.tsv``  at the line **A**. Now the script is ready for your DIA-NN files. 

First change the path at **B** where your DIA-NN output files are stored (*best create a new Folder only containing your files you want to run, on default the output will be saved in the input folder*). If you want to save your report elsewhere you can change it at line **C**. You also can give your experiment a name at line **F**.

You can change the peptides detection threshold at line **D**. Default is **3** (*best use between 2-4*). If using the Full CARD the output can be messy, by using the whitelist option at line **E** (1 = on) only whitelisted hits are shown (``KPC,OXA,SHV,TEM,VAN/AB,MCR,CMY,AAC,CTX-M,NDM,VIM``), turn this off if  you used the Whitelisted CARD.

### Run the script
When everthing is setup you can run the script in RStudio by 
> Run Source or Ctrl+Shift+S

or the inside the bash by 
```
sudo Rscript path/to/your/file/rawDIAtect.R
```
## Output
![Main-Page](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/86-09_main_page.png )
**A**=AMR gene families, **B**=expected Drug classs spectra, **C**= On top possibly protein isoform and bottom the count of unique peptides.

![Details-1](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/86-09_overview_page.png)
**A**=AMR gene families, **B**= All possible protein isoforms in relation.

![Details-2](https://github.com/CptChiler/rawDIAtect/blob/master/readme_png/86-09_PQ_top3.png)
**A**=TOP-3 preqursor quantity per AMR gene family, **B**= Drug spectra proportion of all peptides found.

## Release History
* v1 (18.05.2020)

## Meta
Christian Blumenscheit – [@ChrisMiBiFlower](https://twitter.com/chrismibiflower)
Distributed under the GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007. See ``LICENSE`` for more information.
[https://github.com/CptChiler](https://github.com/CptChiler)

