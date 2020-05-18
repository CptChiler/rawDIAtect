# rawDIAtect

> Short blurb about what your product does.

One to two paragraph statement about your product and what it does.

## Installation

Copy the file ``rawDIAtect.R`` anywhere on your computer an run it with RStudio 
> Run Source or Ctrl+Shift+S
 
 or the inside the bash by 
```
sudo Rscript path/to/your/file/rawDIAtect.R
```
The script will install all dependecies upon the first run.

Dependecies:
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

![Options](https://imgur.com/8STvKgY)

### Setup script and options
Before the first run setup up the path to the ``aro_index_DIA.tsv``  at the line **A**. Now the script is ready for your DIA-NN files. 

First change the path and **B** where your DIA-NN output files are stored (*best create a new Folder only containing your files you want to run, on default the ourput will be saved in the input folder*). If you want to save your report elsewhere you can change it at line **C**. You can give your experiment a name at line **F**.

You can change the peptides detection threshold at line **D**. Default is **3** (*best use between 2-4*). If using the Dull CARD the output can be messy, by using the whitelist option at line **E** (1 = one) only whitelisted hits are shown (``KPC,OXA,SHV,TEM,VAN/AB,MCR,CMY,AAC,CTX-M,NDM,VIM``), turn this off if used the Whitelisted CARD.

### Run the script
When everthing is setup you can run the script in RStudio by 
> Run Source or Ctrl+Shift+S

or the inside the bash by 
```
sudo Rscript path/to/your/file/rawDIAtect.R
```
## Output

![Main-Page](https://imgur.com/8TnjcVx)
**A**=AMR gene families, **B**=expected Drug classs spectra, **C**= On top possibly protein isoform and bottom the count of unique peptides.

![Details-1](https://imgur.com/KSg2NP8)
**A**=AMR gene families, **B**= All possible protein isoforms in relation.

![Details-2](https://imgur.com/8v6JMsO)
**A**=TOP-3 preqursor quantity per AMR gene family, **B**= Drug spectra proportion of all peptides found.


## Release History

* v1

## Meta

Your Name – [@YourTwitter](https://twitter.com/dbader_org) – YourEmail@example.com

Distributed under the XYZ license. See ``LICENSE`` for more information.

[https://github.com/yourname/github-link](https://github.com/dbader/)
