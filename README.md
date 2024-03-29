# Leavitt et al. (2023) supplemental

This repository holds data and supplementary source code for *W.D. Leavitt, S.H. Kopf, Y. Weber, B. Chiu, J.M. McFarlin, F.J. Elling, S. Hoeft-McCann, A. Pearson. Controls on the hydrogen isotope composition of tetraether lipids in an autotrophic ammonia-oxidizing marine archaeon. Geochimica et Cosmochimica Acta, Volume 352 (194-210), 2023*.

## What can I do with this code? <a href="https://creativecommons.org/licenses/by/4.0/"><img src="https://mirrors.creativecommons.org/presskit/buttons/88x31/png/by.png" align = "right" width = "100"/></a>

We hope that this code, or any part of it, might prove useful to other members of the scientific community interested in the subject matter. This repository is released under a [Creative Commons BY (CC-BY)](https://creativecommons.org/licenses/by/4.0/) license, which means all code can be shared and adapted for any purpose as long as appropriate credit is given. See [Attribution section](https://creativecommons.org/licenses/by/4.0/) for details. Most of the code is provided in [R Markdown](http://rmarkdown.rstudio.com/) (see below).

## What is R Markdown?

[R Markdown](http://rmarkdown.rstudio.com/) is a so-called "literate programming" format that enables easy creation of dynamic documents with the open-source [R](http://www.r-project.org/) language. HTML and PDF reports can be generated from R Markdown files using [knitr](http://yihui.name/knitr/) and [pandoc](http://johnmacfarlane.net/pandoc/), which can be installed automatically with [RStudio](http://www.rstudio.com/), and are fully integrated into this cross-platform IDE. All software used for these RMarkdown reports (R, RStudio, etc.) is freely available and completely open-source. 

## How can I run this code?

The quickest and easiest way is to use RStudio.

 1. Download and install [R](http://cran.rstudio.com/) for your operating system
 1. Download and install [RStudio](http://www.rstudio.com/products/rstudio/download/) for your operating system
 1. Download a [zip file of this repository](https://github.com/KopfLab/2022_leavitt_et_al/archive/master.zip) and unpack it in an easy to find directory on your computer
 1. Navigate to the directory and double-click the `project.Rproj` file to start RStudio and load this project.
 1. Install the required libraries by running the following command in the Console in RStudio: `install.packages(c("rlang", "tidyverse", "latex2exp", "cowplot", "readxl", "openxlsx"))` or by installing them manually in RStudio's Packages manager.
 1. Open the `.Rmd` notebooks in the file browser (e.g. `01_calculations.Rmd`)
 1. To generate an HTML report ("knit HTML"), select File --> Knit from the menu. The HTML report will be displayed upon successful completion and is saved as a standalone file in the same directory. All generated data figures are saved as PDF and PNG in the `figures` sub-directory. All generated data tables are saved as XLSX in the `tables` sub-directory.
 
## Troubleshooting notes

The R Markdown files in this repository make use of various R modules for data processing, plotting and modelling. All of these should be installed automatically when the first R Markdown file is knitted (if the knitting fails because of a missing package, please install it manually, an error will indicate which package could not be installed). 
