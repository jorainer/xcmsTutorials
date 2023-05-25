# Exploring and analyzing LC-MS data with *Spectra* and *xcms*

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)


![xcms](man/figures/xcms.png)
![MsExperiment](man/figures/MsExperiment.png)
![Spectra](man/figures/Spectra-rainbow.png)
![MetaboCoreUtils](man/figures/MetaboCoreUtils.png)
![MsCoreUtils](man/figures/MsCoreUtils.png)

This workshop provides an overview of recent developments in Bioconductor to
work with mass spectrometry
([MsExperiment](https://github.com/RforMassSpectrometry/MsExperiment),
[Spectra](https://github.com/RforMassSpectrometry/Spectra)) and specifically
LC-MS data ([xcms](https://github.com/sneumann/xcms)) and walks through the
preprocessing of a small data set emphasizing on selection of data-dependent
settings for the individual pre-processing steps. The present workshop
represents an updated version of the workshop given at the Metabolomics Society
conference 2018 in Seattle (http://metabolomics2018.org).

Covered topics are:
- Data import and representation.
- Accessing, subsetting and visualizing data.
- Centroiding of profile MS data.
- Chromatographic peak detection.
- Empirically determine appropriate settings for the analyzed data set.
- Evaluation of identified peaks.
- Alignment (retention time correction).
- Correspondence (grouping of chromatographic peaks across samples).

The full R code of all examples along with comprehensive descriptions is
provided in the [xcms-preprocessing.Rmd](./xcms-preprocessing.Rmd) file. This
file can be opened with e.g. RStudio which allows execution of the individual R
commands (see section below for additionally required R packages). The R command
`rmarkdown::render("xcms-preprocessing.Rmd")` would generate the html file
[xcms-preprocessing.html](https://jorainer.github.io/metabolomics2018/xcms-preprocessing.html).


## Installation

The analysis in this document requires an R version >= 4.3.0 and recent versions
of the `MsExperiment`, `Spectra` and in particular the `xcms` (version >= 3.99.0
is needed) packages.

```r
install.packages("BiocManager")
BiocManager::install("jorainer/xcmsTutorials",
    dependencies = TRUE, ask = FALSE, update = TRUE)
```

Alternatively, the individual packages used in this tutorial can be installed
using the code below.

```r
#' Install the Bioconductor package manager
install.packages("BiocManager")

#' Install the required packages
BiocManager::install(c("msdata",
                       "Spectra",
                       "MsExperiment",
                       "MetaboCoreUtils",
                       "MsCoreUtils",
                       "png"))
BiocManager::install("sneumann/xcms")
```

The source code for this document along with the test data can be downloaded
from the github repository https://github.com/jorainer/xcmsTutorials
with the command (or alternatively downloading the zip archive directly from the
github page).

```
git clone https://github.com/jorainer/xcmsTutorials
```


## Additional documentation resources and tutorials

- Tutorial with additional examples and explanations for MS2-based
  annotations: https://jorainer.github.io/SpectraTutorials/
- Repository of the `MsCoreUtils` package:
  https://rformassspectrometry.github.io/MsCoreUtils/
- Repository of the `MetaboCoreUtils` package:
  https://rformassspectrometry.github.io/MetaboCoreUtils/
- Repository of the `Spectra` package:
  https://rformassspectrometry.github.io/Spectra/
- Repository of the `MetaboAnnotation` package:
  https://rformassspectrometry.github.io/MetaboAnnotation/
- Repository of the `CompoundDb` package:
  https://rformassspectrometry.github.io/CompoundDb/
