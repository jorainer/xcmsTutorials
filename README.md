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
settings for the individual pre-processing steps.

Covered topics are:

- Data import and representation.

- Accessing, subsetting and visualizing data.

- Centroiding of profile mode MS data.

- Chromatographic peak detection.

- Empirically determine appropriate settings for the analyzed data set.

- Evaluation of identified peaks.

- Alignment (retention time correction).

- Correspondence (grouping of chromatographic peaks across samples).

The full R code of all examples along with comprehensive descriptions is
provided in the [xcms-preprocessing.Rmd](./vignettes/xcms-preprocessing.Rmd)
file. This file can be opened with e.g. RStudio which allows execution of the
individual R commands (see section below for additionally required R
packages). The R command `rmarkdown::render("xcms-preprocessing.Rmd")` would
generate the html file
[xcms-preprocessing.html](https://jorainer.github.io/xcmsTutorials/xcms-preprocessing.html).


## Installation

For on-line code evaluation, the workshop can also be run using a self-contained
docker image with all R packages and a server version of RStudio (Posit)
pre-installed:

- If you don't already have, install [docker](https://www.docker.com/). Find
  installation information [here](https://docs.docker.com/desktop/).
- Get the [docker image](https://hub.docker.com/r/jorainer/xcms_tutorials) of
  this tutorial with `docker pull jorainer/xcms_tutorials:latest`.
- Start docker using
```
  docker run \
      -e PASSWORD=bioc \
      -p 8787:8787 \
      jorainer/xcms_tutorials:latest
  ```

- Enter `http://localhost:8787` in a web browser and log in with username
  `rstudio` and password `bioc`.
- In the RStudio server version: open any of the R-markdown (*.Rmd*) files in
  the *vignettes* folder and evaluate the R code blocks.


For manual installation, an R version >= 4.3.0 is required as well as recent
versions of the packages `MsExperiment`, `Spectra` and in particular the `xcms`
(version >= 3.99.0 is needed). These can be installed using the code below:

```r
install.packages("BiocManager")
BiocManager::install("jorainer/xcmsTutorials",
    dependencies = TRUE, ask = FALSE, update = TRUE)
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
