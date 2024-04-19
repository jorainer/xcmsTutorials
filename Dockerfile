FROM bioconductor/bioconductor_docker:devel

LABEL name="jorainer/xcms_tutorials" \
      url="https://github.com/jorainer/xcmsTutorials" \
      maintainer="johannes.rainer@eurac.edu" \
      description="Docker container to run xcms tutorials." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

## Install the xcmsTutorials package and additional required packages
RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask = FALSE); BiocManager::install('xcms', ask = FALSE, type = 'source')"

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', dependencies = TRUE, type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"
