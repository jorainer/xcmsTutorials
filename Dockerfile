FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL name="jorainer/xcms_tutorials" \
      url="https://github.com/jorainer/xcmsTutorials" \
      maintainer="johannes.rainer@eurac.edu" \
      description="Docker container to run xcms tutorials. This version bases on the Bioconductor release 3.19 docker image." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

## Install the xcmsTutorials package and additional required packages
RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask = FALSE, type = 'source'); BiocManager::install(c('RCurl', 'xcms'), ask = FALSE, dependencies = TRUE, type = 'source')"

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', dependencies = TRUE, type = 'source', build_vignettes = TRUE, repos = BiocManager::repositories())"
