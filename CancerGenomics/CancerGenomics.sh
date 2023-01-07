#! /usr/bin/env bash
# render.sh

# R -e "rmarkdown::render('CancerGenomics.Rmd', clean=TRUE, output_format='html_notebook')"


R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); rmarkdown::render('CancerGenomics.Rmd', output_file= 'CancerGenomics.html')"
