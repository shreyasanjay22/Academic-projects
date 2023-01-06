#!/usr/bin/env bash
# renderDe.sh

R -e "rmarkdown::render('methodsResults.Rmd', output_format='all')"
