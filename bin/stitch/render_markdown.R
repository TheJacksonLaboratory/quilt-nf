#!/usr/bin/env Rscript

# argument information
# 1 - Alignment summary markdown file

# load arguments
args <- commandArgs(trailingOnly = TRUE)

markdown_rmd <- args[1]

rmarkdown::render(markdown_rmd)