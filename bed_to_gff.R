#!/usr/bin/Rscript

###The following is a nice bit of code to prevent installation of already present packages. Code from here: https://gist.github.com/smithdanielle/9913897
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages<-c("rtracklayer")
check.packages(packages)


args<-commandArgs(TRUE)

input <- args[1]
output <- args[2]
species <- args[3]

## import the bed file
bed.ranges <- import.bed(input)

## export as a gff3 file
export.gff3(bed.ranges, source = species, output)
