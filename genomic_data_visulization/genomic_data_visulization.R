###############################################################################
#                                                                             #
#  Intro to visualization of data and genomic data in R                       #
#                                                                             #
#   Script developed to highlight reproducible research methods               #
#     for visulaization creation in R                                         #
#                                                                             #
#    Written in R version 4.3.2 using R Studio 1.3.1073                       #
#                                                                             #
#  Dependencies: TBD                                                          #
#                                                                             #
#                                                                             #
#                   author: Joshua Yukich jyukich@tulane.edu                  #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################

# installation of required packages (Commented) ----
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.18")
#BiocManager::install('ggbio')
#install.packages('rprojroot')
#install.packages('seqinr')
#install.packages('outbreaker2')
#install.packages('rprojroot')
#install.packages('tidyverse')
#install.packages('ape')
#install.packages('phangorn')
#BiocManager::install('DECIPHER')
#install.packages('adephylo')
#install.packages('adegenet')

# loading required libraries ----
library(ggbio)
library(rprojroot)
library(GenomicRanges)

# session info -----

run_date <- date()
run_date

sessionInfo()

# setting working directory and other user preferences -----

if (interactive()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

base_path <- find_root('base_environment.R')

setwd(base_path)
source('base_environment.R')

# preserve graphics settings ----
par_old <- par()


# now for a start ----

# getting citations for tools in R

citation("ggbio")
citation()


p_ideo_1 <- Ideogram(genome = "hg19")
p_ideo_1

p_ideo_Y <- Ideogram(genome = "hg19", subchr = "chrY")
p_ideo_Y
p_ideo_X <- Ideogram(genome = "hg19", subchr = "chrX")
p_ideo_X

help("Ideogram")

p_ideo_1 + xlim(GRanges("chr1", IRanges(1, 100000)))

data("CRC", package = "biovizBase")
head(hg19sub)
autoplot(hg19sub, layout = "circle", fill = "gray70")

gr.crc1 <- crc.gr[values(crc.gr)$individual == "CRC-1"]
gr.crc2 <- crc.gr[values(crc.gr)$individual == "CRC-2"]
gr.crc3 <- crc.gr[values(crc.gr)$individual == "CRC-3"]
gr.crc4 <- crc.gr[values(crc.gr)$individual == "CRC-4"]

ggbio() + #circle(mut.gr, geom = "rect", color = "steelblue") +
  circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3) + 
  #circle(gr.crc1, geom = "point", aes(y = score, size = tumreads),
  #         color = "red", grid = TRUE, radius = 30) + scale_size(range = c(1, 2.5)) + 
  circle(gr.crc1, geom = "link", linked.to = "to.gr", aes(color = rearrangements),
         radius = 35) + 
  circle(gr.crc2, geom = "link", linked.to = "to.gr", aes(color = rearrangements),
         radius = 35) +
  circle(gr.crc3, geom = "link", linked.to = "to.gr", aes(color = rearrangements),
       radius = 35) + 
  circle(gr.crc4, geom = "link", linked.to = "to.gr", aes(color = rearrangements),
       radius = 35)


