###############################################################################
#                                                                             #
#  Intro to Reproducible use of genetic sequences in R and genetic epi        #
#                                                                             #
#   Script developed to highlight reproducible research methods               #
#     and to highlight use of tools realted to sequencing and                 #
#     use and interpretation of genetic sequence data in the R                #
#     programming language and environment.                                   #
#                                                                             #
#    Written in R version 4.0.2 using R Studio 1.3.1073                       #
#                                                                             #
#  Dependndies: tidyverse, bioconductor, seqinr, outbreaker2, rprojroot,      #
#               seqinr, ape, phangorn, DECIPHER, adephylo, adegenet           #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################

# installation of required packages (Commented) ----

#BiocManager::install(version = "3.11")
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
library(rprojroot)
library(tidyverse)
library(seqinr)
library(outbreaker2)
library(ape)
library(phangorn)
library(DECIPHER)
library(adephylo)
library(adegenet)
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

# load sequence data (in FASTA format)

fas <- "./R_example/sequences_from_session_1.fna"
dna <- readDNAStringSet(fas)
dna # the unaligned sequences

#
DNA <- AlignSeqs(dna)
DNA

# write alignments to file. 
writeXStringSet(DNA, file="./R_example/sequences_from_session_1_aligned.fna")

# reload from file (not necessary but worth checking that write worked)
seqs <- read.dna("./R_example/sequences_from_session_1_aligned.fna", format="fasta")
seqs

#create a phyDat object
virus_phyDat <- phyDat(seqs, type = "DNA", levels = NULL)

# test a vareity of nucleotide substitution models
mt <- modelTest(virus_phyDat, model = c("JC", "F81"))

print(mt)

# choose a model and create a distance matrix
dna_dist <- dist.ml(virus_phyDat, model="JC69")
dna_dist


virus_UPGMA <- upgma(dna_dist)

virus_NJ  <- NJ(dna_dist)

plot(virus_UPGMA, main="UPGMA")

plot(virus_NJ, main = "Neighbor Joining")
plot(virus_NJ, type = "unrooted", cex = 0.5)
plot(virus_NJ, type = "cladogram")
plot(virus_NJ, type = "fan")
plot(virus_NJ, type = "radial")
bullseye(virus_NJ, font = 2, cex = 0.5)

# outbreaker2 demo ----
# restore graphics settings 
par(par_old)

fake_outbreak
str(fake_outbreak)

plot(fake_outbreak$w, type = "h", xlim = c(0, 5), 
     lwd = 30, col = "red", lend = 2, 
     xlab = "Days after infection", 
     ylab = "p(new case)", 
     main = "Generation time distribution")

args(outbreaker)

dna <- fake_outbreak$dna
dates <- fake_outbreak$sample
ctd <- fake_outbreak$ctd
w <- fake_outbreak$w

# lets just pretend this is real DNA data

write.dna(dna, "outbreaker_unaligned.fasta", format = "fasta")

fas <- "outbreaker_unaligned.fasta"
dna_2 <- readDNAStringSet(fas)
dna_2 # the unaligned sequences

#
DNA <- AlignSeqs(dna_2)
DNA

writeXStringSet(DNA, file="./R_example/outbreaker_aligned.fna")

seqs <- read.dna("./R_example/outbreaker_aligned.fna", format="fasta")
seqs
outbreaker_phyDat <- phyDat(seqs, type = "DNA", levels = NULL)
mt <- modelTest(outbreaker_phyDat, model = c("JC", "F81"))
print(mt)

# choose a model and create a distance matrix
dna_dist <- dist.ml(outbreaker_phyDat, model="F81")
dna_dist


# contruct tree from matix
outbreaker_UPGMA<- upgma(dna_dist)

plot(outbreaker_UPGMA, main="UPGMA", "unrooted")
plot(outbreaker_UPGMA, main="UPGMA")

bullseye(outbreaker_UPGMA, font = 2, tip.color = any2col(dates, col.pal=seasun)$col)

# contruct the outbreaker model data set

data <- outbreaker_data(dna = dna, dates = dates, ctd = ctd, w_dens = w)

## we set the seed to ensure results won't change
set.seed(1)

# run the analysis with defaults 
res <- outbreaker(data = data)

class(res)

dim(res)

names(res)

plot(res)
plot(res, "mu")
plot(res, "mu", "density", burnin = 2000)
plot(res, type = "alpha", burnin = 2000)
plot(res, type = "t_inf", burnin = 2000)
plot(res, type = "kappa", burnin = 2000)
plot(res, type = "network", burnin = 2000, min_support = 0.01)






