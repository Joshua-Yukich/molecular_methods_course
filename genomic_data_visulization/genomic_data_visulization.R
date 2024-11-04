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

# example installation of required packages (Commented) ----


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.18")
#BiocManager::install('ggbio')
#install.packages('rprojroot')

# loading required libraries ----
library(rprojroot)
library(GenomicRanges)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(treeio)
library(tidytree)
library(ggtreeExtra)
library(phyloseq)
library(adephylo)
library(tidyverse)
library(ggtree)
library(dichromat)
library(ggpubr)
library(pheatmap) ## for heatmap generation
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap
library(airway)
library(ggbio)
library(TDbook)
library(scales)


# session info -----

run_date <- date()
run_date

sessionInfo()

# adapted from: 

#https://bioinformatics.ccr.cancer.gov/docs/data-visualization-with-r/#welcome-to-the-data-visualization-with-r-series
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625
# and https://yulab-smu.top/treedata-book/index.html
# https://r-graph-gallery.com/297-circular-barplot-with-groups.html


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

# most basic genomic data visualization is the phylogenetic tree.

# load sequence data (in FASTA format)

fas <- "./R_example/sequences_from_session_1.fna"
dna <- readDNAStringSet(fas)
dna # the unaligned sequences

#
DNA <- AlignSeqs(dna)
DNA

writeXStringSet(DNA, file="./R_example/sequences_from_session_1_aligned.fna")

# reload from file (not necessary but worth checking that write worked)
seqs <- read.dna("./R_example/sequences_from_session_1_aligned.fna", format="fasta")
seqs

class(seqs)
#create a phyDat object
virus_phyDat <- phyDat(seqs, type = "DNA", levels = NULL)

tree <- NJ(dist.ml(virus_phyDat, "JC"))
# test a vareity of nucleotide substitution models
mt <- modelTest(virus_phyDat, tree = tree, model = c("JC", "F81"))

print(mt)

# choose a model and create a distance matrix
dna_dist <- dist.ml(virus_phyDat, model="JC69")
dna_dist

virus_UPGMA <- upgma(dna_dist)

virus_NJ <- NJ(dna_dist)

plot(virus_UPGMA, main="UPGMA")

plot(virus_NJ, main = "Neighbor Joining")
plot(virus_NJ, type = "unrooted", cex = 0.5)
plot(virus_NJ, type = "cladogram")
plot(virus_NJ, type = "fan")
plot(virus_NJ, type = "radial")
bullseye(virus_NJ, font = 2, cex = 0.5)


class(virus_UPGMA)
virus_UPGMA

virus_UPGMA_tibble <- as_tibble(virus_UPGMA)

virus_UPGMA_tibble$orig <- c("Human", "Human", "Bat", "Bat","Bat","Bat", 
                          "Human", "Human", "Environment", "Bat", "Bat",
                          "Human", "Camel", "Camel", "Human", rep(NA, 14))

virus_UPGMA_full <- as.treedata(virus_UPGMA_tibble) 

virus_UPGMA_full %>% child(16) 
virus_UPGMA_full %>% MRCA(2,3) 

ggtree(virus_UPGMA)
ggtree(virus_UPGMA, layout="roundrect")
ggtree(virus_UPGMA, layout="slanted")
ggtree(virus_UPGMA, layout="ellipse")
ggtree(virus_UPGMA, layout="circular")
ggtree(virus_UPGMA, layout="fan", open.angle=120)
ggtree(virus_UPGMA, layout="equal_angle")
ggtree(virus_UPGMA, layout="daylight")
ggtree(virus_UPGMA, branch.length='none')
ggtree(virus_UPGMA, layout="ellipse", branch.length="none")
ggtree(virus_UPGMA, branch.length='none', layout='circular')
ggtree(virus_UPGMA, layout="daylight", branch.length = 'none')


ggtree(virus_UPGMA) + scale_x_reverse()
ggtree(virus_UPGMA) + coord_flip()
ggtree(virus_UPGMA) + layout_dendrogram()
ggplotify::as.ggplot(ggtree(virus_UPGMA), angle=-30, scale=.7)
ggtree(virus_UPGMA, layout='slanted') + coord_flip()
ggtree(virus_UPGMA, layout='slanted', branch.length='none') + layout_dendrogram()
ggtree(virus_UPGMA, layout='circular') +ylim(0,100)
ggtree(virus_UPGMA) + layout_inward_circular()
ggtree(virus_UPGMA) + layout_inward_circular(xlim=0.85)
ggtree(virus_UPGMA, mrsd = "2024-01-01") + theme_tree2()


beast_file <- system.file("examples/MCC_FluA_H3.tree", 
                          package="ggtree")
beast_tree <- read.beast(beast_file)
ggtree(beast_tree, mrsd = "2024-01-01") + theme_tree2()


ggtree(virus_UPGMA) + geom_treescale()

ggtree(virus_UPGMA) + 
  geom_point(aes(shape=isTip, color=isTip), size=3)

p <- ggtree(virus_UPGMA_full) + 
  geom_nodepoint(color="#b5e521", alpha=1/4, size=10) 
p + geom_tippoint(color="#FDAC4F", shape=8, size=3)

p + geom_tiplab()

ggtree(virus_UPGMA_full, layout = "circular") +
  geom_tiplab(aes(label = orig, angle=angle), color='blue')

ggtree(virus_UPGMA_full, aes(color = orig)) + scale_color_discrete()

ggtree(virus_UPGMA_full, layout='circular', ladderize = TRUE, size=2.8) + 
  geom_tree(aes(color=orig), size=2) +  
  #discreten(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(aes(color=orig)) 



ggtree(virus_UPGMA_full) +geom_label(aes(x=branch, label=orig), fill='lightgreen') +
  geom_cladelab(node=19, label="another clade", align=TRUE, 
                               geom='label', fill='lightblue')

data("GlobalPatterns")
GP <- GlobalPatterns
GP@sam_data
GP <- prune_taxa(taxa_sums(GP) > 600, GP)
sample_data(GP)$human <- get_variable(GP, "SampleType") %in%
  c("Feces", "Skin")
mergedGP <- merge_samples(GP, "SampleType")
mergedGP <- rarefy_even_depth(mergedGP,rngseed=394582)
mergedGP <- tax_glom(mergedGP,"Order")

melt_simple <- psmelt(mergedGP) %>%
  filter(Abundance < 120) %>%
  select(OTU, val=Abundance)

p <- ggtree(mergedGP, layout="fan", open.angle=10) + 
  geom_tippoint(mapping=aes(color=Phylum), 
                size=1.5,
                show.legend=FALSE)
p <- rotate_tree(p, -90)

p <- p +
  geom_fruit(
    data=melt_simple,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=val,
      group=label,
      fill=Phylum,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 

p <- p +
  scale_fill_discrete(
    name="Phyla",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  )


# using dichromat 

q <- p + scale_fill_manual(values = dichromat(scales::hue_pal()(24)))

ggarrange(p, q)

p + scale_fill_viridis_d(option = "rocket")




# Create dataset
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
  value=sample( seq(10,100), 60, replace=T)
)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Make the plot
r <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

r



r + scale_fill_manual(values = dichromat(scales::hue_pal()(4), type = "deutan"))
r + scale_fill_manual(values = dichromat(scales::hue_pal()(4), type = "protan"))
r + scale_fill_manual(values = dichromat(scales::hue_pal()(4), type = "tritan"))

str(r)
show_col(dichromat(unique(ggplot_build(r)$data[[1]]$fill)))

show_col(ggplot_build(r)$data[[1]]$fill)
show_col(dichromat(ggplot_build(r)$data[[1]]$fill))


mat <- read.csv("./genomic_data_visulization/RNAseq_mat_top20.csv",header=TRUE,row.names=1,
              check.names=FALSE)
data(airway)
data(mtcars)
cars <- mtcars
head(cars)

pheatmap(cars)
pheatmap(cars, scale="column")


head(mat)
pheatmap(mat)
pheatmap(mat, scale="row")

pheatmap(mat,scale="row",
         color=dichromat(colorRampPalette(c("navy", "white", "red"))(50),
                         type = "protan"))


dfh <- data.frame(sample=as.character(colnames(mat)),dex="Treatment")%>%
  column_to_rownames("sample")
dfh$dex <- ifelse(rownames(dfh) %in% c("508","512","516","520"),
                "untrt","trt")
dfh

pheatmap(mat,scale="row", annotation_col = dfh,
         annotation_colors=list(dex=c(trt=dichromat("orange"),untrt=dichromat("black"))),
         color=dichromat(colorRampPalette(c("navy", "white", "red"))(50)))



hm_gg<-as.ggplot(pheatmap(mat,scale="row", annotation_col=
            dfh,annotation_colors =list(dex=c(trt="orange",
            untrt="black")),color=colorRampPalette(c("navy", 
                                         "white", "red"))(50)))

hm_gg

heatmaply(mat, scale="row", margins=c(0.5,1,50,1),
          main="Interactively explore airway DE genes",
          col_side_colors=dfh)




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


