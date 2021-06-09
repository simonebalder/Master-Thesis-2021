setwd("C:/Users/sbald/Dropbox/Simone/DTU/10. semester/Speciale/Bioinformatic analysis/R/R scripts")

# Libraries
library(rRDPData)
library("devtools")
library(Rcpp)
library("dada2")
library(plyr)
library("Biostrings")
library("phyloseq")
library(tidyverse)
library(readr)
library(seqinr)
library(ape)
library(vegan)

## Train classifier (rdp) ##

# Load DNA sequences (Qiime2 output)
denoised_seq_qiime = readDNAStringSet("Files/dna-sequences.fasta")
# Initialize random number generator for reproducibility
set.seed(100) 
# Train classifier using the rdp_train_set_16.fa
# https://zenodo.org/record/801828#.YLCVxagzbcc 
rdp_assign_taxonomy <- assignTaxonomy(denoised_seq_qiime, "Files/rdp_train_set_16.fa", multithread=FALSE, tryRC =TRUE)
rdp_add_species_allowMultiple <- assignSpecies(denoised_seq_qiime, "Files/rdp_species_assignment_16.fa.gz", allowMultiple = TRUE, verbose = TRUE)
save.image(file='R workspaces/TrainedRdp05052021_qiime.RData')
# Add species level to the taxonomy
taxonomy_table <- addSpecies(rdp_assign_taxonomy, "Files/rdp_species_assignment_16.fa", verbose=TRUE)
save.image(file = 'R workspaces/TrainedRdp05052021_qiime.RData')

#Loading in the asv table (Qiime2 output)
asv <- read.table("Files/asv_table.tsv", row.names = 1, header = TRUE, sep = "\t")


## Rarefaction curves ##

# Convert from phylo object to dataframe
asv_df <- as.data.frame(asv)
rarecurve(t(asv_df), step=1000, cex=0.6, col = "blue", lwd=2, ylab="ASVs", label=FALSE, xlim=c(0, 25000), ylim=c(0, 3000))


## Make phyloseq object ##

#Convert to phyloseq object
asv_phyloseq <- otu_table(asv, taxa_are_rows = TRUE)

# Import metadata
metadata <- read_tsv('Files/metadata_remove_outliers.tsv') # with positive controls
metadata_pcs <- read_tsv('Files/metadata_pcs.tsv') # without positive controls
# Delete the row with Qiime2 category codes
metadata <- metadata[-c(1),]
metadata_pcs <- metadata_pcs[-c(1),]

# Making phyloseq objects
ASV = otu_table(data.frame(asv_phyloseq), taxa_are_rows = TRUE)
ASV.seq <- phyloseq::phyloseq(ASV, denoised_seq_qiime) # combining asv table with the sequences
rownames(taxonomy_table) <- rownames(otu_table(ASV.seq)) # naming rownames the same as asv.table
TAX = tax_table(taxonomy_table)
META = sample_data(data.frame(metadata, row.names = metadata$`sample-id`))
META_pcs = sample_data(data.frame(metadata_pcs, row.names = metadata_pcs$`sample-id`))

# Make phyloseq object with all three files
ps <- phyloseq(ASV, TAX, META)
ps.pcs <- phyloseq(ASV, TAX, META_pcs)

# Remove the ASVs that stem from P. piscicida B39bio
B39bio_ASV = c("98fc15a699eddfcce377b509fb090405", "6ef0ca040f14accd24943c31630b9e9b", "7cf85f679f8e68e9ec4615ee31774c06", "ef21f351514db2a37c0e9b3b02fcf847")
allTaxa = taxa_names(ps.pcs)
myTaxa <- allTaxa[!(allTaxa %in% B39bio_ASV)]
ps.noB39bio = prune_taxa(myTaxa, ps.pcs)

save.image(file='R workspaces/ps.RData')
