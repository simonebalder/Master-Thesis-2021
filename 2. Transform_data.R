# Libraries
library(microbiome)
library(RColorBrewer)
library(ggpubr)

# Load data from Assign_taxonomy.R
load('R workspaces/ps.RData')


## Richness and evenness (chao1 and shannon) ##

# Richness plot (chao1, Shannon) on the non-transformed data
plot_richness(ps, x = "Setup_2", color="Time", measures=c("Chao1", "Shannon"))

# Estimate richness (chao1, shannon) on the non-transformed data
ps.richness <- estimate_richness(ps, split = TRUE, measures = c("Chao1", "Shannon"))
write.table(ps.richness, file = "Files/ps.richness.csv", sep = ";",
            row.names = TRUE, col.names = NA)


## Transformed data ##

# Transform to even sampling depth
ps.transformed = transform_sample_counts(ps, function(x) 1E6 * x/sum(x)) # with P. piscicida B39bio
ps.transformed.noB39bio = transform_sample_counts(ps.noB39bio, function(x) 1E6 * x/sum(x)) # without P. piscicida B39bio
ps.pcs.transformed = transform_sample_counts(ps.pcs, function(x) 1E6 * x/sum(x))

save.image(file='R workspaces/Transformed.RData')

# NMDS ordination
out.nmds <- ordinate(ps.pcs.transformed, method = "NMDS", distance = "bray")
nmds_plot = plot_ordination(ps.transformed, out.nmds, color = "Setup_2", shape = "Time") + geom_point(size = 3)
nmds_plot


## Relative abundance with B39bio ##

# We need to set Palette
taxic <- as.data.frame(ps.pcs.transformed@tax_table) # this will help in setting large color options
getPalette = colorRampPalette(brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps.pcs.transformed) <- new.tax # incroporate into phyloseq Object

# Edit the unclassified taxa 
tax_table(ps.pcs.transformed)[tax_table(ps.pcs.transformed)[, "Family"] == "", "Family"] <- "Unclassified family"
tax_table(ps.pcs.transformed)[tax_table(ps.pcs.transformed)[, "Order"] == "", "Order"] <- "Unclassified order"
tax_table(ps.pcs.transformed)[tax_table(ps.pcs.transformed)[, "Class"] == "", "Class"] <- "Unclassified class"

# Taxonomic names in italic
guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))

# merge at family level, only the top 12 taxa
ps.transformed.fam <- microbiome::aggregate_top_taxa(ps.pcs.transformed, "Family", top = 12)
# merge at order level, only the top 20 taxa
ps.transformed.order <- microbiome::aggregate_top_taxa(ps.pcs.transformed, "Order", top = 20)
# merge at class level, only the top 12 taxa
ps.transformed.class <- microbiome::aggregate_top_taxa(ps.pcs.transformed, "Class", top = 12)

# Use transform function of microbiome to convert it to relative abundance
ps.transformed.fam.rel <- microbiome::transform(ps.transformed.fam, "compositional")
ps.transformed.order.rel <- microbiome::transform(ps.transformed.order, "compositional")
ps.transformed.class.rel <- microbiome::transform(ps.transformed.class, "compositional")

write.table(otu_table(ps.transformed.fam.rel), file = "Files/rel.abn.fam20.csv", sep = ";",
            row.names = TRUE, col.names = NA)

# Colors
# Colors in R: https://www.r-graph-gallery.com/ggplot2-color.html 
col <- c("Nocardioidaceae" = "darkgoldenrod1", "Micrococcaceae" = "burlywood", "Kordiimonadaceae" = "darkslategray1", "Flavobacteriaceae" = "darksalmon", "Flammeovirgaceae" = "darkgreen", 
         "Comamonadaceae" = "darkolivegreen2", "Chromatiaceae" = "cornflowerblue", "Burkholderiaceae" = "cadetblue3", "Other" = "darkseagreen3", "Pseudoalteromonadaceae" = "coral4",
         "Pseudomonadaceae" = "darkslateblue", "Rhodobacteraceae" = "bisque", "Unknown" = "white", "Zhongshania" = "chocolate")

# Plot relative abundance, family
plot.composition.relAbun.fam <- plot_composition(ps.transformed.fam.rel,
                                             sample.sort = "sample.id",
                                             x.label = "samples") 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + theme(legend.position = "bottom") 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + scale_fill_manual("Family",values = col) + theme_bw() 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun.fam)

# Plot relative abundance, class
plot.composition.relAbun.class <- plot_composition(ps.transformed.class.rel,
                                             sample.sort = "sample.id",
                                             x.label = "samples") 
plot.composition.relAbun.class <- plot.composition.relAbun.class + theme(legend.position = "bottom") 
plot.composition.relAbun.class <- plot.composition.relAbun.class + scale_fill_brewer("Class", palette = "Paired") + theme_bw() 
plot.composition.relAbun.class <- plot.composition.relAbun.class + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun.class <- plot.composition.relAbun.class + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun.class)


## Relative abundance without B39bio ##

# We need to set Palette
taxic <- as.data.frame(ps.transformed.noB39bio@tax_table) # this will help in setting large color options
getPalette = colorRampPalette(brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps.transformed) <- new.tax # incroporate into phyloseq Object

# Edit the unclassified taxa 
tax_table(ps.transformed.noB39bio)[tax_table(ps.transformed.noB39bio)[, "Family"] == "", "Family"] <- "Unclassified family"
tax_table(ps.transformed.noB39bio)[tax_table(ps.transformed.noB39bio)[, "Class"] == "", "Class"] <- "Unclassified class"

# merge at family level, only the top 12 taxa
ps.transformed.fam <- microbiome::aggregate_top_taxa(ps.transformed.noB39bio, "Family", top = 12)
# merge at class level, only the top 12 taxa
ps.transformed.class <- microbiome::aggregate_top_taxa(ps.transformed.noB39bio, "Class", top = 12)

# Use transform function of microbiome to convert it to relative abundance
ps.transformed.fam.rel <- microbiome::transform(ps.transformed.fam, "compositional")
ps.transformed.class.rel <- microbiome::transform(ps.transformed.class, "compositional")

write.table(otu_table(ps.transformed.fam.rel), file = "Files/rel.abn.family.noB39bio.csv", sep = ";",
            row.names = TRUE, col.names = NA)

# Plot relative abundance, family
plot.composition.relAbun.fam <- plot_composition(ps.transformed.fam.rel,
                                                 sample.sort = "sample.id",
                                                 x.label = "samples") 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + theme(legend.position = "bottom") 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + scale_fill_manual("Family", values = col) + theme_bw() 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun.fam <- plot.composition.relAbun.fam + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun.fam)

# Plot relative abundance, class
plot.composition.relAbun.class <- plot_composition(ps.transformed.class.rel,
                                                   sample.sort = "sample.id",
                                                   x.label = "samples") 
plot.composition.relAbun.class <- plot.composition.relAbun.class + theme(legend.position = "bottom") 
plot.composition.relAbun.class <- plot.composition.relAbun.class + scale_fill_brewer("Class", palette = "Paired") + theme_bw() 
plot.composition.relAbun.class <- plot.composition.relAbun.class + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun.class <- plot.composition.relAbun.class + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun.class)



save.image(file='R workspaces/Relative_abundance.RData')



