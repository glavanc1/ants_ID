library(phytools)

# reading the metadata containing the morphological identification
meta.all <- read.table("morpho_ID.txt", h=T)

# reading the phylogeny
tree <- read.tree("Camponotus_aligned.phylip.treefile")

# subsetting the metadata to only individuals included in the tree
meta <- meta.all[match(tree$tip.label, meta.all$CATALOGUENUMBER),]

# defining colours
meta$col.morpho <- "grey"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_aeth")] <- "#BE1E2D"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_fall")] <- "#EF3E23"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_herc")] <- "#27AAE1"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_lign")] <- "#F9ED32"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_pice")] <- "#006838"

legend_Campo <- c("C. ligniperda", "C. herculeanus", "C. aethiops", "C. fallax", "C. piceus")
col_Campo <- c("#F9ED32", "#27AAE1", "#BE1E2D", "#EF3E23", "#006838")

# Plot the tree
plot(tree, cex=0.5, label.offset = 0.01)
# add a square with the colour corresponding to the morphological species identification
tiplabels(pch = 15, offset = 0.005, col=meta$col.morpho)
# add a scale bar
add.scale.bar()

# add a legend
legend("topleft", col=col_Campo, pch=15, legend=legend_Campo, bty="n")

# Used to identify MRCA of individuals belonging to each species
# Not used in the Figure
nodelabels(cex=0.5, frame = "none")

meta$SPECIESID_GENETIC <- meta$SPECIESSUMMARY

# Correcting Camp_lign
# identifying all nodes and tips descendants from 218 (identified visually as the MRCA of all Camp_lign)
Camp_lign.nodes <- tree$tip.label[getDescendants(tree, node = 218)]
# remove nodes (which have indices outside length(tree$tip.label))
Camp_lign.tips <- Camp_lign.nodes[which(is.na(Camp_lign.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Camp_lign.tips)] <- "Camp_lign"

# Correcting Camp_herc
# identifying all nodes and tips descendants from 183 (identified visually as the MRCA of all Camp_herc)
Camp_lign.nodes <- tree$tip.label[getDescendants(tree, node = 183)]
# remove nodes (which have indices outside length(tree$tip.label))
Camp_lign.tips <- Camp_lign.nodes[which(is.na(Camp_lign.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Camp_lign.tips)] <- "Camp_herc"

# Write final table
write.table(meta[,c(1,2,4)], "Camponotus_ID.txt", quote=F, row.name=F)
