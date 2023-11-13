library(phytools)

# reading the metadata containing the morphological identification
meta.all <- read.table("morpho_ID.txt", h=T)

# reading the phylogeny
tree <- read.tree("Tetramorium_aligned.phylip.treefile_rerooted")
# Due to an error with reroot(), the tree was rerooted manually in Figtree

# subsetting the metadata to only individuals included in the tree
meta <- meta.all[match(tree$tip.label, meta.all$CATALOGUENUMBER),]

# defining colours
meta$col.morpho <- "grey"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tetr_caes")] <- "#1B75BC"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tetr_impu")] <- "#C7315D"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tetr_immi")] <- "#86C540"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tetr_semi")] <- "#71489D"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tetr_alpe")] <- "#D4A129"

legend_Tetra <- c("T. caespitum", "T. impurum", "T. immigrans", "T. alpestre", "T. semilaevae")
col_Tetra <- c("#1B75BC", "#C7315D", "#86C540", "#D4A129", "#71489D")

# Plot the tree
plot(tree, cex=0.5, label.offset = 0.01)
# add a square with the colour corresponding to the morphological species identification
tiplabels(pch = 15, offset = 0.005, col=meta$col.morpho[match(tree$tip.label, meta$CATALOGUENUMBER)])
# add a scale bar
add.scale.bar()

# add a legend
legend("topleft", col=col_Tetra, pch=15, legend=legend_Tetra, bty="n")

# Used to identify MRCA of individuals belonging to each species
# Not used in the Figure
nodelabels(cex=0.8, frame = "none")

meta$SPECIESID_GENETIC <- meta$SPECIESSUMMARY

# Correcting Tetr_caes
# identifying all nodes and tips descendants from 566 (identified visually as the MRCA of all Tetr_caes)
Tetr_caes.nodes <- tree$tip.label[getDescendants(tree, node = 566)]
# remove nodes (which have indices outside length(tree$tip.label))
Tetr_caes.tips <- Tetr_caes.nodes[which(is.na(Tetr_caes.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Tetr_caes.tips)] <- "Tetr_caes"

# Correcting Tetr_impu
# identifying all nodes and tips descendants from 465 (identified visually as the MRCA of both Tetr_impu clades as well as Tetr_immi, which are the changed below)
Tetr_impu.nodes <- tree$tip.label[getDescendants(tree, node = 465)]
# remove nodes (which have indices outside length(tree$tip.label))
Tetr_impu.tips <- Tetr_impu.nodes[which(is.na(Tetr_impu.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Tetr_impu.tips)] <- "Tetr_impu"

# Correcting Tetr_immi
# identifying all nodes and tips descendants from 266 (identified visually as the MRCA of all Tetr_immi)
Tetr_immi.nodes <- tree$tip.label[getDescendants(tree, node = 266)]
# remove nodes (which have indices outside length(tree$tip.label))
Tetr_immi.tips <- Tetr_immi.nodes[which(is.na(Tetr_immi.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Tetr_immi.tips)] <- "Tetr_immi"


table(meta$SPECIESSUMMARY, meta$SPECIESID_GENETIC)

# Write final table
write.table(meta[,c(1,2,4)], "Tetramorium_ID.txt", quote=F, row.name=F)
