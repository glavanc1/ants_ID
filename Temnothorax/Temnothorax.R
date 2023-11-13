library(phytools)

# reading the metadata containing the morphological identification
meta.all <- read.table("morpho_ID.txt", h=T)

# reading the phylogeny
tree.unrooted <- read.tree("Temnothorax_aligned.phylip.treefile")

# subsetting the metadata to only individuals included in the tree
meta <- meta.all[match(tree.unrooted$tip.label, meta.all$CATALOGUENUMBER),]

# defining colours
meta$col.morpho <- "grey"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_affi")] <- "#EC008C"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_nigr")] <- "#3E499E"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_nyla")] <- "#F58122"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_unif")] <- "#74C483"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_cort")] <- "#D2D028"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_parv")] <- "#B92025"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_tube")] <- "#6F3686"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_inte")] <- "#81CEC7"

legend_Temno <- c("T. affinis", "T. corticalis", "T. interruptus", "T. nigriceps", "T.nylanderi", "T. parvulus", "T. tuberum", "T. unifasciatus")
col_Temno <- c("#EC008C", "#D2D028", "#81CEC7", "#3E499E", "#F58122", "#B92025", "#6F3686", "#74C483")

# Plot the tree
plot(tree.unrooted, cex=0.5, label.offset = 0.01)
# add a square with the colour corresponding to the morphological species identification
tiplabels(pch = 15, offset = 0.005, col=meta$col.morpho)
# add a scale bar
add.scale.bar()

# add a legend
legend("topright", col=col_Temno, pch=15, legend=legend_Temno, bty="n")

nodelabels(cex=0.5, frame = "none")

# Reroot the tree so that parvulus is an outgroup to other species
tree <- reroot(tree.unrooted, 205)

# Plot the tree
plot(tree, cex=0.5, label.offset = 0.01)
# add a square with the colour corresponding to the morphological species identification
tiplabels(pch = 15, offset = 0.005, col=meta$col.morpho[match(tree$tip.label, meta$CATALOGUENUMBER)])
# add a scale bar
add.scale.bar()

# add a legend
legend("bottomright", col=col_Temno, pch=15, legend=legend_Temno, bty="n")

# Used to identify MRCA of individuals belonging to each species
# Not used in the Figure
nodelabels(cex=0.5, frame = "none")

meta$SPECIESID_GENETIC <- meta$SPECIESSUMMARY

# Correcting Temn_parv
# identifying all nodes and tips descendants from 380 (identified visually as the MRCA of all Temn_parv)
Temn_parv.nodes <- tree$tip.label[getDescendants(tree, node = 380)]
# remove nodes (which have indices outside length(tree$tip.label))
Temn_parv.tips <- Temn_parv.nodes[which(is.na(Temn_parv.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Temn_parv.tips)] <- "Temn_parv"

# Correcting Temn_nyla
# identifying all nodes and tips descendants from 194 (identified visually as the MRCA of all Temn_nyla)
Temn_nyla.nodes <- tree$tip.label[getDescendants(tree, node = 194)]
# remove nodes (which have indices outside length(tree$tip.label))
Temn_nyla.tips <- Temn_nyla.nodes[which(is.na(Temn_nyla.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Temn_nyla.tips)] <- "Temn_nyla"

# Correcting Temn_inte
# identifying all nodes and tips descendants from 377 (identified visually as the MRCA of all Temn_inte)
Temn_inte.nodes <- tree$tip.label[getDescendants(tree, node = 377)]
# remove nodes (which have indices outside length(tree$tip.label))
Temn_inte.tips <- Temn_inte.nodes[which(is.na(Temn_inte.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Temn_inte.tips)] <- "Temn_inte"

# Correcting Temn_unif
# identifying all nodes and tips descendants from 321 (identified visually as the MRCA of all Temn_unif)
Temn_unif.nodes <- tree$tip.label[getDescendants(tree, node = 321)]
# remove nodes (which have indices outside length(tree$tip.label))
Temn_unif.tips <- Temn_unif.nodes[which(is.na(Temn_unif.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Temn_unif.tips)] <- "Temn_unif"

# Correcting Temn_affi
# identifying all nodes and tips descendants from 348 (identified visually as the MRCA of all Temn_affi)
Temn_affi.nodes <- tree$tip.label[getDescendants(tree, node = 348)]
# remove nodes (which have indices outside length(tree$tip.label))
Temn_affi.tips <- Temn_affi.nodes[which(is.na(Temn_affi.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Temn_affi.tips)] <- "Temn_affi"

# Correcting two individuals with a weird haplotype
# identifying all nodes and tips descendants from 347 (identified visually as the MRCA of the two)
Temn.nodes <- tree$tip.label[getDescendants(tree, node = 347)]
# remove nodes (which have indices outside length(tree$tip.label))
Temn.tips <- Temn.nodes[which(is.na(Temn.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Temn.tips)] <- "Temn"

table(meta$SPECIESSUMMARY, meta$SPECIESID_GENETIC)

# Write final table
write.table(meta[,c(1,2,4)], "Temnothorax_ID.txt", quote=F, row.name=F)
