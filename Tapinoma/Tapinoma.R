library(phytools)

# reading the metadata containing the morphological identification
meta.all <- read.table("morpho_ID.txt", h=T)

# reading the phylogeny
tree.unrooted <- read.tree("Tapinoma_aligned.phylip.treefile")

# subsetting the metadata to only individuals included in the tree
meta <- meta.all[match(tree$tip.label, meta.all$CATALOGUENUMBER),]

# defining colours
meta$col.morpho <- "grey"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tapi_erra")] <- "#FCE623"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tapi_nige_gr")] <- "#420053"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Tapi_subb")] <- "#1E8F8A"

legend_Tapi <- c("Tapi_erra", "Tapi_nige_gr", "Tapi_subb")
col_Tapi <- c("#FCE623", "#420053", "#1E8F8A")

# Plot the tree
plot(tree.unrooted, cex=0.5, label.offset = 0.01)
# add a square with the colour corresponding to the morphological species identification
tiplabels(pch = 15, offset = 0.005, col=meta$col.morpho)
# add a scale bar
add.scale.bar()

nodelabels(cex=0.5, frame = "none")

# Reroot the tree so that the nigerrimum group is monophyletic 
# and a sister clade to (erraticum, subboreale)
tree <- reroot(tree.unrooted, 289)

# Plot the tree
plot(tree, cex=0.5, label.offset = 0.01)
# add a square with the colour corresponding to the morphological species identification
tiplabels(pch = 15, offset = 0.005, col=meta$col.morpho[match(tree$tip.label, meta$CATALOGUENUMBER)])
# add a scale bar
add.scale.bar()

# add a legend
legend("topleft", col=col_Tapi, pch=15, legend=legend_Tapi, bty="n")

# Used to identify MRCA of individuals belonging to each species
# Not used in the Figure
nodelabels(cex=0.5, frame = "none")

meta$SPECIESID_GENETIC <- meta$SPECIESSUMMARY

# Correcting Tapi_erra
# identifying all nodes and tips descendants from 298 (identified visually as the MRCA of all Tapi_erra)
Tapi_erra.nodes <- tree$tip.label[getDescendants(tree, node = 298)]
# remove nodes (which have indices outside length(tree$tip.label))
Tapi_erra.tips <- Tapi_erra.nodes[which(is.na(Tapi_erra.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Tapi_erra.tips)] <- "Tapi_erra"

# Correcting Tapi_subb
# identifying all nodes and tips descendants from 472 (identified visually as the MRCA of all Tapi_subb)
Tapi_subb.nodes <- tree$tip.label[getDescendants(tree, node = 472)]
# remove nodes (which have indices outside length(tree$tip.label))
Tapi_subb.tips <- Tapi_subb.nodes[which(is.na(Tapi_subb.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Tapi_subb.tips)] <- "Tapi_subb"

# Correcting Tapi_magn
# identifying all nodes and tips descendants from 266 (identified visually as the MRCA of all Tapi_magn)
Tapi_magn.nodes <- tree$tip.label[getDescendants(tree, node = 266)]
# remove nodes (which have indices outside length(tree$tip.label))
Tapi_magn.tips <- Tapi_magn.nodes[which(is.na(Tapi_magn.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Tapi_magn.tips)] <- "Tapi_magn"

# Correcting Tapi_dari
# identifying all nodes and tips descendants from 257 (identified visually as the MRCA of all Tapi_dari)
Tapi_dari.nodes <- tree$tip.label[getDescendants(tree, node = 257)]
# remove nodes (which have indices outside length(tree$tip.label))
Tapi_dari.tips <- Tapi_dari.nodes[which(is.na(Tapi_dari.nodes)==F)]
# Assign new species ID
meta$SPECIESID_GENETIC[which(meta$CATALOGUENUMBER %in% Tapi_dari.tips)] <- "Tapi_dari"

table(meta$SPECIESSUMMARY, meta$SPECIESID_GENETIC)

# Write final table
write.table(meta[,c(1,2,4)], "Tapinoma_ID.txt", quote=F, row.name=F)
