# Genetic identification of Swiss ants

This repository contains the code used for the genetic identification of ant specimens used in the manuscript

### LEVERAGING CITIZEN SCIENCE TO ASSESS RICHNESS, DIVERSITY, AND ABUNDANCE IN ANT COMMUNITIES
by Tim M. Szewczyk, Guillaume Lavanchy, Anne Freitag, Cleo Bertelsmeier, Tanja Schwander

The COI sequences used in this study are available on GenBank. Accession numbers are given in XXX.

We conducted a separate analysis for each of the four genus: *Camponotus*, *Tapinoma*, *Temnothorax* and *Tetramorium*.
For each genus, the concatenated fasta file including reference individuals is available in `GENUS/$GENUS.fa` (replace "GENUS" by the genus name). The sequences from reference individuals are named with their GenBank accession numbers. There are no reference individuals for *Camponotus* as the purpose of genotyping was to confirm some uncertain indentifications, but there was no suspicion of additional species (this was confirmed *a posteriori*). We added one individual from a different genus as outgroup.

We aligned the sequences using `mafft` v7.475:
```bash
mafft --localpair  --maxiterate 16 --phylipout --inputorder GENUS/GENUS.fa > GENUS/GENUS_aligned.phylip
```

We then reconstructed the phylogenies using `iqtree` v2.0.6

```bash
iqtree2 -T 1 -s GENUS/GENUS_aligned.phylip -m MFP -B 1000
```

The output trees were analysed in R. The scripts are available in `GENUS/GENUS.R`
