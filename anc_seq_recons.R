############################
#### Ancestral Sequence ####
####   reconstruction  #####
############################

rm(list=ls())

library(ape)
library(phangorn)
library(BactDating)

dir_prefix = '/drug_acquisition_results/'

######## FILES NEEDED ##########
mtb_tree_file = paste(dir_prefix, 'mtb.peru.wg.dates.tree', sep = '')
l2_tree_file = paste(dir_prefix, 'l2.assembly.arc.5e7.0.5.wg1_res.RData', sep = '')
l4_tree_file = paste(dir_prefix, 'l4.assembly.arc.1e7.0.5.wg2_res.RData', sep = '')
dna_file = paste(dir_prefix, 'dna_phyDat.assembly.RData', sep = '')
metadata_file = paste(dir_prefix, 'tb_1999_2016_noDups.metadata.csv', sep='')
input_dna = paste(dir_prefix, 'mtb.snps.dels.wg.withDR.assembly.finalSet.snpSites.masked.fa', sep='')
###############################

load(l2_tree_file)
l2_tree = reorder.phylo(res$inputtree, 'cladewise')

load(l4_tree_file)
l4_tree = reorder.phylo(res$tree, 'cladewise')


##### Prepare DNA PhyDat file #####

dna = read.dna(input_dna, format='fasta', as.matrix = T, as.character = T)

# Pure bases: A, C, G, T, -
# Ambiguous: IUPAC ambiguity bases
contrast = matrix(data = c(1,0,0,0,0,
                           0,1,0,0,0,
                           0,0,1,0,0,
                           0,0,0,1,0,
                           0,0,0,0,1,
                           1,1,0,0,0,
                           1,0,1,0,0,
                           1,0,0,1,0,
                           0,1,1,0,0,
                           0,1,0,1,0,
                           0,0,1,1,0,
                           1,1,1,0,0,
                           1,1,0,1,0,
                           1,0,1,1,0,
                           0,1,1,1,0,
                           1,1,1,1,1),
                  ncol = 5, byrow = TRUE)
dimnames(contrast) = list(c("a","c","g","t","-",
                            "m","r","w","s","y",
                            "k","v","h","d","b","n"),
                          c("a", "c", "g", "t", "-"))

dna_phydat = phyDat(dna, type="USER", contrast=contrast)


metadata = read.csv(metadata_file, stringsAsFactors = F, header = T)
metadata = metadata[!is.na(metadata$id),]

# Put dates on DNA names so it matches the phylogeny
for (i in 1:length(names(dna_phydat))) {
  date = 'NA'
  name = names(dna_phydat)[i]
  date = unique(metadata$year_sample[metadata$id == name])
  date = date[which(!is.na(date))]
  new_name = paste(name, '_', date, sep='')
  names(dna_phydat)[i] = gsub('_$','', new_name)
}

save(dna_phydat, file = paste(dir_prefix, 'dna_phyDat.assembly.RData', sep = ''))



##### Ancestral reconstruction #####

load(dna_file)
genome_len = 4411532

#L2 reconstruction
load(l2_tree_file)
ml_tree = reorder.phylo(res$inputtree, 'cladewise')
ml_tree$edge.length = ml_tree$edge.length/genome_len # Branch length in subs/site

# Maximum likelihood
fit = pml(ml_tree, dna_phydat[ml_tree$tip.label])
fit = optim.pml(fit,optEdge=FALSE,optRate=TRUE,optQ=TRUE) # Optimize only the rate and Q matrix
l2.anc.ml = ancestral.pml(fit, "ml") # Ancestral reconstruction
save(l2.anc.ml, file = paste(dir_prefix, 'l2.anc.ml.assemblyTree.RData', sep=''))


#L4 reconstruction
load(l4_tree_file)
ml_tree = reorder.phylo(res$inputtree, 'cladewise')
ml_tree$edge.length = ml_tree$edge.length/genome_len # Branch length in subs/site

# Maximum likelihood
fit = pml(ml_tree, dna_phydat[ml_tree$tip.label])
fit<-optim.pml(fit,optEdge=FALSE,optRate=TRUE,optQ=TRUE) # Optimize only the rate and Q matrix
l4.anc.ml = ancestral.pml(fit, "ml") # Ancestral reconstruction
save(l4.anc.ml, file = paste(dir_prefix, 'l4.anc.ml.assemblyTree.RData', sep=''))

