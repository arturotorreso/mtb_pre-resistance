library(ape)
library(phangorn)
library(BactDating)
library(data.table)
library(seqinr)


old2new.phyDat <- function(obj) {
  att <- attributes(obj)
  l <- length(obj)
  contrast <- attr(obj, "contrast")
  nr <- attr(obj, "nr")
  X <- matrix(rep(rowSums(contrast), each = nr), nrow = nr)
  res <- vector("list", l)
  for (i in 1:l) {
    tmp <- X - tcrossprod(obj[[i]], contrast)
    res[[i]] <- unlist(apply(tmp, 1, function(x) which(x == min(x))[1]))
  }
  attributes(res) <- att
  res
}

syn_nonNyn = function(pos, gen_st, gen_end, ref, alt, ref_genome){
  
  library(Biostrings)
  
  if (alt == '-'){
    return('NonSyn')
  }
  
  if (gen_st == 'intergenic'){
    return('NonGene')
  }
  
  codon = c('a' = 't', 't' = 'a', 'c' = 'g', 'g' = 'c')
  
  seq_gen = ref_genome$Chromosome[gen_st:gen_end]
  seq_gen_rev = rev(ref_genome$Chromosome[gen_st:gen_end])
  prot = translate(DNAString(paste(toupper(codon[seq_gen_rev]), collapse = '')))
  
  ref_genome_mut = ref_genome
  ref_genome_mut$Chromosome[pos] = tolower(alt)
  seq_gen_mut = ref_genome_mut$Chromosome[gen_st:gen_end]
  seq_gen_mut_rev = rev(ref_genome_mut$Chromosome[gen_st:gen_end])
  prot_mut = translate(DNAString(paste(toupper(codon[seq_gen_mut_rev]), collapse = '')))
  
  
  if (as.character(prot) == as.character(prot_mut)){
    return('Syn')
  } else {
    return('NonSyn')
  }
  # split(unname(codon[seq_gen_rev]), ceiling(seq_along(seq_gen_rev)/3))
}


dir_prefix = '/pre_resistance/data/'


######## FILES NEEDED ##########
l2_tree_file = paste(dir_prefix, 'l2.assembly.arc.5e7.0.5.wg1_res.RData', sep = '')
l4_tree_file = paste(dir_prefix, 'l4.assembly.arc.1e7.0.5.wg2_res.RData', sep = '')
l2.anc.seq_file = paste(dir_prefix, 'l2.anc.ml.imp.assemblyTree.RData', sep = '')
l4.anc.seq_file = paste(dir_prefix, 'l4.anc.ml.imp.assemblyTree.RData', sep = '')
l4_dna_file = paste(dir_prefix, 'l4.dna_phyDat.assembly.imputed.RData', sep = '')
l2_dna_file = paste(dir_prefix, 'l2.dna_phyDat.assembly.imputed.RData', sep = '')
genome_idx_file = paste(dir_prefix, 'mtb.snps.dels.wg.withDR.assembly.finalSet.snpSites.masked.idx', sep = '')
tbdr_db_file = paste(dir_prefix, 'tb_profiler_db.complete.good.posStrRef.tsv', sep = '')
metadata_file = paste(dir_prefix, 'tb_1999_2016_noDups.metadata.csv', sep='')
gene_file = paste(dir_prefix, 'h37rv_genes.bed', sep = '')
ref_genome = paste(dir_prefix, 'Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa', sep = '')
###############################


#########################
###### LOAD FILES #######
#########################

# Load phylogenies
load(l2_tree_file)
l2_tree = reorder.phylo(res$tree, 'cladewise')

load(l4_tree_file)
l4_tree = reorder.phylo(res$tree, 'cladewise')

# Load genome-snp index
genome_idx = read.table(genome_idx_file,
                        stringsAsFactors = F, header = T)

# Load metadata
metadata = read.csv(metadata_file, stringsAsFactors = F, header = T)
metadata = metadata[!is.na(metadata$id),]
metadata$id = paste(metadata$id, '_', metadata$year_sample, sep='')
metadata = metadata[metadata$id %in% c(l4_tree$tip.label, l2_tree$tip.label),]

mtb_genes = read.delim(gene_file, sep = '\t', header = T)
ref_genome = read.fasta(ref_genome)

# Load DNA - Alignment and ancestral sequences
load (l2.anc.seq_file)
load (l4.anc.seq_file)

load (l4_dna_file)
l4_phydat = dna_phydat

load(l2_dna_file)
l2_phydat = dna_phydat

# Load the resistance files
l2_res_file = paste(dir_prefix, 'l2_t_to_res_snp.assemblyTree.txt', sep = '')
l4_res_file = paste(dir_prefix, 'l4_t_to_res_snp.assemblyTree.txt', sep = '')

l2_t_to_res_snp = read.table(l2_res_file,stringsAsFactors = F, header = T)
l4_t_to_res_snp = read.table(l4_res_file,stringsAsFactors = F, header = T)

l4_t_to_res_snp$lineage = ifelse (is.na(l4_t_to_res_snp$lineage), 'lineage4',
                                  l4_t_to_res_snp$lineage)

l4_t_to_res_snp = l4_t_to_res_snp[(l4_t_to_res_snp$int_st == 0 & l4_t_to_res_snp$ext_st == 1) | 
                                    (l4_t_to_res_snp$int_st == 0 & l4_t_to_res_snp$ext_st == 0),]

l2_t_to_res_snp = l2_t_to_res_snp[(l2_t_to_res_snp$int_st == 0 & l2_t_to_res_snp$ext_st == 1) | 
                                    (l2_t_to_res_snp$int_st == 0 & l2_t_to_res_snp$ext_st == 0),]

l2_phen_res = l2_t_to_res_snp[,c('int_node', 'int_year', 'ext_year', 'ext_st', 'lineage')]
l4_phen_res = l4_t_to_res_snp[,c('int_node', 'int_year', 'ext_year', 'ext_st', 'lineage')]

colnames(l2_phen_res) = colnames(l4_phen_res) = c('node', 'int_year', 'ext_year', 'status', 'lineage')

l2_phen_res$time = l2_phen_res$ext_year - l2_phen_res$int_year
l4_phen_res$time = l4_phen_res$ext_year - l4_phen_res$int_year



#############################
##  Prepare SNP alignment  ##
#############################

l2_matrix = as.character(l2_phydat)
l4_matrix = as.character(l4_phydat)


###### GET SITES WITH AF > 0.01 #######

pos = 1:ncol(l2_matrix)
l2_pos = pos[apply(l2_matrix, 2, function(x) if (length(unique(x)) > 1) sum(sort(table(x), decreasing = T)[2:length(table(x))])/length(x) >= 0.01 else FALSE)]
l4_pos = pos[apply(l4_matrix, 2, function(x) if (length(unique(x)) > 1) sum(sort(table(x), decreasing = T)[2:length(table(x))])/length(x) >= 0.01 else FALSE)]




###### GET AF IN THE POPULATION #######

bases = c('a', 'c', 'g', 't', '-')
l2_freq = matrix(nrow = ncol(l2_matrix), ncol = (length(bases)), data = 0)
l4_freq = matrix(nrow = ncol(l4_matrix), ncol = (length(bases)), data = 0)

colnames(l2_freq) = colnames(l4_freq) = bases
rownames(l2_freq) = rownames(l4_freq) = 1:ncol(l2_matrix)

for (p in rownames(l2_freq)){
  l2_freq[p,] = round(table(factor(l2_matrix[,as.integer(p)], levels = bases))/length(l2_matrix[,as.integer(p)]),4)
  l4_freq[p,] = round(table(factor(l4_matrix[,as.integer(p)], levels = bases))/length(l4_matrix[,as.integer(p)]),4)
}



###### GET ALIGNMENT #######

l2_seqs = old2new.phyDat(l2.anc.ml[unique(as.character(l2_phen_res$node))])
l4_seqs = old2new.phyDat(l4.anc.ml[unique(as.character(l4_phen_res$node))])

l2_dna = as.character(l2_seqs)
l4_dna = as.character(l4_seqs)

l2_snps = l2_dna[,l2_pos]
l4_snps = l4_dna[,l4_pos]

colnames(l2_snps) = genome_idx$pos[l2_pos]
colnames(l4_snps) = genome_idx$pos[l4_pos]


###### GET BINARY ALIGNMENT #######

snps_matrix = matrix(nrow = nrow(l2_snps), ncol = 0)
row.names(snps_matrix) = rownames(l2_snps)


l2_bin = apply(l2_snps, 2, function(x) 
  matrix(ifelse(x %in% names(sort(table(x), decreasing = T)[-1]), 1, 0)))
rownames(l2_bin) = rownames(l2_snps)
colnames(l2_bin) = sapply(colnames(l2_bin), function(x) 
  paste(x,ifelse(length(unique(l2_snps[,x]))>1,paste(toupper(names(sort(table(l2_snps[,x]),decreasing=T)[-1])),collapse = '_'),
                paste(toupper(names(sort(table(l2_snps[,x]),decreasing=T)[1])),collapse = '_')),sep = '.'))

l4_bin = apply(l4_snps, 2, function(x) 
  matrix(ifelse(x %in% names(sort(table(x), decreasing = T)[-1]), 1, 0)))
rownames(l4_bin) = rownames(l4_snps)
colnames(l4_bin) = sapply(colnames(l4_bin), function(x) 
  paste(x,ifelse(length(unique(l4_snps[,x]))>1,paste(toupper(names(sort(table(l4_snps[,x]),decreasing=T)[-1])),collapse = '_'),
                 paste(toupper(names(sort(table(l4_snps[,x]),decreasing=T)[1])),collapse = '_')),sep = '.'))



###############################
######   Kinship matrix  ######
###############################

### ALL FROM R ###
write.dna(apply(l2_dna, 2, function(x) toupper(x)), format = 'fasta', file = paste(dir_prefix, 'l2_gwas.snps.fa', sep = ''), nbcol = -1, colsep = '')
write.dna(apply(l4_dna, 2, function(x) toupper(x)), format = 'fasta', file = paste(dir_prefix, 'l4_gwas.snps.fa', sep = ''), nbcol = -1, colsep = '')

system(paste("snp-sites -v -o /pre_resistance/data/l2_gwas.snps.vcf /pre_resistance/data/l2_gwas.snps.fa; grep \">\" /pre_resistance/data/l2_gwas.snps.fa | tr -d \">\" > /pre_resistance/data/l2_samples.txt; similarity_pyseer --vcf /pre_resistance/data/l2_gwas.snps.vcf /pre_resistance/data/l2_samples.txt > /pre_resistance/data/l2.gg.snps.txt'", sep = ''),intern=TRUE)
system(paste("snp-sites -v -o /pre_resistance/data/l4_gwas.snps.vcf /pre_resistance/data/l4_gwas.snps.fa; grep \">\" /pre_resistance/data/l4_gwas.snps.fa | tr -d \">\" > /pre_resistance/data/l4_samples.txt; similarity_pyseer --vcf /pre_resistance/data/l4_gwas.snps.vcf /pre_resistance/data/l4_samples.txt > /pre_resistance/data/l4.gg.snps.txt'", sep = ''),intern=TRUE)

l2_gg.snps2 = read.table(paste(dir_prefix, 'l2.gg.snps.txt', sep = ''), header = T, check.names = F)
l4_gg.snps2 = read.table(paste(dir_prefix, 'l4.gg.snps.txt', sep = ''), header = T, check.names = F)


l2_gg.eigen2 = eigen(l2_gg.snps2, symmetric = T)
l4_gg.eigen2 = eigen(l4_gg.snps2, symmetric = T)


row.names(l2_gg.eigen2$vectors) = row.names(l2_gg.snps2)
row.names(l4_gg.eigen2$vectors) = row.names(l4_gg.snps2)

colnames(l2_gg.eigen2$vectors) = paste('Axis', 1:ncol(l2_gg.eigen2$vectors), sep = '')
colnames(l4_gg.eigen2$vectors) = paste('Axis', 1:ncol(l4_gg.eigen2$vectors), sep = '')



##############################
##  Prepare GENE alignment  ##
##############################

colnames(l2_dna) = colnames(l4_dna) = genome_idx$pos

l2_genes = matrix(data = 0, nrow = nrow(l2_dna), ncol = length(mtb_genes$gene),
                  dimnames = list(rownames(l2_dna),mtb_genes$gene))
l4_genes = matrix(data = 0, nrow = nrow(l4_dna), ncol = length(mtb_genes$gene),
                  dimnames = list(rownames(l4_dna),mtb_genes$gene))


l2_dna_snps = l2_dna[,apply(l2_dna, 2, function(x) length(unique(x))>1)]
l4_dna_snps = l4_dna[,apply(l4_dna, 2, function(x) length(unique(x))>1)]


pos = sort(as.integer(colnames(l2_dna_snps)))
snp_df = c()
for (p in pos){
  p = as.integer(p)
  
  ref = ref_genome$Chromosome[p]
  alts = sort(table(l2_dna_snps[,as.character(p)]), decreasing = T)
  alts = names(alts[-which(names(alts) == ref)][1])
  
  snp_df = rbind(snp_df, cbind('Chromosome', p, p, toupper(ref), toupper(alts), '+'))
}

write.table(snp_df, '/pre_resistance/data/snps.txt', col.names = F, row.names = F, quote = F, sep = ' ')

snpEff = read.table(text = system(paste("vep -i snps.txt --gff /pre_resistance/data/reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.37.sorted.gff3.gz --fasta /pre_resistance/data/reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa --regulatory -o stdout --force_overwrite | sed \"/##/d\" | tr -d \"#\" | egrep -v \"downstream_gene_variant|upstream_gene_variant\"'", sep = ''),intern=TRUE), header = T, stringsAsFactors = F)
snpEff2 = snpEff[,c('Uploaded_variation', 'Consequence', 'Protein_position', 'Amino_acids', 'Codons')]
colnames(snpEff2) = c('pos', 'effect', 'aa_pos', 'aa_change', 'codon_change')


pb <- utils::txtProgressBar(min = 1, max = length(colnames(l2_dna_snps)), style = 3)
for (p in colnames(l2_dna_snps)){
  utils::setTxtProgressBar(pb, which(colnames(l2_dna_snps)==p))

  pos = as.integer(p)
  gene_line = mtb_genes[mtb_genes$gene_st <= pos & mtb_genes$gene_end >= pos,]
  gene = as.character(gene_line$gene)
  st = ifelse(length(gene) == 0, 'intergenic', gene_line$gene_st)
  end = ifelse(length(gene) == 0, 'intergenic', gene_line$gene_end)
  ref = ref_genome$Chromosome[pos]
  alts = unique(l2_dna_snps[,p])[!unique(l2_dna_snps[,p]) %in% ref]
  effect = snpEff2[snpEff2$pos == p,'effect']
  
  if (length(effect) == 0){
    l2_genes = cbind(ifelse(l2_dna_snps[,p] == ref, 0, 1), l2_genes)
    colnames(l2_genes)[1] = paste(p, paste(toupper(alts),collapse = '_'),sep = '.')
  } else if (all(effect == "synonymous_variant")){
    next()
  } else {
    l2_genes[names(l2_dna_snps[,p])[l2_dna_snps[,p] != ref],gene] = 1
  }
}




pos = sort(as.integer(colnames(l4_dna_snps)))
snp_df = c()
for (p in pos){
  p = as.integer(p)
  
  ref = ref_genome$Chromosome[p]
  alts = sort(table(l4_dna_snps[,as.character(p)]), decreasing = T)
  alts = names(alts[-which(names(alts) == ref)][1])
  
  snp_df = rbind(snp_df, cbind('Chromosome', p, p, toupper(ref), toupper(alts), '+'))
}

write.table(snp_df, '/pre_resistance/data/snps.txt', col.names = F, row.names = F, quote = F, sep = ' ')

snpEff = read.table(text = system(paste("vep -i snps.txt --gff /pre_resistance/data/reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.37.sorted.gff3.gz --fasta /pre_resistance/data/reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa --regulatory -o stdout --force_overwrite | sed \"/##/d\" | tr -d \"#\" | egrep -v \"downstream_gene_variant|upstream_gene_variant\"'", sep = ''),intern=TRUE), header = T, stringsAsFactors = F)
snpEff2 = snpEff[,c('Uploaded_variation', 'Consequence', 'Protein_position', 'Amino_acids', 'Codons')]
colnames(snpEff2) = c('pos', 'effect', 'aa_pos', 'aa_change', 'codon_change')


pb <- utils::txtProgressBar(min = 1, max = length(colnames(l4_dna_snps)), style = 3)
for (p in colnames(l4_dna_snps)){
  utils::setTxtProgressBar(pb, which(colnames(l4_dna_snps)==p))
  # if (length(unique(l4_dna_snps[,p]))==1){next()}
  pos = as.integer(p)
  gene_line = mtb_genes[mtb_genes$gene_st <= pos & mtb_genes$gene_end >= pos,]
  gene = as.character(gene_line$gene)

  ref = ref_genome$Chromosome[pos]
  alts = unique(l4_dna_snps[,p])[!unique(l4_dna_snps[,p]) %in% ref]
  effect = snpEff2[snpEff2$pos == p,'effect']
  
  if (length(effect) == 0){
    l4_genes = cbind(ifelse(l4_dna_snps[,p] == ref, 0, 1), l4_genes)
    colnames(l4_genes)[1] = paste(p, paste(toupper(alts),collapse = '_'),sep = '.')
  } else if (all(effect == "synonymous_variant")){
    next()
  } else {
    l4_genes[names(l4_dna_snps[,p])[l4_dna_snps[,p] != ref],gene] = 1
  }
}



save(l4_matrix,l2_matrix,
     l4_pos, l2_pos, 
     l4_freq, l2_freq, 
     l2_seqs, l4_seqs, 
     l2_dna, l4_dna,
     l2_snps, l4_snps,
     l4_bin, l2_bin,
     l4_gg.eigen, l2_gg.eigen,
     l4_genes, l2_genes,
     file = paste(dir_prefix, 'l2_l4_data.gwas.RData', sep=''))
 
