
##################################################################################################
# GWAS code inspired by https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6192-1 #
##################################################################################################

library(ape)
library(phangorn)
library(BactDating)
library(data.table)
library(seqinr)


snpEff = function(pos, dna, ref_genome){
  pos = sort(as.integer(pos))
  snp_df = c()
  for (p in pos){
    p = as.integer(p)
    
    ref = ref_genome$Chromosome[p]
    alts = sort(table(dna[,as.character(p)]), decreasing = T)
    alts = names(alts[-which(names(alts) == ref)][1])
    
    # alts = unique(dna[,as.character(p)])[!unique(dna[,as.character(p)]) %in% ref]
    
    # effect = syn_nonNyn(pos, st, end, ref, alts[1], ref_genome)
    snp_df = rbind(snp_df, cbind('Chromosome', p, p, toupper(ref), toupper(alts), '+'))
  }

  write.table(snp_df, '/vep/snps.txt', col.names = F, row.names = F, quote = F, sep = ' ')

  snpEff = read.table(text = system(paste("vep -i /vep/snps.txt --gff /reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.37.sorted.gff3.gz --fasta /reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa --regulatory -o stdout --force_overwrite | sed \"/##/d\" | tr -d \"#\" | egrep -v \"downstream_gene_variant|upstream_gene_variant\"'", sep = ''),intern=TRUE), header = T, stringsAsFactors = F)
  snpEff2 = snpEff[,c('Uploaded_variation', 'Consequence', 'Protein_position', 'Amino_acids', 'Codons')]
  colnames(snpEff2) = c('pos', 'effect', 'aa_pos', 'aa_change', 'codon_change')
  return(snpEff2)  
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




load(paste(dir_prefix, 'l2_l4_data.gwas.assembly.RData', sep=''))



##############################
#######   SNP EFFECT   #######
##############################

l4_snpEff = snpEff(colnames(l4_dna), l4_dna, ref_genome)

l4_snpEff2 = l4_snpEff %>% 
  group_by(pos) %>% 
  filter(!(effect=='synonymous_variant' & n() > 1))

l4_snpEff2 = l4_snpEff2[!duplicated(l4_snpEff2$pos),]

save(l4_snpEff, l4_snpEff2, file = paste(dir_prefix, 'snpEffects.RData', sep = ''))

load(file = paste(dir_prefix, 'snpEffects.RData', sep = ''))


###### RUN THE GENE-BASED GWAS #######

l4_genes2 = l4_genes[,as.character(mtb_genes$gene)]
l4_genes2 = l4_genes2[,apply(l4_genes2, 2, function(x) length(unique(x)) > 1)]
l4_genes2 = l4_genes2[,apply(l4_genes2, 2, function(x) table(x)['1']/length(x) >= 0.01)]

l2_genes2 = l2_genes[,as.character(mtb_genes$gene)]
l2_genes2 = l2_genes2[,apply(l2_genes2, 2, function(x) length(unique(x)) > 1)]
l2_genes2 = l2_genes2[,apply(l2_genes2, 2, function(x) table(x)['1']/length(x) >= 0.01)]


l2_gwas = do.call("rbind", pbapply::pbapply(l2_genes2, 2, function(x){
  res_phen = cbind(l2_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l2_phen_res$node)],
                   l2_gg.eigen$vectors[as.character(l2_phen_res$node),1:50])
  
  colnames(res_phen)[4] = 'base'
  if (nrow(res_phen[res_phen$base == 1 & res_phen$status == 1,]) < 5){return(NULL)}
  
  x = res_phen[, 4:10]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
  return(as.data.frame(d))
}))


l2_gwas$gene = rownames(l2_gwas)
l2_gwas = l2_gwas[!is.na(l2_gwas$pval),]

l2_gwas$pos = sapply(l2_gwas$gene, function(x) mtb_genes$gene_st[mtb_genes$gene == x])
l2_gwas = l2_gwas[order(l2_gwas$pval),]


l4_gwas = do.call("rbind", pbapply::pbapply(l4_genes2, 2, function(x){
  res_phen = cbind(l4_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l4_phen_res$node)],
                   l4_gg.eigen$vectors[as.character(l4_phen_res$node),1:50])
  
  colnames(res_phen)[4] = 'base'
  if (nrow(res_phen[res_phen$base == 1 & res_phen$status == 1,]) < 5){return(NULL)}

  x = res_phen[, 4:54]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
  return(as.data.frame(d))
}))


l4_gwas$gene = rownames(l4_gwas)
l4_gwas = l4_gwas[!is.na(l4_gwas$pval),]

l4_gwas$pos = sapply(l4_gwas$gene, function(x) mtb_genes$gene_st[mtb_genes$gene == x])
l4_gwas = l4_gwas[order(l4_gwas$pval),]



l2_gwas_notCorr = as.data.frame(do.call(rbind, pbapply::pbapply(l2_genes2, 2, function(x){
  res_phen = cbind(l2_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l2_phen_res$node)],
                   l2_gg.eigen$vectors[as.character(l2_phen_res$node),1:50])
  colnames(res_phen)[4] = 'base'
  x = res_phen[, 4]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
})))

rownames(l2_gwas_notCorr) = colnames(l2_genes2)
l2_gwas_notCorr$gene = colnames(l2_genes2)
l2_gwas_notCorr = l2_gwas_notCorr[!is.na(l2_gwas_notCorr$pval),]

l2_gwas_notCorr$pos = sapply(l2_gwas_notCorr$gene, function(x) mtb_genes$gene_st[mtb_genes$gene == x])
l2_gwas_notCorr = l2_gwas_notCorr[order(l2_gwas_notCorr$pval),]



l4_gwas_notCorr = as.data.frame(do.call(rbind, pbapply::pbapply(l4_genes2, 2, function(x){
  res_phen = cbind(l4_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l4_phen_res$node)],
                   l4_gg.eigen$vectors[as.character(l4_phen_res$node),1:50])
  colnames(res_phen)[4] = 'base'
  x = res_phen[, 4]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
})))

rownames(l4_gwas_notCorr) = colnames(l4_genes2)
l4_gwas_notCorr$gene = colnames(l4_genes2)
l4_gwas_notCorr = l4_gwas_notCorr[!is.na(l4_gwas_notCorr$pval),]

l4_gwas_notCorr$pos = sapply(l4_gwas_notCorr$gene, function(x) mtb_genes$gene_st[mtb_genes$gene == x])
l4_gwas_notCorr = l4_gwas_notCorr[order(l4_gwas_notCorr$pval),]


##############################
####### SNP-Based GWAS #######
##############################

l4_snps2 = l4_bin[,apply(l4_bin, 2, function(x) length(unique(x)) > 1)]
l4_snps2 = l4_snps2[,apply(l4_snps2, 2, function(x) min(table(x))/length(x) >= 0.01)]


l2_snps2 = l2_bin[,apply(l2_bin, 2, function(x) length(unique(x)) > 1)]
l2_snps2 = l2_snps2[,apply(l2_snps2, 2, function(x) min(table(x))/length(x) >= 0.01)]




l2_snp_gwas_notCorr = as.data.frame(do.call(rbind, pbapply(l2_snps2, 2, function(x){
  res_phen = cbind(l2_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l2_phen_res$node)],
                   l2_gg.eigen$vectors[as.character(l2_phen_res$node),1:50])
  colnames(res_phen)[4] = 'base'
  x = res_phen[, 4]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
})))

rownames(l2_snp_gwas_notCorr) = colnames(l2_snps2)
l2_snp_gwas_notCorr$id = rownames(l2_snp_gwas_notCorr)
l2_snp_gwas_notCorr = l2_snp_gwas_notCorr[!is.na(l2_snp_gwas_notCorr$pval),]

l2_snp_gwas_notCorr$pos = sapply(rownames(l2_snp_gwas_notCorr), function(x) strsplit(x, '[.]')[[1]][1])
l2_snp_gwas_notCorr$snp = sapply(as.integer(l2_snp_gwas_notCorr$pos), function(x) genome_idx$snp[genome_idx$pos == x])
l2_snp_gwas_notCorr$gene = sapply(as.integer(l2_snp_gwas_notCorr$pos), function(x) paste(mtb_genes$gene[mtb_genes$gene_st <= x & mtb_genes$gene_end >= x], collapse = ';'))

l2_snp_gwas_notCorr = merge(x = l2_snp_gwas_notCorr, y = l2_snpEff2, by = 'pos', all.x = T)
l2_snp_gwas_notCorr = l2_snp_gwas_notCorr[order(l2_snp_gwas_notCorr$pval),]


l4_snp_gwas_notCorr = as.data.frame(do.call(rbind, pbapply(l4_snps2, 2, function(x){
  res_phen = cbind(l4_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l4_phen_res$node)],
                   l4_gg.eigen$vectors[as.character(l4_phen_res$node),1:50])
  colnames(res_phen)[4] = 'base'
  x = res_phen[, 4]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
})))

rownames(l4_snp_gwas_notCorr) = colnames(l4_snps2)
l4_snp_gwas_notCorr$id = rownames(l4_snp_gwas_notCorr)
l4_snp_gwas_notCorr = l4_snp_gwas_notCorr[!is.na(l4_snp_gwas_notCorr$pval),]

l4_snp_gwas_notCorr$pos = sapply(rownames(l4_snp_gwas_notCorr), function(x) strsplit(x, '[.]')[[1]][1])
l4_snp_gwas_notCorr$snp = sapply(as.integer(l4_snp_gwas_notCorr$pos), function(x) genome_idx$snp[genome_idx$pos == x])
l4_snp_gwas_notCorr$gene = sapply(as.integer(l4_snp_gwas_notCorr$pos), function(x) paste(mtb_genes$gene[mtb_genes$gene_st <= x & mtb_genes$gene_end >= x], collapse = ';'))
l4_snp_gwas_notCorr = l4_snp_gwas_notCorr[sapply(as.character(l4_snp_gwas_notCorr$pos), function (x) min(table(l4_matrix[,x]))/nrow(l4_matrix) <= 0.3),]

l4_snp_gwas_notCorr = merge(x = l4_snp_gwas_notCorr, y = l4_snpEff2, by = 'pos', all.x = T)
l4_snp_gwas_notCorr = l4_snp_gwas_notCorr[order(l4_snp_gwas_notCorr$pval),]




l2_snp_gwas = as.data.frame(do.call(rbind, pbapply(l2_snps2, 2, function(x){
  res_phen = cbind(l2_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l2_phen_res$node)],
                   l2_gg.eigen$vectors[as.character(l2_phen_res$node),1:50])
  colnames(res_phen)[4] = 'base'
  x = res_phen[, 4:14]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
})))


rownames(l2_snp_gwas) = colnames(l2_snps2)
l2_snp_gwas$id = rownames(l2_snp_gwas)
l2_snp_gwas = l2_snp_gwas[!is.na(l2_snp_gwas$pval),]

l2_snp_gwas$pos = sapply(rownames(l2_snp_gwas), function(x) strsplit(x, '[.]')[[1]][1])
l2_snp_gwas$snp = sapply(as.integer(l2_snp_gwas$pos), function(x) genome_idx$snp[genome_idx$pos == x])
l2_snp_gwas$gene = sapply(as.integer(l2_snp_gwas$pos), function(x) paste(mtb_genes$gene[mtb_genes$gene_st <= x & mtb_genes$gene_end >= x], collapse = ';'))

l2_snp_gwas = merge(x = l2_snp_gwas, y = l2_snpEff2, by = 'pos', all.x = T)
l2_snp_gwas = l2_snp_gwas[order(l2_snp_gwas$pval),]




l4_snp_gwas = as.data.frame(do.call(rbind, pbapply(l4_snps2, 2, function(x){
  res_phen = cbind(l4_phen_res[,c('node', 'status', 'time')],
                   x[as.character(l4_phen_res$node)],
                   l4_gg.eigen$vectors[as.character(l4_phen_res$node),1:50])
  colnames(res_phen)[4] = 'base'
  x = res_phen[, 4:44]
  x = as.matrix(x)
  y = with(res_phen, Surv(time, status))
  agFit = coxph.fit(x, y, strata = NULL, init = NULL, control = coxph.control(iter.max = 100), method = 'efron', resid = FALSE)
  dVec = c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))
  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
})))


rownames(l4_snp_gwas) = colnames(l4_snps2)
l4_snp_gwas$id = rownames(l4_snp_gwas)
l4_snp_gwas = l4_snp_gwas[!is.na(l4_snp_gwas$pval),]

l4_snp_gwas$pos = sapply(rownames(l4_snp_gwas), function(x) strsplit(x, '[.]')[[1]][1])
l4_snp_gwas$snp = sapply(as.integer(l4_snp_gwas$pos), function(x) genome_idx$snp[genome_idx$pos == x])
l4_snp_gwas$gene = sapply(as.integer(l4_snp_gwas$pos), function(x) paste(mtb_genes$gene[mtb_genes$gene_st <= x & mtb_genes$gene_end >= x], collapse = ';'))

l4_snp_gwas = merge(x = l4_snp_gwas, y = l4_snpEff2, by = 'pos', all.x = T)
l4_snp_gwas = l4_snp_gwas[order(l4_snp_gwas$pval),]

######################
### PLOT MANHATTAN ###
######################

res_gwas = l4_snp_gwas
res_gwas$gene[res_gwas$gene == ''] = 'Intergenic'
res_gwas = res_gwas[order(res_gwas$pval),]
res_gwas$pos = as.integer(as.character(res_gwas$pos))
res_gwas$pval = as.numeric(as.character(res_gwas$pval))

rownames(l4_snp_gwas_notCorr) = l4_snp_gwas_notCorr$id
res_gwas$haz = exp(as.numeric(as.vector(l4_snp_gwas_notCorr[res_gwas$id, 'beta'])))



pdf("/l4.gwas.pdf", useDingbats=F,width = 10, height = 8)
range.i <- range(res_gwas$haz)
extr.min <- min(abs(range.i))
extr.max  <- max(abs(range.i))

z <- seq(extr.min, extr.max, length = 101)

cols <- c(colorRampPalette(c("blue", "white"))(8), 
          colorRampPalette(c("yellow", "red"))(length(z) - 9)
)

cols <- c(colorRampPalette(rev(brewer.pal(n = 11, name="RdBu"))[1:6])(length(which(z < 1))), 
          colorRampPalette(rev(brewer.pal(n = 11, name="RdBu"))[8:11])(length(z) - length(which(z < 1)) + 1)
)


par(mar=c(0,0,0,0), omi=(c(0,0,0,0)+0.2))
layout(matrix(c(1,2), nrow=2, byrow=F), heights=c(0.2, 0.8), widths=c(0.2,0.8))



scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}


par(mar=c(5,1,1,30))
image(z = matrix(z), col = cols, 
      xaxt = "n", yaxt = "n", axes=F)

box(lwd=1)
lv <- pretty(z, n = 10)
xv <- scale01(as.numeric(lv),extr.min, extr.max)
axis(1, at = xv, labels = lv, las=1, font=1, lwd = 1, cex=0.8, cex.axis=1, tck= -0.2, padj=-1)
mtext(1, text = 'Hazard ratio', line = 1.5, font = 2)


ii <- cut(res_gwas$haz, breaks = z, 
          include.lowest = TRUE)

haz_cols = cols[ii]

par(mar=c(6,6,4,4))

log.pval <- -log10(res_gwas$pval)
p_thresh = 0.05 /nrow(res_gwas)
plot(res_gwas$pos, log.pval,
     col = haz_cols,
     pch = 19,
     cex = 0.9,
     ylim=c(0, max(log.pval)),
     main="",xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')

box(lwd = 1.2)
axis(1, cex.axis = 1.3, labels = (axis(1, cex.axis = 1.3, labels = F, tick = F))/1e6, at = (axis(1, cex.axis = 1.3, labels = F, tick = F)))
axis(2, las = 2, cex.axis = 1.3)
mtext('Genome Position (Mb)', 1, 3, cex = 1.3, font = 2)
mtext(expression(bold(-log[10](p-value))), 2, 3, cex = 1.3, font = 2)


thresh <- -log10(p_thresh)
abline(h=thresh, col = "red", lwd=2)

library('basicPlotteR')

sig_n = length(res_gwas$pval[res_gwas$pval < p_thresh])
basicPlotteR::addTextLabels(res_gwas$pos[1:sig_n],
                            log.pval[1:sig_n],
                            res_gwas$gene[1:sig_n],
                            cex.label = 1.2, cex.pt = 1.3,
                            col.background=rgb(0,0,0, 0.75),
                            col.label="white")

dev.off()




#####################
###    QQ Plot    ###
#####################


source('/genomic_inflation_lambda.R')

l4_snp_gwas$pval = as.numeric(as.vector(l4_snp_gwas$pval))
l4_snp_gwas_notCorr$pval = as.numeric(as.vector(l4_snp_gwas_notCorr$pval))

scaleBreaksY = 10^(seq(-39, 0, 3))
scaleLabelsY = as.character(-log10(scaleBreaksY))

scaleBreaksX = 10^(seq(-39, 0, 1))
scaleLabelsX = as.character(-log10(scaleBreaksX))


transNegLog10 = scales::trans_new('neglog10', function(x) -log10(x), function(x) 10^(-x))



pvals = data.frame(pvals = c(l4_snp_gwas_notCorr$pval,
                             l4_snp_gwas$pval),
                   group = c(rep('Raw', nrow(l4_snp_gwas_notCorr)),
                             rep('Corrected', nrow(l4_snp_gwas)))
)
pvals$pvals = as.numeric(as.vector(pvals$pvals))
pvals$group = factor(pvals$group, levels = c('Raw', 'Corrected'))

lambda_corr = round(estlambda(pvals$pvals[pvals$group == 'Corrected'])$estimate, 2)
lambda_raw = round(estlambda(pvals$pvals[pvals$group == 'Raw'])$estimate, 2)


lambda_plot = ggplot(pvals, aes(sample = pvals)) + theme_classic() +
  geom_abline(slope = 1, intercept = 0, color = 'tomato3', size = 1.1) +
  stat_qq(size = 1.3, distribution = stats::qunif, aes(color = group)) +
  labs(
    x = expression(bold(-log[10](p)~expected)),
    y = expression(bold(-log[10](p)~observed))) +
  scale_x_continuous(trans = transNegLog10, breaks = scaleBreaksX, labels = scaleLabelsX) +
  scale_y_continuous(trans = transNegLog10, breaks = scaleBreaksY, labels = scaleLabelsY) +
  scale_color_manual(values=c("lightgray", "black")) +
  theme(legend.title = element_blank(),legend.text = element_text(size=17), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.y = element_text(size = 20, face = 2),
        axis.text.x = element_text(size = 20, face = 2),
        axis.title.y = element_text(size=22, face = 2, margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.title.x = element_text(size=22, face = 2, margin = margin(t = 20, r = 20, b = 0, l = 0))) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  annotate("text", cex = 6, x = 0.0005, y=0.5, label = bquote(lambda ~ paste('=')  ~ .(lambda_corr)), color = "black") +
  annotate("text", cex = 6, x = 0.0005, y=0.1, label = bquote(lambda ~ paste('=')  ~ .(lambda_raw)), color = "grey")




pdf("/Users/arturo/Documents/PhD/TB/pre_resistance/figures/supp7_popCorrection.pdf", useDingbats=F,width = 7, height = 5)
lambda_plot
dev.off()
