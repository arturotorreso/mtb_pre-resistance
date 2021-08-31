library(ape)
library(phangorn)
library(BactDating)
library(survival)
library(RColorBrewer)
library(survminer)
library(survutils)
library(ggtree)


source ('~/Desktop/PhD/TB/pre_resistance/drug_resistance_phylogeny/scripts/final_scripts/tree_functions.R')
source ('~/Desktop/PhD/TB/pre_resistance/drug_resistance_phylogeny/scripts/final_scripts/survTree_functions.R')


dir_prefix = '/'

######## FILES NEEDED ##########
l2_tree_file = paste(dir_prefix, 'l2.assembly.arc.5e7.0.5.wg1_res.RData', sep = '')
l4_tree_file = paste(dir_prefix, 'l4.assembly.arc.1e7.0.5.wg2_res.RData', sep = '')
l2.anc.seq_file = paste(dir_prefix, 'l2.anc.ml.assemblyTree.RData', sep = '')
l4.anc.seq_file = paste(dir_prefix, 'l4.anc.ml.assemblyTree.RData', sep = '')
genome_idx_file = paste(dir_prefix, 'mtb.snps.dels.wg.withDR.assembly.finalSet.snpSites.masked.idx', sep = '')
tbdr_db_file = paste(dir_prefix, 'tb_profiler_db.complete.good.posStrRef.tsv', sep = '')
dna_file = paste(dir_prefix, 'dna_phyDat.assembly.RData', sep = '')
metadata_file = paste(dir_prefix, 'tb_1999_2016_noDups.metadata.csv', sep='')
###############################


######## FILES NEEDED (SAMARA) ##########
l2_tree_file = paste(dir_prefix, 'l2.samara.assembly.wg.raxml.tree', sep = '')
l4_tree_file = paste(dir_prefix, 'l4.samara.assembly.wg.raxml.tree', sep = '')
l2.anc.seq_file = paste(dir_prefix, 'l2.anc.ml.samara.assembly.RData', sep = '')
l4.anc.seq_file = paste(dir_prefix, 'l4.anc.ml.samara.assembly.RData', sep = '')
genome_idx_file = paste(dir_prefix, 'mtb.samara.snps.dels.wg.withDR.assembly.finalSet.snpSites.masked.idx', sep = '')
tbdr_db_file = paste(dir_prefix, 'tb_profiler_db.complete.good.posStrRef.tsv', sep = '')
dna_file = paste(dir_prefix, 'dna_phydat.samara.assembly.RData', sep = '')
metadata_file = paste(dir_prefix, 'samara_lineages.txt', sep='')
##########################################


###############################
####      Load Data       ####
##############################


load(l2_tree_file)
l2_tree =  reorder.phylo(res$tree, 'cladewise')
l2_tree$CI = res$CI

load(l4_tree_file)
l4_tree = reorder.phylo(res$tree, 'cladewise')
l4_tree$CI = res$CI

load(dna_file)
dna_phydat = dna_phydat[names(dna_phydat[which(names(dna_phydat) %in% c(l4_tree$tip.label, l2_tree$tip.label))])]

metadata = read.csv(metadata_file, stringsAsFactors = F, header = T)
metadata = metadata[!is.na(metadata$id),]

# Change names in metadata so it matchers the phylogeny
metadata$id = paste(metadata$id, '_', metadata$year_sample, sep='')


metadata = metadata[metadata$id %in% c(l4_tree$tip.label, l2_tree$tip.label),]

l2_edges = date_edges(l2_tree)
l4_edges = date_edges(l4_tree)

rif_phen = setNames(metadata$rif_final, metadata$id)
inh_phen = setNames(metadata$inh_final, metadata$id)

l4_phen_rif = rif_phen[l4_tree$tip.label]
l4_phen_inh = rif_phen[l4_tree$tip.label]

l2_phen_rif = rif_phen[l2_tree$tip.label]
l2_phen_inh = inh_phen[l2_tree$tip.label]


#############################################
####        Defining subclades          ####
############################################

beijing = metadata[grepl('lineage2', metadata$snp_lineage),]


lam = metadata[grepl('lineage4.3', metadata$snp_lineage),]
lam3 = metadata[grepl('lineage4.3.2', metadata$snp_lineage),]
lam9 = metadata[grepl('lineage4.3.3', metadata$snp_lineage),]
lam11 = metadata[grepl('lineage4.3.4', metadata$snp_lineage),]
l4.3 =  metadata[metadata$snp_lineage == 'lineage4.3',]
l4.3.1 =  metadata[metadata$snp_lineage == 'lineage4.3.1',]

lam9_node = getMRCA(l4_tree, lam9$id)
lam9 = rbind(lam9, l4.3[l4.3$id %in% extract.clade(l4_tree, lam9_node)$tip.label,])

lam11_node = getMRCA(l4_tree, lam11$id)
lam11_node = l4_edges$V1[l4_edges$V2 == lam11_node]
lam11 = rbind(lam11, l4.3[l4.3$id %in% extract.clade(l4_tree, lam11_node)$tip.label,])


haarlem = metadata[grepl('lineage4.1.2.1', metadata$snp_lineage),]
l4.1.2 =  metadata[metadata$snp_lineage == 'lineage4.1.2',]

haarlem = rbind (haarlem, l4.1.2)
x_type =  metadata[grepl('lineage4.1.1', metadata$snp_lineage),]
l4.1 = metadata[metadata$snp_lineage == 'lineage4.1',]
x_type = rbind(x_type, l4.1)

l4.4 = metadata[grepl('lineage4.4', metadata$snp_lineage),]

tur = metadata[metadata$snp_lineage == 'lineage4.2.2',]




t_type = metadata[metadata$snp_lineage == 'lineage4' |
                    metadata$snp_lineage == 'lineage4.5' |
                    metadata$snp_lineage == 'lineage4.7' |
                    metadata$snp_lineage == 'lineage4.8' |
                    metadata$snp_lineage == 'lineage4.9'
                  ,]


l2_subclades = list(beijing=beijing$id)


l4_subclades = list(lam3=lam3$id,
                    lam9=lam9$id,
                    lam11=lam11$id,
                    haarlem=haarlem$id,
                    x_type=x_type$id,
                    t_type=t_type$id)

l4_subclades = list(lam3=lam3$id,
                    lam9=lam9$id,
                    lam11=lam11$id,
                    haarlem=haarlem$id,
                    x_type=x_type$id,
                    l4.4=l4.4$id,
                    tur=tur$id,
                    t_type=t_type$id)


l2_subclades_nodes = sapply(names(l2_subclades),function(x) getMRCA(l2_tree, l2_subclades[[x]]))
l4_subclades_nodes = sapply(names(l4_subclades),function(x) getMRCA(l4_tree, l4_subclades[[x]]))


l2_subclades_tips = lapply(l2_subclades, function(x)
  sapply(x, function(y) which(l2_tree$tip.label == y)))

l4_subclades_tips = lapply(l4_subclades, function(x)
  sapply(x, function(y) which(l4_tree$tip.label == y)))





ggtree(tree, size=1, aes(color=group), mrsd = '2016-01-01') +
  scale_color_manual(values=c('0'="grey40", col_vector[names(cls)])) +
  theme_tree2() + coord_cartesian(clip = 'off') +
  xlab('Year') +
  theme(legend.position='left',
        plot.margin = margin(0.5, 3, 0.5, 0.5, 'cm')) +
  geom_cladelabel(node=getMRCA(l4_tree, tip=haarlem$id),
                  offset=30,offset.text=30,
                  label="Haarlem", align=T,barsize=2,
                  color=col_vector['haarlem'],fontsize=5) +
  geom_cladelabel(node=getMRCA(l4_tree, tip=lam$id),
                  offset=30,offset.text=30,
                  label="LAM", align=T,barsize=2,
                  color=col_vector['lam11'],fontsize=5) + 
  geom_cladelabel(node=getMRCA(l4_tree, tip=t_type$id),
                  offset=30, offset.text=30,
                  label="T type", align=T,barsize=2,
                  color=col_vector['t_type'],fontsize=5) + 
  geom_cladelabel(node=getMRCA(l4_tree, tip=x_type$id),
                  offset=30,offset.text=30,
                  label="X type", align=T,barsize=2,
                  color=col_vector['x_type'],fontsize=5)


ggtree(tree, size=1, aes(color=group), mrsd = '2016-01-01') +
  scale_color_manual(values=c('0'="grey40", col_vector[names(cls)])) +
  theme_tree2() + coord_cartesian(clip = 'off') +
  xlab('Year') +
  xlim(1250, 2200) + 
  theme(legend.position='left',
        plot.margin = margin(0.5, 3, 0.5, 0.5, 'cm')) +
  geom_cladelabel(node=getMRCA(l2_tree, tip=beijing$id),
                  offset=20,offset.text=20,
                  label="Beijing", align=T,barsize=2,
                  color=col_vector['beijing'],fontsize=5)




## SAMARA ##
lineage2 = metadata[grepl('lineage2', metadata$snp_lineage),]
lineage4 = metadata[grepl('lineage4', metadata$snp_lineage),]

l2_subclades = list(lineage2=lineage2$id)
l4_subclades = list(lineage4=lineage4$id)

l2_subclades_nodes = sapply(names(l2_subclades),function(x) getMRCA(l2_tree, l2_subclades[[x]]))
l4_subclades_nodes = sapply(names(l4_subclades),function(x) getMRCA(l4_tree, l4_subclades[[x]]))


l2_subclades_tips = lapply(l2_subclades, function(x)
  sapply(x, function(y) which(l2_tree$tip.label == y)))

l4_subclades_tips = lapply(l4_subclades, function(x)
  sapply(x, function(y) which(l4_tree$tip.label == y)))




###################################################
#####    Analysis using drug resistant snps  #####
##################################################

######  1. Ancestral state reconstruction ####

load (l2.anc.seq_file)
load (l4.anc.seq_file)


######  2. DR metadata ######

genome_idx = read.table(genome_idx_file,
                        stringsAsFactors = F, header = T)

# TB-Profiler DataBase with the drug associated mutations
tbprof = read.table(tbdr_db_file,
                  stringsAsFactors = F, header = T)

phylo_snps = c(7539, 1674434, 2518919, 1472337)

tbprof = tbprof[apply(tbprof, 1, function(x) !all(strsplit(x['pos'], ';')[[1]] %in% phylo_snps)),]

# Remove indels and focus in SNPs (maybe change in future...)

tbprof = tbprof[nchar(tbprof$ref) == nchar(tbprof$alt),]

# For eficiency, only use those positions of resistant mutations present in our alignment
tbprof2 = tbprof[apply(tbprof, 1, function(x) all(strsplit(x['pos'], ';')[[1]] %in% genome_idx$pos)),]


# Weird way to do it to preserve the order of the pos1;pos2 and pos2;pos1
tbprof2$snp = sapply(tbprof2$pos, function(p)
  paste(genome_idx$snp[unlist(lapply(strsplit(p, ';')[[1]], function(x)
    which(genome_idx$pos %in% x)))], collapse = ';'))


# Just change rownames. A bit of OCD...
row.names(tbprof2) = 1:nrow(tbprof2)


######  3. Time to resistance  ######

# Skip this to see the resulting files from running this functions
l2_t_to_res_snp = t_to_res_snp (
  edges = l2_edges,
  tree = l2_tree,
  subclades = l2_subclades,
  subclades_tips = l2_subclades_tips,
  subclades_nodes = l2_subclades_nodes,
  year_threshold = FALSE,
  bases = c('A', 'C', 'G', 'T', '-'),
  anc.seq = l2.anc.ml,
  drug_db = tbprof2[tbprof2$confidence %in% c('high', 'moderate', 'low'),]
)


l4_t_to_res_snp = t_to_res_snp (
  edges = l4_edges,
  tree = l4_tree,
  subclades = l4_subclades,
  subclades_tips = l4_subclades_tips,
  subclades_nodes = l4_subclades_nodes,
  year_threshold = FALSE,
  bases = c('A', 'C', 'G', 'T', '-'),
  anc.seq = l4.anc.ml,
  drug_db = tbprof2[tbprof2$confidence %in% c('high', 'moderate', 'low'),]
)



l4_t_to_res_snp = l4_t_to_res_snp[l4_t_to_res_snp$lineage %in% 
                                    names(l4_subclades),]
l4_t_to_res_snp = l4_t_to_res_snp[!is.na(l4_t_to_res_snp$lineage),]

l2_t_to_res_snp$lineage = 'beijing'

mtb_t_to_res = rbind (l2_t_to_res_snp, l4_t_to_res_snp)


mtb_t_to_res$lineage = ifelse (mtb_t_to_res$lineage == 'beijing', 'lineage2', 'lineage4')
mtb_t_to_res$lineage = factor(mtb_t_to_res$lineage,
                              levels = c('lineage4', 'lineage2'))
mtb_t_to_res = mtb_t_to_res[(mtb_t_to_res$int_st == 0 & mtb_t_to_res$ext_st == 1) | 
                (mtb_t_to_res$int_st == 0 & mtb_t_to_res$ext_st == 0),]
mtb_t_to_res = mtb_t_to_res[,c('int_year', 'ext_year', 'ext_st', 'lineage')]
mtb_t_to_res$time = mtb_t_to_res$ext_year - mtb_t_to_res$int_year 
colnames(mtb_t_to_res) = c('int_year', 'ext_year', 'status',  'lineage',
                            'time')
mtb_t_to_res$time = as.numeric(as.vector(mtb_t_to_res$time))
mtb_t_to_res = mtb_t_to_res[mtb_t_to_res$int_year >= 1940,]


sfit <- survfit(Surv(time, status)~lineage, data=mtb_t_to_res)
fit.coxph <- coxph(Surv(time, status) ~ lineage,
                   data = mtb_t_to_res)




##### PLOT WITH TWO LINEAGES #####

cols = brewer.pal(12, "Set3")[c(5,4)]
col.levels = c('lineage4', 'lineage2')

pdf ('/plots/l2_vs_l4.pdf', height = 6, width = 7, useDingbats = F)

prettySurvPlot(sfit = sfit,
               fit.coxph = fit.coxph,
               xlevels = c('lineage4', 'lineage2'),
               kap_meier = TRUE,
               haz = TRUE,
               conf.int = TRUE,
               km.ylim = c(0.60,1.01),
               km.xlim = c(0,40),
               km.xlab = 'Time (years)',
               km.ylab = 'Probability of\nremaining susceptible',
               cols = brewer.pal(12, "Set3")[c(5,4)],
               col.levels = c('lineage4', 'lineage2')
               )
dev.off()


##### PLOT WITH ALL LINEAGES #####

mtb_t_to_res = rbind (l2_t_to_res_snp, l4_t_to_res_snp)

mtb_t_to_res = mtb_t_to_res[(mtb_t_to_res$int_st == 0 & mtb_t_to_res$ext_st == 1) | 
                              (mtb_t_to_res$int_st == 0 & mtb_t_to_res$ext_st == 0),]
mtb_t_to_res = mtb_t_to_res[,c('int_year', 'ext_year', 'ext_st', 'lineage')]

mtb_t_to_res$time = mtb_t_to_res$ext_year - mtb_t_to_res$int_year 
colnames(mtb_t_to_res) = c('int_year', 'ext_year', 'status',  'lineage',
                           'time')
mtb_t_to_res$time = as.numeric(as.vector(mtb_t_to_res$time))
mtb_t_to_res = mtb_t_to_res[mtb_t_to_res$int_year >= 1940,]


mtb_t_to_res$lineage = factor(mtb_t_to_res$lineage,
                              levels = c('lam3', 'lam9',
                                         'lam11',
                                         'haarlem', 't_type','x_type',
                                         'beijing'))


sfit <- survfit(Surv(time, status)~lineage, data=mtb_t_to_res)
fit.coxph <- coxph(Surv(time, status) ~ lineage,
                   data = mtb_t_to_res)


pdf ('/plots/allLin.pdf', height = 6, width = 7, useDingbats = F)

cols = brewer.pal(12, "Set3")[c(1,3,4,5,6,10,12)]
col.levels = c("beijing","haarlem","lam11","lam3",
               "lam9","t_type", "x_type")


prettySurvPlot(sfit = sfit,
               fit.coxph = fit.coxph,
               kap_meier = TRUE,
               haz = TRUE,
               conf.int = F,
               km.ylab = 'Probability of\nremaining susceptible',
               km.ylim = c(0.6,1.01),
               km.xlim = c(0,40),
               cols = cols,haz.size = 0.35,km.size = 0.35,risk.size = 0.3,
               col.levels = col.levels
)

dev.off()


# PLOT INH TO MDR

mtb_t_to_res = rbind(l2_t_to_res_snp, l4_t_to_res_snp)
mtb_t_to_res = mtb_t_to_res[mtb_t_to_res$int_year >= 1940,]


mtb_t_to_mdr = list()
drug1 = 'isoniazid'
drug2 = 'rifampicin'

for (i in 1:nrow(mtb_t_to_res)){
  
  int_res_present = strsplit(as.character(mtb_t_to_res$int_drug[i]), ';')[[1]]
  ext_res_present = strsplit(as.character(mtb_t_to_res$ext_drug[i]), ';')[[1]]
  
  int_time = mtb_t_to_res$int_year[i]
  ext_time = mtb_t_to_res$ext_year[i]
  if (drug2 %in% int_res_present) {
    next()
  } else if (all(is.na(int_res_present)) &
             !drug2 %in% ext_res_present) {
    time = ext_time - int_time
    mtb_t_to_mdr[['time']] = c(mtb_t_to_mdr[['time']], time)
    mtb_t_to_mdr[['status']] = c(mtb_t_to_mdr[['status']], 0)
    mtb_t_to_mdr[['background']] = c(mtb_t_to_mdr[['background']], 'sensitive')
  } else if (all(is.na(int_res_present)) &
             drug2 %in% ext_res_present) {
  
    time = ext_time - int_time
    mtb_t_to_mdr[['time']] = c(mtb_t_to_mdr[['time']], time)
    mtb_t_to_mdr[['status']] = c(mtb_t_to_mdr[['status']], 1)
    mtb_t_to_mdr[['background']] = c(mtb_t_to_mdr[['background']], 'sensitive')
    
  } else if (all(int_res_present %in% drug1) &
             !drug2 %in% ext_res_present) {
    time = ext_time - int_time
    mtb_t_to_mdr[['time']] = c(mtb_t_to_mdr[['time']], time)
    mtb_t_to_mdr[['status']] = c(mtb_t_to_mdr[['status']], 0)
    mtb_t_to_mdr[['background']] = c(mtb_t_to_mdr[['background']],
                                     drug1)
    
    
  } else if (all(int_res_present %in% drug1) &
             drug2 %in% ext_res_present) {
    time = ext_time - int_time
    mtb_t_to_mdr[['time']] = c(mtb_t_to_mdr[['time']], time)
    mtb_t_to_mdr[['status']] = c(mtb_t_to_mdr[['status']], 1)
    mtb_t_to_mdr[['background']] = c(mtb_t_to_mdr[['background']],
                                     drug1)
    
  }
  
  
}

mtb_t_to_mdr = as.data.frame(mtb_t_to_mdr)
mtb_t_to_mdr[,c("time")] = as.numeric(as.vector(unlist(mtb_t_to_mdr[,c("time")])))

mtb_t_to_mdr$background = factor(mtb_t_to_mdr$background,
                                 levels = c('sensitive', 'isoniazid'))


sfit <- survfit(Surv(time, status)~background, data=mtb_t_to_mdr)
fit.coxph <- coxph(Surv(time, status)~background, data=mtb_t_to_mdr)


cols = brewer.pal(12, "Set3")[c(5,4)]
col.levels = c('sensitive', 'isoniazid')

pdf ('/plots/inh2mdr.samara.pdf', height = 6, width = 7, useDingbats = F)


prettySurvPlot(sfit = sfit,
               fit.coxph = fit.coxph,
               xlevels = c('sensitive', 'isoniazid'),
               kap_meier = TRUE,
               haz = T,
               risk = T,
               conf.int = TRUE,
               km.ylim = c(0,1.01),
               km.xlim = c(0,40),
               km.xlab = 'Time (years)',
               km.ylab = 'Probability of\nremaining monoresistant',
               cols = brewer.pal(12, "Set3")[c(5,4)],
               col.levels = c('sensitive', 'isoniazid')
)

dev.off()

