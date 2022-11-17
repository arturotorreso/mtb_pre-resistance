# Genomic Signatures of Pre-Resistance in Mycobacterium tuberculosis

Code for Genomic Signatures of Pre-Resistance in Mycobacterium tuberculosis.
https://www.nature.com/articles/s41467-021-27616-7

The files repeat the analysis specified in Materials and Mehods:<br/>
1) Assembly, variant calling and pseudosequence:<br/>
>Main script: pseudoseq_pipeline.sh<br/>
>VCF filtering and annotation: addFT.py<br/>
>Pseudosequence creation: vcf2pseudoseq.py<br/>
2) Phylogenetic inference<br/>
> Tree inference: raxml.sh<br/>
> Tree dating: runBactDat.R, run_bactDat.sh<br/>
4) Phylogenetic analysis<br/>
>Ancestral sequence reconstruction: anc_seq_recons.R<br/>
>Survival analysis using phylogenetic tree: survTree_functions.R, tree_functions.R, survival_analysis.R<br/>
>TB profiler DB file: tb_profiler_db.complete.good.posStrRef.tsv<br/>
>Genome index file: mtb.snps.dels.wg.withDR.assembly.finalSet.snpSites.masked.idx<br/>
5) Genome-wide association study<br/>
> Alignment preparation: gwas_alignment.R
> GWAS: gwas2_analysis.R

This code has been tested on R version 4.1.1 (2021-08-10).
Dependencies required are specified within the R scripts.

