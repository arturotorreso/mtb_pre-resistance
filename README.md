# Genomic Signatures of Pre-Resistance in Mycobacterium tuberculosis

Code for Genomic Signatures of Pre-Resistance in Mycobacterium tuberculosis.
https://www.researchsquare.com/article/rs-364747/v1

The files repeat the analysis specified in Materials and Mehods:<br/>
1) Assembly, variant calling and pseudosequence:<br/>
>Main script: pseudoseq_pipeline.sh<br/>
>VCF filtering and annotation: addFT.py<br/>
>Pseudosequence creation: vcf2pseudoseq.py<br/>
2) Phylogenetic inference<br/>
> Tree inference: raxml.sh
> Tree dating: runBactDat.R, run_bactDat.sh
4) Phylogenetic analysis<br/>
>Ancestral sequence reconstruction: anc_seq_recons.R<br/>
>Survival analysis using phylogenetic tree: survTree_functions.R<br/>
5) Genome-wide association study<br/>



