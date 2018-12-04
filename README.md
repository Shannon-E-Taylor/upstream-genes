# upstream-genes
Filters RNAseq data by ovary activation status and asks NCBI for the 10kb upstream of each gene, outputs to file in fasta format. 


NB there is a bug with this- when there are multiple genomes for one species, it just takes gene location from the first one, which may not be what you want.  
