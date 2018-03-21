library("rentrez")
library("Biostrings")
library("dplyr")
library("data.table")

#################
###DESCRIPTION###
#################

#Aim: to find upstream gene regulatory regions to feed into CLARE software. 
#Takes an input set of NCBI gene identifiers, and their species genome downloaded from NCBI.
#then takes the CRE as the 1000bp immediately upstream from the first mRNA start site.  
#TODO: get Tom to check it! 


#############
###GLOBALS###
#############

DEGs <- "Honours_data/DEGs AvsW.txt"
apis_genome <- "Honours_data/GCF_000002195.4_Amel_4.5_genomic.fna"
file <- "Honours_data/QWA.csv"

#open fasta genome file 
scaffolds <- readDNAStringSet(apis_genome, format = "fasta")

# #open DEGs file
# tb <- scan(DEGs)
# gene_ids <- c(tb)

#open and filter RNAseq data file 
f <- fread(file)

gene_ids <- f[(WvA_FDR <= 0.05) & (QvA_FDR <= 0.05), `Gene ID`]



#variables for scaffold names
tmp <- strsplit(names(scaffolds), " ")
formatted_scaffold_names <- NULL

#rename scaffolds to same thing as NCBI- to uppercase to catch chromosome/Chromosome typo
#BUG: This does not deal with the GroupUn scaffolds correctly. 
#Am ignoring for now, as I'm assuming no DEGs will fall in these scaffolds. 
i <- 1
while (i <= length(tmp))
{
  val <- tmp[[i]][8]
  tmp_name <- capture.output(cat("Chromosome ", substr(val, 1, nchar(val)-1), " Reference Amel_4.5 Primary Assembly", sep=""))
  formatted_scaffold_names <- c(formatted_scaffold_names, tmp_name)
  i <- i + 1
  
}

names(scaffolds) <- toupper(formatted_scaffold_names)

#get info about each gene from NCBI
n <- 1
CDS_output <- DNAStringSet()
CDS_output_5000 <- DNAStringSet()

while (n <= length(gene_ids)){     #START while loop 
nuc_links <- entrez_fetch(db='gene', id=gene_ids[n], rettype = "xml")
nuc_lists <- XML::xmlToList(nuc_links)

#find start, strand, first RNA position of genes
start <- nuc_lists$Entrezgene$Entrezgene_locus$`Gene-commentary`$`Gene-commentary_seqs`$`Seq-loc`$`Seq-loc_int`$`Seq-interval`$`Seq-interval_from`
stop <- nuc_lists$Entrezgene$Entrezgene_locus$`Gene-commentary`$`Gene-commentary_seqs`$`Seq-loc`$`Seq-loc_int`$`Seq-interval`$`Seq-interval_to`
strand <- nuc_lists$Entrezgene$Entrezgene_locus$`Gene-commentary`$`Gene-commentary_seqs`$`Seq-loc`$`Seq-loc_int`$`Seq-interval`$`Seq-interval_strand`
RNA1 <- nuc_lists$Entrezgene$Entrezgene_locus$`Gene-commentary`$`Gene-commentary_products`$`Gene-commentary`$`Gene-commentary_genomic-coords`$`Seq-loc`$`Seq-loc_mix`$`Seq-loc-mix`$`Seq-loc`$`Seq-loc_int`$`Seq-interval`$`Seq-interval_from`
scaffold_name <- toupper(nuc_lists$Entrezgene$Entrezgene_locus$`Gene-commentary`$`Gene-commentary_label`)

# #obtain fasta of region between first RNA start and start of gene
# 
# if (strand == "minus"){
#   print ("minus")
#   fasta <- reverseComplement(scaffolds[scaffold_name][[1]][RNA1:stop])
# } else if (strand == "plus") {
#   print ("plus")
#   fasta <- scaffolds[scaffold_name][[1]][start:RNA1]
# } else {
#   print ("Error: geneID " + " " + "has no strand! ")
# }



# #grab whole gene
# if (strand == "minus"){
#   fasta <- reverseComplement(scaffolds[scaffold_name][[1]][start:stop])
# } else if (strand == "plus") {
#   fasta <- scaffolds[scaffold_name][[1]][stop:start]
# } else {
#   print ("Error: geneID " + " " + "has no strand! ")
# }

#catch errors in scaffold name
try({

#grab 1000bp upstream from first RNA start 
if (strand == "minus"){
  fasta <- reverseComplement(scaffolds[scaffold_name][[1]][RNA1:(as.integer(RNA1)+1000)])
  fasta_5000 <- reverseComplement(scaffolds[scaffold_name][[1]][RNA1:(as.integer(RNA1)+5000)])
} else if (strand == "plus") {
  fasta <- scaffolds[scaffold_name][[1]][(as.integer(RNA1)-1000):RNA1]
  fasta5000 <- scaffolds[scaffold_name][[1]][(as.integer(RNA1)-5000):RNA1]
} else {
  print ("Error: geneID " + " " + "has no strand! ")
}
  CDS_output <- append(CDS_output, DNAStringSet(fasta))
  CDS_output_5000 <- append(CDS_output_5000, DNAStringSet(fasta_5000))
  }, silent = FALSE, outFile = "Honours_data/scaffold_ass_errors.txt")  



#sending request every 5 seconds to be nice to NCBI. 
Sys.sleep(runif(1, min=0, max=.8))
Sys.sleep(runif(1, min=0, max=.8))
Sys.sleep(runif(1, min=0, max=.8))
Sys.sleep(runif(1, min=0, max=.8))
Sys.sleep(runif(1, min=0, max=.8))
print (n)

#increment n to avoid infinite loop
n <- n + 1

}# END of while loop through gene_ids 

#write to file!
writeXStringSet(CDS_output, "Honours_data/output_fasta_2", format = "fasta")
writeXStringSet(CDS_output_5000, "Honours_data/output_fasta_5000", format = "fasta")

