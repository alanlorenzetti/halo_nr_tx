# 20201228 alorenzetti

# description ####
# this script will download 
# Halobacterium salinarum annotation from NCBI RefSeq ftp
# it is going to read the custom sequence file 
# derived from Ng et al. 2000

# downloading files from web resources ####

# downloading Halobacterium salinarum genome from
# NCBI assembly ref seq https://www.ncbi.nlm.nih.gov/assembly/GCF_000006805.1/
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz",
              destfile = "data/Hsalinarum.fa.gz")
system2(command = "gunzip", args = "-f data/Hsalinarum.fa.gz")

# downloading Halobacterium salinarum annotation from refseq
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.gff.gz",
              destfile = "data/Hsalinarum.gff.gz")
system2(command = "gunzip", args = "-f data/Hsalinarum.gff.gz")

# parsing annotation ####
NCBIAnnot = rtracklayer::import("data/Hsalinarum.gff") %>% 
  as_tibble() %>% 
  filter(type != "region")

# creating annotation objx for different types of features
# NCBIAnnot$type %>%
#   unlist(use.names = F) %>%
#   as.character() %>%
#   unique()

ncbiStableRNA = list()
ncbiStableRNA[["annot"]] = NCBIAnnot %>% 
  filter(type == "rRNA" | type == "tRNA" | type == "SRP_RNA" | type == "RNase_P_RNA")
ncbiStableRNA[["gr"]] = GRanges(seqnames = ncbiStableRNA[["annot"]]$seqnames,
                                ranges = IRanges(start = ncbiStableRNA[["annot"]]$start,
                                                 end = ncbiStableRNA[["annot"]]$end),
                                strand = ncbiStableRNA[["annot"]]$strand)
ncbiStableRNA[["gr"]]$locus_tag = ncbiStableRNA[["annot"]]$locus_tag

ncbiCDS = list()
ncbiCDS[["annot"]] = NCBIAnnot %>% 
  filter(type == "CDS")
ncbiCDS[["gr"]] = GRanges(seqnames = ncbiCDS[["annot"]]$seqnames,
                          ranges = IRanges(start = ncbiCDS[["annot"]]$start,
                                           end = ncbiCDS[["annot"]]$end),
                          strand = ncbiCDS[["annot"]]$strand)
ncbiCDS[["gr"]]$locus_tag = ncbiCDS[["annot"]]$locus_tag

# loading and parsing custom protein sequence file
custCDS = readAAStringSet(filepath = "data/Hsalinarum_pi2s.fasta",
                           format = "fasta")

# removing non halo seqs and parsing names
custCDS = custCDS[1:2646]
names(custCDS) = str_replace(names(custCDS), " .*$", "")
