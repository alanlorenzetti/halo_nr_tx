# 20201228 alorenzetti

# description ####
# this script will download 
# Halobacterium salinarum annotation from NCBI RefSeq ftp;
# also it is going to read the custom sequence file 
# derived from Ng et al. 2000;
# furthermore it will download Hsal R1 Genbank annotation

# downloading files from NCBI for NRC-1 ####

# downloading Halobacterium salinarum genome from
# NCBI assembly ref seq https://www.ncbi.nlm.nih.gov/assembly/GCF_000006805.1/
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz",
              destfile = "data/Hsalinarum.fa.gz")
system2(command = "gunzip", args = "-f data/Hsalinarum.fa.gz")

# downloading Halobacterium salinarum annotation from refseq
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.gff.gz",
              destfile = "data/Hsalinarum.gff.gz")
system2(command = "gunzip", args = "-f data/Hsalinarum.gff.gz")

# parsing annotation #
NCBIAnnot = rtracklayer::import("data/Hsalinarum.gff") %>% 
  as_tibble() %>% 
  filter(type != "region")

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

# downloading files from NCBI for R1 ####
# downloading Halobacterium salinarum genome from
# NCBI assembly ref seq https://www.ncbi.nlm.nih.gov/assembly/GCF_000006805.1/
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/069/025/GCF_000069025.1_ASM6902v1/GCF_000069025.1_ASM6902v1_genomic.fna.gz",
              destfile = "data/HsalinarumR1.fa.gz")
system2(command = "gunzip", args = "-f data/HsalinarumR1.fa.gz")

# downloading Halobacterium salinarum annotation from refseq
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/069/025/GCF_000069025.1_ASM6902v1/GCF_000069025.1_ASM6902v1_genomic.gff.gz",
              destfile = "data/HsalinarumR1.gff.gz")
system2(command = "gunzip", args = "-f data/HsalinarumR1.gff.gz")

# parsing annotation #
R1Annot = rtracklayer::import("data/HsalinarumR1.gff") %>% 
  as_tibble() %>% 
  filter(type != "region")

R1StableRNA = list()
R1StableRNA[["annot"]] = R1Annot %>% 
  filter((!is.na(old_locus_tag)) & type == "gene" & (gene_biotype == "tRNA" | gene_biotype == "SRP_RNA" | gene_biotype == "RNase_P_RNA"))
R1StableRNA[["gr"]] = GRanges(seqnames = R1StableRNA[["annot"]]$seqnames,
                                ranges = IRanges(start = R1StableRNA[["annot"]]$start,
                                                 end = R1StableRNA[["annot"]]$end),
                                strand = R1StableRNA[["annot"]]$strand)
R1StableRNA[["gr"]]$locus_tag = R1StableRNA[["annot"]]$old_locus_tag %>% 
  str_replace(pattern = ",.*$", "")

R1CDS = list()
R1CDS[["annot"]] = R1Annot %>% 
  filter((!is.na(old_locus_tag)) & type == "gene" & gene_biotype == "protein_coding")
R1CDS[["gr"]] = GRanges(seqnames = R1CDS[["annot"]]$seqnames,
                          ranges = IRanges(start = R1CDS[["annot"]]$start,
                                           end = R1CDS[["annot"]]$end),
                          strand = R1CDS[["annot"]]$strand)
R1CDS[["gr"]]$locus_tag = R1CDS[["annot"]]$old_locus_tag %>% 
  str_replace(pattern = ",.*$", "") %>% 
  str_replace(pattern = "_", "")

# loading and parsing custom protein sequence file ####
custCDS = readAAStringSet(filepath = "data/Hsalinarum_pi2s.fasta",
                           format = "fasta")

# removing non halo seqs and parsing names
custCDS = custCDS[1:2646]
names(custCDS) = str_replace(names(custCDS), " .*$", "")
