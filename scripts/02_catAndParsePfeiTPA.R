# 20201228 alorenzetti

# description ####
# this script will concatenate
# and parse manually downloaded files
# of Pfeiffer et al 2019
# third party annotation

# concatenating files ####
system2("cat", "data/chr.gff3 data/pnrc100.gff3 data/pnrc200.gff3",
        stdout = "data/pfeiTPA.gff3.tmp")

# parsing ####
pfeiAnnot = rtracklayer::import("data/pfeiTPA.gff3.tmp", format = "gff3") %>% 
  as_tibble() %>% 
  mutate(seqnames = case_when(seqnames == "BK010829.1" ~ "NC_002607.1",
                              seqnames == "BK010830.1" ~ "NC_001869.1",
                              seqnames == "BK010831.1" ~ "NC_002608.1",
                              TRUE ~ as.character(seqnames))) %>% 
  filter(type != "region")

# creating annotation objx for different types of features
# pfeiAnnot$type %>% 
#   unlist(use.names = F) %>% 
#   as.character() %>% 
#   unique()
pfeiStableRNA = list()
pfeiStableRNA[["annot"]] = pfeiAnnot %>% 
  filter(type == "rRNA" | type == "tRNA" | type == "SRP_RNA" | type == "RNase_P_RNA")
pfeiStableRNA[["gr"]] = GRanges(seqnames = pfeiStableRNA[["annot"]]$seqnames,
                                ranges = IRanges(start = pfeiStableRNA[["annot"]]$start,
                                        end = pfeiStableRNA[["annot"]]$end),
                                strand = pfeiStableRNA[["annot"]]$strand)
pfeiStableRNA[["gr"]]$locus_tag = pfeiStableRNA[["annot"]]$locus_tag

pfeiCDS = list()
pfeiCDS[["annot"]] = pfeiAnnot %>% 
  filter(type == "CDS")
pfeiCDS[["gr"]] = GRanges(seqnames = pfeiCDS[["annot"]]$seqnames,
                                ranges = IRanges(start = pfeiCDS[["annot"]]$start,
                                                 end = pfeiCDS[["annot"]]$end),
                                strand = pfeiCDS[["annot"]]$strand)
pfeiCDS[["gr"]]$locus_tag = pfeiCDS[["annot"]]$locus_tag

# removing tmp files
system2("rm", "data/*.tmp")
