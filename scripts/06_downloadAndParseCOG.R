#20210303

# description ####
# this script will download
# COG 2020 data to perform
# functional caracterization
# of halo non redundant transcriptome

# downloading files directly from COG ftp ####
cog = list()

# description of each file is available at
# ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/Readme.2020-11-25.txt

# cog obj
# it is going to throw a warning, since number of
# cols is different depending on the entry
cog$cog = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv",
                   col_names = c("gene_id",
                                 "refseq_id",
                                 "protein_id",
                                 "protein_length",
                                 "cog_coord_prot",
                                 "cog_length_prot",
                                 "cog_id",
                                 "reserved",
                                 "cog_memb_class",
                                 "PSI_BLAST_bit_score",
                                 "PSI_BLAST_evalue",
                                 "cog_prof_length",
                                 "prot_coord_cog"))

# cog memb class definition:
# 0: footprint covers most of the protein and most of the COG profile;
# 1: footprint covers most of the COG profile and part of the protein;
# 2: footprint covers most of the protein and part of the COG profile;
# 3: partial match on both protein and COG profile)

# cog definitions
cog$def = read_tsv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab",
                   col_names = c("cog_id",
                                 "cog_funcat",
                                 "cog_name",
                                 "gene_symbol",
                                 "functional_pathway",
                                 "pubmed_id",
                                 "pdb_id")) #%>% 
  # mutate(cog_funcat = str_replace_all(string = cog_funcat,
  #                                     pattern = "([A-Z])",
  #                                     replacement = "\\1,"),
  #        cog_funcat = str_replace(string = cog_funcat,
  #                                 pattern = ",$",
  #                                 replacement = "")) %>% 
  # separate_rows(cog_funcat, sep = ",")

# some cogs have two categories, but that makes it
# complex since a gene can have more than one COG
# a simple solution is to keep only the first
# category, and that is going to be applied here
cog$def = cog$def %>% 
  mutate(cog_funcat = str_replace(cog_funcat,
                                  "([A-Z]).*$",
                                  "\\1"))

# cog patt
# it is going to throw a warning, since number of
# cols is different depending on the entry
# cog$patt = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.patt.txt",
#                    col_names = F)
# 
# # cog tax
# cog$tax = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.tax.csv",
#                     col_names = F)

# cog function
cog$fun = read_tsv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab",
                   col_names = c("cog_funcat",
                                 "cog_color_hex",
                                 "cog_category"))

# list of organisms
cog$org = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.org.csv",
                   col_names = c("refseq_id", "name", "tx_id", "taxon")) %>% 
  filter(name == "Halobacterium_salinarum_NRC-1_ATCC_700922")

# parsing ####
# getting Halobacterium salinarum NRC-1
# refseq assembly id
rsid = cog$org$refseq_id

# filtering cog object to include
# only halo
cog$cog = cog$cog %>%
  filter(refseq_id == rsid)

# unifying datasets
# keeping only member class 0 or 1
# according to what have been described above
cog$final = cog$cog %>%
  left_join(x = ., y = cog$def,
            by = c("cog_id" = "cog_id")) %>% 
  left_join(x = ., y = cog$fun,
            by = c("cog_funcat" = "cog_funcat")) %>% 
  filter(cog_memb_class == 0 | cog_memb_class == 1) %>% 
  select(gene_id,
         cog_id,
         cog_funcat,
         cog_name,
         gene_symbol,
         functional_pathway,
         cog_category,
         cog_memb_class) %>% 
  mutate(functional_pathway = case_when(is.na(functional_pathway) ~ "Undefined",
                                        TRUE ~ as.character(functional_pathway)),
         gene_symbol = case_when(is.na(gene_symbol) ~ "Undefined",
                                 TRUE ~ as.character(gene_symbol))) %>% 
  group_by(gene_id) %>% 
  summarise(across(.cols = everything(),
                   .fns = ~ paste0(.x, collapse = "|"))) %>% 
  ungroup()

# removing redundancy in case
# a gene has two different COGs
# of the same category
removRed = function(x){
  noRed = str_split(string = x,
                    pattern = "\\|") %>%
    unlist() %>% 
    unique() %>% 
    paste0(collapse = "|")
  
  return(noRed)
}

# performing operations
cog$final = cog$final %>% 
  rowwise() %>% 
  mutate(across(.cols = everything(),
                .fns = ~ removRed(.x)))

# wrangling dict obj and merging with cog obj
cogdict = dict %>% 
  separate_rows(locus_tag, sep=",") %>%
  left_join(x = .,
            y = cog$final,
            by = c("locus_tag" = "gene_id")) %>% 
  select(representative,
         gene_symbol,
         cog_id,
         cog_name,
         cog_category,
         functional_pathway) %>% 
  drop_na() %>% 
  distinct()

# writing cog file
write_tsv(x = cogdict,
          file = "data/cog.tsv")
