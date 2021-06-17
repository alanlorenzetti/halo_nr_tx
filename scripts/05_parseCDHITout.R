# 20201229 alorenzetti

# description ####
# this script will take
# cd-hit output files
# and parse the clusters
# identifying a representative
# finally, files will be written

# parsing ####
# declaring parsing function
parseCDHIT = function(cdhitoutfile){
  # code taken from
  # https://rpubs.com/rmurdoch/cdhit_to_mapping_file
  clstr = read.csv(cdhitoutfile, sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
  
  clstr2 = clstr
  n = nrow(clstr)
  x = 0
  numbers_only = function(x) !grepl("\\D", x)
  for (row in c(1:n)) {
    if (numbers_only(clstr2[row,1]) == TRUE) {
      clstr2[row,1] = x}
    else {NULL}
    x = clstr2[row,1]
  }
  
  clstr.sums = data.frame(dplyr::count(clstr2,V1))
  switch = clstr.sums[1,2]
  clstr3 = cbind(clstr2[1], clstr)
  
  clstr4 = clstr2[-which(clstr2$V2 == ""), ]
  
  clstr5 = clstr4
  clstr5[] = lapply(clstr5, gsub, pattern='>', replacement='')
  clstr5.2 = data.frame(str_split_fixed(clstr5$V2, "aa, ", 2))
  clstr5.3 = data.frame(str_split_fixed(clstr5.2$X2, "... ", 2))
  clstr6 = cbind(clstr5[1],clstr5.2[1],clstr5.3[1:2])
  colnames(clstr6) = c("cluster","aa","locus_tag","stat")
  
  # summarizing clusters
  clusters = clstr6 %>% 
    as_tibble() %>% 
    select(cluster,locus_tag) %>% 
    mutate(cluster = str_replace(cluster, " ", "_")) %>% 
    group_by(cluster) %>% 
    summarise(locus_tag = paste0(locus_tag, collapse = ","))
  
  clusters = clusters[mixedorder(x = clusters$cluster),]
  
  # removing clusters of a single entry coming only from NCBI or custom annot
  clusters = clusters %>% 
    filter(!str_detect(locus_tag, pattern = "^VNG_RS")) %>% 
    filter(str_detect(locus_tag, pattern = "VNG_"))
  
  # a representative sequence will be elected
  # as a simplifying rule, the first entry from locus_tag col
  # will be taken as a representative locus
  clusters = clusters %>% 
    mutate(representative = str_replace(locus_tag, ",.*$", "")) %>% 
    select(cluster, representative, locus_tag)
  
  return(clusters)
}

# running for rna and cds datasets
clusters = list()

clusters[["rna"]] = parseCDHIT("data/stableRNAs/cdhitOut.fa.clstr")
clusters[["cds"]] = parseCDHIT("data/cds/cdhitOut.fa.clstr")

# writing a non-redundant transcriptome fasta file for rnas and cds
nrtx = c(seqs$pfei$rna[names(seqs$pfei$rna) %in% sort(clusters$rna$representative)],
         seqs$pfei$cds[names(seqs$pfei$cds) %in% sort(clusters$cds$representative)])
nrtx = nrtx[mixedorder(names(nrtx))]

writeXStringSet(x = nrtx,
                filepath = "data/Hsalinarum_nrtx.fa",
                format = "fasta")

# writing a dictionary for rnas and cds
dict = left_join(x = bind_rows(clusters[["rna"]], clusters[["cds"]]) %>% 
                   select(-cluster),
                 y = bind_rows(pfeiStableRNA$annot, pfeiCDS$annot) %>% select(locus_tag, product),
                 by = c("representative" = "locus_tag")) %>% 
  select(representative, product, locus_tag)
dict = dict[mixedorder(dict$representative),]

write_tsv(x = dict,
          file = "data/dictionary.tsv",
          col_names = T)

## adendum II
# getting tables of custCDS and pfeiCDS
# custCDSTbl = tibble(name = names(custCDS),
#                     seqCust = as.character(custCDS, use.names = F))
# 
# pfeiCDSAA = seqs[["pfei"]][["cds"]] %>%
#   translate(x = ., genetic.code = arc)
# pfeiCDSTbl = tibble(name = names(pfeiCDSAA),
#                     seqPfei = as.character(pfeiCDSAA, use.names = F))
# pfeiCDSTbl$seqPfei = str_replace(string = pfeiCDSTbl$seqPfei,
#                                  pattern = "\\*$",
#                                  replacement = "")
# 
# dictWithSeqs = dict %>%
#   separate_rows(sep = ",", locus_tag) %>% 
#   left_join(x = ., y = custCDSTbl,
#             by = c("locus_tag" = "name")) %>% 
#   left_join(x = ., y = pfeiCDSTbl,
#             by = c("locus_tag" = "name")) %>% 
#   group_by(representative) %>% 
#   summarise(product = c(product)[1],
#             locus_tag = str_c(locus_tag, collapse = ","),
#             seqCust = na.omit(seqCust) %>% as.character() %>% unique() %>% .[1],
#             seqPfei = na.omit(seqPfei) %>% as.character() %>% unique() %>% .[1]) %>% 
#   separate_rows(sep = ",", locus_tag) %>% 
#   filter(str_detect(string = locus_tag, pattern = "VNG_|OE", negate = T)) %>% 
#   group_by(representative) %>% 
#   summarise(representative = c(locus_tag)[1],
#             product = c(product)[1],
#             locus_tag = str_c(locus_tag, collapse = ","),
#             seqCust = c(seqCust)[1],
#             seqPfei = c(seqPfei)[1])
# 
# write_tsv(x = dictWithSeqs,
#           file = "data/custCDSTbl.tsv",
#           col_names = T)
# 
# write_tsv(x = dictWithSeqs %>% 
#             filter(str_detect(string = product,
#                               pattern = "transposase|ISH2")) %>% 
#             filter(str_detect(string = product,
#                               pattern = "nonfunc",
#                               negate = T)),
#           file = "data/custCDSTbl_onlyTransposase.tsv")
#           