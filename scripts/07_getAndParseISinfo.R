#20210617

# description ####
# this script will get and
# parse insertion sequence
# information. this does rely
# on Table 1 of my thesis

# starting ####
# copying data from Table 1
isname = c("Hasal_ISNpe8",
           "HsIRS06",
           "HsIRS12",
           "HsIRS44",
           "HsIRS45",
           "ISH1",
           "ISH2",
           "ISH3B",
           "ISH3C",
           "ISH3D",
           "ISH4",
           "ISH6",
           "ISH7",
           "ISH8A",
           "ISH8B",
           "ISH8C",
           "ISH8D",
           "ISH8E",
           "ISH9",
           "ISH10",
           "ISH11",
           "ISH12",
           "ISH22",
           "ISH29",
           "ISH30",
           "ISH32",
           "ISH34",
           "ISH35",
           "ISH37",
           "ISH38",
           "ISH39",
           "ISH40")

isfamily = c("IS66",
             "IS200/IS605",
             "IS630",
             "IS6",
             "IS200/IS605",
             "IS5",
             "IS4",
             "ISH3",
             "ISH3",
             "ISH3",
             "IS1595",
             "ISH6",
             "ISNCY",
             "IS4",
             "IS4",
             "IS4",
             "IS4",
             "IS4",
             "IS5",
             "IS66",
             "IS5",
             "IS200/IS605",
             "IS200/IS605",
             "IS6",
             "IS4",
             "IS4",
             "IS200/IS605",
             "IS200/IS605",
             "IS200/IS605",
             "IS200/IS605",
             "IS200/IS605",
             "IS200/IS605")

isubgroup = c("ISBst12",
              "IS605",
              NA,
              NA,
              "IS1341",
              "ISH1",
              NA,
              NA,
              NA,
              NA,
              "ISH4",
              NA,
              NA,
              "ISH8",
              "ISH8",
              "ISH8",
              "ISH8",
              "ISH8",
              "ISH1",
              "ISBst12",
              NA,
              "IS605",
              "IS605",
              NA,
              "ISH8",
              "ISH8",
              "IS605",
              "IS605",
              "IS1341",
              "IS1341",
              "IS1341",
              "IS605")

isclass = tibble(ISName = isname,
                 ISFamily = isfamily,
                 ISSubgroup = isubgroup)

# getting the 32 non-redundant potentially functional
# insertion sequence CDS from dict
isnrcds = dict %>% 
  filter(str_detect(string = product,
                    pattern = "transposase|ISH")) %>% 
  filter(str_detect(string = product,
                    pattern = "nonfunc",
                    negate = T)) %>% 
  pull(representative)

# getting the names of IS harboring those CDS based
# on Pfeiffer et al. 2019 annotation
isnrcdsgr = pfeiCDS$gr[pfeiCDS$gr$locus_tag %in% isnrcds]

ishits = GenomicRanges::findOverlaps(query = isnrcdsgr,
                                     subject = pfeiMob$gr,
                                     ignore.strand = T,
                                     minoverlap = 50) %>% 
  as_tibble()

results = tibble(IS = pfeiMob$gr$locus_tag[ishits$subjectHits],
                 locus_tag = isnrcdsgr$locus_tag[ishits$queryHits])

# merging with classification data
ISinfo = left_join(x = results,
                   y = isclass,
                   by = c("IS" = "ISName"))

ISinfo = ISinfo[mixedorder(ISinfo$IS),]

ISinfo = ISinfo %>% 
  select(representative = locus_tag,
         ISName = IS,
         ISFamily,
         ISSubgroup)

# writing IS info file
write_tsv(x = ISinfo,
          file = "data/is_info.tsv")
