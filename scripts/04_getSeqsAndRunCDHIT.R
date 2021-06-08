# 20201228 alorenzetti

# description ####
# this script will read the genome file
# and extract sequences prior 
# to the execution of
# similarity clustering by CD-HIT

# loading genome file ####
genome = readDNAStringSet("data/Hsalinarum.fa")
names(genome) = str_replace(names(genome), " .*$", "")

# loading R1 genome file
genomeR1 = readDNAStringSet("data/HsalinarumR1.fa")
names(genomeR1) = str_replace(names(genomeR1), " .*$", "")

# getting seqs ####
arc = getGeneticCode(id_or_name2 = "11")
seqs = list()

seqs[["pfei"]][["rna"]] = BSgenome::getSeq(genome, pfeiStableRNA$gr)
names(seqs[["pfei"]][["rna"]]) = pfeiStableRNA$gr$locus_tag

seqs[["pfei"]][["cds"]] = BSgenome::getSeq(genome, pfeiCDS$gr)
names(seqs[["pfei"]][["cds"]]) = pfeiCDS$gr$locus_tag

seqs[["ncbi"]][["rna"]] = BSgenome::getSeq(genome, ncbiStableRNA$gr)
names(seqs[["ncbi"]][["rna"]]) = ncbiStableRNA$gr$locus_tag

seqs[["ncbi"]][["cds"]] = BSgenome::getSeq(genome, ncbiCDS$gr)
names(seqs[["ncbi"]][["cds"]]) = ncbiCDS$gr$locus_tag

seqs[["R1"]][["rna"]] = BSgenome::getSeq(genomeR1, R1StableRNA$gr)
names(seqs[["R1"]][["rna"]]) = R1StableRNA$gr$locus_tag

seqs[["R1"]][["cds"]] = BSgenome::getSeq(genomeR1, R1CDS$gr)
names(seqs[["R1"]][["cds"]]) = R1CDS$gr$locus_tag

# unifying and saving stable rnas
c(seqs$pfei$rna, seqs$ncbi$rna, seqs$R1$rna) %>% 
  writeXStringSet(filepath = "data/stableRNAs.fa",
                  format = "fasta")

# unifying and saving cds
seqvec1 = c(seqs$pfei$cds, seqs$ncbi$cds) %>% translate(x = ., genetic.code = arc)
seqvec2 = seqs$R1$cds %>% translate(x = ., genetic.code = arc)
c(seqvec1, custCDS, seqvec2) %>% 
  writeXStringSet(filepath = "data/cds.fa",
                  format = "fasta")

# running cd-hit using stable RNA files
if(!dir.exists("data/stableRNAs")){dir.create("data/stableRNAs")}
cdhit = c("/Users/alorenzetti/CompiledPrograms/cdhit-4.8.1/cd-hit",
          "-i data/stableRNAs.fa",
          "-o data/stableRNAs/cdhitOut.fa",
          "-c 0.99",
          "-M 1000",
          "-T 8",
          "-sc 1")

system2(command = cdhit[1],
        args = paste(cdhit[-1], collapse = " "),
        stdout = "data/stableRNAs/cdhit.log",
        stderr = "data/stableRNAs/cdhit.err")

# running cd-hit using cds files
if(!dir.exists("data/cds")){dir.create("data/cds")}
cdhit = c("/Users/alorenzetti/CompiledPrograms/cdhit-4.8.1/cd-hit",
          "-i data/cds.fa",
          "-o data/cds/cdhitOut.fa",
          "-c 0.95",
          "-M 1000",
          "-T 8",
          "-sc 1")

system2(command = cdhit[1],
        args = paste(cdhit[-1], collapse = " "),
        stdout = "data/cds/cdhit.log",
        stderr = "data/cds/cdhit.err")

# # addendum
# # extracting sequences and aligning transposase CDS
# is4xt = list()
# is4xt$is200 = pfeiCDS$annot %>% filter(str_detect(string = product, pattern = "IS200") & is.na(pseudo))
# is4xt$is1341 = pfeiCDS$annot %>% filter(str_detect(string = product, pattern = "IS1341") & is.na(pseudo))
# is4xt$ish3 = pfeiCDS$annot %>% filter(str_detect(string = product, pattern = "ISH3$") & is.na(pseudo))
# is4xt$ish8 = pfeiCDS$annot %>% filter(str_detect(string = product, pattern = "ISH8$") & is.na(pseudo))
# is4xt$ish2 = pfeiCDS$annot %>% filter(str_detect(string = product, pattern = "ISH2\\)$") & is.na(pseudo))
# is4xt$ish4 = pfeiCDS$annot %>% filter(str_detect(string = product, pattern = "ISH4$") & is.na(pseudo))
# is4xt$ish7 = pfeiCDS$annot %>% filter(str_detect(string = product, pattern = "ISH7$") & is.na(pseudo))
# 
# # function to write fasta files using a
# # list containing several IS orfs
# xtfasta = function(lt){
#   for(i in names(lt)){
#     gr = GRanges(seqnames = lt[[i]]$seqnames,
#                  ranges = IRanges(start = lt[[i]]$start,
#                                   end = lt[[i]]$end),
#                  strand = lt[[i]]$strand)
#     names(gr) = paste0(lt[[i]]$locus_tag, "_", gsub(" ", "_", lt[[i]]$product))
#     
#     seqs = BSgenome::getSeq(genome, gr)
#     
#     trl = translate(x = seqs,
#                     genetic.code = getGeneticCode(id_or_name2 = "11"))
#     
#     writeXStringSet(x = trl,
#                     filepath = paste0("~/Downloads/", i, ".fa"),
#                     format = "fasta")
#   }
# }
# 
# xtfasta(is4xt)

