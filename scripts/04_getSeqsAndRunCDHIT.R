# 20201228 alorenzetti

# description ####
# this script will read the genome file
# and extract sequences prior 
# to the execution of
# similarity clustering by CD-HIT

# loading genome file ####
genome = readDNAStringSet("data/Hsalinarum.fa")
names(genome) = str_replace(names(genome), " .*$", "")

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

# unifying and saving stable rnas
c(seqs$pfei$rna, seqs$ncbi$rna) %>% 
  writeXStringSet(filepath = "data/stableRNAs.fa",
                  format = "fasta")

# unifying and saving cds
c(seqs$pfei$cds, seqs$ncbi$cds) %>% 
  translate(x = ., genetic.code = arc) %>% 
  c(., custCDS) %>% 
  writeXStringSet(filepath = "data/cds.fa",
                  format = "fasta")

# running cd-hit using stable RNA files
if(!dir.exists("data/stableRNAs")){dir.create("data/stableRNAs")}
cdhit = c("cd-hit",
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
cdhit = c("cd-hit",
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
