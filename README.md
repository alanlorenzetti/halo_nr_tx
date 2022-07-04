# Generating a non-redundant transcriptome for _Halobacterium salinarum_

Some of _Halobacterium salinarum_ NRC-1 genes are repeated. Reducing it to a non-redundant set of genes can help simplifying analyses. Moreover, there have been several annotation efforts, causing a lot of pain to integrate findings from distinct studies due to non-standardized locus tags.

This is an attempt to create a non-redundant (pseudo)transcriptome for _Halobacterium salinarum_ NRC-1 using as source the annotation released by [Pfeiffer et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6624760/). A dictionary containing alternative locus_tag was also created using annotation provided by NCBI's RefSeq prokaryotic automatic annotation pipeline and a custom in-house file derived from the original annotation provided by [Ng et al. 2000](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC17314/).

Briefly:  

1.1) Pfeiffer et al. 2019 third party annotation was manually downloaded from:  
  * https://www.ncbi.nlm.nih.gov/nuccore/BK010829 (data/chr.gff3);
  * https://www.ncbi.nlm.nih.gov/nuccore/BK010830 (data/pnrc100.gff3);
  * https://www.ncbi.nlm.nih.gov/nuccore/BK010831 (data/pnrc200.gff3).  
  
1.2) Files were renamed and concatenated;  
1.3) Accession names where replaced in GFF3 file to match those of NCBI RefSeq;  

2.1) NCBI RefSeq annotations and genomes ([ASM680v1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000006805.1/ and [ASM6902v1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000069025.1/)) were automatically downloaded and minimally parsed from the ftp resource;  

3.1) Ng et al. 2000 derived annotation containing protein sequences was loaded and parsed;  

4.1) For CDS, nucleotide sequences were obtained using the genome and both annotations, and converted to amino acid using the translation table #11. The aa sequences were further clustered using CD-HIT (95% identity) to obtain the non-redundant set;  
4.2) For RNAs, nucleotide sequences were obtained using the genome and both annotations. The nt sequences were further clustered using CD-HIT (99% identity) to obtain the non-redundant set;  
4.3) A cluster of a single locus_tag, coming from NCBI will be automatically discarded;  

5.1) Finally, it is going to knit a Rmd file containing the dictionary table;  

Please, check the final report [here](https://alanlorenzetti.github.io/halo_nr_tx/index.html).  