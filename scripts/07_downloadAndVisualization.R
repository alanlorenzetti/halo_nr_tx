#' ---
#' title: "_Halobacterium salinarum_ NRC-1 non-redundant transcriptome"
#' author: "Alan Lorenzetti"
#' date: "March 3rd, 2021"
#' ---
#'
#' Download the non-redundant (pseudo)transcriptome of _Halobacterium salinarum_ NRC-1 [here](https://alanlorenzetti.github.io/halo_nr_tx/data/Hsalinarum_nrtx.fa).  
#'
#' Download the dictionary file [here](https://alanlorenzetti.github.io/halo_nr_tx/data/dictionary.tsv).  
#'
#' Download the [COG (2020)](https://www.ncbi.nlm.nih.gov/research/cog) data for this non-redundant transcriptome [here](https://alanlorenzetti.github.io/halo_nr_tx/data/cog.tsv).  
#'
#' ## Browsable tables {.tabset}
#' ### Locus tag dictionary
#' <div style = "width:100%; height: auto; margin: auto;">
#+ echo=FALSE
datatable(dict, options = list(scrollX = "800px"))
#' </div>
#'
#' ### COG Info
#' <div style = "width:100%; height: auto; margin: auto;">
#+ echo=FALSE
datatable(cogdict, options = list(scrollX = "800px"))
#' </div>
