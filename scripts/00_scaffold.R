# alorenzetti 20201228

# description ####
# this is the scaffold script
# it is going to call all 
# scripts required by this project

source("scripts/01_loadingLibs.R")
source("scripts/02_catAndParsePfeiTPA.R")
source("scripts/03_downloadAndParseNCBI.R")
source("scripts/04_getSeqsAndRunCDHIT.R")
source("scripts/05_parseCDHITout.R")

# loading and rendering using knitr
library(knitr)
rmarkdown::render(input = "scripts/06_downloadAndVisualization.R",
                  output_dir = ".",
                  output_file = "index.html",
                  output_format = "html_document")
