# 20201228 alorenzetti

# description ####
# this script will load all the libs required
# in this project

# if(!require(pacman)){install.packages("pacman")}
library(pacman)

# additional packages
packs = c("tidyverse",
          "rtracklayer",
          "Biostrings",
          "BSgenome",
          "gtools",
          "knitr",
          "DT")

# loading libs
p_load(char = packs)
