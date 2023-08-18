#run all code in this file to start the TMT analysis pipeline

#check that all packages are installed, and load them using pacman
if (!require('pacman', character.only=T, quietly=T)) {
  install.packages('pacman')
  library('pacman', character.only=T)
}else{
  library('pacman',character.only=T)
}

p_load(ggplot2)
p_load(tibble)
p_load(devtools)
p_load(openxlsx)
p_load(plyr)
p_load(dplyr)
p_load(shiny)
p_load(shinyscreenshot)

#install Bioconductor if needed
if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
BiocManager::install(update = TRUE, ask = FALSE)
}

#install PoissonSeq from Github if needed
if (!require('PoissonSeq',quietly=T)) {
  devtools::install_github("cran/PoissonSeq", upgrade = "never")
  library('PoissonSeq')
}else{
  library('PoissonSeq')
}
#install ggbiplot from Github if needed
if (!require('ggbiplot',quietly=T)) {
  devtools::install_github("vqv/ggbiplot", upgrade = "never")
  library('ggbiplot')
}else{
  library('ggbiplot')
}
#install shinyDirectrInput from Github if needed
if (!require('shinyDirectoryInput',quietly=T)) {
  devtools::install_github('wleepang/shiny-directory-input',upgrade = "never")
  library('shinyDirectoryInput')
}else{
  library('shinyDirectoryInput')
}

source("TMT_pseq_pipeline.R")
runApp('app.R')
