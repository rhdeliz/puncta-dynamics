# Load libraries
if("pacman" %in% rownames(installed.packages()) == FALSE)
{install.packages("pacman")}

# Library list
pacman::p_load(
  dplyr,
  stringr,
  parallel,
  data.table,
  tidyr,
  ff,
  igraph,
  ggplot2,
  ggdark,
  XML,
  fitdistrplus,
  scales,
  changepoint,
  ijtiff
) 
