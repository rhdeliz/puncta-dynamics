#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
args = strsplit(args, " ")

# Import variables ----
if(grepl("Darwin", as.character(Sys.info()['sysname']))){
  # args = strsplit("/Users/u_deliz/Desktop/grid /Users/u_deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /Users/u_deliz/Desktop/grid 5 1 199 299 0.01 0.95 0", " ")
  args = strsplit("/Users/u_deliz/Desktop/PhasePortraitAnalysis /Users/u_deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /Users/u_deliz/Desktop/PhasePortraitAnalysis 5 3 199 299 0.01 0.95 0", " ")
  # args = strsplit("/Users/u_deliz/Desktop/20221110 /Users/u_deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /Users/u_deliz/Desktop/20221110/pp 5 5 199 299 0.01 0.95 0", " ")
} else{
  args = strsplit('/raven/u/deliz/new_pipeline /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /raven/u/deliz/new_pipeline/phase_portraits 5 3 199 299 0.01 0.95 0', " ")
}
args = unlist(args)

# Table location
TABLE_PATH = args[1]
print(paste("TABLE_PATH =", TABLE_PATH))

# Script location
SCRIPTS_DIRECTORY = args[2]
print(paste("SCRIPTS_DIRECTORY =", SCRIPTS_DIRECTORY))

# Output folder
OUTPUT_DIRECTORY = args[3]
print(paste("OUTPUT_DIRECTORY =", OUTPUT_DIRECTORY))

# Lead/Lag factor
LEAD_LAG = args[4] #frame(s) before and __ frame(s) after
LEAD_LAG = as.numeric(LEAD_LAG)
print(paste("LEAD_LAG =", LEAD_LAG))

# Grid size for binning and plots
STEP_SIZE = args[5] # molecules
STEP_SIZE = as.numeric(STEP_SIZE)
print(paste("STEP_SIZE =", STEP_SIZE))

MAX_FRAMES = args[6] #FALSE for linear scale
MAX_FRAMES = as.numeric(MAX_FRAMES)
print(paste("MAX_FRAMES =", MAX_FRAMES))

FRAMES_SINCE_LANDING_THRESHOLD = args[7] #FALSE for linear scale
FRAMES_SINCE_LANDING_THRESHOLD = as.numeric(FRAMES_SINCE_LANDING_THRESHOLD)
print(paste("FRAMES_SINCE_LANDING_THRESHOLD =", FRAMES_SINCE_LANDING_THRESHOLD))

LOWER_BOUNDARY = args[8] #FALSE for linear scale
LOWER_BOUNDARY = as.numeric(LOWER_BOUNDARY)
print(paste("LOWER_BOUNDARY =", LOWER_BOUNDARY))

UPPER_BOUNDARY = args[9] #FALSE for linear scale
UPPER_BOUNDARY = as.numeric(UPPER_BOUNDARY)
print(paste("UPPER_BOUNDARY =", UPPER_BOUNDARY))

MOVING_AVERAGE_WINDOW = args[10] #FALSE for linear scale
MOVING_AVERAGE_WINDOW = as.numeric(MOVING_AVERAGE_WINDOW)
print(paste("MOVING_AVERAGE_WINDOW =", MOVING_AVERAGE_WINDOW))

# Import libraries
pacman::p_load(lemon, ggquiver, ggplot2, ggdark, scales, arrow, R.utils, metR, parallel,
               ggforce, data.table, viridis, RcppRoll, tidyr, dplyr, compiler, dtplyr, kit)

if(grepl("macOS", osVersion)){
  pacman::p_load(metR)
}

filter <- dplyr::filter
select <- dplyr::select

if(
  !file.exists(OUTPUT_DIRECTORY)){
  dir.create(OUTPUT_DIRECTORY)
}

# OTHER
# Reference protein
USE_REFERENCE_PROTEIN = "MyD88"
# USE_REFERENCE_PROTEIN = "TRAF6"

# Order of protein appearance
PROTEIN_ORDER = c("MyD88", "IRAK4", "IRAK1", "TRAF6", "TNIP", "TNIP1",
                  "Pellino", "HOIL1", "NEMO", "TAB2", "A20", "RelA")

# Scripts to execute----
ScriptList <- c(
  # # Need to modify tables later so that step size is not included
  "Setup.R",
  "Analysis.R",

  "Norm_All_Summarize.R",
  "Norm_All_Plot.R"
  # "1D Phase Portrait Cuts All.R",
  # 
  # "Norm_Summarize.R"
  # "Norm_Plot.R"
  # 
  # "Summarize.R",
  # "Plot.R"
  # 
  # "Stretch_Summarize.R",
  # "Stretch_Plot.R"
)

# Command
ExecuteScript <- function(x){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::"))
    setwd(SCRIPTS_DIRECTORY)
    source(x, local = T)
  }, error = function(e) {print(paste("Error loading", x))})
}
lapply(ScriptList, ExecuteScript)
