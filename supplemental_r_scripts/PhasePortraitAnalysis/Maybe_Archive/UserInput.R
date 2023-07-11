#!/usr/bin/env Rscript

remove(list = ls())
gc(reset = TRUE)
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Import variables
args = commandArgs(trailingOnly=TRUE)
# 
# args = strsplit("/raven/u/deliz/new_pipeline/Output4 /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /raven/u/deliz/new_pipeline/Output4/PhasePortraitAnalysis 5 2 100 200 F", " ")
# args = unlist(args)

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

# Log scale for plots
USE_MIN_FRAMES_SINCE_LANDING = args[8] #FALSE for linear scale
USE_MIN_FRAMES_SINCE_LANDING = as.logical(USE_MIN_FRAMES_SINCE_LANDING)
print(paste("USE_MIN_FRAMES_SINCE_LANDING =", USE_MIN_FRAMES_SINCE_LANDING))

# Import libraries
pacman::p_load(
  dtplyr, lemon, ggquiver, ggplot2, ggdark, scales, parallel,
  ggforce, data.table, viridis, RcppRoll, tidyr, dplyr,
  arrow)

filter <- dplyr::filter

if(grepl("Darwin", as.character(Sys.info()['sysname']))){
  pacman::p_load(metR)
  
  TABLE_PATH = "/Users/u_deliz/Desktop/PhasePortraitAnalysis"
  SCRIPTS_DIRECTORY = "/Users/u_deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis"
  OUTPUT_DIRECTORY = "/Users/u_deliz/Desktop/PhasePortraitAnalysis"
  LEAD_LAG = 5
  STEP_SIZE = 2
  # TIME PARAMETERS
  # Max number of frames per track
  MAX_FRAMES = 100 # frames
  # Frames since landing
  FRAMES_SINCE_LANDING_THRESHOLD = 200 #frames
  # Use actual frames since landing or min value?
  USE_MIN_FRAMES_SINCE_LANDING = F
}

if(
  !file.exists(OUTPUT_DIRECTORY)){
  dir.create(OUTPUT_DIRECTORY)
}

# OTHER
# Reference protein
USE_REFERENCE_PROTEIN = "MyD88"
# Order of protein appearance
PROTEIN_ORDER = c("MyD88", "IRAK4", "IRAK1", "TRAF6", "Pellino", "HOIL1", "NEMO", "A20")
MOVING_AVERAGE_WINDOW = 0
# Run Setup----
tryCatch({
    print(":::::::::::::::::::: Setup.R ::::::::::::::::::::")
    setwd(SCRIPTS_DIRECTORY)
    print("Running script Setup.R")
    source("Setup.R", local = T)
}, error = function(e) {print("Error with Setup.R")})

# Run Analysis----
tryCatch({
  print(":::::::::::::::::::: Analysis.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Analysis.R")
  source("Analysis.R", local = T)
}, error = function(e) {print("Error with Analysis.R")})

# Run Summary----
tryCatch({
  print(":::::::::::::::::::: Summarize.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Summarize.R")
  source("Summarize.R", local = T)
}, error = function(e) {print("Error with Summarize.R")})

# Plot Phase Portrait Linear----
tryCatch({
  if(LOG_SCALE == F){
    print(":::::::::::::::::::: PlotPhasePortraitLinear.R ::::::::::::::::::::")
    setwd(SCRIPTS_DIRECTORY)
    print("Running script PlotPhasePortraitLinear.R")
    source("PlotPhasePortraitLinear.R", local = T)
  }
}, error = function(e) {print("Error with PlotPhasePortraitLinear.R")})
