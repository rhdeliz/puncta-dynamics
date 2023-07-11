# FOLDERS
# Scripts path
SCRIPTS_DIRECTORY = "/Users/u_deliz/dynamics_pipeline/supplemental_r_scripts/Montage"
# Images path
INPUT_DIRECTORY = "/Volumes/taylor-lab/Rafael/NewPipeline"
# Output path
OUTPUT_DIRECTORY = "/Users/u_deliz/Desktop/Montage"

# Track parameters
COHORT = "MyD88 IRAK4 Kinase_Dead"
SELECT_UNIVERSAL_TRACK_ID = "20210525 5nM 185_B8_MyD88_IRAK4-KO_IRAK-KD 001...7...MyD88...0"
USE_REFERENCE_PROTEIN = "MyD88"
# Spots to plot
SELECT_FRAMES_ADJSUTED <- c(10,15,20) #or NULL
# SELECT_FRAMES_ADJSUTED = NULL
# Plot other protein
PLOT_OTHER_PROTEIN = TRUE
# Pixel size
PIXEL_SIZE = 0.146666

# Image spot size for montage
BOX_SIZE = 15 #px

# Run Setup----
tryCatch({
  print(":::::::::::::::::::: Setup.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Setup.R")
  source("Setup.R", local = T)
}, error = function(e) {print("Error with Setup.R")})

tryCatch({
  print(":::::::::::::::::::: ImageTable.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Setup.R")
  source("ImageTable.R", local = T)
}, error = function(e) {print("Error with ImageTable.R")})


tryCatch({
  print(":::::::::::::::::::: ImagePlot.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script ImagePlot.R")
  source("ImagePlot.R", local = T)
}, error = function(e) {print("Error with ImagePlot.R")})


tryCatch({
  print(":::::::::::::::::::: TrackTable.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script TrackTable.R")
  source("TrackTable.R", local = T)
}, error = function(e) {print("Error with TrackTable.R")})


tryCatch({
  print(":::::::::::::::::::: TrackPlot.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script TrackPlot.R")
  source("TrackPlot.R", local = T)
}, error = function(e) {print("Error with TrackPlot.R")})
