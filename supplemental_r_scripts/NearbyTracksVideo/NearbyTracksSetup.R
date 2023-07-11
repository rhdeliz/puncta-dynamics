# CLEAN ENVIRONMENT----
remove(list = ls())
gc(reset = TRUE)
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Folder containing all scripts
SCRIPTS_DIRECTORY = "/Users/u_deliz/image-analysis/SingleMolecule/08_Analysis"
# Folder containing all data
TOP_DIRECTORY = "/Volumes/taylor-lab/Rafael/SingleMoleculePipeline"
# Folder containing calibration and input data
SETUP_FOLDER = "00_Setup"
# Folder containing all images
ANALYSIS_FOLDER = "08_Analysis"
# List of images to analyze
IMAGE_INPUT_DATA = "2020-03-26 R_Input.csv"
# List of tracks
TRACK_INPUT_DATA = "TracksList.csv"
# Pixel size
PIXEL_SIZE = 0.1467
# Frame rate for output
VIDEO_FRAME_RATE = 3

if("RANN" %in% rownames(installed.packages()) == FALSE)
{install.packages("RANN")}

if("ijtiff" %in% rownames(installed.packages()) == FALSE)
{install.packages("ijtiff")}

if("gridExtra" %in% rownames(installed.packages()) == FALSE)
{install.packages("gridExtra")}
library(gridExtra)

if("viridis" %in% rownames(installed.packages()) == FALSE)
{install.packages("viridis")}
library(viridis)

if("scales" %in% rownames(installed.packages()) == FALSE)
{install.packages("scales")}
library(scales)

if("ggdark" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggdark")}
library(ggdark)

if("parallel" %in% rownames(installed.packages()) == FALSE)
{install.packages("parallel")}
library(parallel)

if("ggplot2" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggplot2")}
library(ggplot2)

if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr")}
library(dplyr)

if("foreach" %in% rownames(installed.packages()) == FALSE)
{install.packages("foreach")}
library(foreach)


if("data.table" %in% rownames(installed.packages()) == FALSE)
{install.packages("data.table")}
library(data.table)


# Import tables to create videos for
setwd(file.path(TOP_DIRECTORY, SETUP_FOLDER, "Input"))
InputImages <- fread(IMAGE_INPUT_DATA)
InputImages <- InputImages %>% select(IMAGE, COHORT) %>% distinct()
InputTracks <- fread(TRACK_INPUT_DATA)
Tracks <- InputTracks$UNIVERSAL_TRACK_ID
InputTracks <- strsplit(Tracks, '...', fixed = T)
InputTracks <- do.call("rbind", InputTracks)
InputTracks <- as_tibble(InputTracks)
names(InputTracks) <- c("IMAGE", "CELL", "PROTEIN", "TRACK_ID")
InputTracks$UNIVERSAL_TRACK_ID <- Tracks
Tracks <- merge(InputTracks, InputImages, by = "IMAGE")
remove(InputTracks, InputImages)

# Images to analyze
Tracks <-
  Tracks %>%
  mutate(
    FOLDER = file.path(TOP_DIRECTORY, ANALYSIS_FOLDER, COHORT, IMAGE, paste0("Cell_", CELL))
  ) %>%
  distinct()

# Create output folder
OUTPUT_FOLDER = file.path(TOP_DIRECTORY, ANALYSIS_FOLDER, "Tracks Graph")
dir.create(OUTPUT_FOLDER)

# List scripts folder
NEARBY_TRACKS_SCRIPT = file.path(SCRIPTS_DIRECTORY, "NearbyTracksVideo")

# Images to analyze
Images <-
  Tracks %>%
  select(
    FOLDER,
    PROTEIN
  ) %>%
  distinct()

ImagesFx <- function(ImageX) {
  tryCatch({
    # Loop parameters
    FOLDER = Images$FOLDER[ImageX]
    print(paste("Working on", FOLDER))
    lp.PROTEIN = Images$PROTEIN[ImageX]
    # Change folder
    setwd(FOLDER)
    # Load protein
    filename = paste(lp.PROTEIN, "Analysis", sep = "_")
    
    # Load tracks
    ExpTracks <-
      data.table::fread(paste(filename, ".csv.gz", sep = ""), integer64="character")
    
    # Get unique tracks
    TracksList <- unique(ExpTracks$UNIVERSAL_TRACK_ID)
    TracksList <-
      Tracks %>%
      filter(
        FOLDER == FOLDER
      )
    
    SubTracks <-
      Tracks %>%
      filter(
        FOLDER %in% Images$FOLDER[ImageX]
      ) %>%
      arrange(
        UNIVERSAL_TRACK_ID
      )
    
    SubTracks <- SubTracks$UNIVERSAL_TRACK_ID
    
    RETURN_TO_DIRECTORY <- getwd()
    setwd(NEARBY_TRACKS_SCRIPT)
    print("Running script NearbyTracksTable.R")
    source("NearbyTracksTable.R", local = T)
    
    # Combine track data with nearest spot data
    ExpTracks <- merge(ExpTracks, Distances, by = "UNIVERSAL_TRACK_ID", all = T)
    
    ExpTracks <-
      ExpTracks %>%
      arrange(
        -CLOSEST_DISTANCE,
        TRACK_ID,
        FRAME
      )
    
    START <- min(ExpTracks$FRAME) - 3
    END <- max(ExpTracks$FRAME) + 3
    
    RETURN_TO_DIRECTORY <- getwd()
    setwd(NEARBY_TRACKS_SCRIPT)
    print("Running script PullTIFFImages.R")
    source("PullTIFFImages.R", local = T)
    
    ObjectList <- ls(envir=environment())
    RETURN_TO_DIRECTORY <- getwd()
    setwd(NEARBY_TRACKS_SCRIPT)
    print("Running script MakeVideoFrames.R")
    source("MakeVideoFrames.R", local = T)
  }, error = function(e){print(paste("     ERROR with ImagesFx ImageX =", ImageX))})
}
lapply(1:NROW(Images), ImagesFx)

ggdark::invert_geom_defaults()
