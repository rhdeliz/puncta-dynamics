# Load libraries
pacman::p_load(ggplot2, ggdark, parallel, data.table, signal, dplyr, ijtiff, ggfx, ggforce, ggquiver, scales)
filter <- dplyr::filter

# Get track parameters
IMAGE = strsplit(SELECT_UNIVERSAL_TRACK_ID, '...', fixed = T)[[1]][1]
CELL = strsplit(SELECT_UNIVERSAL_TRACK_ID, '...', fixed = T)[[1]][2]
REFERENCE_PROTEIN = strsplit(SELECT_UNIVERSAL_TRACK_ID, '...', fixed = T)[[1]][3]
TRACK_ID = strsplit(SELECT_UNIVERSAL_TRACK_ID, '...', fixed = T)[[1]][4]

# Track paths
CELL_PATH = file.path(INPUT_DIRECTORY, COHORT, IMAGE, paste0("Cell_", CELL))
# Imports
# Import table
Table <- fread(file.path(CELL_PATH, "Analysis.csv.gz"))
Table <- Table %>% filter(UNIVERSAL_TRACK_ID == SELECT_UNIVERSAL_TRACK_ID)

# Import image
# Time parameters
if(is.null(SELECT_FRAMES_ADJSUTED)){
  IMPORT_FRAMES = min(Table$FRAME):max(Table$FRAME)
  # SELECT_FRAMES_ADJSUTED = min(Table$FRAMES_ADJUSTED):max(Table$FRAMES_ADJUSTED)
} else{
  IMPORT_FRAMES = SELECT_FRAMES_ADJSUTED + min(Table$FRAME)-1
}

OUTPUT_PATH <- file.path(OUTPUT_DIRECTORY, SELECT_UNIVERSAL_TRACK_ID)
dir.create(OUTPUT_PATH)
