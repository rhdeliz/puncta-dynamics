# args = strsplit("/raven/u/deliz/new_pipeline/pending_processing /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /raven/u/deliz/new_pipeline/phase_portraits 5 2 100 200 0.01 0.95 0", " ")
if(grepl("Darwin", as.character(Sys.info()['sysname']))){
  args = strsplit("/Users/u_deliz/Desktop/figure_drafts /Users/u_deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /Users/u_deliz/Desktop/figure_drafts 5 3 199 299 0.01 0.95 0", " ")
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
pacman::p_load(dtplyr, lemon, ggquiver, ggplot2, ggdark, scales, arrow,
               parallel, ggforce, data.table, viridis, RcppRoll, tidyr, dplyr, compiler)

if(grepl("macOS", osVersion)){
  pacman::p_load(metR)
}
filter <- dplyr::filter

if(
  !file.exists(OUTPUT_DIRECTORY)){
  dir.create(OUTPUT_DIRECTORY)
}

# OTHER
# Reference protein
USE_REFERENCE_PROTEIN = "MyD88"
# Order of protein appearance
# No Pellino!
PROTEIN_ORDER = c("MyD88", "IRAK4", "IRAK1", "TRAF6", "TNIP", "TNIP1",
                  "Pellino", "HOIL1", "NEMO", "TAB2", "A20", "RelA")

pacman::p_load(data.table, dplyr, dtplyr, arrow, ggplot2, viridis, tidyr, zoo, signal)

filter <- dplyr::filter

StatTable <- "/Users/u_deliz/Desktop/PhasePortraitAnalysis/Normalized StatTable - LeadLag 5.gz.parquet"

OUTPUT_DIRECTORY = "/Users/u_deliz/Desktop/figure_drafts"
setwd(OUTPUT_DIRECTORY)
# 
# StatTable <- lapply(StatTable, fread)
# StatTable <- rbindlist(StatTable)
# 
# 
# Table <- read_parquet(
#   "Essential.gz.parquet",
#   col_select = c(
#     "COHORT",
#     "LIGAND_DENSITY_CAT",
#     "IMAGE",
#     "CELL",
#     "UNIVERSAL_TRACK_ID",
#     "PROTEIN",
#     "NORMALIZED_INTENSITY",
#     "COMPLEMENTARY_NORMALIZED_INTENSITY_1",
#     "FRAME",
#     "TIME_ADJUSTED"
#   ))

SelectData <- read.csv(file.path(OUTPUT_DIRECTORY, "plot_selection.csv"))

SelectData <-
  SelectData %>%
  mutate(TEST = paste(COHORT, LIGAND_DENSITY_CAT, sep = "...")) %>%
  as.data.table()

StatTable <- read_parquet(StatTable)

StatTable <-
  StatTable %>%
  mutate(TEST = paste(COHORT, LIGAND_DENSITY_CAT, sep = "...")) %>%
  filter(
    TEST %in% SelectData$TEST,
    COHORT %in% paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER),
    FPS == 0.25
  ) %>% 
  select(-c(TEST)) %>% 
  as.data.table()

remove(SelectData)

SplitNormalization <-
  StatTable %>% 
  mutate(
    NORMALIZATION_GROUP = paste(LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, COHORT)
  ) %>% 
  as_tibble() %>% 
  group_split(
    NORMALIZATION_GROUP
  )

SplitNormalization <- SplitNormalization[order(sapply(SplitNormalization,nrow))]
SplitNormalization <- SplitNormalization[c(NROW(SplitNormalization):1)]

# Normalize intensities
NormalizeFx <- function(TableX){
  
  TableX <- as.data.table(TableX)
  
  ReferenceNormalization <-
    TableX %>% 
    group_by(
      COHORT,
      QUERY_PROTEIN,
      IMAGE
    ) %>% 
    summarize(
      RANGE = MAX_REFERENCE[1] - MIN_REFERENCE[1],
      MIN_REFERENCE = MIN_REFERENCE[1],
      MAX_REFERENCE = MAX_REFERENCE[1]
    ) %>%
    ungroup() %>% 
    mutate(
      MAX_REFERENCE_DISTANCE = abs(MAX_REFERENCE-median(MAX_REFERENCE))/median(MAX_REFERENCE)
    ) %>% 
    filter(
      MAX_REFERENCE_DISTANCE <= 1,
    ) %>% 
    group_by(
      IMAGE
    ) %>% 
    summarize(
      RANGE = MAX_REFERENCE[1] - MIN_REFERENCE[1],
      MIN_REFERENCE = MIN_REFERENCE[1],
      MAX_REFERENCE = MAX_REFERENCE[1]
    ) %>% 
    ungroup() %>% 
    mutate(
      FINAL_RANGE = median(RANGE)
    ) %>% 
    mutate(
      TEST = which.min(abs(RANGE - FINAL_RANGE)),
      ID = 1:n()
    ) %>% 
    filter(
      ID == TEST
    ) %>% 
    as.data.table()
  
  Normalization <-
    TableX %>% 
    group_by(
      COHORT,
      QUERY_PROTEIN,
      IMAGE
    ) %>% 
    summarize(
      RANGE = MAX_QUERY[1] - MIN_QUERY[1],
      MIN_QUERY = MIN_QUERY[1],
      MAX_QUERY = MAX_QUERY[1]
    ) %>%
    group_by(
      COHORT,
      QUERY_PROTEIN
    ) %>% 
    mutate(
      MAX_QUERY_DISTANCE = abs(MAX_QUERY-median(MAX_QUERY))/median(MAX_QUERY)
    ) %>% 
    filter(
      MAX_QUERY_DISTANCE <= 1
    ) %>% 
    group_by(
      COHORT,
      QUERY_PROTEIN,
      IMAGE
    ) %>% 
    summarize(
      RANGE = MAX_QUERY[1] - MIN_QUERY[1],
      MIN_REFERENCE = ReferenceNormalization$MIN_REFERENCE[1],
      MAX_REFERENCE = ReferenceNormalization$MAX_REFERENCE[1],
      MIN_QUERY = MIN_QUERY[1],
      MAX_QUERY = MAX_QUERY[1],
    ) %>% 
    group_by(
      COHORT,
      QUERY_PROTEIN
    ) %>% 
    mutate(
      FINAL_RANGE = median(RANGE)
      # FINAL_ASPECT_RATIO = median(ASPECT_RATIO),
    ) %>% 
    mutate(
      TEST = which.min(abs(RANGE - FINAL_RANGE)),
      ID = 1:n()
    ) %>% 
    filter(
      ID == TEST
    ) %>% 
    as.data.table()
  
  NormalizedTable <-
    TableX %>%
    as_tibble() %>% 
    mutate(
      MIN_REFERENCE = Normalization$MIN_REFERENCE[1],
      MAX_REFERENCE = Normalization$MAX_REFERENCE[1],
      MIN_QUERY = Normalization$MIN_QUERY[1],
      MAX_QUERY = Normalization$MAX_QUERY[1],
    ) %>%
    mutate(
      REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY*(MAX_REFERENCE-MIN_REFERENCE)+MIN_REFERENCE,
      ROUNDED_REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY*(MAX_REFERENCE-MIN_REFERENCE),
      ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY*(MAX_REFERENCE-MIN_REFERENCE),
      
      QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY*(MAX_QUERY - MIN_QUERY)+MIN_QUERY,
      ROUNDED_QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY*(MAX_QUERY - MIN_QUERY),
      ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY*(MAX_QUERY - MIN_QUERY)
    ) %>%
    as.data.table()
  
  NormalizedTable$NORMALIZATION_GROUP = NULL
  
  return(NormalizedTable)
}
NormStatTable <- lapply(SplitNormalization, NormalizeFx)
NormStatTable <- rbindlist(NormStatTable)
remove(SplitNormalization)
remove(StatTable)

# Set colors
PROTEIN_ORDER <- c("MyD88", PROTEIN_ORDER[PROTEIN_ORDER %in% unique(NormStatTable$QUERY_PROTEIN)])

ColorTable <- c("#000000", turbo(NROW(PROTEIN_ORDER), direction = -1))#[1:NROW(PROTEIN_ORDER)]
names(ColorTable) <- PROTEIN_ORDER

