setwd(OUTPUT_DIRECTORY)

# Table Name
SaveName <-
  paste0(
    "Normalized StatTable - ",
    "LeadLag ", LEAD_LAG,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)
StatTable <- read_parquet(SaveName)

CohortList <- StatTable[!kit::fduplicated(StatTable$COHORT), ]
CohortList <- CohortList$COHORT
CohortList <- CohortList[grep("rid", CohortList)]
CohortList <- c(CohortList, paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER))

SelectData <- read.csv("/Users/u_deliz/Desktop/figure_drafts/plot_selection.csv")
Images <- read.csv("/Users/u_deliz/Desktop/new_ligand.csv")
Images <- Images %>% filter(KEEP != "NO") %>% as.data.table()

SelectData <-
  SelectData %>%
  mutate(TEST = paste(COHORT, LIGAND_DENSITY_CAT, sep = "...")) %>%
  as.data.table()

StatTable <-
  StatTable %>%
  mutate(TEST = paste(COHORT, LIGAND_DENSITY_CAT, sep = "...")) %>%
  filter(
    TEST %in% SelectData$TEST,
    FPS == 0.25,
    IMAGE %in% Images$IMAGE
  ) %>%
  as.data.table()
# 
# Images <- c("20221108 10nM 231-G2_TRAF6_NEMO 001",
#             "20221108 10nM 231-G5_TRAF6_NEMO 001")

StatTable <-
  StatTable %>%
  filter(
    # IMAGE %in% Images,
    COHORT %in% CohortList
  ) %>%
  # mutate(
  #   COHORT = "TRAF6 NEMO"
  # ) %>%
  group_by(
    IMAGE
  ) %>%
  mutate(
    N = n()
  ) %>%
  filter(
    N >= 175
  ) %>%
  mutate(
    NORMALIZATION_GROUP = paste(LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, COHORT)
  ) %>%
  as.data.table()

# Separate by landing time
LandingTimeFx <- function(ThresholdX){
  TempTable <-
    StatTable %>%
    filter(
      FRAMES_SINCE_LANDING_CAT <= ThresholdX
    ) %>%
    mutate(
      FRAMES_SINCE_LANDING_CAT = ThresholdX
    ) %>%
    as.data.table()

  return(TempTable)
}
LandingTimes <- unique(StatTable$FRAMES_SINCE_LANDING_CAT)
StatTable <- lapply(LandingTimes, LandingTimeFx)
StatTable <- rbindlist(StatTable)

SplitNormalization <-
  StatTable %>% 
  as_tibble() %>% 
  group_split(
    NORMALIZATION_GROUP
  )

SplitNormalization <- SplitNormalization[order(sapply(SplitNormalization,nrow))]
SplitNormalization <- SplitNormalization[c(NROW(SplitNormalization):1)]

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

# Table Name
SaveName <-
  paste0(
    "NormStatTable - ",
    "LeadLag ", LEAD_LAG,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)
write_parquet(NormStatTable, file.path(OUTPUT_DIRECTORY, SaveName))

SplitStatTable <-
  NormStatTable %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    FPS = round(FPS, 2)
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    FACET = paste(COHORT, LIGAND_DENSITY_CAT, FPS, FRAMES_SINCE_LANDING_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN),
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN)
  ) %>% 
  as_tibble() %>% 
  group_split(
    FACET
  )

SplitStatTable <- SplitStatTable[order(sapply(SplitStatTable,nrow))]
SplitStatTable <- SplitStatTable[c(NROW(SplitStatTable):1)]

SplitSummaryFx <- function(TableX){
  print(TableX$FACET[1])
  # Summarize results dealing with either reference OR query increasing
  PhasePortrait <-
    TableX %>% 
    filter(
      !is.infinite(DELTA_REFERENCE_TOTAL_INTENSITY),
      !is.infinite(DELTA_QUERY_TOTAL_INTENSITY),
      
      !is.infinite(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
      !is.infinite(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY)
    ) %>% 
    as_tibble() %>%
    group_by(
      COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN,
      FACET, PLOT_FACETS,
      ROUNDED_REFERENCE_TOTAL_INTENSITY, ROUNDED_QUERY_TOTAL_INTENSITY,
      FRAMES_SINCE_LANDING_CAT
    ) %>% 
    summarize(
      # Get number of spots in bin so that bins with few spots can be filtered out
      N = n(),
      FRAMES_ADJUSTED = median(FRAMES_ADJUSTED),
      TIME_ADJUSTED = median(TIME_ADJUSTED),
      
      # Deviation from median
      MAD_DELTA_REFERENCE_TOTAL_INTENSITY = mad(DELTA_REFERENCE_TOTAL_INTENSITY),
      MAD_DELTA_QUERY_TOTAL_INTENSITY = mad(DELTA_QUERY_TOTAL_INTENSITY),
      MAD_ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = mad(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
      MAD_ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = mad(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      MAD_DELTA_RATIO = mad(DELTA_REFERENCE_TOTAL_INTENSITY/DELTA_QUERY_TOTAL_INTENSITY),
      MAD_ADJUSTED_DELTA_RATIO = mad(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY/ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),

      DELTA_REFERENCE_Q1 = quantile(DELTA_REFERENCE_TOTAL_INTENSITY, .25, na.rm = T),
      DELTA_REFERENCE_Q3 = quantile(DELTA_REFERENCE_TOTAL_INTENSITY, .75, na.rm = T),
      DELTA_QUERY_Q1 = quantile(DELTA_QUERY_TOTAL_INTENSITY, .25, na.rm = T),
      DELTA_QUERY_Q3 = quantile(DELTA_QUERY_TOTAL_INTENSITY, .75, na.rm = T),
      ADJUSTED_DELTA_REFERENCE_Q1 = quantile(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, .25, na.rm = T),
      ADJUSTED_DELTA_REFERENCE_Q3 = quantile(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, .75, na.rm = T),
      ADJUSTED_DELTA_QUERY_Q1 = quantile(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, .25, na.rm = T),
      ADJUSTED_DELTA_QUERY_Q3 = quantile(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, .75, na.rm = T),

      # Get median of intensity changes. mean could skew the data, thus not used
      DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
      DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T),
      
      ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = median(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
      ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = median(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, na.rm = T)
    ) %>% 
    group_by(
      COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN,
      FACET, PLOT_FACETS,
      ROUNDED_REFERENCE_TOTAL_INTENSITY, ROUNDED_QUERY_TOTAL_INTENSITY,
      FRAMES_SINCE_LANDING_CAT
    ) %>% 
    mutate(
      THETA = atan2(DELTA_QUERY_TOTAL_INTENSITY, DELTA_REFERENCE_TOTAL_INTENSITY),
      ADJUSTED_THETA = atan2(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
      
      # Get angle to account for negative magnitudes
      THETA = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
      ADJUSTED_THETA = atan2(-ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, -ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      
      DELTA_RATIO = DELTA_REFERENCE_TOTAL_INTENSITY/DELTA_QUERY_TOTAL_INTENSITY,
      ADJUSTED_DELTA_RATIO = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY/ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY,
    ) %>% 
    mutate(
      # Make angle
      THETA = (THETA*180/pi+180)/360,
      ADJUSTED_THETA = (ADJUSTED_THETA*180/pi+180)/360,
      
      DELTA_RATIO = ifelse(DELTA_RATIO<1, -1/DELTA_RATIO, DELTA_RATIO),
      ADJUSTED_DELTA_RATIO = ifelse(ADJUSTED_DELTA_RATIO<1, -1/ADJUSTED_DELTA_RATIO, ADJUSTED_DELTA_RATIO)
    ) %>% 
    filter(
      N >= 5
    ) %>% 
    group_by(
      COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN,
      FRAMES_SINCE_LANDING_CAT,
      FACET, PLOT_FACETS
    ) %>% 
    mutate(
      N_TEST = ifelse(N >= 175, T, F)
    ) %>%
    as.data.table()
    
  PhasePortrait$QUERY_PROTEIN = factor(PhasePortrait$QUERY_PROTEIN, levels = PROTEIN_ORDER)
  return(PhasePortrait)
}
PhasePortraitAll <- mclapply(SplitStatTable, SplitSummaryFx)
PhasePortraitAll <- rbindlist(PhasePortraitAll)

# Write table
SaveName <-
  paste0(
    "Normalized All PhasePortrait - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", STEP_SIZE,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)

file.remove(file.path(OUTPUT_DIRECTORY, SaveName))
write_parquet(PhasePortraitAll, file.path(OUTPUT_DIRECTORY, SaveName))
