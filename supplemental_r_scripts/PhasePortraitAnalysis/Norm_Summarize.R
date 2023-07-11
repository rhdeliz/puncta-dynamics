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

StatTable <-
  StatTable %>%
  filter(COHORT %in% paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER)) %>%
  group_by(
    IMAGE
  ) %>% 
  mutate(
    N = n(),
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY*(MAX_REFERENCE-MIN_REFERENCE)+MIN_REFERENCE,
    ROUNDED_REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY*(MAX_REFERENCE-MIN_REFERENCE),
    ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY*(MAX_REFERENCE-MIN_REFERENCE),
    
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY*(MAX_QUERY - MIN_QUERY)+MIN_QUERY,
    ROUNDED_QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY*(MAX_QUERY - MIN_QUERY),
    ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY*(MAX_QUERY - MIN_QUERY)
  ) %>%
  filter(
    N >= 250
  ) %>% 
  as.data.table()

# Separate by landing time
LandingTimeFx <- function(ThresholdX){
  TempTable <- 
    StatTable %>%
    filter(
      FRAMES_SINCE_LANDING_CAT == ThresholdX
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

StatTable <-
  StatTable %>% 
  mutate(
    FPS = round(FPS, 2)
  ) %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT, IMAGE
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, IMAGE
  ) %>%
  mutate(
    IMAGENUMBER = cur_group_id()
  ) %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT
  ) %>%
  mutate(
    IMAGENUMBER = IMAGENUMBER - min(IMAGENUMBER) + 1
  ) %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    FACET = paste(COHORT, LIGAND_DENSITY_CAT, FPS, FRAMES_SINCE_LANDING_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER),
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER)
  ) %>% 
  as.data.table()

SplitStatTable <-
  StatTable %>% 
  as_tibble() %>% 
  group_split(
    FACET
  )

SplitStatTable <- SplitStatTable[order(sapply(SplitStatTable,nrow))]
SplitStatTable <- SplitStatTable[c(NROW(SplitStatTable):1)]

SplitSummaryFx <- function(TableX){
  
  # Summarize results dealing with either reference OR query increasing
  PhasePortrait <-
    as.data.table(TableX) %>% 
    filter(
      !is.infinite(DELTA_REFERENCE_TOTAL_INTENSITY),
      !is.infinite(DELTA_QUERY_TOTAL_INTENSITY),
      
      !is.infinite(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
      !is.infinite(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY)
    ) %>% 
    # as_tibble() %>%
    group_by(
      IMAGE, IMAGENUMBER,
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
      IMAGE, IMAGENUMBER,
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
      IMAGE, IMAGENUMBER,
      COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN,
      FRAMES_SINCE_LANDING_CAT,
      FACET, PLOT_FACETS
    ) %>% 
    mutate(
      N_TEST = ifelse(N >= 175, T, F)
      
      # # Get absolute magnitude so that the log can be taken later
      # MAD_MAGNITUDE = abs(MAD_DELTA_REFERENCE_TOTAL_INTENSITY) + abs(MAD_DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # MAD_RAD_ANGLE = atan2(-MAD_DELTA_REFERENCE_TOTAL_INTENSITY, -MAD_DELTA_QUERY_TOTAL_INTENSITY),
      # # Make anglez
      # MAD_DEG_ANGLE = MAD_RAD_ANGLE*180/pi+180,
      # MAD_DEG_ANGLE = MAD_DEG_ANGLE/360,
      # 
      # # Get absolute magnitude so that the log can be taken later
      # MAD_ADJUSTED_MAGNITUDE = abs(MAD_ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY) + abs(MAD_ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # MAD_ADJUSTED_RAD_ANGLE = atan2(-MAD_ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, -MAD_ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Make angle
      # MAD_ADJUSTED_DEG_ANGLE = MAD_ADJUSTED_RAD_ANGLE*180/pi+180,
      # MAD_ADJUSTED_DEG_ANGLE = MAD_ADJUSTED_DEG_ANGLE/360,
      # 
      # # # Get absolute magnitude so that the log can be taken later
      # MAGNITUDE = abs(DELTA_REFERENCE_TOTAL_INTENSITY) + abs(DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # RAD_ANGLE = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
      # # Make angle
      # DEG_ANGLE = RAD_ANGLE*180/pi+180,
      # DEG_ANGLE = DEG_ANGLE/360,
      # 
      # # Get absolute magnitude so that the log can be taken later
      # ADJUSTED_MAGNITUDE = abs(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY) + abs(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # ADJUSTED_RAD_ANGLE = atan2(-ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, -ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Make angle
      # ADJUSTED_DEG_ANGLE = ADJUSTED_RAD_ANGLE*180/pi+180,
      # ADJUSTED_DEG_ANGLE = ADJUSTED_DEG_ANGLE/360,
    ) %>%
    as.data.table()
  
  PhasePortrait$QUERY_PROTEIN = factor(PhasePortrait$QUERY_PROTEIN, levels = PROTEIN_ORDER)
  return(PhasePortrait)
}
PhasePortrait <- mclapply(SplitStatTable, SplitSummaryFx, mc.cores = detectCores(logical = F))
PhasePortrait <- rbindlist(PhasePortrait)
remove(SplitStatTable)

print(paste("PhasePortrait Images=", NROW(unique(PhasePortrait$IMAGE))))

SaveName <-
  paste0(
    "Normalized PhasePortrait - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", STEP_SIZE,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)

file.remove(file.path(OUTPUT_DIRECTORY, SaveName))
write_parquet(PhasePortrait, file.path(OUTPUT_DIRECTORY, SaveName))
