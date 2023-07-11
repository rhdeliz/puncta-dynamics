setwd(OUTPUT_DIRECTORY)

StatTable <- fread("/Users/u_deliz/Desktop/PhasePortraitAnalysis/Normalized StatTable - LeadLag 5 - StepSize 2.csv.gz")

FRAMES_SINCE_LANDING_CAT_BIN = 200

setwd(OUTPUT_DIRECTORY)

# Table Name
SaveName <-
  paste0(
    "Normalized StatTable - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", round(STEP_SIZE, 2),
    ".csv.gz"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)

StatTable <- fread(file.path(OUTPUT_DIRECTORY, SaveName))

# StatTable <- 
#   StatTable %>% 
#   filter(
#     QUERY_PROTEIN == "TRAF6"
#   ) %>% 
#   as.data.table()

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
StatTable <- lapply(FRAMES_SINCE_LANDING_CAT_BIN, LandingTimeFx)
StatTable <- rbindlist(StatTable)

SplitNormalization <-
  StatTable %>% 
  mutate(
    NORMALIZATION_GROUP = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN)
  ) %>% 
  as_tibble() %>% 
  group_split(
    NORMALIZATION_GROUP
  )

SplitNormalization <- SplitNormalization[order(sapply(SplitNormalization,nrow))]
SplitNormalization <- SplitNormalization[c(NROW(SplitNormalization):1)]

NormalizeFx <- function(TableX){
  
  Normalization <-
    TableX %>%
    group_by(
      IMAGE
    ) %>% 
    summarize(
      RATIO = mean(MAX_REFERENCE/MAX_QUERY),
      MIN_REFERENCE = mean(MIN_REFERENCE),
      MIN_QUERY = mean(MIN_QUERY),
      MAX_REFERENCE = mean(MAX_REFERENCE),
      MAX_QUERY = mean(MAX_QUERY)
    ) %>% 
    ungroup() %>% 
    mutate(
      FINAL_RATIO = median(RATIO),
    ) %>% 
    mutate(
      TEST = which.min(abs(RATIO - FINAL_RATIO)),
      ID = 1:n()
    ) %>% 
    filter(
      ID == TEST
    ) %>% 
    as_tibble()
  
  NormalizedTable <-
    TableX %>%
    mutate(
      MIN_REFERENCE = Normalization$MIN_REFERENCE[1],
      MAX_REFERENCE = Normalization$MAX_REFERENCE[1],
      MIN_QUERY = Normalization$MIN_QUERY[1],
      MAX_QUERY = Normalization$MAX_QUERY[1],
    ) %>%
    mutate(
      REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE+MIN_REFERENCE,
      ROUNDED_REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE,
      ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE,
      
      QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY*MAX_QUERY+MIN_QUERY,
      ROUNDED_QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY*MAX_QUERY,
      ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY*MAX_QUERY
    ) %>%
    as.data.table()
  
  NormalizedTable$NORMALIZATION_GROUP = NULL
  
  return(NormalizedTable)
}
NormStatTable <- mclapply(SplitNormalization, NormalizeFx)
NormStatTable <- rbindlist(NormStatTable)
remove(SplitNormalization)

SplitStatTable <-
  NormStatTable %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    FPS = round(FPS, 2),
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
  
  LifetimeTable <-
    TableX %>% 
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>% 
    mutate(
      REFERENCE_TOTAL_INTENSITY = rollmedian(REFERENCE_TOTAL_INTENSITY, 11, align = c("right"), fill = NA),
    ) %>% 
    tidyr::drop_na(
      REFERENCE_TOTAL_INTENSITY
    ) %>% 
    mutate(
      STOICHIOMETRY = REFERENCE_TOTAL_INTENSITY/(REFERENCE_TOTAL_INTENSITY+QUERY_TOTAL_INTENSITY),
      MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY)
    ) %>% 
    mutate(
      STARTING_REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY[1:4], na.rm = T),
      STARTING_STOICHIOMETRY = median(STOICHIOMETRY[1:4], na.rm = T)
    ) %>% 
    filter(
      REFERENCE_TOTAL_INTENSITY == MAX_REFERENCE_TOTAL_INTENSITY
    ) %>% 
    mutate(
      # ID = 1:n(),
      ROUNDED_MAX_REFERENCE_TOTAL_INTENSITY = round(MAX_REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      ROUNDED_STARTING_REFERENCE_TOTAL_INTENSITY = round(STARTING_REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      ROUNDED_STARTING_STOICHIOMETRY = round(STARTING_STOICHIOMETRY*10)/10
    ) %>% 
    # filter(
      # ID == 1
    # ) %>%
    group_by(
      PLOT_FACETS,
      ROUNDED_STARTING_STOICHIOMETRY,
      ROUNDED_STARTING_REFERENCE_TOTAL_INTENSITY
    ) %>% 
    summarize(
      N = n(),
      TIME_ADJUSTED = median(TIME_ADJUSTED, na.rm = T)
    ) %>% 
    filter(
      N >=50
    ) %>% 
    as.data.table()
  
  
  ggplot(
    LifetimeTable,
    aes(
      x = ROUNDED_STARTING_STOICHIOMETRY,
      y = TIME_ADJUSTED,
      color = ROUNDED_STARTING_REFERENCE_TOTAL_INTENSITY,
      group = ROUNDED_STARTING_REFERENCE_TOTAL_INTENSITY
    )
  ) +
    geom_path() +
    geom_point(
      aes(
        size = N
      )
    ) +
    labs(
      x = "Starting Stoichiometry",
      y = "Time (s)",
      title = LifetimeTable$PLOT_FACETS[1]
    ) +
    scale_y_continuous(
      limits = c(0, max(LifetimeTable$TIME_ADJUSTED))
    ) +
    theme_classic()
  
  ggsave(
    file = paste0("StartStoichiometry vs Lifetime ", LifetimeTable$PLOT_FACETS[1], ".pdf"),
    height = 3,
    width = 4
  )
  
}
LinearPhasePortrait <- mclapply(SplitStatTable, SplitSummaryFx, mc.cores = detectCores(logical = F))
LinearPhasePortrait <- rbindlist(LinearPhasePortrait)
