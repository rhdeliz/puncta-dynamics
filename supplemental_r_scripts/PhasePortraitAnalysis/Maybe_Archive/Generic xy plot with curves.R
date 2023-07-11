# time vs total intensity
# color = initial stoichiometry

X_VARIABLE = "TOTAL_INTENSITY"
X_TITLE = "Total Intensity of Cluster"
X_BIN_SIZE = 2
# 
# X_VARIABLE = "TIME_ADJUSTED"
# X_TITLE = "Cluster Time (s)"
# X_BIN_SIZE = 8
# 
Y_VARIABLE = "RATIO"
Y_TITLE = "Stoichiometry"
# 
# Y_VARIABLE = "TOTAL_INTENSITY"
# Y_TITLE = "Total Intensity of Cluster"

CURVE = "STARTING_RATIO"
CURVE_TITLE = "Starting\nStoichiometry"
CURVE_BIN_SIZE = 1/6

# Frames since landing
FRAMES_SINCE_LANDING_BIN = 100

# Table Name
SaveName <-
  paste0(
    "Normalized StatTable - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", round(STEP_SIZE, 2),
    ".csv.gz"
  )

SaveName <- gsub(" - \\.", "\\.", SaveName)

# Import table
# StatTable <- fread(file = file.path(OUTPUT_DIRECTORY, SaveName))

# Pre-processing
TempTable <-
  StatTable %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT <= FRAMES_SINCE_LANDING_BIN,
    # COHORT == "MyD88 TRAF6",
    LIGAND_DENSITY_CAT == 32,
    # To mirror phase portrait
    # Filter out unusually large clusters
    REFERENCE_TOTAL_INTENSITY <= 2,
    REFERENCE_TOTAL_INTENSITY >= 0,
    QUERY_TOTAL_INTENSITY <= 2,
    QUERY_TOTAL_INTENSITY >= 0,
    FRAMES_ADJUSTED > min(FRAMES_ADJUSTED),
    FRAMES_ADJUSTED < max(FRAMES_ADJUSTED)
  ) %>% 
  arrange(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    IMAGE,
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY*max(MAX_REFERENCE),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY*max(MAX_QUERY)
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED),
    # Get track parameters
    LIFETIME = max(FRAMES_ADJUSTED) - min(FRAMES_ADJUSTED) + 1,
    TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY + QUERY_TOTAL_INTENSITY,
  ) %>%
  mutate(
    RATIO = REFERENCE_TOTAL_INTENSITY/(TOTAL_INTENSITY)
  ) %>% 
  mutate(
    STARTING_RATIO = median(RATIO[1:9])
  ) %>% 
  filter(
    # To mirror phase portrait
    LIFETIME >= (LEAD_LAG*2+1)
  ) %>%
  as.data.frame()

TempTable <- 
  TempTable %>% 
  mutate(
    X = !!as.name(X_VARIABLE),
    Y = !!as.name(Y_VARIABLE),
    CURVE = !!as.name(CURVE),
    SAMPLE_ID = UNIVERSAL_TRACK_ID,
    # R Facet
    GROUP = paste(COHORT, LIGAND_DENSITY_CAT, FPS)
  ) %>% 
  arrange(
    GROUP,
    CURVE,
    SAMPLE_ID,
    X,
    Y
  ) %>% 
  group_by(
    SAMPLE_ID
  ) %>% 
  mutate(
    ROUNDED_X = round(X/X_BIN_SIZE)*X_BIN_SIZE,
    CURVE = round(CURVE/CURVE_BIN_SIZE)*CURVE_BIN_SIZE
  ) %>% 
  group_by(
    GROUP,
    CURVE,
    ROUNDED_X
  ) %>% 
  summarize(
    Y = median(Y),
    N = n()
  ) %>% 
  filter(
    N >= 25
  ) %>% 
  group_by(
    GROUP,
    ROUNDED_X
  ) %>% 
  mutate(
    N = (N - min(N) + 1)/(max(N) - min(N) + 1)
  ) %>% 
  distinct() %>% 
  drop_na() %>% 
  as.data.table()

# Plot
ggplot(
  data = TempTable,
  aes(
    x = ROUNDED_X,
    y = Y,
    group = CURVE,
    color = as.factor(round(CURVE, 2))
  )
) +
  geom_point(
    aes(
      size = N
    )
  ) +
  geom_path() +
  scale_color_viridis(
    option = "plasma",
    discrete = T
  )  +
  labs(
    x = X_TITLE,
    y = Y_TITLE,
    color = CURVE_TITLE,
    size = "Fractional Clusters\nby x"
  ) +
  facet_wrap(
    ~GROUP
  ) +
  theme_classic() 