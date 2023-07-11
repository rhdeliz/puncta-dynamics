# Set scaling of plots
TOTAL_INTENSITY_SCALE = 1/4
RATIO_SCALE = 10

# Frames since landing
FRAMES_SINCE_LANDING_BIN = 200

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
StatTable <- fread(file = file.path(OUTPUT_DIRECTORY, SaveName))


# Import table
TracksTable <-
  StatTable %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT <= 100,
    # COHORT == "MyD88 TRAF6",
    LIGAND_DENSITY_CAT == 32,
    # To mirror phase portrait
    # Filter out unusually large clusters
    REFERENCE_TOTAL_INTENSITY <= 2,
    REFERENCE_TOTAL_INTENSITY >= 0,
    QUERY_TOTAL_INTENSITY <= 2,
    QUERY_TOTAL_INTENSITY >= 0
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
    # Get track parameters
    LIFETIME = max(FRAMES_ADJUSTED) - min(FRAMES_ADJUSTED) + 1,
    TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY+REFERENCE_TOTAL_INTENSITY,
  ) %>%
  mutate(
    RATIO = REFERENCE_TOTAL_INTENSITY/(TOTAL_INTENSITY)
  ) %>% 
  filter(
    # To mirror phase portrait
    LIFETIME >= (LEAD_LAG*2+1)
  ) %>%
  mutate(
    # Get ± ratio
    LEAD_RATIO = lead(RATIO, LEAD_LAG),
    LAG_RATIO = lag(RATIO, LEAD_LAG),
    
    LEAD_TOTAL_INTENSITY = lead(TOTAL_INTENSITY, LEAD_LAG),
    LAG_TOTAL_INTENSITY = lag(TOTAL_INTENSITY, LEAD_LAG),
    
    ROUNDED_TOTAL_INTENSITY = round(TOTAL_INTENSITY*TOTAL_INTENSITY_SCALE)/TOTAL_INTENSITY_SCALE,
    ROUNDED_RATIO = round(RATIO*RATIO_SCALE)/RATIO_SCALE
  ) %>% 
  mutate(
    DELTA_RATIO = LEAD_RATIO - LAG_RATIO,
    DELTA_TOTAL_INTENSITY = LEAD_TOTAL_INTENSITY - LAG_TOTAL_INTENSITY
  ) %>% 
  drop_na(
    DELTA_RATIO,
    DELTA_TOTAL_INTENSITY
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    ROUNDED_TOTAL_INTENSITY,
    ROUNDED_RATIO
  ) %>% 
  summarize(
    DELTA_RATIO = median(DELTA_RATIO),
    DELTA_TOTAL_INTENSITY = median(DELTA_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    N >= 50
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS
  ) %>% 
  mutate(
    DELTA_RATIO = DELTA_RATIO/max(abs(DELTA_RATIO)),
    DELTA_TOTAL_INTENSITY = DELTA_TOTAL_INTENSITY/max(abs(DELTA_TOTAL_INTENSITY))
  ) %>% 
  distinct() %>% 
  as.data.table()


StreamTempTable <-
  TracksTable %>%
  as_tibble() %>%
  arrange(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    ROUNDED_TOTAL_INTENSITY,
    ROUNDED_RATIO
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS
  ) %>%
  mutate(
    ROUNDED_TOTAL_INTENSITY = ROUNDED_TOTAL_INTENSITY*TOTAL_INTENSITY_SCALE,
    ROUNDED_RATIO = ROUNDED_RATIO*RATIO_SCALE
  ) %>%
  complete(
    ROUNDED_TOTAL_INTENSITY = full_seq(c(min(ROUNDED_TOTAL_INTENSITY), max(ROUNDED_TOTAL_INTENSITY)), period = 1),
    ROUNDED_RATIO = full_seq(c(min(ROUNDED_RATIO), max(ROUNDED_RATIO)), period = 1)
  ) %>% 
  mutate(
    ROUNDED_TOTAL_INTENSITY = ROUNDED_TOTAL_INTENSITY/TOTAL_INTENSITY_SCALE,
    ROUNDED_RATIO = ROUNDED_RATIO/RATIO_SCALE
  )

# Add derivatives
StreamTempTable$DELTA_TOTAL_INTENSITY[is.na(StreamTempTable$DELTA_TOTAL_INTENSITY)] = 0
StreamTempTable$DELTA_RATIO[is.na(StreamTempTable$DELTA_RATIO)] = 0
StreamTempTable$N[is.na(StreamTempTable$N)] = 0



ggplot(
) +
geom_streamline(
  data = StreamTempTable,
  aes(
    x = ROUNDED_TOTAL_INTENSITY,
    y = ROUNDED_RATIO,
    dx = DELTA_TOTAL_INTENSITY,
    dy = DELTA_RATIO,
    color = sqrt(..dx..^2 + ..dy..^2),
    size = ..step..
    # alpha = ..step..
  ),
  arrow = NULL,
  n = 10,
  # arrow.length = 0.3,
  # jitter = 4,
  L = 10, res = 10, lineend = "round"
) +
geom_arrow(
  data = TracksTable,
  aes(
    x = ROUNDED_TOTAL_INTENSITY,
    y = ROUNDED_RATIO,
    mag = 0.5,
    angle = atan2(DELTA_RATIO, DELTA_TOTAL_INTENSITY)*180/pi,
    color = sqrt(DELTA_RATIO^2 + DELTA_TOTAL_INTENSITY^2),
    # color = FRAMES_ADJUSTED
  ),
  # size = .01, # Small arrow head
  # arrow.length = .5, # Small arrow head

  size = 1, # Big arrow head
  arrow.length = 1,# Big arrow head
  lineend = "square"
) +
  scale_size(range = c(.01, 1), guide = "none") +
  scale_alpha(guide = "none") +
  scale_mag(
    # max = 3, # Small arrow head
    max = 1, # Big arrow head
    guide = 'none'
  ) +
  scale_color_viridis(
    option = "plasma",
    guide = "none"
    # direction = -1,
    # trans = "log10",
    # breaks = trans_breaks("log10", function(x) 10^x),
    # labels = trans_format("log10", math_format(10^.x))
  )  +
  labs(
    x = "Total Intensity of Cluster",
    y = "Stiochiometry\nMyD88 / (MyD88 + TRAF6)"
  ) +
  facet_wrap(
    ~COHORT
  ) +
  theme_classic() #+
  # coord_fixed()
