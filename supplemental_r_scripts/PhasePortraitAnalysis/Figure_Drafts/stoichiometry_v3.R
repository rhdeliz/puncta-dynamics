# Needs to be incorporated into heatmap stoichiometries paths.R

# Set scaling of plots
SCALE = 5
# Set number of facets
STARTING_RATIO_CAT_BINS = 10
# Frames since landing
FRAMES_SINCE_LANDING_BIN = 200
# 
# # Table Name
# SaveName <-
#   paste0(
#     "Normalized StatTable - ",
#     "LeadLag ", LEAD_LAG, " - ",
#     "StepSize ", round(STEP_SIZE, 2),
#     ".csv.gz"
#   )
# 
# SaveName <- gsub(" - \\.", "\\.", SaveName)
# 
# # Import table
# StatTableSeparate <- fread(file.path(OUTPUT_DIRECTORY, SaveName))
# 
# # Make combined table
# StatTableCombined <- NormStatTable
# StatTableCombined$IMAGE = "All"
# # Merge individual and combined tables
# StatTable <- rbindlist(list(StatTableSeparate, StatTableCombined))
# remove(StatTableSeparate, StatTableCombined)

# Import table
TracksTable <-
  NormStatTable %>% 
  filter(
    # To mirror phase portrait
    FRAMES_SINCE_LANDING_CAT %in% FRAMES_SINCE_LANDING_BIN,
    # COHORT == "MyD88 TRAF6"
    # Filter out unusually large clusters
    # REFERENCE_TOTAL_INTENSITY <= 2,
    # REFERENCE_TOTAL_INTENSITY >= 0,
    # QUERY_TOTAL_INTENSITY <= 2,
    # QUERY_TOTAL_INTENSITY >= 0
  ) %>% 
  arrange(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    # IMAGE,
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  # filter(row_number()==1 | row_number()==n()) %>%
  mutate(
    # Get track parameters
    # TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED),
    # LIFETIME = max(FRAMES_ADJUSTED) - min(FRAMES_ADJUSTED) + 1,
    RATIO = REFERENCE_TOTAL_INTENSITY/(QUERY_TOTAL_INTENSITY+REFERENCE_TOTAL_INTENSITY)
  ) %>%
  # filter(
  #   # To mirror phase portrait
  #   LIFETIME >= (LEAD_LAG*2+1)
  # ) %>%
  mutate(
    # Average ratios
    STARTING_RATIO = median(RATIO[1:SCALE], na.rm = T),
    ENDING_RATIO =  median(RATIO[n()-5:n()], na.rm = T),
    # Get ± ratio
    LEAD_RATIO = lead(RATIO, LEAD_LAG),
    LAG_RATIO = lag(RATIO, LEAD_LAG)
  ) %>% 
  mutate(
    # Bin data
    TIME_ADJUSTED = round(TIME_ADJUSTED/(4*SCALE))*(4*SCALE),
    # STARTING_RATIO_CAT = ifelse(STARTING_RATIO_CAT <= 0.45, 0.4, ifelse(STARTING_RATIO_CAT>=0.75, 0.8, 0.6)),
    STARTING_RATIO_CAT = round(STARTING_RATIO*STARTING_RATIO_CAT_BINS)/STARTING_RATIO_CAT_BINS,
    # RATIO = round(RATIO*(50/SCALE))/(50/SCALE),
    # Calculate ratio difference
    DELTA_RATIO = LEAD_RATIO - LAG_RATIO
  ) %>% 
  drop_na(
    DELTA_RATIO,
    STARTING_RATIO_CAT
  ) %>% 
  as.data.table()

PlotData <- 
  TracksTable %>% 
  filter(
    # QUERY_PROTEIN == "MyD88 TRAF6"
    # FRAMES_ADJUSTED == LEAD_LAG
    # QUERY_PROTEIN %in% c("IRAK1", "TRAF6", "TAB2", "NEMO", "RelA")
  ) %>% 
  as.data.table()

# ThresholdTable <- NULL
# ThresholdTable$QUERY_PROTEIN <- c("IRAK1", "TRAF6", "TAB2", "NEMO", "RelA")
# ThresholdTable$THRESHOLD <- c(.5, .65, .6, .55, .65)
# ThresholdTable <- as.data.table(ThresholdTable)
# ThresholdTable$QUERY_PROTEIN = factor(ThresholdTable$QUERY_PROTEIN, levels = PROTEIN_ORDER)
# 
# PlotData <- merge(PlotData, ThresholdTable, by = "QUERY_PROTEIN")
# PlotData$QUERY_PROTEIN <- factor(PlotData$QUERY_PROTEIN , levels = PROTEIN_ORDER)

PlotData <-
  PlotData %>% 
  mutate(
    # STARTING_RATIO_CAT = 0
    # STARTING_RATIO_CAT = floor(STARTING_RATIO_CAT*10)/10
    # STARTING_RATIO_CAT = ifelse(STARTING_RATIO >= THRESHOLD+.15, "MyD88-Rich",
    #                             ifelse(STARTING_RATIO <= THRESHOLD-.15, "Downstream-Rich", "Intermediate"))
  ) %>%
  mutate(
    STARTING_RATIO_CAT = 0
    # STARTING_RATIO_CAT = factor(STARTING_RATIO_CAT, levels = c("MyD88-Rich", "Intermediate", "Downstream-Rich"))
  ) %>% 
  group_by(
    QUERY_PROTEIN,
    STARTING_RATIO_CAT
  ) %>% 
  # mutate(
  #   N =n(),
  #   STARTING_RATIO_CAT = round(STARTING_RATIO*5)/5
  # ) %>%
  # filter(
  #   N >= 50
  # ) %>% 
  mutate(
    # Bin data
    # STARTING_RATIO_CAT = floor(STARTING_RATIO*STARTING_RATIO_CAT_BINS)/STARTING_RATIO_CAT_BINS,
    # Calculate ratio difference
    DELTA_RATIO = LEAD_RATIO - LAG_RATIO
  ) %>% 
  mutate(
    # STARTING_RATIO_CAT = 0,
    FRAMES_ADJUSTED = FRAMES_ADJUSTED/FPS
  ) %>% 
  drop_na(
    DELTA_RATIO,
    STARTING_RATIO_CAT
  ) %>%
  as.data.table()

PlotArrows <-
  PlotData %>% 
  mutate(
    FRAMES_ADJUSTED = round(FRAMES_ADJUSTED/80)*80,
    BIN_RATIO = round(RATIO*5)/5
  ) %>% 
  group_by(
    QUERY_PROTEIN,
    STARTING_RATIO_CAT,
    FRAMES_ADJUSTED,
    BIN_RATIO
  ) %>% 
  summarize(
    DELTA_RATIO = median(DELTA_RATIO),
    N = n()
  ) %>% 
  filter(
    N >= 25
  ) %>% 
  group_by(
    QUERY_PROTEIN
  ) %>% 
  mutate(
    DELTA_RATIO = DELTA_RATIO/max(abs(DELTA_RATIO))
  ) %>% 
  mutate(
    DELTA_RATIO_ANGLE = atan2(DELTA_RATIO, 1)*180/pi
  ) %>% 
  mutate(
    DELTA_RATIO_MAGNITUDE = abs(DELTA_RATIO)/quantile(abs(DELTA_RATIO), 0.75)
  ) %>% 
  mutate(
    DELTA_RATIO_MAGNITUDE = ifelse(DELTA_RATIO_MAGNITUDE >= 1, 1, DELTA_RATIO_MAGNITUDE)
  ) %>% 
  as.data.table()


  ggplot(
    PlotData,# %>% filter(STARTING_RATIO_CAT == "Intermediate") %>% as.data.table(),
    aes(
      FRAMES_ADJUSTED,
      RATIO
    )
  ) +
  stat_density2d(
    geom="tile", aes(fill = ..ndensity..), adjust = 3, contour = FALSE
  ) +
    geom_point(
      data = ThresholdTable,
      aes(
        x = -20,
        y = THRESHOLD
      ),
      size = 3,
      shape = 13,
      color = "blue"
    ) +
  geom_arrow(
    data = PlotArrows, #%>% filter(STARTING_RATIO_CAT == "Intermediate") %>% as.data.table(),
    aes(
      x = FRAMES_ADJUSTED,
      y = BIN_RATIO,
      mag = DELTA_RATIO_MAGNITUDE,
      angle = DELTA_RATIO_ANGLE
    ),
    # size  = .5,
    color = "black",
    # size = .01, # Small arrow head
    # arrow.length = .5, # Small arrow head
    lineend = "square"
  ) +
  scale_mag(
      max = 2, # Big arrow head
      guide = "none"
    ) +
  scale_fill_viridis(
    # limits = c(0, .75),
    option = "plasma",
    # na.value = "#0D1687",
    direction = -1
  ) + 
  scale_y_continuous(
    limits = c(-.1, 1.1)
  ) +
  scale_x_continuous(
    limits = c(-40, 800)
  ) +
  labs(
    x = "Cluster Time (s)",
    y = "% MyD88 in cluster",
    fill = "Scaled Density"
  ) +
  facet_grid(
    QUERY_PROTEIN~STARTING_RATIO_CAT
    
    # STARTING_RATIO_CAT~QUERY_PROTEIN
  ) +
  # facet_wrap(
  #   ~QUERY_PROTEIN,
  #   nrow = 1,
  #   scales = "free"
  # ) +
  theme_classic(
    base_size = 15
  ) +
    theme(
      strip.background = element_blank(),
      legend.position = "none"
    )
  
  ggsave(
    "stoichiometry_threshold_split.pdf",
    height = 2.25,
    width = 15 
  )
  
  ggsave(
    "stoichiometry_time_mag.pdf",
    height = 4,
    width = 16
  )
  
  ggsave(
    "stoichiometry_threshold2.pdf",
    height = 5,
    width = 5*1.5
  )
  
Correlations <-
  PlotData %>% 
  filter(
    FRAMES_ADJUSTED == 40
  ) %>% 
    mutate(
      STARTING_RATIO_CAT = round(STARTING_RATIO*20)/20
    ) %>% 
    group_by(
      QUERY_PROTEIN,
      STARTING_RATIO_CAT
    ) %>% 
    summarize(
      N = n(),
      START_END = cor(STARTING_RATIO, ENDING_RATIO),
      SD_ENDING_RATIO = sd(ENDING_RATIO),
    ) %>% 
  filter(
    N >= 25
  ) %>% 
  group_by(
    QUERY_PROTEIN
  ) %>% 
  mutate(
    SMOOTH_START_END = signal::sgolayfilt(START_END, p = 1, n = 3),
    SMOTH_SD_ENDING_RATIO = signal::sgolayfilt(SD_ENDING_RATIO, p = 1, n = 3),
  ) %>%
  mutate(
    DELTA_SMOOTH_START_END = lead(SMOOTH_START_END) - lag(SMOOTH_START_END),
    DELTA_SMOTH_SD_ENDING_RATIO = lead(SMOTH_SD_ENDING_RATIO) - lag(SMOTH_SD_ENDING_RATIO)
  #   SMOOTH_START_END = signal::sgolayfilt(SMOOTH_START_END, p = 1, n = 5),
  #   SMOTH_SD_ENDING_RATIO = signal::sgolayfilt(SMOTH_SD_ENDING_RATIO, p = 1, n = 5),
  ) %>%
    as.data.table()
  
ggplot(
  Correlations,
) +
  geom_hline(
    yintercept = 0,
    color = "Red"
  ) +
  geom_vline(
    data = ThresholdTable,
    aes(
      xintercept = THRESHOLD
    ),
    color = "blue"
  ) +
  geom_line(
    aes(
      STARTING_RATIO_CAT,
      START_END
    ),
    color = "grey"
  ) +
  geom_line(
    aes(
      STARTING_RATIO_CAT,
      SMOOTH_START_END
    )
  ) +
  labs(
    x = "Starting Stoichiometric Ratio (% MyD88 in cluster)",
    y = "R\n(Starting vs Ending\nStoichiometric Ratio)",
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  scale_y_continuous(
    limits = c(min(Correlations$START_END), max(Correlations$START_END))
  ) +
  facet_wrap(
    ~QUERY_PROTEIN,
    nrow = 1,
    # scales = "free"
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    strip.background = element_blank()
  )

ggsave(
  "stoichiometry_corr.pdf",
  height = 2,
  width = 15
)

ggplot(
  Correlations,
) +
  geom_vline(
    data = ThresholdTable,
    aes(
      xintercept = THRESHOLD
    ),
    color = "blue"
  ) +
  geom_line(
    aes(
      STARTING_RATIO_CAT,
      SD_ENDING_RATIO
    ),
    color = "grey"
  ) +
  # geom_line(
  #   aes(
  #     STARTING_RATIO_CAT,
  #     SMOTH_SD_ENDING_RATIO
  #   )
  # ) +
  labs(
    x = "Starting Stoichiometric Ratio (% MyD88 in cluster)",
    y = "S.Dev., Ending\nStoichiometric Ratio\n(% MyD88 in cluster)",
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  scale_y_continuous(
    limits = c(min(Correlations$SD_ENDING_RATIO), max(Correlations$SD_ENDING_RATIO))
  ) +
  facet_wrap(
    ~QUERY_PROTEIN,
    nrow = 1,
    scales = "free"
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    strip.background = element_blank()
  )  

ggsave(
  "stoichiometry_deviation.pdf",
  height = 2,
  width = 16
)

PlotData %>% 
filter(
  FRAMES_ADJUSTED == 40
) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      STARTING_RATIO,
      ENDING_RATIO
    )
  ) +
    stat_density2d(geom="tile", aes(fill = ..ndensity..), adjust = 2, contour = FALSE) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "black"
    ) +
    geom_point(
      data = ThresholdTable,
      aes(
        x = THRESHOLD,
        y = THRESHOLD
      ),
      color = "blue",
      size = 3,
      shape = 13
    ) +
    scale_fill_viridis(
      # trans = "log1p",
      option = "plasma",
      direction = -1
    ) + 
    labs(
      x = "Starting Stoichiometric Ratio\n(% MyD88 in cluster)",
      y = "Ending Stoichiometric Ratio\n(% MyD88 in cluster)",
      # fill = "Scaled Density\n(log1p)"
      fill = "Scaled Density"
    ) +
    scale_x_continuous(
      limits = c(-.1, 1.1)
    ) +
    scale_y_continuous(
      limits = c(-.1, 1.1)
    ) +
    facet_wrap(
      ~QUERY_PROTEIN,
      nrow = 1,
      scales = "free"
    ) +
    theme_classic(
      base_size = 15
    ) +
    theme(
      strip.background = element_blank(),
      aspect.ratio = 1
    )
  
  ggsave(
    "stoichiometry_start_end.png",
    height = 3,
    width = 16
  )
  
  ggsave(
    "stoichiometry_start_end.pdf",
    height = 3,
    width = 16
  )
  
  
  
  
  