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
    STARTING_RATIO_CAT = median(RATIO[1:SCALE], na.rm = T),
    # Get ± ratio
    LEAD_RATIO = lead(RATIO, LEAD_LAG),
    LAG_RATIO = lag(RATIO, LEAD_LAG)
  ) %>% 
  mutate(
    # Bin data
    TIME_ADJUSTED = round(TIME_ADJUSTED/(4*SCALE))*(4*SCALE),
    STARTING_RATIO_CAT = ifelse(STARTING_RATIO_CAT <= 0.45, 0.4, ifelse(STARTING_RATIO_CAT>=0.75, 0.8, 0.6)),
    # STARTING_RATIO_CAT = round(STARTING_RATIO_CAT*STARTING_RATIO_CAT_BINS)/STARTING_RATIO_CAT_BINS,
    RATIO = round(RATIO*(50/SCALE))/(50/SCALE),
    # Calculate ratio difference
    DELTA_RATIO = LEAD_RATIO - LAG_RATIO
  ) %>% 
  drop_na(
    DELTA_RATIO,
    STARTING_RATIO_CAT
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    # IMAGE,
    STARTING_RATIO_CAT,
    TIME_ADJUSTED,
    RATIO
  ) %>% 
  summarize(
    # Get typical ratio
    DELTA_RATIO = median(DELTA_RATIO, na.rm = T),
    N = n()
  ) %>% 
  filter(
    # Filter sparce data
    N >= 25
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    # IMAGE,
    STARTING_RATIO_CAT
  ) %>%
  mutate(
    DELTA_RATIO = DELTA_RATIO/max(abs(DELTA_RATIO)),
    STARTING_RATIO_CAT = formatC(STARTING_RATIO_CAT, digits = 2, format = "f")
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    # IMAGE,
    STARTING_RATIO_CAT,
    TIME_ADJUSTED
  ) %>%
  mutate(
    # Use sum to get column distribution
    # Use max to scale by column max
    N = N/sum(N)*100,
    # For making phase portrait of stoichiometry by time
    DELTA_RATIO_MAGNITUDE = DELTA_RATIO,
    DELTA_RATIO = -90*DELTA_RATIO+180
  ) %>%
  mutate(
    # For making phase portrait of stoichiometry by time
    N_SCALED = (N-min(N)) / (max(N) - min(N))*100,
    DELTA_RATIO = -DELTA_RATIO+180
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    # IMAGE
  ) %>%
  # mutate(
  #   IMAGENUMBER = ifelse(IMAGE== "All", "All Combined", paste("Replicate", cur_group_id() ))
  # ) %>%
  mutate(
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS),
    FACET = paste("Init. Fract", STARTING_RATIO_CAT)
    # FACET = paste(IMAGENUMBER, "\nInit. Fract", STARTING_RATIO_CAT)
  ) %>% 
  filter(
    # Focus starting groups
    # STARTING_RATIO_CAT >= .4,
    # STARTING_RATIO_CAT <= .8
    # !STARTING_RATIO_CAT %in% c(0, 1)
  ) %>%
  as.data.table()

TracksTable$STARTING_RATIO_CAT = paste("Init. Fract.", TracksTable$STARTING_RATIO_CAT)

Plots <- unique(TracksTable$PLOT_FACETS)

HeatmapFx <- function(PlotX){
  
  # Get table
  TempTracksTable <-
    TracksTable %>% 
    filter(
      PLOT_FACETS == Plots[PlotX],
      # IMAGE == "All"
    ) %>% 
    as.data.table()
  
  COHORT_FOLDER = file.path(OUTPUT_DIRECTORY,
                            paste0(
                              "Cell Line ", TempTracksTable$COHORT[[1]], " - ",
                              TempTracksTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
                              TempTracksTable$FPS[[1]], " Hz"
                            ))
  
  if(!file.exists(COHORT_FOLDER)){
    dir.create(COHORT_FOLDER)
  }
  
  ggplot(
    TempTracksTable
  ) +
    geom_tile(
      aes(
        x = TIME_ADJUSTED,
        y = RATIO,
        fill = N
      )
    ) +
    geom_arrow(
      aes(
        x = TIME_ADJUSTED,
        y = RATIO,
        # dx = abs(DELTA_RATIO),
        # dy = DELTA_RATIO,
        mag = N_SCALED,
        # size = N_SCALED,
        angle = DELTA_RATIO,
        color = N_SCALED
      ),
      color = "black",
      # size = .01, # Small arrow head
      # arrow.length = .5, # Small arrow head
      lineend = "square"
    ) +
    # scale_size(range = c(.01, 1), guide = "none") +
    scale_mag(
      max = 300, # Small arrow head
      guide = 'none'
    ) +
    scale_color_distiller(
      palette = "Greys",
      direction = 1,
      # limits = c(-1, 1)
    ) +
    scale_fill_viridis(
      option = "plasma",
      limits = c(0, quantile(TracksTable$N, 0.8)[[1]]),
      na.value = "yellow"
    ) +
    labs(
      y = "MyD88 Fraction\n(MyD88/Cluster Size)",
      x = "Cluster Time (s)",
      color  = "Cluster Fraction\nby Time",
      fill  = "Cluster Fraction\nby Time"
    ) +
    scale_y_continuous(
      labels = percent,
      limits = c(-0.1, 1.1)
    ) +
    scale_x_continuous(
      limits = c(-20, max(TracksTable$TIME_ADJUSTED)+20)
    ) +
    facet_grid(
      ~STARTING_RATIO_CAT,
      # ~FACET,
      scales = "free"
      # nrow = NROW(unique(TracksTable$IMAGE))
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom"
    )
  
  SaveName <-
    paste0(
      "StoichiometryTime ",
      "HeatBins_",
      SCALE,
      " StartBins_",
      STARTING_RATIO_CAT_BINS,
      ".pdf"
    )
  
  ggsave(
    file.path(COHORT_FOLDER, SaveName),
    height = 2.75,
    width = 2.75*3
  )
  
}
lapply(1:NROW(Plots), HeatmapFx)






# 
# MedianTracks <-
#   StatTable %>% 
#   arrange(
#     FRAMES_ADJUSTED
#   ) %>% 
#   group_by(
#     UNIVERSAL_TRACK_ID
#   ) %>% 
#   mutate(
#     LIFETIME = max(FRAMES_ADJUSTED) - min(FRAMES_ADJUSTED) + 1
#   ) %>%
#   filter(
#     LIFETIME >= 11
#   ) %>%
#   mutate(
#     RATIO = REFERENCE_TOTAL_INTENSITY/(REFERENCE_TOTAL_INTENSITY+QUERY_TOTAL_INTENSITY)
#   ) %>% 
#   group_by(
#     IMAGE,
#     UNIVERSAL_TRACK_ID
#   ) %>% 
#   summarize(
#     RATIO = median(RATIO),
#     LIFETIME = max(LIFETIME)
#   ) %>% 
#   group_by(
#     IMAGE
#   ) %>% 
#   filter(
#     RATIO >= quantile(RATIO, .01)[[1]],
#     RATIO <= quantile(RATIO, .95)[[1]]
#   ) %>%
#   as.data.table()
# 
# ggplot(
#   MedianTracks,
#   aes(
#     x = LIFETIME,
#     y = RATIO
#   )
# ) +
#   geom_hex(
#     aes(
#       fill = ..ncount..
#     ),
#     bins = 5
#   )+
#   geom_jitter(
#     alpha = .25,
#     size = .25
#   ) +
#   scale_fill_viridis(
#     option = "plasma"
#   ) +
#   facet_wrap(
#     ~IMAGE
#   ) +
#   theme_classic()




PortraitTracks <-
  StatTable %>%
  filter(
    # To mirror phase portrait
    FRAMES_SINCE_LANDING <= FRAMES_SINCE_LANDING_BIN,
    # Filter out unusually large clusters
    REFERENCE_TOTAL_INTENSITY <= 2,
    REFERENCE_TOTAL_INTENSITY >= 0,
    QUERY_TOTAL_INTENSITY <= 2,
    QUERY_TOTAL_INTENSITY >= 0,
    COHORT == "MyD88 TRAF6",
    LIGAND_DENSITY_CAT == 32,
    IMAGE == "All"
  ) %>% 
  arrange(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = max(FRAMES_ADJUSTED) - min(FRAMES_ADJUSTED) + 1,
    RATIO = REFERENCE_TOTAL_INTENSITY/(REFERENCE_TOTAL_INTENSITY + QUERY_TOTAL_INTENSITY)
  ) %>% 
  filter(
    LIFETIME >= (LEAD_LAG*2+1)
  ) %>% 
  mutate(
    STARTING_RATIO = median(RATIO[1:5], na.rm = T),
    ENDING_RATIO =  median(RATIO[n()-5:n()], na.rm = T),
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  filter(
    # FRAMES_ADJUSTED == max(FRAMES_ADJUSTED)
    FRAMES_ADJUSTED == 75
  ) %>%
  mutate(
    STARTING_RATIO_CAT = round(STARTING_RATIO*STARTING_RATIO_CAT_BINS)/STARTING_RATIO_CAT_BINS
  ) %>% 
  arrange(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    IMAGE
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    IMAGE
  ) %>%
  mutate(
    IMAGENUMBER = ifelse(IMAGE== "All", 999, cur_group_id())
  ) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS
  ) %>%
  mutate(
    IMAGENUMBER = IMAGENUMBER-min(IMAGENUMBER) +1
  ) %>% 
  mutate(
    STARTING_RATIO_CAT_FACET = ifelse(STARTING_RATIO_CAT<=0.4, "TRAF6\nRich", ifelse(STARTING_RATIO_CAT>=0.8, "MyD88\nRich", "Goldilock\nZone")),
    # STARTING_RATIO_CAT = ifelse(STARTING_RATIO_CAT<=0.4, .4, ifelse(STARTING_RATIO_CAT>=0.8, 0.8, 0.6))
  ) %>%
  mutate(
    STARTING_RATIO_CAT_FACET = as.factor(STARTING_RATIO_CAT_FACET)
  ) %>% 
  as_tibble()

PortraitTracks$IMAGENUMBER = ifelse(PortraitTracks$IMAGENUMBER == 999, "All Combined", paste("Replicate", PortraitTracks$IMAGENUMBER))

X_MAX = max(PortraitTracks$REFERENCE_TOTAL_INTENSITY)
Y_MAX = max(PortraitTracks$QUERY_TOTAL_INTENSITY)


STARTING_RATIO_CAT_BINS = 10
# plot(Test$STARTING_RATIO, Test$RATIO)
SamplingLoop <- function(LoopX){
  Test <-
    PortraitTracks %>% 
    drop_na(
      ENDING_RATIO
    ) %>% 
    mutate(
      STARTING_RATIO_CAT = round(STARTING_RATIO*STARTING_RATIO_CAT_BINS)/STARTING_RATIO_CAT_BINS,
      FRAMES_ADJUSTED = round(FRAMES_ADJUSTED/50)*50
    ) %>% 
    group_by(
      STARTING_RATIO_CAT,
      FRAMES_ADJUSTED
    ) %>% 
    mutate(
      N = n()
    ) %>% 
    filter(
      N >=LoopX
    ) %>% 
    group_by(
      STARTING_RATIO_CAT,
      FRAMES_ADJUSTED
    ) %>% 
    sample_n(
      LoopX
    ) %>% 
    summarize(
      Q1 = quantile(ENDING_RATIO, .25),
      Q3 = quantile(ENDING_RATIO, .75),
      SPEARMAN_CORRELATION = cor(STARTING_RATIO, ENDING_RATIO, method = "spearman"),
      PEARSON_CORRELATION = cor(STARTING_RATIO, ENDING_RATIO, method = "pearson"),
      RATIO_VARIANCE = var(ENDING_RATIO),
      N = n()
    ) %>% 
    as_tibble()
  return(Test)
}
Data <- mclapply(rep(50, 50), SamplingLoop)
Data <- rbindlist(Data)

Data <-
  Data %>% 
  group_by(
    FRAMES_ADJUSTED,
    STARTING_RATIO_CAT
  ) %>% 
  summarise_all(
    median
  ) %>%
  as.data.frame()

ggplot() +
  geom_path(
    data = Data %>% mutate(FACET = "Variance") %>% as_tibble(),
    aes(
      x = STARTING_RATIO_CAT,
      y = RATIO_VARIANCE,
      # color = as.factor(FRAMES_ADJUSTED),
      # group = as.factor(FRAMES_ADJUSTED)
    )
  ) +
  geom_path(
    data = Data %>% mutate(FACET = "Spearman\nCorrelation") %>% as_tibble(),
    aes(
      x = STARTING_RATIO_CAT,
      y = SPEARMAN_CORRELATION,
      # y = RATIO_VARIANCE,
      # color = as.factor(FRAMES_ADJUSTED),
      # group = as.factor(FRAMES_ADJUSTED)
    )
  ) +
  geom_path(
    data = Data %>% mutate(FACET = "Pearson\nCorrelaton") %>% as_tibble(),
    aes(
      x = STARTING_RATIO_CAT,
      y = PEARSON_CORRELATION,
      # y = RATIO_VARIANCE,
      # color = as.factor(FRAMES_ADJUSTED),
      # group = as.factor(FRAMES_ADJUSTED)
    )
  ) +
  scale_color_brewer(
    palette = "Set1"
  ) +
  facet_wrap(
    ~FACET,
    scales = "free_y"
  ) +
  labs(
    x = "Starting Stoichiometry",
    y = "Measurement",
    color = "Lifetime Bin",
    title = "50 random samples with N=50 per bin"
  ) +
  theme_classic()


PortraitTracksSummary <-
  PortraitTracks %>% 
  mutate(
    STARTING_RATIO_CAT = ifelse(STARTING_RATIO<=0.55, "TRAF6\nRich", ifelse(STARTING_RATIO>=0.7, "MyD88\nRich", "Goldilock\nZone"))
  ) %>% 
  mutate(
    ENDING_RATIO_FACET = ifelse(RATIO<=0.55, "TRAF6\nRich", ifelse(RATIO>=0.7, "MyD88\nRich", "Goldilock\nZone")),
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    IMAGE,
    STARTING_RATIO_CAT_FACET,
    ENDING_RATIO_FACET
  ) %>% 
  summarize(
    # ERROR = mad(RATIO),
    # RATIO = median(RATIO),
    N = n()
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    IMAGE,
    STARTING_RATIO_CAT_FACET
  ) %>% 
  mutate(
    SCALED_N = N/sum(N)
  ) %>% 
  mutate(
    SCALED_N = round(SCALED_N, 2)
  )

ggplot(
  data = PortraitTracksSummary,
  aes(
    x = STARTING_RATIO_CAT_FACET,
    y = ENDING_RATIO_FACET,
    fill = N,
    label = SCALED_N
  )
) +
  geom_tile() +
  geom_label() +
  scale_fill_viridis(
    option = "plasma"
  ) +
  labs(
    x = "Start",
    y = "End",
    fill = "N"
  ) +
  theme_classic()


PortraitTracks %>% select(COHORT, IMAGE, STARTING_RATIO_CAT_FACET) %>% distinct() %>% as_tibble() %>% 
ggplot() +
  # geom_hex(
  #   data = PortraitTracks,
  #   aes(
  #     x = REFERENCE_TOTAL_INTENSITY,
  #     y = QUERY_TOTAL_INTENSITY,
  #     fill = ..ncount..
  #   ),
  #   bins = 10
  # ) +
  geom_abline(
    intercept = 0, slope = (10-6)/6,
    color = "black"
  ) +
  geom_abline(
    intercept = 0, slope = 1,
    color = "black"
    # linetype = "dashed"
  ) +
  geom_label(
    aes(
      x = X_MAX*0.8,
      y = Y_MAX*0.8,
      label = "Equal\nMix"
    ),
    color = "white",
    fill = "black",
    fontface = "bold"
  ) +
  geom_label(
    aes(
      x = X_MAX*0.8,
      y = Y_MAX*0.5,
      label = "Goldilock\nZone",
      color = "Goldilock\nZone"
    ),
    fill = "darkgrey",
    fontface = "bold"
  ) +
  geom_label(
    aes(
      x = X_MAX*0.8,
      y = Y_MAX*0.2,
      label = "MyD88\nRich",
      color = "MyD88\nRich"
    ),
    fill = "darkgrey",
    fontface = "bold"
  ) +
  geom_label(
    aes(
      x = X_MAX*0.2,
      y = Y_MAX*0.8,
      label = "TRAF6\nRich",
      color = "TRAF6\nRich"
    ),
    fill = "darkgrey",
    fontface = "bold"
  ) +
  geom_jitter(
    data = PortraitTracks,
    aes(
      x = REFERENCE_TOTAL_INTENSITY,
      y = QUERY_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID,
      color = STARTING_RATIO_CAT_FACET,
      size = FRAMES_ADJUSTED
    ),
    alpha = 0.5,
    size = 2
  ) +
  scale_fill_distiller(
    palette = "Greys",
    trans = "log1p"
  ) +
  scale_color_brewer(
    palette = "Set1"
  ) +
  labs(
    x = "MyD88 Scaled Int.",
    y = "TRAF6 Scaled Int.",
    color = "Starting\nStoichiometry",
    title = "Complex Final Size"
  ) +
  facet_grid(
    COHORT~IMAGE
  ) +
  theme_classic() +
  coord_fixed()


