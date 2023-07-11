# Set scaling of plots
SCALE = 5
# Set number of facets
STARTING_RATIO_CAT_BINS = 3


setwd(OUTPUT_DIRECTORY)

# Table Name
SaveName <-
  paste0(
    "NormStatTable - ",
    "LeadLag ", LEAD_LAG,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)
NormStatTable <- read_parquet(file.path(OUTPUT_DIRECTORY, SaveName))


# Import table
TracksTable <-
  NormStatTable %>% 
  mutate(
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS)
  ) %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT == 200,
    QUERY_TOTAL_INTENSITY >= 0,
    REFERENCE_TOTAL_INTENSITY >=0
    # COHORT == "MyD88 TRAF6"
  ) %>%
  arrange(
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    # Get track parameters
    FRAMES_SINCE_LANDING_CAT = FRAMES_SINCE_LANDING_CAT/FPS,
    TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED),
    RATIO = REFERENCE_TOTAL_INTENSITY/(QUERY_TOTAL_INTENSITY+REFERENCE_TOTAL_INTENSITY)
  ) %>%
  mutate(
    # Average ratios
    STARTING_RATIO = median(RATIO[1:SCALE], na.rm = T),
    ENDING_RATIO =  median(RATIO[n()-5:n()], na.rm = T),
    # Get ± ratio
    LEAD_RATIO = lead(RATIO, LEAD_LAG),
    LAG_RATIO = lag(RATIO, LEAD_LAG)
  ) %>% 
  mutate(
    STARTING_RATIO_CAT = round(STARTING_RATIO, 1)
  ) %>% 
  as.data.table()


TrackBoundaries <-
  TracksTable %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  filter(
    FRAMES_ADJUSTED == 5
  ) %>% 
  drop_na(
    STARTING_RATIO, ENDING_RATIO
  ) %>% 
  group_by(
    PLOT_FACETS,
    STARTING_RATIO_CAT
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  filter(
    N >= 50
  ) %>% 
  as.data.table()

SampleFx <- function(x){
  
  TempBoundaries <-
    TrackBoundaries %>% 
    group_by(
      PLOT_FACETS,
      STARTING_RATIO_CAT
    ) %>% 
    sample_n(
      50
    ) %>% 
    summarize(
      COR = cor(STARTING_RATIO, ENDING_RATIO),
      VAR = sd(ENDING_RATIO),
      N = n()
    ) %>% 
    as.data.table()
  
  return(TempBoundaries)
  
}
Boundaries <- mclapply(1:101, SampleFx)
Boundaries <- rbindlist(Boundaries)

Boundaries <-
  Boundaries %>% 
  group_by(
    PLOT_FACETS,
    STARTING_RATIO_CAT
  ) %>% 
  summarize(
    COR = median(COR, na.rm = T),
    VAR = median(VAR)
  ) %>% 
  drop_na(
    COR, VAR
  ) %>% 
  group_by(
    PLOT_FACETS
  ) %>% 
  mutate(
    SCALE_COR = COR-min(COR),
    SCALE_VAR = VAR-min(VAR)
  ) %>% 
  mutate(
    SCALE_COR = SCALE_COR/max(SCALE_COR),
    SCALE_VAR = SCALE_VAR/max(SCALE_VAR)
  ) %>% 
  mutate(
    SCALE_COR = SCALE_COR
  ) %>% 
  mutate(
    SPECIAL = (SCALE_VAR + 1- SCALE_COR)/2
  ) %>% 
  mutate(
    SCALE_COR = signal::sgolayfilt(SCALE_COR, p = 1, n = 3),
    SCALE_VAR = signal::sgolayfilt(SCALE_VAR, p = 1, n = 3),
    SPECIAL = signal::sgolayfilt(SPECIAL, p = 1, n = 3)
  ) %>%
  drop_na(
    SPECIAL
  ) %>%
  group_by(
    PLOT_FACETS
  ) %>% 
  mutate(
    THRESHOLD = min(case_when(SPECIAL == min(SPECIAL) ~ STARTING_RATIO_CAT), na.rm = T)
  ) %>% 
  as.data.table()


ggplot(
  Boundaries
) +
  geom_line(
    aes(
      STARTING_RATIO_CAT, SCALE_COR
    ),
    color = "red"
  ) +
  geom_line(
    aes(
      STARTING_RATIO_CAT, SCALE_VAR
    ),
    color = "blue"
  ) +
  geom_line(
    aes(
      STARTING_RATIO_CAT, SPECIAL
    ),
    color = "green"
  ) +
  geom_vline(
    aes(
      xintercept = THRESHOLD
    )
  ) +
  # geom_point(
  #   aes(
  #     STARTING_RATIO_CAT, SCALE_VAR
  #   ),
  #   color = "blue"
  # ) +
  facet_wrap(
    ~PLOT_FACETS
  ) +
  theme_classic()
  





TracksTable <-
  TracksTable %>% 
  mutate(
    # Bin data
    TIME_ADJUSTED = round(TIME_ADJUSTED/(4*SCALE))*(4*SCALE),
    STARTING_RATIO_CAT = round(STARTING_RATIO_CAT*STARTING_RATIO_CAT_BINS)/STARTING_RATIO_CAT_BINS,
    RATIO = round(RATIO*(50/SCALE))/(50/SCALE),
    # Calculate ratio difference
    DELTA_RATIO = LEAD_RATIO - LAG_RATIO
  ) %>% 
  drop_na(
    DELTA_RATIO,
    STARTING_RATIO_CAT
  ) %>% 
  mutate(
    # STARTING_RATIO_CAT = 0
  ) %>%
  filter(
    FRAMES_SINCE_LANDING_CAT == max(FRAMES_SINCE_LANDING_CAT)
  ) %>% 
  group_by(
    PLOT_FACETS,
    
    FRAMES_SINCE_LANDING_CAT,
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
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
  as.data.table()

TracksTable$STARTING_RATIO_CAT = paste("Init. Fract.", TracksTable$STARTING_RATIO_CAT)

Plots <- unique(TracksTable$PLOT_FACET)

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
    geom_hline(
      yintercept = .65,
      color = "green"
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
        # mag = N_SCALED,
        mag = .5,
        # size = N_SCALED,
        angle = DELTA_RATIO,
        color = N_SCALED
      ),
      color = "black",
      # size = .01, # Small arrow head
      # arrow.length = .5, # Small arrow head
      lineend = "square"
    ) +
    # scale_size(range = c(.1, 1), guide = "none") +
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
      y = "Fractional Reference\n(Reference/Cluster Size)",
      x = "Spot Time (s)",
      color  = "Scaled Magnitude",
      fill  = "Fractional Clusters\nby Spot Time"
    ) +
    scale_y_continuous(
      limits = c(-0.1, 1.1)
    ) +
    scale_x_continuous(
      limits = c(-20, max(TracksTable$TIME_ADJUSTED)+20)
    ) +
    facet_wrap(
    # facet_grid(
      ~STARTING_RATIO_CAT,
      # STARTING_RATIO_CAT~FRAMES_SINCE_LANDING_CAT,
      # ~FACET,
      # scales = "free"
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
    height = 9,
    width = 16
  )
  
}
mclapply(1:NROW(Plots), HeatmapFx)
