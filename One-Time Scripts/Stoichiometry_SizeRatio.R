# Set scaling of plots
TOTAL_INTENSITY_SCALE = 5
RATIO_SCALE = .05

# Import table
TOTAL_INTENSITY_SCALE = 1/TOTAL_INTENSITY_SCALE
RATIO_SCALE = 1/RATIO_SCALE

TracksTable <-
  NormStatTable %>% 
  filter(
    # FRAMES_SINCE_LANDING_CAT == max(FRAMES_SINCE_LANDING_CAT)
  ) %>% 
  mutate(
    # FRAMES_SINCE_LANDING_CAT = FRAMES_SINCE_LANDING_CAT/FPS,
    # Get track parameters
    TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY+REFERENCE_TOTAL_INTENSITY,
  ) %>%
  filter(
    TOTAL_INTENSITY >= 2,
    # REFERENCE_TOTAL_INTENSITY >= 1 |
    # QUERY_TOTAL_INTENSITY >= 1
  ) %>% 
  mutate(
    RATIO = REFERENCE_TOTAL_INTENSITY/(TOTAL_INTENSITY)
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
    FRAMES_SINCE_LANDING_CAT,
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
    N >= 100
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
  drop_na(
    FRAMES_SINCE_LANDING_CAT
  ) %>% 
  group_by(
    FRAMES_SINCE_LANDING_CAT,
    COHORT,
    LIGAND_DENSITY_CAT,
    ROUNDED_TOTAL_INTENSITY
  ) %>%
  mutate(
    # N = N/max(N)*100
  ) %>%
  as.data.table()


StreamTempTable <-
  TracksTable %>%
  drop_na(
    FRAMES_SINCE_LANDING_CAT
  ) %>%
  as_tibble() %>%
  arrange(
    FRAMES_SINCE_LANDING_CAT,
    COHORT,
    LIGAND_DENSITY_CAT,
    FPS,
    ROUNDED_TOTAL_INTENSITY,
    ROUNDED_RATIO
  ) %>%
  group_by(
    FRAMES_SINCE_LANDING_CAT,
    COHORT,
    LIGAND_DENSITY_CAT
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
  ) %>%
  as.data.table()

# Add derivatives
StreamTempTable$DELTA_TOTAL_INTENSITY[is.na(StreamTempTable$DELTA_TOTAL_INTENSITY)] = 0
StreamTempTable$DELTA_RATIO[is.na(StreamTempTable$DELTA_RATIO)] = 0
StreamTempTable$N[is.na(StreamTempTable$N)] = 100


ggplot(
  data = TracksTable
) +
  geom_tile(
    aes(
      x = ROUNDED_TOTAL_INTENSITY,
      y = ROUNDED_RATIO,
      fill = DELTA_TOTAL_INTENSITY
    )
  ) +
  geom_hline(
    yintercept = .2,
    color = "green"
  ) +
  # geom_streamline(
  #   data = StreamTempTable,
  #   aes(
  #     x = ROUNDED_TOTAL_INTENSITY,
  #     y = ROUNDED_RATIO,
  #     dx = DELTA_TOTAL_INTENSITY/max(DELTA_TOTAL_INTENSITY),
  #     dy = DELTA_RATIO,
  #     # color = sqrt(..dx..^2 + ..dy..^2),
  #     # size = ..step..
  #     # alpha = ..step..
  #   ),
  #   color = "black",
  #   arrow = NULL,
  #   n = 8,
  #   size = .1,
  #   # arrow.length = 0.3,
  #   jitter = 4,
  #   L = 2, res = 1, lineend = "round"
  # ) +
  geom_arrow(
    aes(
      x = ROUNDED_TOTAL_INTENSITY,
      y = ROUNDED_RATIO,
      dx = DELTA_TOTAL_INTENSITY,
      dy = DELTA_RATIO
      
      # mag = 0.5,
      # angle = atan2(DELTA_RATIO, DELTA_TOTAL_INTENSITY)*180/pi,
      # color = sqrt(DELTA_RATIO^2 + DELTA_TOTAL_INTENSITY^2),
      # color = FRAMES_ADJUSTED
    ),
    color = "black",
    # size = .01, # Small arrow head
    # arrow.length = .5, # Small arrow head
    # size = 1, # Big arrow head
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
  scale_color_distiller(
    palette = "Greys"
  ) +
  scale_fill_viridis(
    option = "plasma",
    # guide = "none"
    # direction = -1,
    # trans = "log10",
    # breaks = trans_breaks("log10", function(x) 10^x),
    # labels = trans_format("log10", math_format(10^.x))
  )  +
  labs(
    x = "Cluster Size (TRAF6 + NEMO)",
    y = "Stoichometric Ratio\nTRAF6 / (TRAF6 + NEMO)",
    color = "Magnitude",
    # fill = "Fractional\nN"
  ) +
  facet_wrap(
    ~FRAMES_SINCE_LANDING_CAT
  ) +
  theme_classic() #+
# coord_fixed()
