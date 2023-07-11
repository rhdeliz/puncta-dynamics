NormStatTable %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    MAX_REFERENCE = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  filter(
    MAX_REFERENCE >= STEP_SIZE*2,
    MAX_QUERY >= STEP_SIZE
  ) %>% 
  group_by(
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    TIME_ADJUSTED = median(TIME_ADJUSTED),
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY),
    
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    N >= 175
  ) %>% 
  group_by(
    QUERY_PROTEIN
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 11),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 11),
    
    DELTA_REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(DELTA_REFERENCE_TOTAL_INTENSITY, p = 1, n = 31),
    DELTA_QUERY_TOTAL_INTENSITY = signal::sgolayfilt(DELTA_QUERY_TOTAL_INTENSITY, p = 1, n = 31)
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(abs(DELTA_QUERY_TOTAL_INTENSITY))
  ) %>% 
  as.data.table() %>% 
  ggplot() +
  geom_path(
    aes(
      TIME_ADJUSTED,
      DELTA_REFERENCE_TOTAL_INTENSITY,
      # DELTA_QUERY_TOTAL_INTENSITY,
      
      # REFERENCE_TOTAL_INTENSITY,
      # QUERY_TOTAL_INTENSITY,
      # color = TIME_ADJUSTED
    ),
    color = "black"
  ) +
  geom_path(
    aes(
      TIME_ADJUSTED,
      # DELTA_REFERENCE_TOTAL_INTENSITY,
      DELTA_QUERY_TOTAL_INTENSITY,
      
      # REFERENCE_TOTAL_INTENSITY,
      # QUERY_TOTAL_INTENSITY,
      # color = TIME_ADJUSTED
      color = QUERY_PROTEIN
    )
  ) +
  geom_hline(
    yintercept = 0
  ) +
  geom_vline(
    xintercept = 0
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    # guide = "none"
  ) +
  # scale_y_continuous(
  #   limits = c(-0.005, 0.0025)
  # ) +
  # scale_color_viridis() +
  facet_wrap(
    ~QUERY_PROTEIN,
    scales = "free"
  ) +
  theme_classic()
