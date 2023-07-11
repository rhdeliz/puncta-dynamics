NormStatTable %>% 
  group_by(
    QUERY_PROTEIN
  ) %>% 
  filter(
    REFERENCE_TOTAL_INTENSITY >= STEP_SIZE,
    QUERY_TOTAL_INTENSITY >= STEP_SIZE
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = round(DELTA_REFERENCE_TOTAL_INTENSITY, 2),
    # DELTA_QUERY_TOTAL_INTENSITY = round(DELTA_QUERY_TOTAL_INTENSITY, 2),
    # QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY)
  ) %>% 
  group_by(
    QUERY_PROTEIN,
    # QUERY_TOTAL_INTENSITY,
    # DELTA_QUERY_TOTAL_INTENSITY,
    DELTA_REFERENCE_TOTAL_INTENSITY,
    FRAMES_SINCE_LANDING_CAT
  ) %>% 
  summarize(
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY),
    # DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    N >= 100,
    # DELTA_QUERY_TOTAL_INTENSITY > 0,
    # FRAMES_SINCE_LANDING_CAT > 100
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      # QUERY_TOTAL_INTENSITY,
      DELTA_REFERENCE_TOTAL_INTENSITY,
      DELTA_QUERY_TOTAL_INTENSITY,
      # DELTA_REFERENCE_TOTAL_INTENSITY,
      group = FRAMES_SINCE_LANDING_CAT,
      color = FRAMES_SINCE_LANDING_CAT
    )
  ) +
  geom_path() +
  geom_hline(
    yintercept = 0
  ) +
  scale_color_viridis() +
  facet_wrap(
    ~QUERY_PROTEIN,
    scales = "free"
  ) +
  theme_classic()
