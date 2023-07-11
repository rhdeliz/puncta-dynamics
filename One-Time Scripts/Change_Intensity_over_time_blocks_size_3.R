
ScaledDeltaRefTime <-
  NormStatTable %>% 
  group_by(
    COHORT,
    ROUNDED_REFERENCE_TOTAL_INTENSITY,
    ROUNDED_QUERY_TOTAL_INTENSITY
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  filter(
    N >= 100,
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  mutate(
    ROUNDED_REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY*10)/10,
    ROUNDED_QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY*10)/10,
  ) %>% 
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_SINCE_LANDING_CAT,
    ROUNDED_QUERY_TOTAL_INTENSITY,
    ROUNDED_REFERENCE_TOTAL_INTENSITY
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    ROUNDED_QUERY_TOTAL_INTENSITY >= 0.5,
    ROUNDED_REFERENCE_TOTAL_INTENSITY >= 0.5
  ) %>% 
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_SINCE_LANDING_CAT
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    N = sum(N)
  ) %>% 
  as.data.table()



ggplot(
  ScaledDeltaRefTime,
  aes(
    FRAMES_SINCE_LANDING_CAT,
    DELTA_REFERENCE_TOTAL_INTENSITY,
    color = QUERY_PROTEIN,
    group = QUERY_PROTEIN
  )
) +
  geom_point(
    aes(
      size = N
    )
  ) +
  geom_path() +
  geom_hline(
    yintercept = 0
  ) +
  # geom_vline(
  #   xintercept = 150
  # ) +
  labs(
    x = "Cell Time (s)",
    y = "d-MyD88",
    color ="Cell Line"
  ) +
  scale_x_continuous(
    # limits = c(100, 200)
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    # guide = "none"
  ) +
  # facet_wrap(
  #   ~ROUNDED_QUERY_TOTAL_INTENSITY,
  #   scales = "free_y"
  # ) +
  theme_classic()
