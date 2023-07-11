IntensityByTime <-
  StatTable %>% 
  filter(
    MAX_REFERENCE_TOTAL_INTENSITY_CAT == "≥4.5x MyD88"
  ) %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT <= 100,
    ROUNDED_REFERENCE_TOTAL_INTENSITY >= 1
  ) %>% 
  mutate(
    TIME_ADJUSTED = round(TIME_ADJUSTED),
  ) %>% 
  group_by(
    IMAGE,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    ROUNDED_REFERENCE_TOTAL_INTENSITY
  ) %>% 
  summarize(
    ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = median(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
    ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = median(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    # TIME_ADJUSTED <= 300,
    N >= 5
  ) %>% 
  group_by(
    IMAGE,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  filter(
    N >= 5
  )

ggplot(
  IntensityByTime
) +
  geom_hline(
    yintercept = 0
  ) +
  geom_path(
    aes(
      x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
      y = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    aes(
      x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
      y = ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  scale_color_manual(
    values = c("magenta", "green", "pink")
  ) +
  labs(
    x = "MyD88 Size",
    y = "Size Change",
    color = "Protein"
  ) +
  dark_theme_classic(
    base_size = 20
  ) +
  facet_grid(
    ~QUERY_PROTEIN
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(OUTPUT_DIRECTORY, "AdjDeltaIntensity.pdf"),
  height = 4.76,
  width = 11.5
)
