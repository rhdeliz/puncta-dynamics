IntensityByTime <-
  StatTable %>% 
  filter(
    MAX_REFERENCE_TOTAL_INTENSITY_CAT == "≥4.5x MyD88"
  ) %>% 
  mutate(
    TIME_ADJUSTED = round(TIME_ADJUSTED),
    FRAMES_SINCE_LANDING_CAT = round(FRAMES_SINCE_LANDING_CAT/100)*100
  ) %>% 
  group_by(
    IMAGE,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    TIME_ADJUSTED,
    FRAMES_SINCE_LANDING_CAT
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    # TIME_ADJUSTED <= 300,
    N >= 5
  )

ggplot(
  IntensityByTime
) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  scale_color_manual(
    values = c("magenta", "green", "pink", "purple")
  ) +
  labs(
    x = "Time (s)",
    y = "Size",
    color = "Protein"
  ) +
  dark_theme_classic(
    base_size = 20
  ) +
  facet_grid(
    QUERY_PROTEIN~FRAMES_SINCE_LANDING_CAT,
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(OUTPUT_DIRECTORY, "IntensityTimeLanding.pdf"),
  height = 4.76,
  width = 11.5
)
