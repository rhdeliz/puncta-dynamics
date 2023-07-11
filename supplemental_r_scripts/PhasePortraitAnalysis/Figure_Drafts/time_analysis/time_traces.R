setwd(OUTPUT_DIRECTORY)

IRAK4 <-
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 IRAK4",
    FRAMES_SINCE_LANDING_CAT == 150
    # !grepl("rid", COHORT)
    # REFERENCE_TOTAL_INTENSITY > 0,
    # QUERY_TOTAL_INTENSITY > 0,
    # LIGAND_DENSITY_CAT == 100
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    START_REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY[1:3]),
    START_QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY[1:3]),
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  as.data.table()

IRAK4_Summary <- 
  IRAK4 %>% 
  filter(
    # START_REFERENCE_TOTAL_INTENSITY <= 2,
    # START_QUERY_TOTAL_INTENSITY <= 2,
    # MAX_REFERENCE_TOTAL_INTENSITY >= 1,
    # MAX_QUERY_TOTAL_INTENSITY >= 3
  ) %>% 
  group_by(
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    TIME_ADJUSTED = median(TIME_ADJUSTED),
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY)
  ) %>% 
  ungroup() %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7)
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 15)
  ) %>% 
  mutate(
    REFERENCE_PROTEIN = factor(REFERENCE_PROTEIN, levels = PROTEIN_ORDER),
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = PROTEIN_ORDER)
  ) %>% 
  as.data.table()


ggplot(
  data = IRAK4_Summary
) +
  geom_path(
    aes(
      TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    aes(
      TIME_ADJUSTED,
      QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  scale_color_manual(
    values = ColorTable,
    breaks = names(ColorTable),
    limits = force
  ) +
  labs(
    x = "Cluster Time (s)",
    y =  "Norm. Int. (x)",
    color = "Protein"
  ) +
  theme_classic(
    base_size = 24
  )  +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "time_trace.pdf",
  height = 2*3,
  width = 2.5*2
)

ggplot(
  data = IRAK4_Summary
) +
  geom_path(
    aes(
      REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY,
      color = TIME_ADJUSTED
    )
  ) +
  scale_color_viridis(
    breaks = c(1:6*150)
  ) +
  labs(
    x = "(x) MyD88",
    y =  "(x) IRAK4",
    color = "Cluster\nTime (s)"
  ) +
  # coord_fixed() +
  theme_classic(
    base_size = 24
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1,"cm"),
    legend.text = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(
  "phase_space_sample_track.pdf",
  height = 2*3,
  width = 2.5*2
)

IRAK4 %>% 
  mutate(
    # REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    # QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7)
  ) %>% 
  as.data.table() %>% 
  ggplot() +
  geom_path(
    aes(
      REFERENCE_TOTAL_INTENSITY,
      QUERY_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID
    ),
    size = .1,
    alpha = .1
  ) +
  scale_x_continuous(
    limits = c(0, STEP_SIZE*10)
  ) +
  scale_y_continuous(
    limits = c(0, STEP_SIZE*2)
  ) +
  labs(
    x = "(x) MyD88",
    y =  "(x) IRAK4",
  ) +
  theme_classic(
    base_size = 18
  ) +
  coord_fixed()


ggsave(
  "phase_space_all_tracks.pdf",
  height = 2*1.25,
  width = 2.5*3
)

ggplot(
  IRAK4,
  aes(
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY
  )
) +
  stat_density2d(geom="raster", aes(fill = ..count../max(..count..)), adjust = 1, contour = FALSE) +
  scale_fill_viridis(
    limit = c(0, .05),
    na.value = viridis(10)[10],
    labels = percent,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, STEP_SIZE*10)
  ) +
  scale_y_continuous(
    limits = c(0, STEP_SIZE*2)
  ) +
  labs(
    x = "(x) MyD88",
    y =  "(x) IRAK4",
    fill = "Scaled Density"
  ) +
  theme_classic(
    base_size = 18
  ) +
  coord_fixed()


ggsave(
  "phase_space_density.pdf",
  height = 2*1.25,
  width = 2.5*3
)
