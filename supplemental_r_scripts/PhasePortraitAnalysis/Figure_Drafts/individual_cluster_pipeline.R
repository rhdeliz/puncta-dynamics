          setwd(OUTPUT_DIRECTORY)

THRESHOLD = 0.5

ClusterTime <- 
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25,
    CELL == 5
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    N = 1:n()
  ) %>% 
  filter(
    N == 1
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    N = n()
  ) %>% 
  filter(
    N > 16
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY =
      REFERENCE_TOTAL_INTENSITY*median(MAX_REFERENCE) + median(MIN_REFERENCE),
    QUERY_TOTAL_INTENSITY =
      QUERY_TOTAL_INTENSITY*median(MAX_QUERY) + median(MIN_QUERY),
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    # REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    # QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>%
  mutate(
    # REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    # QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  mutate(
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  as.data.table()

ClusterTime %>% 
  arrange(
    TIME_ADJUSTED
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID
    )
  ) +
  geom_path(
    # size = .1,
    alpha = .25
  ) +
  theme_void()

ggsave(
  "all_tracks.pdf",
  height = 0.5,
  width = 1.5
)


ClusterTime <- 
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25,
    CELL == 5
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    N = 1:n()
  ) %>% 
  filter(
    N == 1
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    N = n()
  ) %>% 
  filter(
    N > 16
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY =
      REFERENCE_TOTAL_INTENSITY*median(MAX_REFERENCE) + median(MIN_REFERENCE),
    QUERY_TOTAL_INTENSITY =
      QUERY_TOTAL_INTENSITY*median(MAX_QUERY) + median(MIN_QUERY),
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  mutate(
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  as.data.table()

ClusterTime %>% 
  arrange(
    TIME_ADJUSTED
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID
    )
  ) +
  geom_path(
    # size = .1,
    alpha = .25
  ) +
  theme_void()

ggsave(
  "smooth_tracks.pdf",
  height = 0.5,
  width = 1.5
)


ClusterTime <- 
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25,
    CELL == 5
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    N = 1:n()
  ) %>% 
  filter(
    N == 1
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    N = n()
  ) %>% 
  filter(
    N > 16
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY =
      REFERENCE_TOTAL_INTENSITY*median(MAX_REFERENCE) + median(MIN_REFERENCE),
    QUERY_TOTAL_INTENSITY =
      QUERY_TOTAL_INTENSITY*median(MAX_QUERY) + median(MIN_QUERY),
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  mutate(
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  filter(
    MAX_REFERENCE_TOTAL_INTENSITY >= 3,
    MAX_QUERY_TOTAL_INTENSITY >= 3
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  as.data.table()

ClusterTime %>% 
  arrange(
    TIME_ADJUSTED
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID
    )
  ) +
  geom_hline(
    yintercept = 3,
    color = "red"
  ) +
  geom_path(
    alpha = .25
  ) +
  theme_void()

ggsave(
  "grow_tracks.pdf",
  height = 0.5,
  width = 1.5
)


ClusterTime <- 
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25,
    CELL == 5
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    N = 1:n()
  ) %>% 
  filter(
    N == 1
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    N = n()
  ) %>% 
  filter(
    N > 16
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY =
      REFERENCE_TOTAL_INTENSITY*median(MAX_REFERENCE) + median(MIN_REFERENCE),
    QUERY_TOTAL_INTENSITY =
      QUERY_TOTAL_INTENSITY*median(MAX_QUERY) + median(MIN_QUERY),
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  mutate(
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  filter(
    MAX_REFERENCE_TOTAL_INTENSITY >= 3,
    MAX_QUERY_TOTAL_INTENSITY >= 3
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/MAX_REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/MAX_QUERY_TOTAL_INTENSITY
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  mutate(
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY)
  ) %>% 
  as.data.table()

ClusterTime %>% 
  arrange(
    TIME_ADJUSTED
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID
    )
  ) +
  geom_path(
    alpha = .25
  ) +
  theme_void()

ggsave(
  "scaled.pdf",
  height = 0.5,
  width = 1.5
)

ClusterTime <- 
  NormStatTable %>% 
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25,
    CELL == 5
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  mutate(
    N = 1:n()
  ) %>% 
  filter(
    N == 1
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    N = n()
  ) %>% 
  filter(
    N > 16
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY =
      REFERENCE_TOTAL_INTENSITY*median(MAX_REFERENCE) + median(MIN_REFERENCE),
    QUERY_TOTAL_INTENSITY =
      QUERY_TOTAL_INTENSITY*median(MAX_QUERY) + median(MIN_QUERY),
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  mutate(
    MAX_REFERENCE_TOTAL_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    MAX_QUERY_TOTAL_INTENSITY = max(QUERY_TOTAL_INTENSITY),
  ) %>% 
  filter(
    MAX_REFERENCE_TOTAL_INTENSITY >= 3,
    MAX_QUERY_TOTAL_INTENSITY >= 3
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/MAX_REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/MAX_QUERY_TOTAL_INTENSITY
  ) %>%
  # Determine if it crosses threshold
  mutate(
    REFERENCE_TOTAL_INTENSITY_50 = ifelse(REFERENCE_TOTAL_INTENSITY >= THRESHOLD, 1, 0),
    QUERY_TOTAL_INTENSITY_50 = ifelse(QUERY_TOTAL_INTENSITY >= THRESHOLD, 1, 0)
  ) %>% 
  # Look at previous one to see if it's closer to threhsold
  group_by(
    COHORT,
    QUERY_PROTEIN,
    REFERENCE_PROTEIN,
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY_MID =  min(case_when(REFERENCE_TOTAL_INTENSITY_50 == 1 ~ REFERENCE_TOTAL_INTENSITY), na.rm = T),
    QUERY_TOTAL_INTENSITY_MID =  min(case_when(QUERY_TOTAL_INTENSITY_50 == 1 ~ QUERY_TOTAL_INTENSITY), na.rm = T),
    REFERENCE_TIME_ADJUSTED = min(case_when(REFERENCE_TOTAL_INTENSITY_50 == 1 ~ TIME_ADJUSTED), na.rm = T),
    QUERY_TIME_ADJUSTED = min(case_when(QUERY_TOTAL_INTENSITY_50 == 1 ~ TIME_ADJUSTED), na.rm = T),
  ) %>%
  mutate(
    REFERENCE_FLAG = ifelse(REFERENCE_TOTAL_INTENSITY_MID == REFERENCE_TOTAL_INTENSITY, 1, 0)
  ) %>% 
  as.data.table()

ClusterTime %>% 
  arrange(
    TIME_ADJUSTED
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID
    )
  ) +
  geom_hline(
    yintercept = 1,
    color = "red"
  ) +
  geom_point(
    aes(
      REFERENCE_TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY_MID,
      group = UNIVERSAL_TRACK_ID
    )
  ) + 
  geom_path(
    # size = .1,
    alpha = .25
  ) +
  theme_void()

ggsave(
  "mid_tracks.pdf",
  height = 0.5,
  width = 1.5
)


ClusterTime %>% 
  mutate(
    FRAMES_ADJUSTED=  FRAMES_ADJUSTED -min(FRAMES_ADJUSTED)
  ) %>% 
  filter(
    FRAMES_ADJUSTED == 0
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    DELTA_TIME = QUERY_TIME_ADJUSTED - REFERENCE_TIME_ADJUSTED,
    MIN_TIME = min(REFERENCE_TIME_ADJUSTED, QUERY_TIME_ADJUSTED)
  ) %>% 
  ungroup() %>% 
  arrange(
    MIN_TIME,
    DELTA_TIME
  ) %>% 
  mutate(
    ID = 1:n()
  ) %>% 
  as.data.table() %>% 
  ggplot(
    aes(
      TIME_ADJUSTED,
      ID,
      group = UNIVERSAL_TRACK_ID
    )
  ) +
  geom_segment(
    aes(
      x = REFERENCE_TIME_ADJUSTED,
      xend = QUERY_TIME_ADJUSTED,
      y = ID,
      yend = ID
    )
  ) +
  geom_point(
    aes(
      REFERENCE_TIME_ADJUSTED,
      ID,
      group = UNIVERSAL_TRACK_ID,
      color = REFERENCE_PROTEIN
    )
  ) + 
  geom_point(
    aes(
      QUERY_TIME_ADJUSTED,
      ID,
      group = UNIVERSAL_TRACK_ID,
      color = QUERY_PROTEIN
    )
  ) + 
  geom_path(
    alpha = .25
  ) +
  scale_color_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  theme_void()

ggsave(
  "mid_point.pdf",
  height = 0.5,
  width = 1.5
)


ClusterTime %>% 
  as.data.table() %>% 
  ggplot() +
  geom_path(
    aes(
      TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    aes(
      TIME_ADJUSTED,
      QUERY_TOTAL_INTENSITY,
      group = UNIVERSAL_TRACK_ID,
      color = QUERY_PROTEIN
    )
  ) +
  geom_segment(
    aes(
      x = REFERENCE_TIME_ADJUSTED,
      xend = QUERY_TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY_MID,
      yend = QUERY_TOTAL_INTENSITY_MID,
      color = QUERY_PROTEIN
    ),
    size = 1.5
  ) +
  geom_point(
    aes(
      REFERENCE_TIME_ADJUSTED,
      REFERENCE_TOTAL_INTENSITY_MID,
      group = UNIVERSAL_TRACK_ID,
      color = REFERENCE_PROTEIN
    )
  ) + 
  geom_point(
    aes(
      QUERY_TIME_ADJUSTED,
      QUERY_TOTAL_INTENSITY_MID,
      group = UNIVERSAL_TRACK_ID,
      color = QUERY_PROTEIN
    )
  ) + 
  scale_color_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  facet_wrap(
    ~UNIVERSAL_TRACK_ID,
    ncol = 3
  ) +
  theme_void() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(
  "all_deltas.pdf",
  height = 0.5*3,
  width = 1.5*1.5
)





