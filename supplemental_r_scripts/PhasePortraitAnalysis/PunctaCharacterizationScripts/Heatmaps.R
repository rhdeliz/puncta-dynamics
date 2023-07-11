# Heatmap
# Table <- fread(file.path(OUTPUT_DIRECTORY, "Essential.csv.gz"))

FilteredTable <-
  Table %>% 
  group_by(
    IMAGE,
    CELL,
    PROTEIN
  ) %>% 
  mutate(
    FRAME_LAND = FRAME - min(FRAME)
  ) %>% 
  ungroup() %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME_FRAMES = max(FRAMES_ADJUSTED),
    FIRST_FRAME = min(FRAME_LAND)
  ) %>% 
  filter(
    LIFETIME_FRAMES >= 50,
    FRAMES_ADJUSTED <= 100,
    # FIRST_FRAME <= 1000
    # PROTEIN == "TRAF6"
  ) %>% 
  ungroup() %>% 
  arrange(
    FIRST_FRAME
  ) %>% 
  mutate(
    ID = group_indices(., IMAGE, PROTEIN, FIRST_FRAME, LIFETIME_FRAMES, UNIVERSAL_TRACK_ID)
    # ID = group_indices(., IMAGE, CELL, PROTEIN, FIRST_FRAME, LIFETIME_FRAMES, UNIVERSAL_TRACK_ID)
  ) %>%
  group_by(
    IMAGE, PROTEIN,
    # CELL
  ) %>% 
  mutate(
    ID = ID - min(ID)
  ) %>% 
  group_by(
    IMAGE,
    CELL,
    PROTEIN
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/max(NORMALIZED_INTENSITY)
  ) %>% 
  group_by(
    IMAGE, PROTEIN,
    # CELL
  ) %>% 
  mutate(
    N = NROW(unique(ID)),
    PROTEIN = factor(PROTEIN, levels = c("MyD88", PROTEIN_ORDER)),
    COHORT = factor(COHORT, levels = paste("MyD88", PROTEIN_ORDER))
  ) %>% 
  filter(
    N >= 10
  )


ggplot(
  FilteredTable# %>% filter(COHORT == "MyD88 TRAF6", PROTEIN == "MyD88")
) +
  geom_tile(
    aes(
      x = FRAMES_ADJUSTED,
      y = ID,
      fill = NORMALIZED_INTENSITY
    )
  ) +
  scale_fill_viridis(
   option = "inferno"
  ) +
  facet_wrap(
    ~COHORT+PROTEIN, scales = "free_y",
    ncol = 2
  ) +
  labs(
    x = "Time (Frames)",
    y = "Tracks",
    fill = "Scaled\nNorm. Int."
  ) +
  dark_theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  "heatmaps.png",
  height = 4.56,
  width = 11.5
)
