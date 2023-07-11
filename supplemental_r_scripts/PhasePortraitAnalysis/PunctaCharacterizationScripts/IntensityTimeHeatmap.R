SummaryFlow2 <-
  AnalysisStatTable %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT %in% c(50, 100, 200),
    !is.infinite(DELTA_REFERENCE_TOTAL_INTENSITY),
    !is.infinite(DELTA_QUERY_TOTAL_INTENSITY),
    
    !is.infinite(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
    !is.infinite(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY)
  ) %>% 
  ungroup() %>%
  mutate(
    IMAGENUMBER = group_indices(., COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGE),
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN),
    COHORT = factor(COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER))
  ) %>%
  group_by(
    PLOT_FACETS, COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    FPS = round(FPS, 2),
    IMAGENUMBER = IMAGENUMBER - min(IMAGENUMBER) + 1
  ) %>% 
  group_by(
    PLOT_FACETS, COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, FRAMES_SINCE_LANDING_CAT, IMAGENUMBER,
    ROUNDED_REFERENCE_TOTAL_INTENSITY, ROUNDED_QUERY_TOTAL_INTENSITY,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    N = n(),
    FRAMES_ADJUSTED = median(FRAMES_ADJUSTED)
  ) %>% 
  group_by(
    PLOT_FACETS,
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, FRAMES_SINCE_LANDING_CAT, IMAGENUMBER,
    ROUNDED_REFERENCE_TOTAL_INTENSITY, ROUNDED_QUERY_TOTAL_INTENSITY
  ) %>% 
  mutate(
    N_TEST = ifelse(N >= quantile(N, 0.50), T, F)
  ) %>% 
  filter(
    N_TEST == T
  )

ggplot(
) +
  geom_tile(
    data = SummaryFlow2,
    aes(
      x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
      y = ROUNDED_QUERY_TOTAL_INTENSITY,
      fill = FRAMES_ADJUSTED
    )
  ) +
  scale_fill_viridis(
  ) +
  labs(
    x = "MyD88 Norm. Intensity",
    y = "mScarlet Norm. Intensity",
    fill = "Time\n(Frames)"
  ) +
  facet_grid(
    QUERY_PROTEIN~FRAMES_SINCE_LANDING_CAT,
    # scales = "free"
  ) +
  dark_theme_classic(
    base_size = 20
  ) +
  coord_fixed()

ggsave(
  # Save vector image
  "TimeHeatmap.pdf",
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)

