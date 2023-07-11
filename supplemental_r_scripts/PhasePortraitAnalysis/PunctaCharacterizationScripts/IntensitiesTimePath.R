SummaryFlow <-
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
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    N = n(),
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY)
  ) %>% 
  group_by(
    PLOT_FACETS,
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, FRAMES_SINCE_LANDING_CAT, IMAGENUMBER
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = roll_median(REFERENCE_TOTAL_INTENSITY, n = 5, fill = NA, align = "center"),
    QUERY_TOTAL_INTENSITY = roll_median(QUERY_TOTAL_INTENSITY, n = 5, fill = NA, align = "center"),
    N_TEST = ifelse(N >= quantile(N, 0.50), T, F)
  ) %>% 
  filter(
    N_TEST == T
  )

ggplot(
) +
  geom_path(
    data = SummaryFlow,
    aes(
      x = REFERENCE_TOTAL_INTENSITY,
      y = QUERY_TOTAL_INTENSITY,
      color = FRAMES_ADJUSTED
    )
  ) +
  scale_color_viridis(
  ) +
  labs(
    x = "MyD88 Norm. Intensity",
    y = "mScarlet Norm. Intensity",
    color = "Time\n(Frames)"
  ) +
  facet_grid(
    QUERY_PROTEIN~FRAMES_SINCE_LANDING_CAT,
    scales = "free"
  ) +
  dark_theme_classic(
    base_size = 20
  ) #+
# coord_fixed()

ggsave(
  # Save vector image
  "TimePath.pdf",
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)


