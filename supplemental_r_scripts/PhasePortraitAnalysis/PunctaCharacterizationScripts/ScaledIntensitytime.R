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
    MAD_REFERENCE_TOTAL_INTENSITY = mad(REFERENCE_TOTAL_INTENSITY),
    MAD_QUERY_TOTAL_INTENSITY = mad(QUERY_TOTAL_INTENSITY),
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY)
  ) %>% 
  group_by(
    PLOT_FACETS,
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, FRAMES_SINCE_LANDING_CAT, IMAGENUMBER
  ) %>% 
  mutate(
    # REFERENCE_TOTAL_INTENSITY = roll_median(REFERENCE_TOTAL_INTENSITY, n = 5, fill = NA, align = "center"),
    # QUERY_TOTAL_INTENSITY = roll_median(QUERY_TOTAL_INTENSITY, n = 5, fill = NA, align = "center"),
    N_TEST = ifelse(N >= quantile(N, 0.50), T, F)
  ) %>% 
  filter(
    # N_TEST == T
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER
  ) %>% 
  mutate(
    SCALED_REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    SCALED_QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY),
    SCALED_MAD_REFERENCE_TOTAL_INTENSITY = MAD_REFERENCE_TOTAL_INTENSITY/max(MAD_REFERENCE_TOTAL_INTENSITY),
    SCALED_MAD_QUERY_TOTAL_INTENSITY = MAD_QUERY_TOTAL_INTENSITY/max(MAD_QUERY_TOTAL_INTENSITY)
  )


ggplot(
) +
  geom_ribbon(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      ymin = REFERENCE_TOTAL_INTENSITY - MAD_REFERENCE_TOTAL_INTENSITY,
      ymax =  REFERENCE_TOTAL_INTENSITY + MAD_REFERENCE_TOTAL_INTENSITY,
      fill = REFERENCE_PROTEIN,
      group = REFERENCE_PROTEIN
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_ribbon(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      ymin = QUERY_TOTAL_INTENSITY - MAD_QUERY_TOTAL_INTENSITY,
      ymax =  QUERY_TOTAL_INTENSITY + MAD_QUERY_TOTAL_INTENSITY,
      fill = QUERY_PROTEIN,
      group = QUERY_PROTEIN
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_path(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  scale_color_manual(
    values = c("magenta",  "magenta", "green", "magenta")
  ) +
  scale_fill_manual(
    values = c("magenta",  "magenta", "green", "magenta")
  ) +
  labs(
    x = "Time (Frames)",
    y = "Norm. Intensity",
    color = "Protein",
    fill = "Protein"
  ) +
  facet_grid(
    QUERY_PROTEIN~FRAMES_SINCE_LANDING_CAT,
    scales = "free"
  ) +
  dark_theme_classic(
    base_size = 20
  ) 


ggsave(
  # Save vector image
  "IntensityTimeLine.pdf",
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)



ggplot(
) +
  geom_ribbon(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      ymin = SCALED_REFERENCE_TOTAL_INTENSITY - SCALED_MAD_REFERENCE_TOTAL_INTENSITY,
      ymax =  SCALED_REFERENCE_TOTAL_INTENSITY + SCALED_MAD_REFERENCE_TOTAL_INTENSITY,
      fill = REFERENCE_PROTEIN,
      group = REFERENCE_PROTEIN
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_ribbon(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      ymin = SCALED_QUERY_TOTAL_INTENSITY - SCALED_MAD_QUERY_TOTAL_INTENSITY,
      ymax =  SCALED_QUERY_TOTAL_INTENSITY + SCALED_MAD_QUERY_TOTAL_INTENSITY,
      fill = QUERY_PROTEIN,
      group = QUERY_PROTEIN
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_path(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      y = SCALED_REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = SummaryFlow,
    aes(
      x = FRAMES_ADJUSTED,
      y = SCALED_QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  scale_color_manual(
    values = c("magenta",  "magenta", "green", "magenta")
  ) +
  scale_fill_manual(
    values = c("magenta",  "magenta", "green", "magenta")
  ) +
  labs(
    x = "Time (Frames)",
    y = "Norm. Intensity",
    color = "Protein",
    fill = "Protein"
  ) +
  facet_grid(
    QUERY_PROTEIN~FRAMES_SINCE_LANDING_CAT,
    scales = "free"
  ) +
  dark_theme_classic(
    base_size = 20
  ) 


ggsave(
  # Save vector image
  "ScaledIntensityTimeLine.pdf",
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)