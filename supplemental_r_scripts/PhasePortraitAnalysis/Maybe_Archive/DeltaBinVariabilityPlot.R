DeltaTable <-
  AnalysisStatTable %>% 
  filter(
    COHORT == "MyD88 IRAK4",
    FPS == 0.5,
    LIGAND_DENSITY_CAT == 100,
    MAX_REFERENCE_TOTAL_INTENSITY_CAT == "All",
    DELTA_CAT == "Increase",
    SPOTS_WITHIN_RADIUS_CAT == "All"
  ) %>% 
  group_by(
    IMAGE
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  ungroup() %>% 
  filter(
    N == max(N)
  ) %>%
  group_by(
    ROUNDED_REFERENCE_TOTAL_INTENSITY,
    ROUNDED_QUERY_TOTAL_INTENSITY
  ) %>%
  mutate(
    N = n()
  ) %>%
  filter(
    N >= 100
  ) %>% 
  arrange(
    -N
  )

PlotTitle <-
  paste0(
    "Cell Line ", DeltaTable$COHORT[[1]], "\n",
    DeltaTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2\n",
    DeltaTable$FPS[[1]], " FPS"
  )

ggplot(
  DeltaTable,
  aes(
    # x = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY,
    x = DELTA_REFERENCE_TOTAL_INTENSITY,
    y = ..ndensity..,
    color = ROUNDED_REFERENCE_TOTAL_INTENSITY,
    group = ROUNDED_REFERENCE_TOTAL_INTENSITY
  )
) +
  geom_density() +
  scale_color_distiller(
    palette = "Greys",
    direction = 1
  ) +
  labs(
    title = PlotTitle,
    x = "Scaled MyD88 Size Change",
    y = "Scaled Denstity",
    color = "MyD88 Size"
  ) +
  # facet_wrap(
    # ~ROUNDED_QUERY_TOTAL_INTENSITY
  # )+
  dark_theme_classic() +
  theme(
    legend.position = "bottom"
  )

setwd(OUTPUT_DIRECTORY)

ggsave(
  # Save vector image
  "ScaledDeltaVariance.pdf",
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)