RatioSize <-
  NormStatTable %>% 
  mutate(
    PLOT_FACET = paste(COHORT, LIGAND_DENSITY_CAT, FPS)
  ) %>% 
  filter(
    # COHORT == "MyD88 TRAF6",
    # ROUNDED_REFERENCE_TOTAL_INTENSITY > 1 & ROUNDED_QUERY_TOTAL_INTENSITY > 1
  ) %>% 
  mutate(
    RATIO = REFERENCE_TOTAL_INTENSITY/(QUERY_TOTAL_INTENSITY),
  ) %>%
  mutate(
    RATIO = ifelse(RATIO <1, -1/RATIO, RATIO)
  ) %>%
  mutate(
    RATIO = round(RATIO/2)*2
  ) %>%
# 
#   mutate(
#     RATIO = REFERENCE_TOTAL_INTENSITY/(QUERY_TOTAL_INTENSITY+REFERENCE_TOTAL_INTENSITY),
#   ) %>%
#   mutate(
#     RATIO = round(RATIO, 1),
#   ) %>%
  
  mutate(
    TOTAL_SIZE = REFERENCE_TOTAL_INTENSITY + QUERY_TOTAL_INTENSITY
  ) %>% 
  mutate(
    TOTAL_SIZE = round(TOTAL_SIZE/2)*2
  ) %>% 
  group_by(
    PLOT_FACET,
    FPS,
    LIGAND_DENSITY_CAT,
    
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_SINCE_LANDING_CAT,
    RATIO,
    TOTAL_SIZE
  ) %>% 
  summarize(
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY),
    N = n()
  ) %>% 
  filter(
    # TOTAL_SIZE > 0,
    N >= 100,
    # TOTAL_SIZE <= (50+30)
  ) %>% 
  group_by(
    PLOT_FACET,
    RATIO
  ) %>% 
  mutate(
    N = n()
  ) %>% 
  filter(
    N >= 3
  ) %>% 
  as.data.table()

if(min(RatioSize$RATIO) < 0 &
   max(RatioSize$RATIO) > 0){
  
  plot_values <-
    c(
      min(RatioSize$RATIO)/200*(200:1),
      0,
      max(RatioSize$RATIO)/200*(1:200)
    )
  plot_values <- plot_values - min(plot_values)
  plot_values <- plot_values/max(plot_values)
  
  plot_colors <-
    c(
      colorRampPalette(c("navy","royalblue","lightskyblue"))(200),
      # "#FFFFFF",
      colorRampPalette(c("mistyrose", "red2","darkred"))(200)
    )
  
} else{
  plot_values <- 0:100/100
  plot_colors <- viridis(101, option = "turbo")
}

RatioSize1 <- RatioSize %>% mutate(PROTEIN = REFERENCE_PROTEIN) %>% as.data.table()
RatioSize2 <- RatioSize %>% mutate(
  PROTEIN = QUERY_PROTEIN,
  DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY
  ) %>% as.data.table()

RatioSize <- rbindlist(list(RatioSize1, RatioSize2))

Plots <- unique(RatioSize$PLOT_FACET)

PlotFx <- function(PlotX){
  
  
  TempRatioSize <-
    RatioSize %>% 
    filter(
      PLOT_FACET == Plots[PlotX]
    ) %>% 
    as.data.table()
  
  
  COHORT_FOLDER = file.path(OUTPUT_DIRECTORY,
                            paste0(
                              "Cell Line ", TempRatioSize$COHORT[[1]], " - ",
                              TempRatioSize$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
                              TempRatioSize$FPS[[1]], " Hz"
                            ))
  
  
  ggplot(
    TempRatioSize,
    aes(
      TOTAL_SIZE,
      DELTA_REFERENCE_TOTAL_INTENSITY*60,
      color = RATIO,
      group = RATIO
    )
  ) +
    geom_hline(
      yintercept = 0,
      size = 2
    ) +
    geom_point() +
    geom_path() +
    scale_color_gradientn(
      colors = plot_colors, na.value = "black",
      values = plot_values
    ) +
    labs(
      x = "Cluster Size",
      y = "Rate (mol./min)",
      color = "Ratio"
    ) +
    facet_grid(
      PROTEIN~FRAMES_SINCE_LANDING_CAT,
      scales = "free_y"
    ) +
    theme_classic()
  
  SaveName <- "1D_Stoichiometry_Rates.pdf"
  
  ggsave(
    file.path(COHORT_FOLDER, SaveName),
    height = 9,
    width = 16
  )
  
}
PlotList <- mclapply(1:NROW(Plots), PlotFx)
PlotList <- unlist(PlotList)
