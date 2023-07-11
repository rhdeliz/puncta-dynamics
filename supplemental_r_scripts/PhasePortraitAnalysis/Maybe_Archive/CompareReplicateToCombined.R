# Import tables
SaveName <-
  paste0(
    "Normalized PhasePortrait - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", STEP_SIZE,
    ".csv.gz"
  )

PhasePortaitCombined <-
  fread(file.path(OUTPUT_DIRECTORY, SaveName))

SaveName <-
  paste0(
    "Normalized All PhasePortrait - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", STEP_SIZE,
    ".csv.gz"
  )

PhasePortaitSeparate <-
  fread(file.path(OUTPUT_DIRECTORY, SaveName))

# Label combined table
PhasePortaitCombined$IMAGENUMBER = "All"

LinearPhasePortrait <-
  rbindlist(list(PhasePortaitSeparate, PhasePortaitCombined), fill = T)

remove(PhasePortaitSeparate, PhasePortaitCombined)

LinearPhasePortrait <-
  LinearPhasePortrait %>% 
  filter(
    N_TEST == T
  ) %>% 
  mutate(
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, FRAMES_SINCE_LANDING_CAT)
  ) %>% 
  group_by(
    PLOT_FACETS,
    IMAGENUMBER,
    ROUNDED_REFERENCE_TOTAL_INTENSITY,
    ROUNDED_QUERY_TOTAL_INTENSITY
  ) %>% 
  mutate(
    SEM_DELTA_REFERENCE = MAD_DELTA_REFERENCE_TOTAL_INTENSITY/sqrt(N),
    SEM_DELTA_QUERY = MAD_DELTA_QUERY_TOTAL_INTENSITY/sqrt(N),
  ) %>% 
  as.data.table()

PlotFx <- function(PlotX){
  
  TempLinearPhasePortrait <-
    LinearPhasePortrait %>% 
    filter(
      PLOT_FACETS == Plots[PlotX]
    ) %>% 
    as.data.table()
  
  COHORT_FOLDER = file.path(OUTPUT_DIRECTORY,
                            paste0(
                              "Cell Line ", TempLinearPhasePortrait$COHORT[[1]], " - ",
                              TempLinearPhasePortrait$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
                              TempLinearPhasePortrait$FPS[[1]], " Hz"
                            ))
  
  if(!file.exists(COHORT_FOLDER)){
    dir.create(COHORT_FOLDER)
  }
  
  setwd(COHORT_FOLDER)
  
  ggplot() +
    geom_path(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER == "All") %>% as_tibble(),
      aes(
        x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        y = DELTA_REFERENCE_TOTAL_INTENSITY,
        group = IMAGENUMBER
      ),
      color = "black"
    ) +
    geom_ribbon(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER == "All") %>% as_tibble(),
      aes(
        x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        y = DELTA_REFERENCE_TOTAL_INTENSITY,
        ymin = DELTA_REFERENCE_TOTAL_INTENSITY - SEM_DELTA_REFERENCE,
        ymax = DELTA_REFERENCE_TOTAL_INTENSITY + SEM_DELTA_REFERENCE
      ),
      fill = "black",
      alpha = .5
    ) +
    geom_ribbon(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER != "All") %>% as_tibble(),
      aes(
        x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        y = DELTA_REFERENCE_TOTAL_INTENSITY,
        ymin = DELTA_REFERENCE_TOTAL_INTENSITY - SEM_DELTA_REFERENCE,
        ymax = DELTA_REFERENCE_TOTAL_INTENSITY + SEM_DELTA_REFERENCE,
        color = IMAGENUMBER,
        fill = IMAGENUMBER
      ),
      alpha = .25
    ) +
    geom_path(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER != "All") %>% as_tibble(),
      aes(
        x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        y = DELTA_REFERENCE_TOTAL_INTENSITY,
        color = IMAGENUMBER,
        group = IMAGENUMBER
      )
    ) +
    scale_fill_brewer(
      palette = "Set1"
    ) +
    scale_color_brewer(
      palette = "Set1"
    ) +
    labs(
      x = "MyD88 Scaled Int.",
      y = "dMyD88",
      color = "Replicate",
      fill = "Replicate"
    ) +
    facet_wrap(
      ~ROUNDED_QUERY_TOTAL_INTENSITY,
      scales = "free"
    ) +
    theme_classic()
  
  ggsave(
    paste0("SEM_Reference FramesSinceLanding_",TempLinearPhasePortrait$FRAMES_SINCE_LANDING_CAT[1],".pdf"),
    height = 9/1.5,
    width = 16/1.5
  )
  
  
  ggplot() +
    geom_path(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER == "All") %>% as_tibble(),
      aes(
        x = ROUNDED_QUERY_TOTAL_INTENSITY,
        y = DELTA_QUERY_TOTAL_INTENSITY,
        group = IMAGENUMBER
      ),
      color = "black"
    ) +
    geom_ribbon(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER == "All") %>% as_tibble(),
      aes(
        x = ROUNDED_QUERY_TOTAL_INTENSITY,
        y = DELTA_QUERY_TOTAL_INTENSITY,
        ymin = DELTA_QUERY_TOTAL_INTENSITY - SEM_DELTA_QUERY,
        ymax = DELTA_QUERY_TOTAL_INTENSITY + SEM_DELTA_QUERY
      ),
      fill = "black",
      alpha = .5
    ) +
    geom_ribbon(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER != "All") %>% as_tibble(),
      aes(
        x = ROUNDED_QUERY_TOTAL_INTENSITY,
        y = DELTA_QUERY_TOTAL_INTENSITY,
        ymin = DELTA_QUERY_TOTAL_INTENSITY - SEM_DELTA_QUERY,
        ymax = DELTA_QUERY_TOTAL_INTENSITY + SEM_DELTA_QUERY,
        color = IMAGENUMBER,
        fill = IMAGENUMBER
      ),
      alpha = .25
    ) +
    geom_path(
      data = TempLinearPhasePortrait %>% filter(IMAGENUMBER != "All") %>% as_tibble(),
      aes(
        x = ROUNDED_QUERY_TOTAL_INTENSITY,
        y = DELTA_QUERY_TOTAL_INTENSITY,
        color = IMAGENUMBER,
        group = IMAGENUMBER
      )
    ) +
    scale_fill_brewer(
      palette = "Set1"
    ) +
    scale_color_brewer(
      palette = "Set1"
    ) +
    labs(
      x = "Query Scaled Int.",
      y = "d Query",
      color = "Replicate",
      fill = "Replicate"
    ) +
    facet_wrap(
      ~ROUNDED_REFERENCE_TOTAL_INTENSITY,
      scales = "free"
    ) +
    theme_classic()
  
  ggsave(
    paste0("SEM_Query FramesSinceLanding_",TempLinearPhasePortrait$FRAMES_SINCE_LANDING_CAT[1],".pdf"),
    height = 9/1.5,
    width = 16/1.5
  )
  
}
Plots <- unique(LinearPhasePortrait$PLOT_FACETS)
mclapply(1:NROW(Plots), PlotFx)
