
library(ggExtra)
library(gridExtra)

StatTableStoichiometry <-
  NormStatTable %>% 
  filter(
    COHORT %in% paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER),
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25
  ) %>%
  arrange(
    LIGAND_DENSITY_CAT, COHORT, IMAGE
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT, COHORT, IMAGE
  ) %>% 
  mutate(
    IMAGENUMBER = cur_group_id()
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT, COHORT, IMAGE
  ) %>% 
  mutate(
    IMAGENUMBER = IMAGENUMBER - min(IMAGENUMBER) + 1
  ) %>% 
  filter(
    IMAGENUMBER == 1
  ) %>% 
  arrange(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  # filter(
  #   FRAMES_ADJUSTED <= 97
  # ) %>% 
  mutate(
    N = n(),
    TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED),
    # REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/GFP_SummaryDF$P[FRAMES_ADJUSTED-1],
    # QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/mScarlet_SummaryDF$P[FRAMES_ADJUSTED-1],
    # REFERENCE_TOTAL_INTENSITY = exp(0.1357531)*exp(0.8871245*FRAMES_ADJUSTED)/REFERENCE_TOTAL_INTENSITY,
    # QUERY_TOTAL_INTENSITY = exp(0.02411673)*exp(0.9536553*FRAMES_ADJUSTED)/QUERY_TOTAL_INTENSITY,
  ) %>% 
  filter(
    N >= 75
  ) %>% 
  mutate(
    MAX_TIME_ADJUSTED = max(TIME_ADJUSTED, na.rm =T),
    # SMOOTH_REFERENCE = REFERENCE_TOTAL_INTENSITY,
    # SMOOTH_QUERY = QUERY_TOTAL_INTENSITY,
    SMOOTH_REFERENCE = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    SMOOTH_QUERY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
    COMPLEX_SIZE = REFERENCE_TOTAL_INTENSITY+QUERY_TOTAL_INTENSITY
  ) %>% 
  mutate(
    SMOOTH_REFERENCE = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    SMOOTH_QUERY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  mutate(
    MAX_REFERENCE = max(SMOOTH_REFERENCE, na.rm = T),
    MAX_QUERY = max(SMOOTH_QUERY, na.rm = T),
    STOICHIOMETRY = REFERENCE_TOTAL_INTENSITY/(COMPLEX_SIZE),
  ) %>% 
  mutate(
    MAX_COMPLEX_SIZE = MAX_REFERENCE + MAX_QUERY,
    MAX_REFERENCE_TIME = min(case_when(SMOOTH_REFERENCE == MAX_REFERENCE ~ TIME_ADJUSTED), na.rm = T),
    MAX_QUERY_TIME = min(case_when(SMOOTH_QUERY == MAX_QUERY ~ TIME_ADJUSTED), na.rm = T),
    STARTING_STOICHIOMETRY = median(STOICHIOMETRY[1:5], na.rm = T),
    ENDING_STOICHIOMETRY = median(STOICHIOMETRY[(n()-5):n()], na.rm = T)
  ) %>% 
  mutate(
    MAX_COMPLEX_TIME =  max(c(MAX_REFERENCE_TIME, MAX_QUERY_TIME))
  ) %>% 
  as_tibble()

Subset <-
  StatTableStoichiometry %>% 
  mutate(
    COHORT = factor(COHORT, levels = paste("MyD88", PROTEIN_ORDER))
  ) %>% 
  mutate(
    COHORT = forcats::as_factor(COHORT)
  ) %>% 
  arrange(
    COHORT
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  filter(
    FRAMES_ADJUSTED==min(FRAMES_ADJUSTED)
  ) %>%
  as_tibble() %>% 
  group_split(
    COHORT
  )

TrackVariables <-
  c("MAX_COMPLEX_SIZE", "MAX_TIME_ADJUSTED",
    "STARTING_STOICHIOMETRY", "ENDING_STOICHIOMETRY",
    "MAX_REFERENCE", "MAX_REFERENCE_TIME",
    "MAX_QUERY", "MAX_QUERY_TIME"
    )

TrackVariables <- combn(TrackVariables, 2)
TrackVariables <- t(TrackVariables)
TrackVariables <- as.data.table(TrackVariables)
names(TrackVariables) <- c("x", "y")

COMPARISON_DIRECTORY <- file.path(OUTPUT_DIRECTORY, "variable_comparison")
if(!file.exists(COMPARISON_DIRECTORY)){
  dir.create(COMPARISON_DIRECTORY)
}

TracksVariableFx <- function(CombinationX){
  
  x = TrackVariables$x[CombinationX]
  y = TrackVariables$y[CombinationX]
  
  R.devices::suppressGraphics({
    MakePanels <- function(df){
      
      df <-
        df %>% 
        mutate(
          x = !!as.name(x),
          y = !!as.name(y)
        ) %>% 
        as.data.table()
      
      x_min = ifelse(grepl("STOICHIOMETRY", x) | grepl("TIME", x), 0, as.numeric(quantile(df$x, 0.05)))
      y_min = ifelse(grepl("STOICHIOMETRY", y) | grepl("TIME", y), 0, as.numeric(quantile(df$y, 0.05)))
      
      x_max = ifelse(grepl("STOICHIOMETRY", x), 1, as.numeric(quantile(df$x, 0.95)))
      y_max = ifelse(grepl("STOICHIOMETRY", y), 1, as.numeric(quantile(df$y, 0.95)))
      
      p <-
        ggplot(
          df,
          aes(
            x, y
          )
        ) +
        geom_point(
          size = 0,
          color = "white"
        ) +
        stat_density2d(geom="tile", aes(fill = ..ndensity..), adjust = 2, contour = FALSE) +
        # stat_density2d(geom="tile", aes(fill = log(..ndensity..*1000+1, 10)), adjust = 2, contour = FALSE) +
        # geom_hline(yintercept = 0.5) +
        # geom_abline(slope = 1) +
        scale_x_continuous(limits = c(x_min, x_max))+
        scale_y_continuous(limits = c(y_min, y_max))+
        scale_fill_viridis() +
        labs(
          title = df$COHORT[1],
          x = x,
          y = y
          # x = "Maximum Query Size",
          # y = "Time to Max (s)"
        ) +
        theme_classic()  +
        theme(
          legend.position = "bottom"
        )
      
      p <- ggMarginal(p, type = "density", adjust = 2)
      
      return(p)
    }
    Panels <- lapply(Subset, MakePanels)
    n_row = sqrt(NROW(Panels))
    n_row = floor(n_row)
    
    # CombinedPanels <- do.call("grid.arrange", c(Panels, nrow=n_row))
    CombinedPanels <- do.call("grid.arrange", c(Panels, nrow=1))
    
    SaveName <-
      paste0(x, " ", y, ".png")
    
    ggsave(
      file = file.path(COMPARISON_DIRECTORY, SaveName),
      plot = CombinedPanels,
      height = 4,
      width = 16
    )
    
    SaveName <-
      paste0(x, " ", y, ".svg")
    
    ggsave(
      file = file.path(COMPARISON_DIRECTORY, SaveName),
      plot = CombinedPanels,
      height = 4,
      width = 16
    )
  })
  return(SaveName)
}
lapply(1:NROW(TrackVariables), TracksVariableFx)
# 
# PointVariables <-
#   c("STOICHIOMETRY", "TIME_ADJUSTED")
