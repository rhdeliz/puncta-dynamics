# Filter values
# FRAMES_SINCE_LANDING_CAT_BIN <- 150
FRAMES_SINCE_LANDING_CAT_BIN <- 300
N_THRESHOLD = 100
MAX_X = 15

# OUTPUT_DIRECTORY = "/Users/u_deliz/Desktop/PhasePortraitAnalysis"

setwd(OUTPUT_DIRECTORY)

# Table Name
SaveName <-
  paste0(
    "NormStatTable - ",
    "LeadLag ", LEAD_LAG,
    ".gz.parquet"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)
NormStatTable <- read_parquet(file.path(OUTPUT_DIRECTORY, SaveName))

SplitStatTable <-
  NormStatTable %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT == FRAMES_SINCE_LANDING_CAT_BIN
  ) %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    FPS = round(FPS, 2)
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    FACET = paste(COHORT, LIGAND_DENSITY_CAT, FPS, FRAMES_SINCE_LANDING_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN),
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN),
    # ROUNDED_REFERENCE_TOTAL_INTENSITY = ifelse(ROUNDED_REFERENCE_TOTAL_INTENSITY <= STEP_SIZE*2, paste("<", STEP_SIZE*2), paste("≥", STEP_SIZE*2))
  ) %>% 
  arrange(
    PLOT_FACETS
  ) %>% 
  as_tibble()

# Variables list
Facets <- unique(SplitStatTable$PLOT_FACETS)
x_variables <- c("ROUNDED_REFERENCE_TOTAL_INTENSITY", "ROUNDED_QUERY_TOTAL_INTENSITY")
y_variables <-c(
  "DELTA_REFERENCE_TOTAL_INTENSITY", "DELTA_QUERY_TOTAL_INTENSITY"
)
z_variables <- rev(x_variables)

PlotList <- NULL
# x-axis
rep_x_variables = rep(x_variables, each = NROW(y_variables))
rep_x_variables = rep(rep_x_variables, NROW(Facets))
PlotList$x_axis = rep_x_variables
# y-axis
rep_y_variables = rep(y_variables, NROW(x_variables))
PlotList$y_axis = rep(rep_y_variables, NROW(Facets))
# z-axis
rep_z_variables = rep(z_variables, each = NROW(y_variables))
rep_z_variables = rep(rep_z_variables, NROW(Facets))
PlotList$z_axis = rep_z_variables
# category
Facets <- sort(Facets)
PlotList$Facets = rep(Facets, each = NROW(x_variables)*NROW(y_variables))
# Combine
PlotList <- as_tibble(PlotList)

SplitStatTable <-
  SplitStatTable %>%
  arrange(
    FACET
  )# %>% 
  # group_split(
  #   FACET
  # )

PlotTable <- function(FacetX){

  TempTable <-
    SplitStatTable %>%
    as_tibble() %>% 
    mutate(
      # Define plot axis replacing variables with generic names
      x_axis = (!!as.name(PlotList$x_axis[FacetX])),
      y_axis = (!!as.name(PlotList$y_axis[FacetX]))*60,
      z_axis = (!!as.name(PlotList$z_axis[FacetX])),
      y_name = PlotList$y_axis[FacetX]
    ) %>%
    filter(
      x_axis <= MAX_X
    ) %>%
    mutate(
      z_axis = ifelse(z_axis <= STEP_SIZE, paste0("<", STEP_SIZE, "x"), paste0("≥", STEP_SIZE, "x"))
    ) %>% 
    group_by(
      PLOT_FACETS,
      COHORT, LIGAND_DENSITY_CAT, FPS,
      FRAMES_SINCE_LANDING_CAT,
      REFERENCE_PROTEIN, QUERY_PROTEIN,
      x_axis,
      z_axis,
      y_name
    ) %>% 
    summarize(
      N = n(),
      FRAMES_ADJUSTED = median(FRAMES_ADJUSTED),
      TIME_ADJUSTED = median(TIME_ADJUSTED),
      y_axis_sem = mad(y_axis),
      y_axis = median(y_axis)
    ) %>%
    mutate(
      y_axis_sem = y_axis_sem/sqrt(N),
      N_TEST = ifelse(N >= N_THRESHOLD, T, F)
    ) %>%
    filter(
      N_TEST == T
    ) %>%
    ungroup() %>% 
    mutate(
      MIN_DELTA = min(y_axis - y_axis_sem),
      MAX_DELTA = max(y_axis + y_axis_sem)
    ) %>% 
    filter(
      PLOT_FACETS == PlotList$Facets[FacetX]
    ) %>%
    as.data.table()
  
  return(TempTable)
  
}
PhasePortraitAll <- mclapply(1:NROW(PlotList), PlotTable)
# 
# LimitsTable <- rbindlist(PhasePortraitAll)
# 
# LimitsTable <- 
#   LimitsTable %>% 
#   group_by(
#     y_name
#   ) %>% 
#   summarize(
#     MIN_DELTA = min(y_axis - y_axis_sem)*1.2,
#     MAX_DELTA = max(y_axis + y_axis_sem)*1.2
#   ) %>% 
#   as.data.table()

PlotFx <- function(FacetX){
  tryCatch({
    
    # PlotLims <-
    #   LimitsTable %>% 
    #   filter(
    #     y_name == PlotList$y_axis[FacetX]
    #   ) %>% 
    #   as.data.table()
    
    PlotLims <-
      PhasePortraitAll[[FacetX]]
    
    TempTable <- PhasePortraitAll[[FacetX]]
    
    # Labels
    x_title <- PlotList$x_axis[FacetX]
    x_title <- gsub("_TOTAL_INTENSITY", "", x_title)
    x_title <- gsub("ROUNDED_", "", x_title)
    x_title <- gsub("REFERENCE", paste(TempTable$REFERENCE_PROTEIN[1], "(mol.)"), x_title)
    x_title <- gsub("QUERY", paste(TempTable$QUERY_PROTEIN[1], "(mol.)"), x_title)
    
    y_title <- PlotList$y_axis[FacetX]
    y_title <- gsub("_TOTAL_INTENSITY", " Rate (mol./min)", y_title)
    y_title <- gsub("ROUNDED_", "", y_title)
    y_title <- gsub("REFERENCE", TempTable$REFERENCE_PROTEIN[1], y_title)
    y_title <- gsub("QUERY", TempTable$QUERY_PROTEIN[1], y_title)
    y_title <- gsub("DELTA_", "", y_title)
    y_title <- gsub("ADJUSTED_", "Scaled ", y_title)
    y_title <- gsub("FRAMES_ADJUSTED", "Age", y_title)
    y_title <- gsub("RATIO", "Stoichiometry", y_title)
    
    z_title <- PlotList$z_axis[FacetX]
    z_title <- gsub("_TOTAL_INTENSITY", "", z_title)
    z_title <- gsub("ROUNDED_", "", z_title)
    z_title <- gsub("REFERENCE", paste0("", TempTable$REFERENCE_PROTEIN[1]), z_title)
    z_title <- gsub("QUERY", paste0("", TempTable$QUERY_PROTEIN[1]), z_title)
    z_title <- gsub("DELTA_", "Change ", z_title)
    z_title <- gsub("ADJUSTED_", "Scaled ", z_title)
    z_title <- gsub("FRAMES_ADJUSTED", "Age ", z_title)
    
    
    InterpolationTable <-
      TempTable %>% 
      select(
        x_axis,
        y_axis,
        y_axis_sem,
        z_axis,
        FRAMES_SINCE_LANDING_CAT
      ) %>% 
      as_tibble() %>% 
      group_by(
        z_axis,
        FRAMES_SINCE_LANDING_CAT
      ) %>% 
      complete(
        x_axis = full_seq(c(min(x_axis), max(x_axis)+STEP_SIZE), period = STEP_SIZE/21),
      ) %>% 
      fill(
        one_of(
          "y_axis",
          "y_axis_sem"
        ),
        .direction = "downup"
      ) %>% 
      group_by(
        z_axis,
        FRAMES_SINCE_LANDING_CAT
      ) %>% 
      mutate(
        y_axis_old = y_axis,
        y_axis_sem_old = y_axis_sem,
        y_axis = roll_mean(y_axis, 11, fill = NA, align = "left"),
        y_axis_sem = roll_mean(y_axis_sem, 11, fill = NA, align = "left"),
        y_axis = roll_mean(y_axis, 11, fill = NA, align = "left"),
        y_axis_sem = roll_mean(y_axis_sem, 11, fill = NA, align = "left")
      ) %>% 
      tidyr::drop_na()
    
    # Save folder
    COHORT_FOLDER = file.path(OUTPUT_DIRECTORY,
                              paste0(
                                "Cell Line ", TempTable$COHORT[1], " - ",
                                TempTable$LIGAND_DENSITY_CAT[1],  " mol. µm^2 - ",
                                TempTable$FPS[1], " Hz"
                              ))
    
    if(!file.exists(COHORT_FOLDER)){
      dir.create(COHORT_FOLDER)
    }
    
    # File name
    SaveName <-
      paste0(
        gsub(" \\(mol.\\)|Change |\\(x\\) ", "", x_title), " - ",
        gsub(" Norm. Int.|Change | Rate \\(mol.\\/min\\)|\\(x\\) ", "", y_title), " - ",
        gsub(" Norm. Int.|Change |\\(x\\) ", "", z_title), " - ",
        "StepSize ", round(STEP_SIZE, 2), " - ",
        "paper",
        ".pdf"
      )
    
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    
    ggplot() +
      geom_hline(
        yintercept = 0,
        size = 2,
        color = "black"
      ) +
      geom_ribbon(
        data = InterpolationTable,
        aes(
          x = x_axis,
          y = y_axis,
          color = z_axis,
          fill = z_axis,
          group = z_axis,
          ymin = y_axis - y_axis_sem,
          ymax =  y_axis + y_axis_sem
        ),
        alpha = .25
      ) +
      geom_point(
        data = TempTable,
        aes(
          x = x_axis,
          y = y_axis,
          color = z_axis,
          fill = z_axis,
          group = z_axis,
          size = N
        )
      ) +
      scale_color_brewer(
        palette = "Set1"
      ) +
      scale_size_continuous(
        guide = "none"
      ) +
      scale_x_continuous(
        limits = c(0, MAX_X)
      ) +
      scale_y_continuous(
        limits = c(PlotLims$MIN_DELTA[1], PlotLims$MAX_DELTA[1])
      ) +
      labs(
        x = x_title,
        # y = y_title,
        y = gsub(" \\(mol", "\n\\(mol", y_title),
        color = gsub(" Norm. Int.", "\nNorm. Int.", z_title),
        fill = gsub(" Norm. Int.", "\nNorm. Int.", z_title)
      ) +
      theme_classic(
        base_size = 20
      ) +
      theme(
        # legend.key.height= unit(1.25, 'cm'),
        legend.position = "right",
        rect = element_rect(fill = "transparent"),
        strip.background = element_blank()
      )
    
    ggsave(
      # Save vector image
      file.path(COHORT_FOLDER, SaveName),
      height = 3,
      width = 3*sqrt(2)*1.25,
      device=cairo_pdf
    )
    
    return(file.path(COHORT_FOLDER, SaveName))
    
  }, error = function(e) {print(paste("Error with PlotFx. FacetX =",FacetX))})
  
}
PlotPaths <- mclapply(1:NROW(PlotList), PlotFx, mc.cores = detectCores(logical = F))
PlotPaths <- unlist(PlotPaths)
