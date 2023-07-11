library(gganimate)

StatTable <- fread("/Users/u_deliz/Desktop/PhasePortraitAnalysis/StatTable - LeadLag 5 - StepSize 2.csv.gz")
AnimationStatTable <- StatTable# %>% filter(LIGAND_DENSITY_CAT == 32, COHORT %in% c("MyD88 TRAF6", "MyD88 NEMO", "MyD88 HOIL1")) %>% as.data.table()

# ProcessTable
AnimationStatTable <-
  AnimationStatTable %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT, IMAGE
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, IMAGE
  ) %>%
  mutate(
    IMAGENUMBER = cur_group_id()
  ) %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT
  ) %>%
  mutate(
    IMAGENUMBER = IMAGENUMBER - min(IMAGENUMBER) + 1
  ) %>%
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    FPS = round(FPS, 2),
    FACET = paste(COHORT, LIGAND_DENSITY_CAT, FPS, FRAMES_SINCE_LANDING_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER),
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER),
    
    RATIO = REFERENCE_TOTAL_INTENSITY/QUERY_TOTAL_INTENSITY,
    RATIO = ifelse(RATIO<1, -1/RATIO, RATIO),
    
    THETA = atan2(DELTA_QUERY_TOTAL_INTENSITY, DELTA_REFERENCE_TOTAL_INTENSITY),
    THETA = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
    THETA = THETA*180/pi+180,
    THETA = THETA/360,
    
    ADJUSTED_THETA = atan2(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
    ADJUSTED_THETA = atan2(-ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, -ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
    ADJUSTED_THETA = ADJUSTED_THETA*180/pi+180,
    ADJUSTED_THETA = ADJUSTED_THETA/360
  ) %>% 
  as.data.table()

# Variables list
Facets <- unique(AnimationStatTable$PLOT_FACETS)
z_variables <-c(
  "TIME_ADJUSTED", "N",
  "THETA","ADJUSTED_THETA",
  "RATIO"
)

PlotList <- NULL
# z-axis
rep_z_variables = rep(z_variables, NROW(Facets))
PlotList$z_axis = rep_z_variables
# category
Facets <- sort(Facets)
PlotList$Facets = rep(Facets, each = NROW(z_variables))
# Combine
PlotList <- as_tibble(PlotList)

PlotFx <- function(FacetX){
  tryCatch({
    
    if(PlotList$z_axis[FacetX] == "N"){
      mult_format <- function(l) {
        l <- format(l, digits = 0 , scientific = F)
        l <- paste0("2^'", l, "'")
        parse(text=l)
      }
    } else{
      mult_format <- function(l) {l}
    }
    
    # Pull table
    TempTable <-
      AnimationStatTable %>%
      as_tibble() %>% 
      filter(
        PLOT_FACETS == PlotList$Facets[FacetX],
      ) %>% 
      as.data.table()
    
    TempTable <-
      TempTable %>% 
      group_by(
        ROUNDED_REFERENCE_TOTAL_INTENSITY,
        ROUNDED_QUERY_TOTAL_INTENSITY
      ) %>% 
      mutate(
        FRAMES_SINCE_LANDING_CAT = round(FRAMES_SINCE_LANDING/10)*10
      ) %>% 
      group_by(
        ROUNDED_REFERENCE_TOTAL_INTENSITY,
        ROUNDED_QUERY_TOTAL_INTENSITY,
        FRAMES_SINCE_LANDING_CAT
      ) %>% 
      mutate(
        N = n(),
        N = log(N+1, 2),
        z_axis = (!!as.name(PlotList$z_axis[FacetX]))
      ) %>% 
      as.data.table()

    # Axis names
    if(PlotList$z_axis[FacetX] == "TIME_ADJUSTED"){
      z_title = "Age (s)"
    }
    if(PlotList$z_axis[FacetX] == "N"){
      z_title = "N"
    }
    if(PlotList$z_axis[FacetX] == "RATIO"){
      z_title = "Ratio"
    }
    if(PlotList$z_axis[FacetX] == "ADJUSTED_RATIO"){
      z_title = "Scaled Ratio"
    }
    if(PlotList$z_axis[FacetX] == "THETA"){
      z_title = "Scaled Angle"
    }
    if(PlotList$z_axis[FacetX] == "ADJUSTED_THETA"){
      z_title = "Scaled Angle"
    }
    
    PlotTitle <-
      bquote(
        "Cell Line " *
          .(TempTable$COHORT[1]) *
          " ~ " *
          .( TempTable$LIGAND_DENSITY_CAT[1]) *
          " mol. " *
          µm^2 *
          " ~ " *
          .(TempTable$FPS[1]) *
          " ~ Rep. "*
          .(TempTable$IMAGENUMBER[1])
      )

    # Generate tables of frames
    ThresholdLoop <- function(FrameX){
      TempFrameTable <-
        TempTable %>%
        filter(
          FRAMES_SINCE_LANDING >= (FrameX-10),
          FRAMES_SINCE_LANDING <= FrameX
        ) %>%
        as_tibble() %>%
        group_by(
          REFERENCE_PROTEIN,
          QUERY_PROTEIN,
          ROUNDED_REFERENCE_TOTAL_INTENSITY,
          ROUNDED_QUERY_TOTAL_INTENSITY
        ) %>%
        summarize(
          z_axis = median(z_axis),
          N = n(),
          THRESHOLD = FrameX,
          
          REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
          QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY)
        ) %>%
        filter(
          N >= 5
        ) %>%
        ungroup() %>%
        filter(
          N >= quantile(N, .5, na.rm = T)
        ) %>%
        as.data.table()
      return(TempFrameTable)
    }
    FrameTable <- mclapply((1:11)*10, ThresholdLoop)
    FrameTable <- rbindlist(FrameTable)
    FrameTable <- FrameTable %>% arrange(THRESHOLD) %>% as.data.table()

    Plot <-
      ggplot(
        FrameTable,
        aes(
          x = ROUNDED_REFERENCE_TOTAL_INTENSITY,
          y = ROUNDED_QUERY_TOTAL_INTENSITY,
          fill = z_axis
        )
      ) +
      geom_tile() +
      geom_abline(
        intercept = 0,
        slope = 1
      ) +
      scale_fill_viridis(
        option = "turbo",
        labels = mult_format
      ) +
      transition_manual(
        frames = THRESHOLD
      ) +
      labs(
        x = paste(TempTable$REFERENCE_PROTEIN[1], "Norm. Int."),
        y = paste(TempTable$QUERY_PROTEIN[1], "Norm. Int."),
        fill = z_title,
        title = as.expression(PlotTitle),
        subtitle = 'Time: {frame*10-10} s'
      ) +
      theme_classic() +
      theme(
        rect = element_rect(fill = "transparent")
      ) +
      coord_fixed()

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

    SaveName <-
      paste0(
        z_title, " - ",
        "Replicate ", TempTable$IMAGENUMBER[1],
        # " - Cummulative.gif"
        ".gif"
      )

    SaveName <- gsub(" - \\.", "\\.", SaveName)
    SaveName <- file.path(COHORT_FOLDER, SaveName)
    # SaveName <- if(file.exists(SaveName)){file.remove(SaveName)}

    anim_save(
      SaveName,
      Plot,
      height = 4.76*300,
      width = 11.5*300,
      res = 300, fps = 10, end_pause = 20
    )

    print(SaveName)

  }, error = function(e) {print(paste("Error with PlotFx. FacetX =",FacetX))})
}
PlotPaths <- lapply(1:NROW(PlotList), PlotFx)
PlotPaths <- unlist(PlotPaths)
