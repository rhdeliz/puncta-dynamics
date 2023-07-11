setwd(OUTPUT_DIRECTORY)

Facets <- unique(LinearPhasePortrait$LIGAND_DENSITY_CAT)

Variables <- 
  c(
    "MAD_ADJUSTED_",
    "SD_ADJUSTED_",
    "ADJUSTED_",
    ""
  )

PlotList <- NULL
PlotList$Facets = rep(Facets, each = NROW(Variables))
PlotList$Variables = rep(Variables, NROW(Facets))
PlotList <- as.data.table(PlotList)
remove(Facets, Variables)

PlotFx <- function(FacetX){
  tryCatch({
    
    TempVariables <-
      paste0(
        PlotList$Variables[FacetX],
        c(
          "DELTA_REFERENCE_TOTAL_INTENSITY",
          "DELTA_QUERY_TOTAL_INTENSITY",
          "MAGNITUDE",
          "RAD_ANGLE",
          "DEG_ANGLE"
        )
      )
    
    TempTable <-
      LinearPhasePortrait %>% 
      group_by(
        QUERY_PROTEIN, COHORT
      ) %>% 
      mutate(
        TEST_COHORT = grepl(QUERY_PROTEIN, COHORT)
      ) %>% 
      filter(
        TEST_COHORT == TRUE
      ) %>% 
      select(-c(
        TEST_COHORT
      )) %>% 
      filter(
        # Plot filters
        ROUNDED_REFERENCE_TOTAL_INTENSITY >= 1,
        # ROUNDED_REFERENCE_TOTAL_INTENSITY <= 12,
        ROUNDED_QUERY_TOTAL_INTENSITY <= 8,
        LIGAND_DENSITY_CAT == PlotList$Facets[FacetX],
        N_TEST == TRUE,
        
        # Category filters
        MAX_REFERENCE_TOTAL_INTENSITY_CAT == "All",
        # DELTA_CAT == "All",
        SPOTS_WITHIN_RADIUS_CAT == "All"
      ) %>%
      group_by(
        COHORT, IMAGENUMBER
      ) %>% 
      mutate(
        # To select image with most spots
        IMAGE_N = n()
      ) %>% 
      group_by(
        COHORT
      ) %>% 
      filter(
        # Select image with most bins
        IMAGE_N == max(IMAGE_N)
      ) %>% 
      filter(
        IMAGENUMBER == min(IMAGENUMBER)
      ) %>% 
      group_by(
        COHORT, LIGAND_DENSITY_CAT, FPS, MAX_REFERENCE_TOTAL_INTENSITY_CAT, SPOTS_WITHIN_RADIUS_CAT, DELTA_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER, IMAGE,
        FACET, PLOT_FACETS
      ) %>% 
      mutate(
        FPS = paste(FPS, "Hz"),
        DELTA_CAT = gsub("ease", ".", DELTA_CAT),
        DELTA_CAT = factor(DELTA_CAT, levels = c("All", "Decr.", "Incr.")),
        QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = PROTEIN_ORDER),
        SD_ADJUSTED_MAGNITUDE = SD_ADJUSTED_MAGNITUDE/max(SD_ADJUSTED_MAGNITUDE, na.rm = T),
        MAD_MAGNITUDE = MAD_MAGNITUDE/max(MAD_MAGNITUDE, na.rm = T),
        MAD_ADJUSTED_MAGNITUDE = MAD_ADJUSTED_MAGNITUDE/max(MAD_ADJUSTED_MAGNITUDE, na.rm = T),
        MAGNITUDE = MAGNITUDE/max(MAGNITUDE, na.rm = T),
        ADJUSTED_MAGNITUDE = ADJUSTED_MAGNITUDE/max(ADJUSTED_MAGNITUDE, na.rm = T)
      ) %>% 
      group_by(
        PLOT_FACETS, IMAGENUMBER
      ) %>% 
      mutate(
        # Define plot axis replacing variables with generic names
        ref_int = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        qry_int = ROUNDED_QUERY_TOTAL_INTENSITY,
        
        d_ref_int = (!!as.name(TempVariables[1])),
        d_qry_int = (!!as.name(TempVariables[2])),
        magnitude = (!!as.name(TempVariables[3])),
        rad_angle = (!!as.name(TempVariables[4])),
        deg_angle = (!!as.name(TempVariables[5]))
      )
    
    PlotTitle <-
      paste0(
        TempTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2"
      )
    
    COHORT_FOLDER = PlotTitle
    COHORT_FOLDER = file.path(OUTPUT_DIRECTORY, COHORT_FOLDER)
    
    if(!file.exists(COHORT_FOLDER)){
      dir.create(COHORT_FOLDER)
    }
    setwd(COHORT_FOLDER)
    
    ggplot(
      TempTable# %>% filter(deg_angle <= 0.25)
    ) +
      # geom_raster(
      #   aes(
      #     x = ref_int,
      #     y = qry_int,
      #     fill = hsv(deg_angle*4)
      #   )
      # ) +
      # scale_fill_identity(guide = "none") +
      geom_vline(
        xintercept = c(6)
      ) +
      geom_arrow(
        aes(
          x = ref_int+1,
          y = qry_int+1,
          mag = 0.1,
          angle = atan2(d_qry_int, d_ref_int)*180/pi,
          color = hsv(deg_angle)
          # color = N
        ),
        # color = "white",
        size = .5,
        arrow.length = .33
      ) +
      scale_mag(
        max = 3,
        guide = 'none'
      ) +
      # Color of magnitude. Goes from gray to white
      # scale_color_gradient(
      #   low = "black",
      #   high = "white",
      #   trans = "log2",
      #   breaks = trans_breaks("log2", function(x) 2^x),
      #   labels = trans_format("log2", math_format(2^.x))
      # ) +
      # scale_color_distiller(
      #   palette = "RdPu",
      #   trans = "log2",
    #   breaks = trans_breaks("log2", function(x) 2^x),
    #   labels = trans_format("log2", math_format(2^.x))
    # ) +
    # Use angle of arrows to color the triangle. Full spectrum
    scale_color_identity(guide = "none") +
      # scale_x_continuous(
      #   trans = "log2",
      #   breaks = trans_breaks("log2", function(x) round(2^x, 1)),
      #   sec.axis = sec_axis(
      #     trans = ~.,
      #     breaks = trans_breaks("log2", function(x) 2^x),
      #     labels = trans_format("log2", math_format(2^.x))
      #   )
      # ) + 
      # scale_y_continuous(
      #   trans = "log2",
    #   breaks = trans_breaks("log2", function(x) round(2^x, 1)),
    #   sec.axis = sec_axis(
    #     trans = ~.,
    #     breaks = trans_breaks("log2", function(x) 2^x),
    #     labels = trans_format("log2", math_format(2^.x))
    #   )
    # ) +
    labs(
      title = PlotTitle,
      x = "MyD88 Size",
      y = "Query Size"
    ) +
      facet_grid(
        DELTA_CAT~QUERY_PROTEIN+FPS
        # nrow = 2
      ) +
      # Make dark to make the colors pop in a monitor
      # Use theme classic for printing
      dark_theme_classic(
        base_size = 20
      ) +
      theme(
        # legend.key.width= unit(1.5, 'cm'),
        legend.position = "none"
      ) +
      coord_fixed()
    
    # File name
    SaveName <-
      paste0(
        PlotList$Variables[FacetX],
        "PhasePortrait - ",
        "LigandCat ", PlotList$Facets[FacetX], " - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(1/STEP_SIZE, 2), " - ",
        ifelse(MOVING_AVERAGE_WINDOW > 0, paste("MvAvg", MOVING_AVERAGE_WINDOW, "- "), ""),
        ifelse(LOG_SCALE == T, "Log - ", ""),
        "Arrow Direction.pdf"
      )
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    
    ggsave(
      # Save vector image
      SaveName,
      # height = 4.76*1.25,
      # width = 11.5*1.25
      height = 4.76,
      width = 11.5
    )
    
    ggplot(
      TempTable
    ) +
      geom_vline(
        xintercept = c(6)
      ) +
      geom_arrow(
        aes(
          x = ref_int+1,
          y = qry_int+1,
          mag = 0.1,
          angle = atan2(d_qry_int, d_ref_int)*180/pi,
          color = magnitude
        ),
        size = .5,
        arrow.length = .33
      ) +
      scale_mag(
        max = 3,
        guide = 'none'
      ) +
      scale_color_viridis(
        option = "inferno"
      ) +
      labs(
        title = PlotTitle,
        x = "MyD88 Size",
        y = "Query Size",
        color = "Scaled\nMagnitude"
      ) +
      facet_grid(
        DELTA_CAT~QUERY_PROTEIN+FPS
        # nrow = 2
      ) +
      # Make dark to make the colors pop in a monitor
      # Use theme classic for printing
      dark_theme_classic(
        base_size = 20
      ) +
      theme(
        legend.key.height= unit(1.25, 'cm'),
        legend.position = "right"
      ) +
      coord_fixed()
    
    # File name
    SaveName <-
      paste0(
        PlotList$Variables[FacetX],
        "PhasePortrait - ",
        "LigandCat ", PlotList$Facets[FacetX], " - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(1/STEP_SIZE, 2), " - ",
        ifelse(MOVING_AVERAGE_WINDOW > 0, paste("MvAvg", MOVING_AVERAGE_WINDOW, "- "), ""),
        ifelse(LOG_SCALE == T, "Log - ", ""),
        "Arrow Magnitude.pdf"
      )
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    
    ggsave(
      # Save vector image
      SaveName,
      # height = 4.76*1.25,
      # width = 11.5*1.25
      height = 4.76,
      width = 11.5
    )
    
  }, error = function(e) {print(paste("Error with PlotFx", PlotList$Variables[FacetX], PlotList$Facets[FacetX]))})}
mclapply(1:NROW(PlotList), PlotFx)
