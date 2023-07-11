setwd(OUTPUT_DIRECTORY)




Facets <- unique(LinearPhasePortrait$PLOT_FACETS)

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
      filter(
        # Plot filters
        # ROUNDED_REFERENCE_TOTAL_INTENSITY <= 8,
        PLOT_FACETS == PlotList$Facets[FacetX],
        N_TEST == TRUE,
        
        # Category filters
        MAX_REFERENCE_TOTAL_INTENSITY_CAT == "All",
        # DELTA_CAT == "All",
        SPOTS_WITHIN_RADIUS_CAT == "All"
      ) %>%
      group_by(
        COHORT, LIGAND_DENSITY_CAT, FPS, MAX_REFERENCE_TOTAL_INTENSITY_CAT, SPOTS_WITHIN_RADIUS_CAT, DELTA_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER, IMAGE,
        FACET, PLOT_FACETS
      ) %>% 
      mutate(
        # 95th percentile to make sure large triangles don't alter results
        # SD_MAGNITUDE = SD_MAGNITUDE/quantile(SD_MAGNITUDE, .95, na.rm = T),
        SD_ADJUSTED_MAGNITUDE = SD_ADJUSTED_MAGNITUDE/quantile(SD_ADJUSTED_MAGNITUDE, .95, na.rm = T),
        MAD_MAGNITUDE = MAD_MAGNITUDE/quantile(MAD_MAGNITUDE, .95, na.rm = T),
        MAD_ADJUSTED_MAGNITUDE = MAD_ADJUSTED_MAGNITUDE/quantile(MAD_ADJUSTED_MAGNITUDE, .95, na.rm = T),
        MAGNITUDE = MAGNITUDE/quantile(MAGNITUDE, .95, na.rm = T),
        ADJUSTED_MAGNITUDE = ADJUSTED_MAGNITUDE/quantile(ADJUSTED_MAGNITUDE, .95, na.rm = T)
      ) %>% 
      group_by(
        PLOT_FACETS, IMAGENUMBER
      ) %>% 
      mutate(
        # To select image with most spots
        IMAGE_N = n(),
        
        # Define plot axis replacing variables with generic names
        ref_int = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        qry_int = ROUNDED_QUERY_TOTAL_INTENSITY,
        
        d_ref_int = (!!as.name(TempVariables[1])),
        d_qry_int = (!!as.name(TempVariables[2])),
        magnitude = (!!as.name(TempVariables[3])),
        rad_angle = (!!as.name(TempVariables[4])),
        deg_angle = (!!as.name(TempVariables[5]))
      ) %>% 
      ungroup() %>% 
      filter(
        # Select image with most bins
        IMAGE_N == max(IMAGE_N)
      ) %>% 
      filter(
        IMAGENUMBER == min(IMAGENUMBER)
      )
    
    PlotTitle <-
      paste0(
        "Cell Line ", TempTable$COHORT[[1]], "\n",
        TempTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2\n",
        TempTable$FPS[[1]], " FPS"
      )
    
    # Plot quiver
    ggplot(
      TempTable
    ) +
      # Make triangles of the arrows
      geom_regon(
        aes(
          x0 = ref_int,
          y0 = qry_int,
          # Size of triangle is the magnitude
          r = (magnitude)/4, # log(N/max(N)+1), # (MAGNITUDE)/10
          # Angle of triangle is the derivative
          angle = -(rad_angle+pi)+2*pi,
          # 3 sides for a triangle
          sides = 3,
          # Angle of the derivatives is the color
          fill = hsv(deg_angle),
          # Redundancy so that the outline of the triangle also tells you the magnitude
          # Also helps pair the arrow with the triangle
          color = N
        ),
        # Make transparent to see what's happening behind
        alpha = 0.75,
        size = 0.1
      ) +
      # Puts arrows inside triangles so that we know which tip of the triangle to look at
      geom_quiver(
        aes(
          x = ref_int,
          y = qry_int,
          u = d_ref_int,
          v = d_qry_int,
          # Arrow color is the magnitude
          # Same as the triangle so that they can be visually paired
          color = N,
          # vecsize = 1/magnitude
        ),
        # Arrow size blowup
        vecsize = 2,
        # Center the arrow
        center = T,
        size = 0.25
      ) +
      # Color of magnitude. Goes from gray to white
      scale_color_gradient(
        low = "white",
        high = "white",
        trans = "log2",
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))
      ) +
      # Use angle of arrows to color the triangle. Full spectrum
      scale_fill_identity(guide = "none") +
      # Log transform axis. The secondary axis reinforces that it's in log as it'll say 2^x
      scale_x_continuous(
        trans = "log2",
        sec.axis = sec_axis(
          trans = ~.,
          breaks = trans_breaks("log2", function(x) 2^x),
          labels = trans_format("log2", math_format(2^.x))
        )
      ) +
      scale_y_continuous(
        trans = "log2",
        sec.axis = sec_axis(
          trans = ~.,
          breaks = trans_breaks("log2", function(x) 2^x),
          labels = trans_format("log2", math_format(2^.x))
        )
      ) +
      labs(
        title = PlotTitle,
        x = "MyD88 Size",
        y = paste(TempTable$QUERY_PROTEIN[1], "Size")
      ) +
      facet_grid(
        IMAGENUMBER+MAX_REFERENCE_TOTAL_INTENSITY_CAT~DELTA_CAT
      ) +
      # Make dark to make the colors pop in a monitor
      # Use theme classic for printing
      dark_theme_classic(
        base_size = 18
      ) +
      theme(
        legend.key.width= unit(1.5, 'cm'),
        # legend.key.width= unit(1.5, 'cm'),
        legend.position = "none"
      ) +
      coord_fixed()
    
    COHORT_FOLDER = gsub("\n", " - ", PlotTitle)
    COHORT_FOLDER = file.path(OUTPUT_DIRECTORY, COHORT_FOLDER)
    
    if(!file.exists(COHORT_FOLDER)){
      dir.create(COHORT_FOLDER)
    }
    
    # File name
    SaveName <-
      paste0(
        PlotList$Variables[FacetX],
        "PhasePortrait - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(1/STEP_SIZE, 2), " - ",
        ifelse(MOVING_AVERAGE_WINDOW > 0, paste("MvAvg", MOVING_AVERAGE_WINDOW, "- "), ""),
        ifelse(LOG_SCALE == T, "Log - ", ""),
        "Triangle.pdf"
      )
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    # Save
    ggsave(
      # Save vector image
      file.path(COHORT_FOLDER, SaveName),
      # height = 4.76*1.25,
      # width = 11.5*1.25
      height = 4.76*2,
      width = 11.5*2
    )
    
    ggplot(
      TempTable
    ) +
      geom_arrow(
        aes(
          x = ref_int,
          y = qry_int,
          mag = .75,
          angle = atan2(d_qry_int, d_ref_int)*180/pi,
          # color = N
          color = hsv(deg_angle)
        ),
        # color = "white",
        size = 1,
        arrow.length = 1
      ) +
      scale_color_identity(guide = "none") +
      scale_mag(
        max = .75,
        guide = 'none'
      ) +
      # Log transform axis. The secondary axis reinforces that it's in log as it'll say 2^x
      scale_x_continuous(
        trans = "log2",
        sec.axis = sec_axis(
          trans = ~.,
          breaks = trans_breaks("log2", function(x) 2^x),
          labels = trans_format("log2", math_format(2^.x))
        )
      ) +
      scale_y_continuous(
        trans = "log2",
        sec.axis = sec_axis(
          trans = ~.,
          breaks = trans_breaks("log2", function(x) 2^x),
          labels = trans_format("log2", math_format(2^.x))
        )
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
    labs(
      title = PlotTitle,
      x = "MyD88 Size",
      y = paste(TempTable$QUERY_PROTEIN[1], "Size")
    ) +
      facet_grid(
        IMAGENUMBER+MAX_REFERENCE_TOTAL_INTENSITY_CAT~DELTA_CAT
      ) +
      # Make dark to make the colors pop in a monitor
      # Use theme classic for printing
      dark_theme_classic(
        base_size = 18
      ) +
      theme(
        legend.key.width= unit(1.5, 'cm'),
        legend.position = "none"
      ) +
      coord_fixed()
    
    
    COHORT_FOLDER = gsub("\n", " - ", PlotTitle)
    COHORT_FOLDER = file.path(OUTPUT_DIRECTORY, COHORT_FOLDER)
    
    if(!file.exists(COHORT_FOLDER)){
      dir.create(COHORT_FOLDER)
    }
    
    # File name
    SaveName <-
      paste0(
        PlotList$Variables[FacetX],
        "PhasePortrait - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(1/STEP_SIZE, 2), " - ",
        ifelse(MOVING_AVERAGE_WINDOW > 0, paste("MvAvg", MOVING_AVERAGE_WINDOW, "- "), ""),
        ifelse(LOG_SCALE == T, "Log - ", ""),
        "Arrow.pdf"
      )
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    
    ggsave(
      # Save vector image
      file.path(COHORT_FOLDER, SaveName),
      # height = 4.76*1.25,
      # width = 11.5*1.25
      height = 4.76*2,
      width = 11.5*2
    )
  }, error = function(e) {print(paste("Error with PlotFx", PlotList$Variables[FacetX], PlotList$Facets[FacetX]))})}
mclapply(1:NROW(PlotList), PlotFx)
