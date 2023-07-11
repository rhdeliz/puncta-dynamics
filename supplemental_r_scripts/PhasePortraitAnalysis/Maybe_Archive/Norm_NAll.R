setwd(OUTPUT_DIRECTORY)

Facets <- unique(LinearPhasePortrait$PLOT_FACETS)

FRAMES_SINCE_LANDING_CAT_BINS <- c(50, 100)

Variables <-
  c(
    # "MAD_ADJUSTED_",
    "ADJUSTED_",
    ""
  )

PlotList <- NULL
PlotList$Facets = rep(Facets, each = NROW(Variables))
PlotList$Variables = rep(Variables, NROW(Facets))
PlotList <- as.data.table(PlotList)
PlotList <- PlotList %>% arrange(Facets) %>% as_tibble()
remove(Facets, Variables)

PlotFx <- function(FacetX){
  tryCatch({
    
    TempVariables <-
      paste0(
        PlotList$Variables[FacetX],
        c(
          "DELTA_REFERENCE_TOTAL_INTENSITY",
          "DELTA_QUERY_TOTAL_INTENSITY"
          # "MAGNITUDE",
          # "RAD_ANGLE",
          # "DEG_ANGLE"
        )
      )
    
    TempTable <-
      LinearPhasePortrait %>%
      select(-c(
        FPS
      )) %>% 
      filter(
        FRAMES_SINCE_LANDING_CAT %in% FRAMES_SINCE_LANDING_CAT_BINS,
        PLOT_FACETS == PlotList$Facets[FacetX]
      ) %>%
      # group_by(
      #   COHORT, LIGAND_DENSITY_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER, IMAGE,
      #   FACET, PLOT_FACETS
      # ) %>%
      # mutate(
      #   MAD_MAGNITUDE = MAD_MAGNITUDE/max(MAD_MAGNITUDE, na.rm = T),
      #   MAD_ADJUSTED_MAGNITUDE = MAD_ADJUSTED_MAGNITUDE/max(MAD_ADJUSTED_MAGNITUDE, na.rm = T),
      #   MAGNITUDE = MAGNITUDE/max(MAGNITUDE, na.rm = T),
      # ) %>%
      group_by(
        PLOT_FACETS
      ) %>%
      mutate(
        # Define plot axis replacing variables with generic names
        ref_int = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        qry_int = ROUNDED_QUERY_TOTAL_INTENSITY,
        
        d_ref_int = (!!as.name(TempVariables[1])),
        d_qry_int = (!!as.name(TempVariables[2])),
        # magnitude = (!!as.name(TempVariables[3])),
        # rad_angle = (!!as.name(TempVariables[4])),
        # deg_angle = (!!as.name(TempVariables[5])),
        
        FPS = paste(as.character(LinearPhasePortrait$FPS[[1]]), "Hz")
      ) %>% 
      as.data.table()
    
    PlotTitle <-
      paste0(
        "Cell Line ", TempTable$COHORT[[1]], " - ",
        TempTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
        TempTable$FPS[[1]]
      )
    
    PlotTitle <-
      bquote(
        "Cell Line " *
          .(TempTable$COHORT[[1]]) *
          " ~ " *
          .( TempTable$LIGAND_DENSITY_CAT[[1]]) *
          " mol. " *
          µm^2 *
          " ~ " *
          .(TempTable$FPS[[1]])
      )
    
    ggplot(
    ) +
      geom_tile(
        data = TempTable,
        aes(
          x = ref_int,
          y = qry_int,
          mag = 0.5,
          fill = log(N+1),
          # color = FRAMES_ADJUSTED
        ),
        # size = .01, # Small arrow head
        # arrow.length = .5, # Small arrow head
        
        size = 1, # Big arrow head
        arrow.length = 1,# Big arrow head
        lineend = "square"
      ) +
      scale_size(range = c(.01, 1), guide = "none") +
      scale_alpha(guide = "none") +
      scale_mag(
        # max = 3, # Small arrow head
        max = 1, # Big arrow head
        guide = 'none'
      ) +
      scale_fill_viridis(
        option = "plasma",
        # direction = -1,
        # trans = "log10",
        # breaks = trans_breaks("log10", function(x) 10^x),
        # labels = trans_format("log10", math_format(10^.x))
      ) +
      labs(
        # title = PlotTitle,
        x = paste(TempTable$REFERENCE_PROTEIN[1], "Norm. Int."),
        y = paste(TempTable$QUERY_PROTEIN[1], "Norm. Int."),
        color = "Magnitude"
        # color = "Scaled\nMagnitude",
        # color = "Time\n(Frames)"
      ) +
      facet_grid(
        ~FRAMES_SINCE_LANDING_CAT
      ) +
      # Make dark to make the colors pop in a monitor
      # Use theme classic for printing
      theme_classic(
        base_size = 20
      ) +
      theme(
        legend.key.height= unit(1.25, 'cm'),
        legend.position = "right",
        rect = element_rect(fill = "transparent")
      ) +
      coord_fixed(
        # xlim = c(0, max(TempTable$ref_int, GreyTempTable$ref_int)),
        # ylim = c(0, max(TempTable$qry_int, GreyTempTable$qry_int))
      )
    
    
    # Save folder
    COHORT_FOLDER = file.path(OUTPUT_DIRECTORY,
                              paste0(
                                "Cell Line ", TempTable$COHORT[[1]], " - ",
                                TempTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
                                TempTable$FPS[[1]]
                              ))
    
    if(!file.exists(COHORT_FOLDER)){
      dir.create(COHORT_FOLDER)
    }
    
    # File name
    SaveName <-
      paste0(
        "Normalized All ",
        PlotList$Variables[FacetX],
        "PhasePortrait - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(STEP_SIZE, 2),
        ".pdf"
      )
    
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    
    ggsave(
      # Save vector image
      file.path(COHORT_FOLDER, SaveName),
      # height = 4.76*1.25,
      # width = 11.5*1.25
      height = 4.76*1.5,
      width = 11.5*1.5
    )
    
  }, error = function(e) {print(paste("Error with PlotFx", PlotList$Variables[FacetX], PlotList$Facets[FacetX]))})}
mclapply(1:NROW(PlotList), PlotFx, mc.cores = detectCores(logical = F))
