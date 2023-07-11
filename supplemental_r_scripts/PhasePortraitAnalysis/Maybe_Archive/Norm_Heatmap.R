setwd(OUTPUT_DIRECTORY)

Facets <- unique(LinearPhasePortrait$PLOT_FACETS)

FRAMES_SINCE_LANDING_CAT_BINS <- c(50, 100)

Variables <-
  c(
    "N",
    "TIME_ADJUSTED",
    # "ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY",
    # "ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY",
    "DELTA_QUERY_TOTAL_INTENSITY",
    "DELTA_REFERENCE_TOTAL_INTENSITY"
  )

PlotList <- NULL
PlotList$Facets = rep(Facets, each = NROW(Variables))
PlotList$Variables = rep(Variables, NROW(Facets))
PlotList <- as.data.table(PlotList)
PlotList <- PlotList %>% arrange(Facets) %>% as_tibble()
remove(Facets, Variables)

PlotFx <- function(FacetX){
  tryCatch({
    
    TempVariables <- PlotList$Variables[FacetX]
    
    TempTable <-
      LinearPhasePortrait %>%
      filter(
        N >= 50
      ) %>% 
      as_tibble() %>% 
      select(-c(
        FPS
      )) %>% 
      filter(
        FRAMES_SINCE_LANDING_CAT %in% FRAMES_SINCE_LANDING_CAT_BINS,
        PLOT_FACETS == PlotList$Facets[FacetX]
      ) %>%
      group_by(
        PLOT_FACETS, IMAGENUMBER
      ) %>%
      mutate(
        # Define plot axis replacing variables with generic names
        ref_int = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        qry_int = ROUNDED_QUERY_TOTAL_INTENSITY,
        fill = (!!as.name(TempVariables)),
        FPS = paste(as.character(LinearPhasePortrait$FPS[[1]]), "Hz")
      )
    
    if(TempVariables == "N"){
      TempTable <-
        TempTable %>%
        as_tibble() %>% 
        mutate(
          fill = log(fill+1)
        ) %>% 
        as.data.table()
    }
    
    PlotTitle <-
      paste0(
        "Cell Line ", TempTable$COHORT[[1]], " - ",
        TempTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
        TempTable$FPS[[1]],
        " - Replicate ",TempTable$IMAGENUMBER[[1]]
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
          .(TempTable$FPS[[1]]) *
          " ~ Rep. "*
          .(TempTable$IMAGENUMBER[[1]])
      )
    
      ggplot(
      ) +
        geom_tile(
          data = TempTable,
          aes(
            x = ref_int,
            y = qry_int,
            fill = fill
          )
        ) +
        # scale_fill_gradient2(midpoint = 0, low = "#2c7bb6", mid = "#ffffbf",
        #                       high = "#d7191c", space = "Lab" )+
        # 
        scale_fill_viridis(
          option = "plasma",
          # direction = -1,
          # trans = "log10",
          # breaks = trans_breaks("log10", function(x) 10^x),
          # labels = trans_format("log10", math_format(10^.x))
        ) +
        labs(
          title = PlotTitle,
          x = paste(TempTable$REFERENCE_PROTEIN[1], "Norm. Int."),
          y = paste(TempTable$QUERY_PROTEIN[1], "Norm. Int."),
          fill = TempVariables
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
        TempVariables, " - ",
        "PhasePortrait - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(STEP_SIZE, 2), " - ",
        "Replicate ",TempTable$IMAGENUMBER[[1]],
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
