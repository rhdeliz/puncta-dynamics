FRAMES_SINCE_LANDING_CAT_BIN <- c(1:8)*50
# FRAMES_SINCE_LANDING_CAT_BIN <- 100
N_THRESHOLD = 100

setwd(OUTPUT_DIRECTORY)

if(exists("PhasePortrait") == FALSE){
  SaveName <-
    paste0(
      "Normalized PhasePortrait - ",
      "LeadLag ", LEAD_LAG, " - ",
      "StepSize ", STEP_SIZE,
      ".gz.parquet"
    )
  SaveName <- gsub(" - \\.", "\\.", SaveName)
  PhasePortrait <- read_parquet(SaveName)
}

mult_format <- function(l) {
  l <- format(l*1000, digits = 3, scientific = F)
  l <- paste0("'", l, "'%*%10^-3")
  parse(text=l)
}

Facets <- unique(PhasePortrait$PLOT_FACETS)

Variables <-
  c(
    # "MAD_ADJUSTED_",
    # "ADJUSTED_",
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
      PhasePortrait %>%
      filter(
        FRAMES_SINCE_LANDING_CAT %in% FRAMES_SINCE_LANDING_CAT_BIN,
        PLOT_FACETS == PlotList$Facets[FacetX],
        # ROUNDED_REFERENCE_TOTAL_INTENSITY <= 30,
        # ROUNDED_QUERY_TOTAL_INTENSITY <= 50
      ) %>%
      mutate(
        FRAMES_SINCE_LANDING_CAT = round(FRAMES_SINCE_LANDING_CAT/FPS)
      ) %>% 
      as.data.table()
    
    Limits <-
      PhasePortrait %>% 
      filter(
        COHORT == TempTable$COHORT[1],
        LIGAND_DENSITY_CAT == TempTable$LIGAND_DENSITY_CAT[1],
        FPS == TempTable$FPS[1],
        REFERENCE_PROTEIN == TempTable$REFERENCE_PROTEIN[1],
        QUERY_PROTEIN == TempTable$QUERY_PROTEIN[1],
        FRAMES_SINCE_LANDING_CAT %in% FRAMES_SINCE_LANDING_CAT_BIN
      ) %>%
      filter(N_TEST == T) %>% 
      ungroup() %>% 
      mutate(
        DELTA = sqrt(DELTA_QUERY_TOTAL_INTENSITY^2 + DELTA_REFERENCE_TOTAL_INTENSITY^2)
      ) %>% 
      summarise(
        MIN_DELTA = min(DELTA),
        MAX_DELTA = max(DELTA),
        MAX_REF = max(ROUNDED_REFERENCE_TOTAL_INTENSITY) + STEP_SIZE/2,
        MAX_QRY = max(ROUNDED_QUERY_TOTAL_INTENSITY) + STEP_SIZE/2
      ) %>% 
      as.data.table()
    
    TIME_LEVELS = sort(TempTable$FRAMES_SINCE_LANDING_CAT)
    TIME_LEVELS = unique(TIME_LEVELS)
    TIME_LEVELS = paste(TIME_LEVELS, "s")
    
    TempTable <-
      TempTable %>%  
      mutate(
        FRAMES_SINCE_LANDING_CAT = paste(FRAMES_SINCE_LANDING_CAT, "s")
      ) %>% 
      mutate(
        FRAMES_SINCE_LANDING_CAT = factor(FRAMES_SINCE_LANDING_CAT, levels = TIME_LEVELS)
      ) %>% 
      # mutate(
      #   MIN = floor(FRAMES_SINCE_LANDING_CAT/60)
      # ) %>% 
      # mutate(
      #   SEC = round(FRAMES_SINCE_LANDING_CAT - MIN*60)
      # ) %>% 
      # mutate(
      #   MIN = str_pad(MIN, 2, pad = "0"),
      #   SEC = str_pad(SEC, 2, pad = "0")
      # ) %>% 
      # mutate(
    #   FRAMES_SINCE_LANDING_CAT = paste0("00:", MIN, ":", SEC)
    # ) %>% 
    select(-c(
      FPS
    )) %>% 
      group_by(
        PLOT_FACETS, IMAGENUMBER
      ) %>%
      mutate(
        # Define plot axis replacing variables with generic names
        ref_int = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        qry_int = ROUNDED_QUERY_TOTAL_INTENSITY,
        
        d_ref_int = (!!as.name(TempVariables[1])),
        d_qry_int = (!!as.name(TempVariables[2])),
        
        FPS = paste(as.character(PhasePortrait$FPS[[1]]), "Hz")
      ) %>%
      filter(
        # ref_int<= 10,
        # qry_int <= 10
      ) %>% 
      as.data.table()
    
    GreyTempTable <-
      TempTable %>%
      group_by(
        ROUNDED_REFERENCE_TOTAL_INTENSITY,
        ROUNDED_QUERY_TOTAL_INTENSITY
      ) %>% 
      mutate(
        N_TEST = T
        # N_TEST = ifelse(N >= N_THRESHOLD, T, F)
      ) %>% 
      filter(
        N_TEST == FALSE
      ) %>%
      as.data.table()
    
    TempTable <-
      TempTable %>%
      group_by(
        ROUNDED_REFERENCE_TOTAL_INTENSITY,
        ROUNDED_QUERY_TOTAL_INTENSITY
      ) %>% 
      mutate(
        N_TEST = ifelse(N >= N_THRESHOLD, T, F)
      ) %>% 
      filter(
        N_TEST == TRUE
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
    
    StreamTempTable <-
      TempTable %>%
      # bind_rows(TempTable, GreyTempTable) %>%
      as_tibble() %>%
      group_by(
        FRAMES_SINCE_LANDING_CAT
      ) %>%
      complete(
        ref_int = full_seq(c(min(ref_int), max(ref_int)), period = STEP_SIZE),
        qry_int = full_seq(c(min(qry_int), max(qry_int)), period = STEP_SIZE),
      ) %>%
      select(FRAMES_SINCE_LANDING_CAT, ref_int, qry_int, d_ref_int, d_qry_int) %>%
      mutate(
        d_ref_int = ifelse(is.na(d_ref_int), 0, d_ref_int),
        d_qry_int = ifelse(is.na(d_qry_int), 0, d_qry_int)
      ) %>% 
      arrange(
        ref_int,
        qry_int
      ) %>% 
      distinct() %>%
      as.data.table()
    
    ASPECT_RATIO = Limits$MAX_REF/Limits$MAX_QRY
    
    if(NROW(FRAMES_SINCE_LANDING_CAT_BIN) > 1){
      PLOT_NCOL = ifelse(ASPECT_RATIO >= 10, 1, ifelse(ASPECT_RATIO <= 0.75, NROW(FRAMES_SINCE_LANDING_CAT_BIN), round(sqrt(NROW(FRAMES_SINCE_LANDING_CAT_BIN)))))
    } else{
      PLOT_NCOL = 1
    }
    
    ARROW_LENGTH = STEP_SIZE/max(Limits$MAX_REF/2.415966, Limits$MAX_QRY)*8/PLOT_NCOL
    
    if(PLOT_NCOL == NROW(FRAMES_SINCE_LANDING_CAT_BIN)){
      ARROW_LENGTH = STEP_SIZE/max(Limits$MAX_REF, Limits$MAX_QRY)*10
    }
    
    HeatmapTable <-
      bind_rows(TempTable, GreyTempTable) %>% 
      group_by(
        qry_int
      ) %>% 
      mutate(
        QUERY_N = N/sum(N)*100
      ) %>% 
      group_by(
        ref_int
      ) %>% 
      mutate(
        REFERENCE_N = N/sum(N)*100
      ) %>% 
      mutate(
        QUERY_N = ifelse(QUERY_N >= 25, 25, QUERY_N*4),
        REFERENCE_N = ifelse(REFERENCE_N >= 25, 25, REFERENCE_N*4)
      ) %>% 
      mutate(
        FILL = rgb(QUERY_N/100, REFERENCE_N/100, QUERY_N/100)
      ) %>%
      # ungroup() %>% 
      # mutate(
      #   TIME_ADJUSTED = 100-TIME_ADJUSTED/max(TIME_ADJUSTED)*100
      # ) %>% 
      # group_by(
      #   ref_int
      # ) %>% 
      # mutate(
      #   N = N/max(N)*100
      # ) %>% 
      # mutate(
    #   FILL = rgb(TIME_ADJUSTED/100, N/100, TIME_ADJUSTED/100)
    # ) %>%
    as.data.table()
    
    ggplot(
    ) +
      # geom_tile(
      #   data = HeatmapTable,
      #   aes(
      #     x = ref_int,
      #     y = qry_int,
      #     fill = FILL
      #   ),
      #   color = "black"
      # ) +
      # scale_fill_identity()+
      geom_streamline(
        data = StreamTempTable,
        aes(
          x = ref_int,
          y = qry_int,
          dx = d_ref_int,
          dy = d_qry_int,
          color = sqrt(..dx..^2 + ..dy..^2),
          # size = ..step..,
          # alpha = ..step..
        ),
        arrow = NULL,
        size  = .25,
        # n = 20,
        # arrow.length = 0.3,
        # jitter = 4,
        L = 50, res = 10, lineend = "round"
      ) +
      # geom_arrow(
      #   data = GreyTempTable,
      #   aes(
      #     x = ref_int,
      #     y = qry_int,
      #     # dx = d_ref_int,
      #     # dy = d_qry_int,
      #     mag = .5,
      #     angle = atan2(d_qry_int, d_ref_int)*180/pi
      #   ),
      #   color = "black", # "#3A3B3C",
    #   arrow.length = ARROW_LENGTH,
    #   lineend = "square"
    # ) +
    geom_arrow(
      data = TempTable,
      aes(
        x = ref_int,
        y = qry_int,
        # dx = d_ref_int,
        # dy = d_qry_int,
        mag = .5,
        angle = atan2(d_qry_int, d_ref_int)*180/pi,
        color = sqrt(d_qry_int^2 + d_ref_int^2),
        # color = FRAMES_ADJUSTED
      ),
      arrow.length = ARROW_LENGTH,
      lineend = "square"
    ) +
      scale_size(range = c(ARROW_LENGTH/10, ARROW_LENGTH), guide = "none") +
      scale_alpha(guide = "none") +
      scale_mag(
        max = 1/ARROW_LENGTH, # Big arrow head
        guide = 'none'
      ) +
      scale_x_continuous(
        limits = c(-STEP_SIZE/2, Limits$MAX_REF)
        # limits = c(-STEP_SIZE/2, max(HeatmapTable$ref_int)+STEP_SIZE/2)
      ) +
      scale_y_continuous(
        limits = c(-STEP_SIZE/2, Limits$MAX_QRY)
        # limits = c(-STEP_SIZE/2, max(HeatmapTable$qry_int)+STEP_SIZE/2)
      ) +
      # scale_fill_distiller(
      #   palette = "Greys"
      # ) +
      scale_color_viridis(
        # option = "plasma",
        labels = mult_format,
        limits = c(Limits$MIN_DELTA, Limits$MAX_DELTA)
        # direction = -1,
        # trans = "log10",
        # breaks = trans_breaks("log10", function(x) 10^x),
        # labels = trans_format("log10", math_format(10^.x))
      ) +
      labs(
        # title = PlotTitle,
        x = paste(TempTable$REFERENCE_PROTEIN[1], "Norm. Int."),
        y = paste(TempTable$QUERY_PROTEIN[1], "Norm. Int."),
        color = "Growth\n(mol./s)",
        # fill = "Fraction (by TRAF6)"
        # fill = "Age (s)"
        # color = "Scaled\nMagnitude",
        # color = "Time\n(Frames)"
      ) +
      facet_wrap(
        ~FRAMES_SINCE_LANDING_CAT,
        ncol = PLOT_NCOL,
        scales = "free"
      ) +
      # Make dark to make the colors pop in a monitor
      # Use theme classic for printing
      theme_classic(
        # base_size = 20
      ) +
      theme(
        legend.key.height= unit(1, 'cm'),
        legend.position = "right",
        rect = element_rect(fill = "transparent"),
        # aspect.ratio = 1/ASPECT_RATIO,
        strip.background = element_blank()
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
        "Normalized ",
        PlotList$Variables[FacetX],
        "PhasePortrait - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(STEP_SIZE, 2),
        " Replicate ", TempTable$IMAGENUMBER[1],
        ".pdf"
      )
    
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    
    ggsave(
      # Save vector image
      file.path(COHORT_FOLDER, SaveName),
      # height = 4.76*1.25,
      # width = 11.5*1.25
      # height = 4.76*1.5,
      # width = 11.5*1.5
      height = 3*3,
      width = 3*sqrt(2)*3
    )
    
    return(file.path(COHORT_FOLDER, SaveName))
    
  }, error = function(e) {print(paste("Error with PlotFx", PlotList$Variables[FacetX], PlotList$Facets[FacetX]))})
}
PlotList <- mclapply(1:NROW(PlotList), PlotFx, mc.cores = detectCores(logical = F))
# lapply(1:NROW(PlotList), PlotFx)
