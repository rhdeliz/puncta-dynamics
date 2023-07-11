TIME_BINS = 3

DBSCAN =  F
N_CLUSTERS = 10
D_DISTANCE_QUANTILE = .75

SplitNormStatTable <-
  NormStatTable %>% 
  filter(
    COHORT %in% c("MyD88 TRAF6", "MyD88 NEMO", "MyD88 HOIL1")
    # COHORT == "MyD88 TRAF6",
    # LIGAND_DENSITY_CAT == 32
  ) %>% 
  as_tibble() %>% 
  group_split(
    COHORT,
    LIGAND_DENSITY_CAT
  )

CohortLigandFx <- function(df){
  
  ScaledTracks <-
    df %>% 
    group_by(
      UNIVERSAL_TRACK_ID,
      FRAMES_ADJUSTED
    ) %>% 
    mutate(
      N = n()
    ) %>% 
    filter(
      N == 1
    ) %>% 
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>% 
    mutate(
      MAX_REFERENCE = max(REFERENCE_TOTAL_INTENSITY),
      MAX_QUERY = max(QUERY_TOTAL_INTENSITY),
      N = n(),
      STARTING_INTENSITY = median(REFERENCE_TOTAL_INTENSITY[1:3]) + median(QUERY_TOTAL_INTENSITY[1:3])
    ) %>% 
    filter(
      MAX_REFERENCE >= 2,
      MAX_QUERY >= 2,
      N > 21,
      STARTING_INTENSITY <= 6,
      FRAMES_ADJUSTED <= FRAMES_CUTOFF
    ) %>% 
    filter(
      MAX_REFERENCE >= 6 | MAX_QUERY >= 6
    ) %>% 
    arrange(
      FRAMES_ADJUSTED
    ) %>% 
    mutate(
      DELTA_REFERENCE = lead(REFERENCE_TOTAL_INTENSITY, LEAD_LAG) - lag(REFERENCE_TOTAL_INTENSITY, LEAD_LAG),
      DELTA_QUERY = lead(QUERY_TOTAL_INTENSITY, LEAD_LAG) - lag(QUERY_TOTAL_INTENSITY, LEAD_LAG)
    ) %>%
    drop_na(
      DELTA_REFERENCE,
      DELTA_QUERY
    ) %>%
    arrange(
      FRAMES_ADJUSTED
    ) %>% 
    mutate(
      SCALED_REFERENCE_TOTAL_INTENSITY = (REFERENCE_TOTAL_INTENSITY-min(REFERENCE_TOTAL_INTENSITY))/(max(REFERENCE_TOTAL_INTENSITY)-min(REFERENCE_TOTAL_INTENSITY)),
      SCALED_QUERY_TOTAL_INTENSITY = (QUERY_TOTAL_INTENSITY-min(QUERY_TOTAL_INTENSITY))/(max(QUERY_TOTAL_INTENSITY)-min(QUERY_TOTAL_INTENSITY)),
      SCALED_FRAMES_ADJUSTED = (FRAMES_ADJUSTED-min(FRAMES_ADJUSTED))/(max(FRAMES_ADJUSTED)-min(FRAMES_ADJUSTED)),
      QRY_FUNCTION = cumsum(DELTA_QUERY)
    ) %>% 
    arrange(
      -FRAMES_ADJUSTED
    ) %>% 
    mutate(
      REF_FUNCTION = cumsum(DELTA_REFERENCE)
    ) %>% 
    arrange(
      FRAMES_ADJUSTED
    ) %>% 
    mutate(
      SCALED_FRAMES_ADJUSTED = round(SCALED_FRAMES_ADJUSTED*TIME_BINS)/TIME_BINS
    ) %>% 
    group_by(
      COHORT,
      LIGAND_DENSITY_CAT,
      UNIVERSAL_TRACK_ID,
      SCALED_FRAMES_ADJUSTED
    ) %>% 
    summarize(
      REF_FUNCTION = mean(REF_FUNCTION),
      QRY_FUNCTION = mean(QRY_FUNCTION),
      
      SD_REFERENCE_TOTAL_INTENSITY = sd(REFERENCE_TOTAL_INTENSITY)/sqrt(n()),
      SD_QUERY_TOTAL_INTENSITY = sd(QUERY_TOTAL_INTENSITY)/sqrt(n()),
      
      SCALED_REFERENCE_TOTAL_INTENSITY = mean(SCALED_REFERENCE_TOTAL_INTENSITY),
      SCALED_QUERY_TOTAL_INTENSITY = mean(SCALED_QUERY_TOTAL_INTENSITY),
      
      REFERENCE_TOTAL_INTENSITY = mean(REFERENCE_TOTAL_INTENSITY),
      QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY),
      
      N = NROW(unique(UNIVERSAL_TRACK_ID))
    ) %>% 
    as.data.table()
  
  TrackIDs <- unique(ScaledTracks$UNIVERSAL_TRACK_ID)
  IndexIDs <- 1:NROW(TrackIDs)
  
  TrackIDs <- t(combn(TrackIDs, 2, simplify = T))
  IndexIDs <- t(combn(IndexIDs, 2, simplify = T))
  
  ScaledTracks <-
    ScaledTracks %>% 
    as_tibble() %>% 
    group_split(UNIVERSAL_TRACK_ID)
  
  TrackDistance <- function(PairX){
    Track1 <- ScaledTracks[IndexIDs[PairX, 1]][[1]]
    Track2 <- ScaledTracks[IndexIDs[PairX, 2]][[1]]
    
    # x_dist <- sqrt(sum(c(Track1$REFERENCE_TOTAL_INTENSITY - Track2$REFERENCE_TOTAL_INTENSITY)^2))
    # y_dist <- sqrt(sum(c(Track1$QUERY_TOTAL_INTENSITY - Track2$QUERY_TOTAL_INTENSITY)^2))
    
    # x_dist <- sqrt(sum(c(Track1$SCALED_REFERENCE_TOTAL_INTENSITY - Track2$SCALED_REFERENCE_TOTAL_INTENSITY)^2))
    # y_dist <- sqrt(sum(c(Track1$SCALED_QUERY_TOTAL_INTENSITY - Track2$SCALED_QUERY_TOTAL_INTENSITY)^2))
    
    x_dist <- sqrt(sum(c(Track1$REF_FUNCTION - Track2$REF_FUNCTION)^2))
    y_dist <- sqrt(sum(c(Track1$QRY_FUNCTION - Track2$QRY_FUNCTION)^2))
    # 
    result <- NULL
    result$x <- Track1$UNIVERSAL_TRACK_ID[1]
    result$y <- Track2$UNIVERSAL_TRACK_ID[1]
    result$dist <- x_dist + y_dist
    result <- as.data.table(result)
    
    return(result)
  }
  DistanceTable <- mclapply(1:NROW(IndexIDs), TrackDistance, mc.cores = 8)
  DistanceTable <- rbindlist(DistanceTable)
  
  DistanceTable2 <- NULL
  DistanceTable2$x <- DistanceTable$y
  DistanceTable2$y <- DistanceTable$x
  DistanceTable2$dist <- DistanceTable$dist
  
  DistanceTable <- rbind(DistanceTable, DistanceTable2)
  remove(DistanceTable2)
  
  DistanceTable2 <- NULL
  DistanceTable2$x <- unique(DistanceTable$x)
  DistanceTable2$y <- unique(DistanceTable$x)
  DistanceTable2$dist <- 0
  
  DistanceTable <- rbind(DistanceTable, DistanceTable2)
  remove(DistanceTable2)
  
  DistanceTable <-
    DistanceTable %>%
    arrange(x, y) %>%
    as.data.table()
  
  DistanceTableMatrix <-
    matrix(
      DistanceTable$dist,
      nrow = NROW(unique(DistanceTable$x)),
      dimnames = list(
        unique(DistanceTable$x),
        unique(DistanceTable$y)
      )
    )
  
  ScaledTracks <- rbindlist(ScaledTracks)
  
  D_DISTANCE = as.numeric(quantile(DistanceTableMatrix[DistanceTableMatrix!=0], D_DISTANCE_QUANTILE))
  
  if(DBSCAN == T){
    dbscan_result <- fpc::dbscan(DistanceTableMatrix, D_DISTANCE) #method = "dist",
    dbscan_group <- predict(dbscan_result)
    
    ID_Cluster <- NULL
    ID_Cluster$UNIVERSAL_TRACK_ID <- unique(DistanceTable$x)
    ID_Cluster$GROUP_ID <- dbscan_group
    
  } else{
    k_result <- kmeans(DistanceTableMatrix, N_CLUSTERS)
    
    ID_Cluster <- NULL
    ID_Cluster$UNIVERSAL_TRACK_ID <- names(k_result$cluster)
    ID_Cluster$GROUP_ID <- as.numeric(k_result$cluster)
    
  }
  ID_Cluster <- as.data.table(ID_Cluster)
  
  ClusterVisualization2 <- merge(ScaledTracks, ID_Cluster, by = "UNIVERSAL_TRACK_ID")
  
  # 
  # DistanceTable %>%
  #   arrange(dist, x, y) %>%
  #   mutate(
  #     x = factor(x, levels = unique(DistanceTable$x)),
  #     y = factor(y, levels = unique(DistanceTable$x))
  #   ) %>%
  #   as.data.table() %>%
  # ggplot(
  #   aes(
  #     x,y, fill = dist
  #   )
  # ) +
  #   geom_tile() +
  #   scale_fill_viridis(
  #     option = "turbo"
  #   )+
  #   theme_void()
  
  TempClusterVisualization <-
    ClusterVisualization2 %>% 
    group_by(
      COHORT,
      LIGAND_DENSITY_CAT,
      SCALED_FRAMES_ADJUSTED,
      GROUP_ID
    ) %>% 
    summarize(
      REF_FUNCTION = mean(REF_FUNCTION),
      QRY_FUNCTION = mean(QRY_FUNCTION),
      
      SD_REFERENCE_TOTAL_INTENSITY = sd(REFERENCE_TOTAL_INTENSITY)/sqrt(n()),
      SD_QUERY_TOTAL_INTENSITY = sd(QUERY_TOTAL_INTENSITY)/sqrt(n()),
      
      SCALED_REFERENCE_TOTAL_INTENSITY = mean(SCALED_REFERENCE_TOTAL_INTENSITY),
      SCALED_QUERY_TOTAL_INTENSITY = mean(SCALED_QUERY_TOTAL_INTENSITY),
      
      REFERENCE_TOTAL_INTENSITY = mean(REFERENCE_TOTAL_INTENSITY),
      QUERY_TOTAL_INTENSITY = mean(QUERY_TOTAL_INTENSITY),
      N = n()
    ) %>% 
    group_by(
      COHORT,
      LIGAND_DENSITY_CAT,
      GROUP_ID
    ) %>% 
    mutate(
      N = ifelse(SCALED_FRAMES_ADJUSTED == min(SCALED_FRAMES_ADJUSTED), N, 0),
      
      SCALED_X_LABEL = SCALED_REFERENCE_TOTAL_INTENSITY[n()],
      SCALED_Y_LABEL = SCALED_QUERY_TOTAL_INTENSITY[n()],
      
      X_LABEL = REFERENCE_TOTAL_INTENSITY[n()],
      Y_LABEL = QUERY_TOTAL_INTENSITY[n()]
    ) %>% 
    group_by(
      COHORT,
      LIGAND_DENSITY_CAT
    ) %>% 
    mutate(
      N_TOTAL = sum(N)
    ) %>% 
    group_by(
      COHORT,
      LIGAND_DENSITY_CAT,
      GROUP_ID
    ) %>% 
    mutate(
      N = round(sum(N)/N_TOTAL*100, 1),
      COHORT = factor(COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER))
    ) %>% 
    as.data.table()
  
  return(TempClusterVisualization)
}
ClusterVisualization <- lapply(SplitNormStatTable, CohortLigandFx)
ClusterVisualization <- rbindlist(ClusterVisualization)


ggplot(
  data = ClusterVisualization,
  aes(
    # SCALED_REFERENCE_TOTAL_INTENSITY,
    # SCALED_QUERY_TOTAL_INTENSITY,
    
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY,
    group = GROUP_ID, 
    color = GROUP_ID,
    size = N
  )
) +
geom_errorbar(
  aes(
    ymin = QUERY_TOTAL_INTENSITY - SD_QUERY_TOTAL_INTENSITY,
    ymax = QUERY_TOTAL_INTENSITY + SD_QUERY_TOTAL_INTENSITY
  ),
  # color = "black",
  size = 1
) +
geom_errorbar(
  aes(
    xmin = REFERENCE_TOTAL_INTENSITY -SD_REFERENCE_TOTAL_INTENSITY,
    xmax = REFERENCE_TOTAL_INTENSITY +SD_REFERENCE_TOTAL_INTENSITY,
  ),
  # color = "black",
  size = 1
) +
  geom_path() +
  geom_text(
    # aes(x = SCALED_X_LABEL, y = SCALED_Y_LABEL, label = paste(N, "%")),
    aes(x = X_LABEL, y = Y_LABEL, label = paste(N, "%")),
    color = "black"
  ) +
  scale_color_viridis(
    option = "turbo",
    guide = "none"
  ) +
  labs(
    x = "Normalized Reference Intensity",
    y = "Normalized Query Intensity"
  ) +
  facet_wrap(
    ~COHORT,
    scales = "free"
  ) +
  # coord_fixed() +
  theme_classic()

