# Generates 2D interpolation
Bins = 75 #for 2d interpolation
k = 15 #median blur
min_quantile = 0.01
max_quantile = 0.95

TracksTable <-
  Table %>%
  filter(
    !grepl("half", COHORT, ignore.case = T),
    !grepl("DMSO", COHORT, ignore.case = T),
    !grepl("Inhibitor", COHORT, ignore.case = T),
    !grepl("grid", COHORT, ignore.case = T),
    LIGAND_DENSITY_CAT == 32,
    # LIGAND_DENSITY_CAT <= 320,
    COHORT == "MyD88 TRAF6"
  ) %>% 
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    # NORMALIZED_INTENSITY = NORMALIZED_INTENSITY - min(NORMALIZED_INTENSITY),
    # COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1 - min(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
    NORMALIZED_INTENSITY =
      (NORMALIZED_INTENSITY - quantile(NORMALIZED_INTENSITY, min_quantile)) /
      (quantile(COMPLEMENTARY_NORMALIZED_INTENSITY_1, max_quantile) - quantile(NORMALIZED_INTENSITY, min_quantile)),
    
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 =
      (COMPLEMENTARY_NORMALIZED_INTENSITY_1 - quantile(COMPLEMENTARY_NORMALIZED_INTENSITY_1, min_quantile)) /
      (quantile(COMPLEMENTARY_NORMALIZED_INTENSITY_1, max_quantile) - quantile(COMPLEMENTARY_NORMALIZED_INTENSITY_1, min_quantile)),
    
    MAX_TIME_ADJUSTED = max(TIME_ADJUSTED)
  ) %>% 
  filter(
    # MAX_TIME_ADJUSTED >= 200
  ) %>% 
  as.data.table()


NormalizedTracksTable <-
  TracksTable %>% 
  filter(
    NORMALIZED_INTENSITY >= 0,
    NORMALIZED_INTENSITY <= 1,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 0,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 <= 1,
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = round(NORMALIZED_INTENSITY*Bins),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = round(COMPLEMENTARY_NORMALIZED_INTENSITY_1*Bins)
  ) %>% 
  group_by(
    NORMALIZED_INTENSITY,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1
  ) %>% 
  summarize(
    N = n(),
    TIME_ADJUSTED = quantile(TIME_ADJUSTED, .75)
  ) %>% 
  ungroup() %>% 
  as.data.table()
  
ToBlur <- 
  NormalizedTracksTable %>% 
  select(-c(N)) %>%
  complete(
    NORMALIZED_INTENSITY = full_seq(c(min(NORMALIZED_INTENSITY), max(NORMALIZED_INTENSITY)), period = 1),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = full_seq(c(min(COMPLEMENTARY_NORMALIZED_INTENSITY_1), max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)), period = 1),
  ) %>%
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/Bins,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/Bins
  ) %>% 
  as.data.table()

NormalizedTracksTable <-
  NormalizedTracksTable %>% 
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/Bins,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/Bins
  ) %>% 
  as.data.table()

# Applies function circularly and interpolates missing values
# That is, Gaussian blur if using mean
# df must be 3 columns (xy value)
interpolate_2d <- function(df, k, FUN){
  headers <- names(df)
  names(df) <- c("x", "y", "value")
  df <- as.data.table(df)
  df <- df %>% arrange(x,y) %>% mutate(id = paste(x, y, sep = "...")) %>% as.data.table()
  
  delta_x = unique(df$x)
  delta_x = delta_x[2] - delta_x[1]
  delta_y = unique(df$y)
  delta_y = delta_y[2] - delta_y[1]

  circle_mask = (rep(1:k, k) - (k/2+.5))^2 + (rep(1:k, each=k) - (k/2+.5))^2 <= (k/2)^2
  circle_mask = matrix(circle_mask, nrow = k)
  circle_mask = melt(circle_mask)
  names(circle_mask) = c("x", "y", "multiplier")
  circle_mask$x <- circle_mask$x*delta_x - (k/2*delta_x) + delta_x/2
  circle_mask$x <- round(circle_mask$x/delta_x)*delta_x
  circle_mask$y <- circle_mask$y*delta_y - (k/2*delta_y) + delta_y/2
  circle_mask$y <- round(circle_mask$y/delta_y)*delta_y
  circle_mask$multiplier = as.numeric(circle_mask$multiplier)
  
  scan_pixel <- function(row){
    
    temp_circle_mask <-
      circle_mask %>% 
      mutate(
        x = x + df$x[row],
        y = y + df$y[row]
      ) %>% 
      mutate(
        id = paste(x, y, sep = "...")
      ) %>% 
      filter(
        multiplier == 1
      ) %>% 
      as.data.table()
    
    temp_df <-
      df %>% 
      filter(
        id %in% temp_circle_mask$id
      ) %>% 
      as.data.table()
    
    temp_df = temp_df[!is.na(temp_df$value)]
    
    if(NROW(temp_df$value) >0 ){
      new_value <- FUN(temp_df$value)
    } else(
      new_value = NA
    )
    return(new_value)
  }
  new_values <- mclapply(1:NROW(df), scan_pixel)
  new_values <- unlist(new_values)
  df$id = NULL
  df$value = new_values
  names(df) <- headers
  return(df)
}
Blurred <- interpolate_2d(ToBlur, k, median)

ggplot() + 
  geom_tile(
    data = Blurred,
    aes(NORMALIZED_INTENSITY,COMPLEMENTARY_NORMALIZED_INTENSITY_1, fill = TIME_ADJUSTED)
  ) +
  geom_contour(
    data = NormalizedTracksTable %>% filter(N >= 5) %>% as.data.table(),
    aes(NORMALIZED_INTENSITY, COMPLEMENTARY_NORMALIZED_INTENSITY_1, z = log(N, 10)),
    color = "black",
    alpha = .25
  ) +
  # geom_point(
  #   data = IsolationPoint,
  #   aes(.88,.1)
  # ) +
  scale_fill_viridis(
    option = "turbo"
  ) +
  labs(
    title = "Stability Points",
    subtitle = paste(Bins, "bins with", k, "bin median blur"),
    x = "Scaled MyD88 Size",
    y = "Scaled TRAF6 Size",
    fill = "Average\nCluster\nAge (s)"
  ) +
  theme_classic() +
  coord_fixed()




IsolationPoint <- NULL
IsolationPoint$NORMALIZED_INTENSITY = .88 # .15 .85
IsolationPoint$COMPLEMENTARY_NORMALIZED_INTENSITY_1 = .1 # .45 .55
IsolationPoint <- as.data.table(IsolationPoint)

IncludeTracks <- 
  TracksTable %>% 
  filter(
    NORMALIZED_INTENSITY >= IsolationPoint$NORMALIZED_INTENSITY - .075,
    NORMALIZED_INTENSITY <= IsolationPoint$NORMALIZED_INTENSITY + .075,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= IsolationPoint$COMPLEMENTARY_NORMALIZED_INTENSITY_1 - .075,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 <= IsolationPoint$COMPLEMENTARY_NORMALIZED_INTENSITY_1 + .075
  ) %>% 
  select(
    UNIVERSAL_TRACK_ID
  ) %>% 
  distinct() %>% 
  as.data.table()


TracksTable %>% 
  filter(
    # UNIVERSAL_TRACK_ID %in% IncludeTracks$UNIVERSAL_TRACK_ID
  ) %>% 
  as.data.table() %>% 
  ggplot() +
  geom_path(
    aes(NORMALIZED_INTENSITY,COMPLEMENTARY_NORMALIZED_INTENSITY_1, group = UNIVERSAL_TRACK_ID),
    size = .25,
    alpha = .25
  ) +
  # geom_point(
  #   data = IsolationPoint,
  #   aes(NORMALIZED_INTENSITY,COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  # ) +
  # geom_point(
  #   data = IsolationPoint,
  #   aes(NORMALIZED_INTENSITY-.05,COMPLEMENTARY_NORMALIZED_INTENSITY_1-.05)
  # ) +
  labs(
    title = "Tracks Passing Trough Point",
    x = "Scaled MyD88 Size",
    y = "Scaled TRAF6 Size",
    fill = "Average\nCluster\nAge (s)"
  ) +
  scale_x_continuous(
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  theme_classic() +
  coord_fixed()





ggplot(
) + 
  geom_tile(
    data = Blurred,
    aes(NORMALIZED_INTENSITY,COMPLEMENTARY_NORMALIZED_INTENSITY_1, fill = TIME_ADJUSTED)
  ) +
  geom_path(
    data = TracksTable,
    aes(NORMALIZED_INTENSITY,COMPLEMENTARY_NORMALIZED_INTENSITY_1, group = UNIVERSAL_TRACK_ID),
    size = .25,
    alpha = .1
  ) +
  # geom_point(
  #   data = IsolationPoint,
  #   aes(.88,.1)
  # ) +
  scale_x_continuous(
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  scale_fill_viridis(
    option = "viridis"
  ) +
  labs(
    title = "Stability Points",
    subtitle = "50 bins with 5 bin median blur",
    x = "Scaled MyD88 Size",
    y = "Scaled TRAF6 Size",
    fill = "Average\nCluster\nAge (s)"
  ) +
  theme_classic() +
  coord_fixed()
