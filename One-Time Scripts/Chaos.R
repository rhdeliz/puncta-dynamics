pacman::p_load(data.table, dplyr, tseriesChaos, ggplot2, tidyr, parallel, viridis)
output <-lyap_k(lorenz.ts, m=3, d=2, s=200, t=40, ref=1700, k=2, eps=4)

test <-
  as.data.table(lorenz.ts) %>% 
  mutate(
    y = 1:n()
  ) %>% 
  mutate(
    block = ceiling(y/20)
  ) %>% 
  mutate(
    new_y = 20 - block*20 + y
  )

ggplot(
  test%>% filter(block <= 20),
  aes(
    x = y,
    y = x,
    group = block,
    color = block
  )
) +
  geom_path() +
  geom_vline(
    xintercept = 1:20*20
  )+
  geom_vline(
    xintercept = 1:4*100,
    color ="red"
  )



Table <- fread("/Users/u_deliz/Desktop/PhasePortraitAnalysis/aEssential.csv.gz")


sample_n_groups = function(grouped_df, size, replace = FALSE, weight=NULL) {
  grp_var <- grouped_df %>% 
    groups %>%
    unlist %>% 
    as.character
  random_grp <- grouped_df %>% 
    summarise() %>% 
    sample_n(size, replace, weight) %>% 
    mutate(unique_id = 1:NROW(.))
  grouped_df %>% 
    right_join(random_grp, by=grp_var) %>% 
    group_by_(grp_var) 
}

FilteredTable <-
  Table %>% 
  filter(
    PROTEIN == "MyD88",
    COHORT == "MyD88 RelA",
    IMAGE== "20220517 1.2nM 206-C4_RelA_MyD88_IL1 001",
    FRAMES_ADJUSTED <= 99,
    FRAMES_SINCE_LANDING <= 199
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = n()
  ) %>% 
  filter(
    LIFETIME >= 11
  ) %>% 
  mutate(
    STOICHIOMETRY = NORMALIZED_INTENSITY/(NORMALIZED_INTENSITY + COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  mutate(
    STARTING_STOICHIOMETRY = median(STOICHIOMETRY[1:5], na.rm = T),
    STARTING_NORMALIZED_INTENSITY = median(NORMALIZED_INTENSITY[1:5], na.rm = T)
  ) %>% 
  mutate(
    STARTING_STOICHIOMETRY = round(STARTING_STOICHIOMETRY, 1),
    STARTING_NORMALIZED_INTENSITY = round(STARTING_NORMALIZED_INTENSITY)
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    SCALED_TIME_ADJUSTED = TIME_ADJUSTED/max(TIME_ADJUSTED)
  ) %>% 
  group_by(
    STARTING_STOICHIOMETRY
  ) %>% 
  mutate(
    N = NROW(unique(UNIVERSAL_TRACK_ID))
  ) %>% 
  filter(
    N >= 50
  ) %>% 
  group_split()

SampleFx <- function(df){
  
  ids <- unique(df$UNIVERSAL_TRACK_ID)
  
  if(NROW(ids)>=100){
    df <-
      df %>% 
      group_by(
        UNIVERSAL_TRACK_ID
      ) %>% 
      sample_n_groups(100)
  } else{
    df <-
      df %>% 
      group_by(
        UNIVERSAL_TRACK_ID
      ) %>% 
      sample_n_groups(NROW(ids))
  }
  return(df)
}
SampleTable <- mclapply(FilteredTable, SampleFx)
SampleTable <- rbindlist(SampleTable)

SampleTable %>%
  filter(
    # TIME_ADJUSTED <= 200
  ) %>% 
  ggplot(
    aes(
      # TIME_ADJUSTED,
      SCALED_TIME_ADJUSTED,
      STOICHIOMETRY
    )
  ) +
  geom_path(
    aes(
      group = UNIVERSAL_TRACK_ID
    ),
    alpha = .25,
    size = .25
  ) +
  # stat_density2d(geom="tile", aes(fill = log(..ndensity..+1)), adjust = 1.5, contour = FALSE) +
  facet_wrap(
    ~STARTING_STOICHIOMETRY
    # ~STARTING_NORMALIZED_INTENSITY
  ) +
  scale_color_viridis() +
  scale_fill_viridis() +
  theme_classic()





norm_across_y <- function(v, x, y){
  data.frame(v=v, x=x, y=y) %>%
    group_by(x) %>%
    mutate(v=v/((max(y)-min(y))/n()*sum(v))) %>%
    ungroup() %>%
    pull(v)
}


SampleTable %>%
  ggplot(
    aes(
      SCALED_TIME_ADJUSTED,
      STOICHIOMETRY
    )
  ) +
  # geom_path(
  #   aes(
  #     group = UNIVERSAL_TRACK_ID
  #   ),
  #   alpha = .25,
  #   size = .25
  # ) +
  stat_density_2d_filled(aes(fill=after_stat(norm_across_y(density, x, y))), geom="raster", contour=FALSE, n=500) +
  facet_wrap(
    ~STARTING_STOICHIOMETRY
  ) +
  scale_color_viridis() +
  scale_fill_viridis() +
  theme_classic()



library(nonlinearTseries)

ToStrech <-
  rbindlist(FilteredTable) %>% 
  select(
    UNIVERSAL_TRACK_ID,
    FRAMES_ADJUSTED,
    STOICHIOMETRY
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  group_split()


StretchFx <- function(df){
  df <-
    df %>%
    as_tibble() %>%
    arrange(
      FRAMES_ADJUSTED
    ) %>%
    complete(
      FRAMES_ADJUSTED = full_seq(c(min(0), max(99)), period = 1)
    ) %>%
    # arrange(
    #   SCALED_TIME_ADJUSTED
    # ) %>%
    # mutate(
    #   SCALED_TIME_ADJUSTED = SCALED_TIME_ADJUSTED*99
    # ) %>% 
    # mutate(
    #   SCALED_TIME_ADJUSTED = floor(SCALED_TIME_ADJUSTED)
    # ) %>% 
    # complete(
    #   SCALED_TIME_ADJUSTED = full_seq(c(min(0), max(99)), period = 1)
    # ) %>%
    mutate(
      x = FRAMES_ADJUSTED,
      y = STOICHIOMETRY
    ) %>% 
    select(
      x,
      y
    ) %>% 
    as.data.table()
  
  return(df)
}
StretchedTable <- mclapply(ToStrech, StretchFx)
StretchedTable <- rbindlist(StretchedTable)
remove(ToStrech)

StretchedTable$y[is.na(StretchedTable$y)] = 0
TimeSeries <- ts(StretchedTable$y, frequency = 100)


nonlinearTseries::maxLyapunov(TimeSeries, radius = 2, number.boxes = 1)

df <- NULL
df$x <- rep(1:100, 100)
df$y <- rlogis(100*100)#, 100, 25)
# df$y <- rep(sin(1:100*.1), 100)
df <- as.data.table(df)
df <- ts(df$y, frequency = 100, start = 1, end = 100)

output <-lyap_k(lorenz.ts, m=3, d=2, s=200, t=40, ref=1700, k=2, eps=4)

nonlinearTseries::maxLyapunov(TimeSeries, radius = 100)




AnalysisTable <-
  Table %>% 
  filter(
    PROTEIN == "MyD88",
    COHORT == "MyD88 RelA",
    # IMAGE== "20220517 1.2nM 206-C4_RelA_MyD88_IL1 001",
    FRAMES_ADJUSTED <= 99,
    FRAMES_SINCE_LANDING <= 199
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = n()
  ) %>% 
  filter(
    LIFETIME >= 50
  ) %>% 
  filter(
    FRAMES_ADJUSTED <= 49
  ) %>%
  mutate(
    STOICHIOMETRY = NORMALIZED_INTENSITY/(NORMALIZED_INTENSITY + COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) 

Test <-
  AnalysisTable %>% 
  group_by(
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    MEAN = mean(NORMALIZED_INTENSITY),
    MAD = mad(NORMALIZED_INTENSITY)
  )

plot(Test$MEAN, Test$MAD)



TimeSeries <- ts(AnalysisTable$STOICHIOMETRY, frequency = 50)
# TimeSeries <- ts(AnalysisTable$NORMALIZED_INTENSITY, frequency = 50)
# TimeSeries <- ts(AnalysisTable$COMPLEMENTARY_NORMALIZED_INTENSITY_1, frequency = 50, start = 1)
output <- nonlinearTseries::maxLyapunov(TimeSeries, radius = 1)



output <- lyap_k(TimeSeries, m=3, d=2, s=100, t=100, ref=1062, k=2, eps=1)
output
plot(output)

output <- nonlinearTseries::maxLyapunov(lorenz.ts, radius = 1)
output <-lyap_k(lorenz.ts, m=3, d=2, s=200, t=40, ref=1700, k=2, eps=4)





FFT_Table <-
  AnalysisTable %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  group_split()


FFT_batch <- function(df){
  
  minus_mean <- df$NORMALIZED_INTENSITY - mean(AnalysisTable$NORMALIZED_INTENSITY)
  transformation_table <- fft(minus_mean)
  transformation_table <- Re(transformation_table)
  
  real_table <- NULL
  real_table$fft <- transformation_table
  real_table$fft_shift <- waved::fftshift(transformation_table)
  real_table$fft_time <- real_table$fft_shift * 1/(2*nrow(df)*4)
  
  real_table$id <- df$UNIVERSAL_TRACK_ID[1]
  real_table$t <- df$FRAMES_ADJUSTED*4
  real_table$value <- df$NORMALIZED_INTENSITY
  real_table$minus_mean <- minus_mean
  real_table <- as.data.table(real_table)
  
  return(real_table)
}
table <- mclapply(FFT_Table, FFT_batch)
table <- rbindlist(table)

get_density <-
  table %>% 
  group_by(id) %>% 
  group_split()

PeakFinder <- function(df){
  temp_df <- density(df$fft)
  # temp_df <- density(runif(50, min = -10, max = 10))
  # temp_df <- density(sin(1:50))
  peak <-temp_df$x[max(temp_df$y)==temp_df$y]
  temp_df <- NULL
  temp_df$peak <- peak
  temp_df$id <- df$id[1]
  temp_df <- as.data.table(temp_df)
  return(temp_df)
}
density_results <- mclapply(get_density, PeakFinder)
density_results <- rbindlist(density_results)

density_results %>% 
  ggplot(
    aes(
      x = peak,
      y = log(..count.. +1)
    )
  ) +
  geom_density(
    adjust = 1
  ) +
  scale_x_continuous(
    limits = c(-25, 25)
  ) +
  scale_color_viridis(
    discrete = T, guide = "none"
  ) +
  theme_classic()

plot(density(Re(fft(lorenz.ts))), xlim = c(-1000,1000))


