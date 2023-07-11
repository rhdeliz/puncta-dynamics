
library(tictoc)

# Numbers <- rnorm(10^7)
Window <- 5
# 
# tic()
# a <- RcppRoll::roll_median(Numbers, n = Window, fill = NA, align = "center")
# toc()
# 
# tic()
# b <- zoo::rollmedian(Numbers, Window, fill = NA, align = "center")
# toc()


Test <- TrackTable
Test <- Test[, SMOOTH_NORMALIZED_INTENSITY5 :=  RcppRoll::roll_median(NORMALIZED_INTENSITY, n = 5, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]
Test <- Test[, SMOOTH_NORMALIZED_INTENSITY7 :=  RcppRoll::roll_median(NORMALIZED_INTENSITY, n = 7, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]
Test <- Test[, SMOOTH_NORMALIZED_INTENSITY9 :=  RcppRoll::roll_median(NORMALIZED_INTENSITY, n = 9, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]


ggplot() +
  geom_path(
    data = Test,
    aes(
      x = FRAMES_ADJUSTED,
      y = NORMALIZED_INTENSITY
    ),
    color = "#505050"
  ) +
  geom_path(
    data = Test,
    aes(
      x = FRAMES_ADJUSTED,
      y = SMOOTH_NORMALIZED_INTENSITY5
    ),
    color = "red"
  ) +
  geom_path(
    data = Test,
    aes(
      x = FRAMES_ADJUSTED,
      y = SMOOTH_NORMALIZED_INTENSITY7
    ),
    color = "green"
  ) +
  geom_path(
    data = Test,
    aes(
      x = FRAMES_ADJUSTED,
      y = SMOOTH_NORMALIZED_INTENSITY9
    ),
    color = "blue"
  ) +
  # scale_y_log10()+
  facet_wrap(
    ~UNIVERSAL_TRACK_ID
  ) +
  dark_theme_classic()








Window <- 5

Test <- fread("/Users/u_deliz/Desktop/Cluster/Output/Calibrations/20190704 GFP calibration 004/Cell_1/GFP_intensity.csv.gz")
Test <- Test[, SMOOTH_NORMALIZED_INTENSITY_MEDIAN :=  RcppRoll::roll_median(TOTAL_INTENSITY, n = 9, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]
Test <- Test[, SMOOTH_NORMALIZED_INTENSITY_MEAN :=  RcppRoll::roll_mean(TOTAL_INTENSITY, n = 9, fill = NA, align = "center"), by = UNIVERSAL_TRACK_ID]

Test <-
  Test %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    FRAMES_ADJUSTED = FRAME - min(FRAME) +1,
    LIFETIME = max(FRAMES_ADJUSTED)
  ) %>% 
  filter(
    LIFETIME >= 50,
    TRACK_ID <= 50
  )


ggplot() +
  geom_path(
    data = Test,
    aes(
      x = FRAMES_ADJUSTED,
      y = TOTAL_INTENSITY,
      color = "Input"
    )
  ) +
  geom_path(
    data = Test,
    aes(
      x = FRAMES_ADJUSTED,
      y = SMOOTH_NORMALIZED_INTENSITY_MEDIAN,
      color = paste("Mv Median", Window)
    )
  ) +
  geom_path(
    data = Test,
    aes(
      x = FRAMES_ADJUSTED,
      y = SMOOTH_NORMALIZED_INTENSITY_MEAN,
      color = paste("Mv Mean", Window)
    )
  ) +
  scale_color_manual(
    values = c("#505050", "green", "magenta")
  ) +
  facet_wrap(
    ~UNIVERSAL_TRACK_ID
  ) +
  dark_theme_classic() +
  theme(
    legend.position =  "bottom"
  )

