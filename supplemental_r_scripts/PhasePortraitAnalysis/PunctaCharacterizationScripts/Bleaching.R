Table1 <- fread("/Volumes/taylor-lab/Rafael/NewPipeline/Calibrations/20190910 mscarlet calibration007/Cell_1/mScarlet_intensity.csv.gz")
Table2 <- fread("/Volumes/taylor-lab/Rafael/NewPipeline/Calibrations/20190910 gfp calibration005/Cell_1/GFP_intensity.csv.gz")

Table <- rbindlist(list(Table1, Table2))

SpotsPerTime <-
  Table %>% 
  group_by(
    PROTEIN,
    FRAME
  ) %>% 
  summarize(
    N = n()
  ) %>% 
  group_by(
    PROTEIN
  ) %>% 
  mutate(
    N = N/max(N)
  )

ggplot(
  SpotsPerTime,
  aes(
    x = FRAME,
    y = N,
    color = PROTEIN,
    group = PROTEIN
  )
) +
  geom_hline(
    yintercept = 0.5
  ) +
  scale_y_log10()+
  geom_path() +
  labs(
    x = "Frame",
    y = "Spots (%)",
    color = "Protein"
  ) +
  scale_color_manual(
    values = c("green", "magenta")
  ) +
  dark_theme_classic()



ggsave(
  "Bleaching.pdf",
  height = 4.76,
  width = 11.5/2
)

IntensityPerTime <-
  Table %>% 
  group_by(
    PROTEIN,
    FRAME
  ) %>% 
  summarize(
    TOTAL_INTENSITY = median(TOTAL_INTENSITY)
  )

ggplot(
  IntensityPerTime,
  aes(
    x = FRAME,
    y = TOTAL_INTENSITY,
    color = PROTEIN,
    group = PROTEIN
  )
) +
  geom_path() +
  labs(
    x = "Frame",
    y = "Total Intensity (a.u.)",
    color = "Protein"
  ) +
  scale_color_manual(
    values = c("green", "magenta")
  ) +
  dark_theme_classic()



ggsave(
  "Intensity.pdf",
  height = 4.76,
  width = 11.5
)