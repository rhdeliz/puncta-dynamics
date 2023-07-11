pacman::p_load(data.table, dplyr, dtplyr, arrow, ggplot2, viridis, tidyr, zoo, signal)

filter <- dplyr::filter

setwd("/Users/u_deliz/Desktop/time_analysis")

Table <- read_parquet(
  "Essential.gz.parquet",
  col_select = c(
    "COHORT",
    "LIGAND_DENSITY_CAT",
    "IMAGE",
    "CELL",
    "UNIVERSAL_TRACK_ID",
    "PROTEIN",
    "NORMALIZED_INTENSITY",
    "COMPLEMENTARY_NORMALIZED_INTENSITY_1",
    "FRAME",
    "TIME_ADJUSTED"
  ))

Table <-
  Table %>% 
  filter(
    PROTEIN == "MyD88"
  ) %>% 
  group_by(
    IMAGE,
    CELL
  ) %>% 
  mutate(
    FRAMES_SINCE_LANDING = FRAME - min(FRAME),
  ) %>% 
  filter(
    FRAMES_SINCE_LANDING <= 199
  ) %>% 
  select(-c(
    FRAMES_SINCE_LANDING
  )) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    FRAMES_ADJUSTED = FRAME - min(FRAME)
  ) %>% 
  filter(
    FRAMES_ADJUSTED <= 99
  ) %>% 
  mutate(
    LIFETIME = max(FRAMES_ADJUSTED)
  ) %>% 
  filter(
    LIFETIME >= 11
  ) %>% 
  select(-c(
    FRAME
  )) %>% 
  as.data.table()


SummaryTable <-
  Table %>%
  filter(
    !grepl("half", COHORT, ignore.case = T),
    !grepl("DMSO", COHORT, ignore.case = T),
    !grepl("Inhibitor", COHORT, ignore.case = T),
    !grepl("grid", COHORT, ignore.case = T),
    LIGAND_DENSITY_CAT == 32,
    # LIGAND_DENSITY_CAT <= 320,
    # COHORT == "MyD88 HOIL1"
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    COHORT,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    NORMALIZED_INTENSITY = median(NORMALIZED_INTENSITY),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = median(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
    TIME_ADJUSTED = median(TIME_ADJUSTED),
  ) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = sgolayfilt(NORMALIZED_INTENSITY, p = 1, n = 5),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = sgolayfilt(COMPLEMENTARY_NORMALIZED_INTENSITY_1, p = 1, n = 5)
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = signal::sgolayfilt(NORMALIZED_INTENSITY, p = 3, n = 15),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = sgolayfilt(COMPLEMENTARY_NORMALIZED_INTENSITY_1, p = 3, n = 15)
  ) %>% 
  drop_na() %>% 
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY - min(NORMALIZED_INTENSITY),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1 - min(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/max(NORMALIZED_INTENSITY),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/max(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
    MAX_TIME_ADJUSTED = max(TIME_ADJUSTED)
  ) %>% 
  filter(
    MAX_TIME_ADJUSTED >= 200
  ) %>% 
  as.data.table()

ggplot(
  SummaryTable
) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = NORMALIZED_INTENSITY,
      color = "MyD88"
    )
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = COMPLEMENTARY_NORMALIZED_INTENSITY_1,
      color = "Query"
    )
  ) +
  scale_color_viridis(
    option = "turbo",
    # trans = "log10"
    discrete = T
  ) +
  facet_wrap(
    ~COHORT
  ) +
  theme_classic()






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
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY - quantile(NORMALIZED_INTENSITY, 0.05),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1 - quantile(COMPLEMENTARY_NORMALIZED_INTENSITY_1, 0.05)
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/quantile(NORMALIZED_INTENSITY, 0.95),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/quantile(COMPLEMENTARY_NORMALIZED_INTENSITY_1, 0.95),
    MAX_TIME_ADJUSTED = max(TIME_ADJUSTED)
  ) %>% 
  filter(
    MAX_TIME_ADJUSTED >= 200
  ) %>% 
  as.data.table()

TracksTable %>% 
  filter(
    # UNIVERSAL_TRACK_ID %in% unique(TracksTable$UNIVERSAL_TRACK_ID)[1:100]
  ) %>%
  as.data.table() %>% 
  ggplot(
    aes(
      NORMALIZED_INTENSITY,
      COMPLEMENTARY_NORMALIZED_INTENSITY_1,
      # NORMALIZED_INTENSITY/(COMPLEMENTARY_NORMALIZED_INTENSITY_1+NORMALIZED_INTENSITY),
      # group = UNIVERSAL_TRACK_ID,
    )
  ) +
  # geom_tile(
  #   aes(
  #     fill = TIME_ADJUSTED
  #   )
  # ) +
  stat_density2d(geom="tile", aes(fill = ..count..), adjust = 3, contour = FALSE) +
  # scale_fill_viridis(
  #   option = "turbo",
  #   # trans = "log10"
  # ) +
  # geom_path(
  #   size = .1,
  #   alpha = .25
  # ) +
  scale_fill_gradientn(
    limits = c(0, 10^6),
    colors = viridis::turbo(100),
    values = c(0:100/100),
    na.value = viridis::turbo(100)[100]
  ) +
  scale_x_continuous(
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  theme_classic() +
  coord_fixed()

ToBlur <-
  TracksTable %>% 
  filter(
    NORMALIZED_INTENSITY >= 0,
    NORMALIZED_INTENSITY <= 1,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 0,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 <= 1,
  ) %>% 
  mutate(
    NORMALIZED_INTENSITY = round(NORMALIZED_INTENSITY*100),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = round(COMPLEMENTARY_NORMALIZED_INTENSITY_1*100)
  ) %>% 
  group_by(
    NORMALIZED_INTENSITY,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1
  ) %>% 
  summarize(.
    TIME_ADJUSTED = median(TIME_ADJUSTED)
  ) %>% 
  ungroup() %>% 
  complete(
    NORMALIZED_INTENSITY = full_seq(c(min(NORMALIZED_INTENSITY), max(NORMALIZED_INTENSITY)), period = 1),
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = full_seq(c(min(COMPLEMENTARY_NORMALIZED_INTENSITY_1), max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)), period = 1),
  ) %>%
  mutate(
    NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/100,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/100
  ) %>% 
  distinct() %>% 
  as.data.table()

# Matrix <- matrix(ToBlur$TIME_ADJUSTED, nrow = NROW(unique(ToBlur$NORMALIZED_INTENSITY)))
# ijtiff::display(Test)


  ggplot(
    ToBlur,
    aes(
      NORMALIZED_INTENSITY,
      COMPLEMENTARY_NORMALIZED_INTENSITY_1,
      fill = TIME_ADJUSTED
    )
  ) + 
  geom_tile() +
  # stat_density2d(geom = "tile", contour = F) +
  scale_fill_viridis(
    option = "turbo"
  ) +
  scale_x_continuous(
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  theme_classic()



