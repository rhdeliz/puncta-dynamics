# Recruitment time
TimeTable <-
  NormStatTable %>% 
  drop_na(REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  arrange(
    UNIVERSAL_TRACK_ID,
    TIME_ADJUSTED 
  ) %>% 
  mutate(
    START_REF = median(REFERENCE_TOTAL_INTENSITY[1:3]),
    START_QRY = median(QUERY_TOTAL_INTENSITY[1:3]),
    END_REF = median(REFERENCE_TOTAL_INTENSITY[n()-3:n()]),
    END_QRY = median(QUERY_TOTAL_INTENSITY[n()-3:n()])
  ) %>% 
  filter(
    # START_REF <= 3,
    # START_QRY <= 3
    # END_REF >= 3,
    # END_QRY >= 3
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%  
  # mutate(
  #   END_REF > 6,
  #   END_QRY > 3
  # ) %>%
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  drop_na(REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>%
  summarize(
    N = n(),
    MAD_REF = mad(REFERENCE_TOTAL_INTENSITY, na.rm = T),
    MAD_QRY = mad(QUERY_TOTAL_INTENSITY, na.rm = T),
    
    MAD_DELTA_REF = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    MAD_DELTA_QRY = median(DELTA_QUERY_TOTAL_INTENSITY),
    TIME_ADJUSTED = median(TIME_ADJUSTED, na.rm = T),
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY, na.rm = T),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY, na.rm = T),
    
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY)
  ) %>%
  filter(
    N >= 175
  ) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>%
  mutate(
    SEM_DELTA_REF = MAD_DELTA_REF/sqrt(N),
    SEM_DELTA_QRY = MAD_DELTA_QRY/sqrt(N),
    
    SEM_DELTA_REF = signal::sgolayfilt(MAD_DELTA_REF, p = 1, n = 31)/sqrt(N),
    SEM_DELTA_QRY = signal::sgolayfilt(MAD_DELTA_QRY, p = 1, n = 31)/sqrt(N),

    SEM_REF = signal::sgolayfilt(MAD_REF, p = 1, n = 31)/sqrt(N),
    SEM_QRY = signal::sgolayfilt(MAD_QRY, p = 1, n = 31)/sqrt(N),

    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 31),
    QUERY_TOTAL_INTENSITY =  signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 31),
    
    # DELTA_REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(DELTA_REFERENCE_TOTAL_INTENSITY, p = 1, n = 31),
    # DELTA_QUERY_TOTAL_INTENSITY =  signal::sgolayfilt(DELTA_QUERY_TOTAL_INTENSITY, p = 1, n = 31)
  ) %>% 
  drop_na(REFERENCE_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY) %>%
  group_by(
    COHORT,
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY-min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY-min(QUERY_TOTAL_INTENSITY),
    
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY - min(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY - min(abs(DELTA_QUERY_TOTAL_INTENSITY)),
    
    TIME_ADJUSTED = TIME_ADJUSTED - min(TIME_ADJUSTED)
  ) %>% 
  mutate(
    SEM_REF = SEM_REF/max(REFERENCE_TOTAL_INTENSITY),
    SEM_QRY = SEM_QRY/max(QUERY_TOTAL_INTENSITY),

    SEM_DELTA_REF = SEM_DELTA_REF/max(abs(SEM_DELTA_REF)),
    SEM_DELTA_QRY = SEM_DELTA_QRY/max(abs(SEM_DELTA_QRY)),

    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY),
    
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(abs(DELTA_QUERY_TOTAL_INTENSITY)),
  ) %>%
  as.data.table()

StackedOrder <-
  TimeTable %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT,
    TIME_ADJUSTED
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, QUERY_PROTEIN, REFERENCE_PROTEIN
  ) %>% 
  summarize(
    RECRUITMENT_TIME = TIME_ADJUSTED[which.min(QUERY_TOTAL_INTENSITY<=.49)],
    RECRUITMENT_TIME_SEM_LO = TIME_ADJUSTED[which.min((QUERY_TOTAL_INTENSITY - SEM_QRY)<=.49)],
    RECRUITMENT_TIME_SEM_HI = TIME_ADJUSTED[which.min((QUERY_TOTAL_INTENSITY + SEM_QRY)<=.49)],
    
    REF_RECRUITMENT_TIME = TIME_ADJUSTED[which.min(REFERENCE_TOTAL_INTENSITY<=.49)],
    REF_RECRUITMENT_TIME_SEM_LO = TIME_ADJUSTED[which.min((REFERENCE_TOTAL_INTENSITY - SEM_REF)<=.49)],
    REF_RECRUITMENT_TIME_SEM_HI = TIME_ADJUSTED[which.min((REFERENCE_TOTAL_INTENSITY + SEM_REF)<=.49)],
    
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY[which.min(REFERENCE_TOTAL_INTENSITY<=.49)],
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY[which.min(QUERY_TOTAL_INTENSITY<=.49)]
  ) %>% 
  arrange(
    RECRUITMENT_TIME
  ) %>% 
  as.data.table()

TimeTable$COHORT <- factor(TimeTable$COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER))
TimeTable$QUERY_PROTEIN <- factor(TimeTable$QUERY_PROTEIN, levels = PROTEIN_ORDER)
StackedOrder$COHORT <- factor(StackedOrder$COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER))
StackedOrder$QUERY_PROTEIN <- factor(StackedOrder$QUERY_PROTEIN, levels = PROTEIN_ORDER)


TimeTable <- 
  TimeTable %>% 
  mutate(
    YMIN_REF = REFERENCE_TOTAL_INTENSITY - SEM_REF,
    YMAX_REF = REFERENCE_TOTAL_INTENSITY + SEM_REF,
    YMIN_QRY = QUERY_TOTAL_INTENSITY - SEM_QRY,
    YMAX_QRY =  QUERY_TOTAL_INTENSITY + SEM_QRY
  ) %>% 
  mutate(
    YMIN_REF = ifelse(YMIN_REF <0, 0, YMIN_REF),
    YMAX_REF = ifelse(YMAX_REF > 1, 1, YMAX_REF),
    YMIN_QRY = ifelse(YMIN_QRY <0, 0, YMIN_QRY),
    YMAX_QRY = ifelse(YMAX_QRY > 1, 1, YMAX_QRY)
  ) %>% 
  as.data.table()

ggplot(TimeTable) +
  geom_ribbon(
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN,
      fill = REFERENCE_PROTEIN,
      group = COHORT,
      ymin = YMIN_REF,
      ymax = YMAX_REF
    ),
    alpha = .25
  ) +
  geom_ribbon(
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      fill = QUERY_PROTEIN,
      group = COHORT,
      ymin = YMIN_QRY,
      ymax =  YMAX_QRY
    ),
    alpha = .25
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN,
      group = COHORT
    ),
    size = 2
  ) +
  geom_point(
    data = StackedOrder,
    aes(
      x = REF_RECRUITMENT_TIME,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    ),
    size = 6
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      group = COHORT
    ),
    size = 2
  ) +
  geom_point(
    data = StackedOrder,
    aes(
      x = RECRUITMENT_TIME,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    ),
    size = 6
  ) +
  scale_x_continuous(
    breaks = c(0, 400, 800),
    limits = c(0, 800)
  ) +
  scale_y_continuous(
    breaks = c(0, .5, 1)
    # limits = c(0, max(c(TimeTable$SEM_REF, TimeTable$SEM_QRY)+1)
  ) +
  scale_fill_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  labs(
    x = "Cluster Time (s)",
    # y = "Normalized Intensity",
    y = "Scaled Avg. Intensity",
    color = "Protein"
  ) +
  facet_wrap(
    ~QUERY_PROTEIN,
    scales = "free_y",
    nrow = 1
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  )

ggsave(
  # Save vector image
  file.path(OUTPUT_DIRECTORY, "recruitment_time_facet.pdf"),
  height = 2,
  width = 2/.17
)

ggplot(TimeTable) +
  geom_ribbon(
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      # color = QUERY_PROTEIN,
      fill = QUERY_PROTEIN,
      group = COHORT,
      ymin = YMIN_QRY,
      ymax =  YMAX_QRY
    ),
    color = NA,
    alpha = .25
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      group = COHORT
    ),
    size = 2
  ) +
  geom_point(
    data = StackedOrder,
    aes(
      x = RECRUITMENT_TIME,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      group = COHORT
    ),
    size = 6
  ) +
  labs(
    x = "Cluster Time (s)",
    y = "Scaled Avg. Intensity",
    color = "Protein"
  ) +
  scale_fill_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  )

ggsave(
  # Save vector image
  file.path(OUTPUT_DIRECTORY, "recruitment_time.pdf"),
  height = 2,
  width = 9.5/2.25
)


# ggplot(
#   TimeTable
# ) +
#   geom_path(
#     aes(
#       x = TIME_ADJUSTED,
#       y = REFERENCE_TOTAL_INTENSITY,
#       color = QUERY_PROTEIN,
#       group = COHORT
#     ),
#     size = 2
#   ) +
#   geom_point(
#     data = StackedOrder,
#     aes(
#       x = REF_RECRUITMENT_TIME,
#       y = REFERENCE_TOTAL_INTENSITY,
#       color = QUERY_PROTEIN,
#       group = COHORT
#     ),
#     size = 6
#   ) +
#   labs(
#     x = "Cluster Time (s)",
#     y = "Avg. Scaled Int.",
#     color = "Protein"
#   ) +
#   scale_color_manual(
#     values = ColorTable,
#     # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
#     limits = force,
#     # drop = F,
#     guide = "none"
#   ) +
#   theme_classic() +
#   theme(
#     strip.background = element_blank()
#   )
# 
# ggsave(
#   # Save vector image
#   file.path(OUTPUT_DIRECTORY, "recruitment_time_myd88.pdf"),
#   height = 2,
#   width = 9.5/2.25
# )


StackedOrder %>% 
  arrange(
    RECRUITMENT_TIME
  ) %>% 
  mutate(
    # QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = unique(QUERY_PROTEIN))
  ) %>% 
  as.data.table() %>% 
ggplot() +
  geom_col(
    aes(
      x = QUERY_PROTEIN,
      y = RECRUITMENT_TIME,
      fill = QUERY_PROTEIN,
      group = COHORT
    ),
    color = "black"
  ) +
  geom_errorbar(
    aes(
      x = QUERY_PROTEIN,
      ymin = RECRUITMENT_TIME_SEM_LO,
      ymax = RECRUITMENT_TIME_SEM_HI,
    ),
    color = "black",
    width = 0.5
  ) +
  labs(
    # x = "Cell Line",
    y = "Cluster Time at Mid-Intensity (s)"
  ) +
  scale_fill_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  scale_color_manual(
    values = ColorTable,
    limits = force,
    guide = "none"
  ) +
  scale_x_discrete(limits=rev) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_flip()

ggsave(
  # Save vector image
  file.path(OUTPUT_DIRECTORY, "recruitment_time_bar.pdf"),
  height = 2,
  width = 9.5/2.25
)






DeltaTable <-
  TimeTable %>% 
  drop_na() %>% 
  group_by(
    COHORT
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = lead(REFERENCE_TOTAL_INTENSITY) - lag(REFERENCE_TOTAL_INTENSITY),
    DELTA_QUERY_TOTAL_INTENSITY = lead(QUERY_TOTAL_INTENSITY) - lag(QUERY_TOTAL_INTENSITY)
  ) %>% 
  drop_na() %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY =  signal::sgolayfilt(DELTA_REFERENCE_TOTAL_INTENSITY, p = 3, n = 51),
    DELTA_QUERY_TOTAL_INTENSITY =  signal::sgolayfilt(DELTA_QUERY_TOTAL_INTENSITY, p = 3, n = 51)
  ) %>%
  mutate(
    # DELTA_QUERY_TOTAL_INTENSITY =  signal::sgolayfilt(DELTA_QUERY_TOTAL_INTENSITY, p = 1, n = 51)
  ) %>%
  mutate(
    # DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY - min(DELTA_QUERY_TOTAL_INTENSITY)
  ) %>% 
  mutate(
    DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY/max(abs(DELTA_REFERENCE_TOTAL_INTENSITY)),
    DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY/max(abs(DELTA_QUERY_TOTAL_INTENSITY))
  ) %>% 
  mutate(
    RECRUITMENT_TIME = case_when(DELTA_QUERY_TOTAL_INTENSITY == max(DELTA_QUERY_TOTAL_INTENSITY) ~ TIME_ADJUSTED),
    REF_RECRUITMENT_TIME = case_when(DELTA_REFERENCE_TOTAL_INTENSITY == max(DELTA_REFERENCE_TOTAL_INTENSITY) ~ TIME_ADJUSTED),
    MAX_DELTA_QUERY = max(DELTA_QUERY_TOTAL_INTENSITY),
    MAX_DELTA_REFERENCE = max(DELTA_REFERENCE_TOTAL_INTENSITY)
  ) %>% 
  as.data.table()

SummaryDeltaTable <-
  DeltaTable %>% 
  select(
    COHORT,
    RECRUITMENT_TIME
  ) %>% 
  distinct() %>% 
  drop_na() %>% 
  arrange(
    RECRUITMENT_TIME
  ) %>% 
  as.data.table()

MIN_X = min(c(DeltaTable$DELTA_REFERENCE_TOTAL_INTENSITY, DeltaTable$DELTA_QUERY_TOTAL_INTENSITY))*1.1
MIN_X = floor(MIN_X*4)/4

ggplot(
  DeltaTable
) +
  geom_hline(
    yintercept = 0
  ) +
  # geom_ribbon(
  #   aes(
  #     x = TIME_ADJUSTED,
  #     y = DELTA_REFERENCE_TOTAL_INTENSITY,
  #     color = REFERENCE_PROTEIN,
  #     fill = REFERENCE_PROTEIN,
  #     group = COHORT,
  #     ymin = DELTA_REFERENCE_TOTAL_INTENSITY - MAD_DELTA_REF,
  #     ymax = DELTA_REFERENCE_TOTAL_INTENSITY + MAD_DELTA_REF
  #   ),
  #   alpha = .25
# ) +
# geom_ribbon(
#   aes(
#     x = TIME_ADJUSTED,
#     y = DELTA_QUERY_TOTAL_INTENSITY,
#     color = QUERY_PROTEIN,
#     fill = QUERY_PROTEIN,
#     group = COHORT,
#     ymin = DELTA_QUERY_TOTAL_INTENSITY - MAD_DELTA_QRY,
#     ymax =  DELTA_QUERY_TOTAL_INTENSITY + MAD_DELTA_QRY
#   ),
#   alpha = .25
# ) +
geom_path(
  aes(
    x = TIME_ADJUSTED,
    y = DELTA_REFERENCE_TOTAL_INTENSITY,
    color = REFERENCE_PROTEIN,
    group = COHORT
  ),
  size = 2
) +
  geom_point(
    # data = StackedOrder,
    aes(
      x = REF_RECRUITMENT_TIME,
      y = MAX_DELTA_REFERENCE,
      color = REFERENCE_PROTEIN
    ),
    size = 6
  ) +
  geom_path(
    aes(
      x = TIME_ADJUSTED,
      y = DELTA_QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN,
      group = COHORT
    ),
    size = 2
  ) +
  geom_point(
    aes(
      x = RECRUITMENT_TIME,
      y = MAX_DELTA_QUERY,
      color = QUERY_PROTEIN
    ),
    size = 6
  ) +
  scale_x_continuous(
    breaks = c(0, 400, 800),
    limits = c(0, 800)
  ) +
  scale_y_continuous(
    breaks = c(MIN_X, 0, 0.5, 1),
    limits = c(MIN_X, 1.1)
    # limits = c(0, max(c(TimeTable$SEM_REF, TimeTable$SEM_QRY)+1)
  ) +
  scale_fill_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  scale_color_manual(
    values = ColorTable,
    # labels = PROTEIN_ORDER[2:NROW(PROTEIN_ORDER)],
    limits = force,
    # drop = F,
    guide = "none"
  ) +
  labs(
    x = "Cluster Time (s)",
    # y = "Normalized Intensity",
    y = "Scaled Avg. Int. Change",
    color = "Protein"
  ) +
  facet_wrap(
    ~QUERY_PROTEIN,
    scales = "free_y",
    nrow = 1
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  )


ggsave(
  # Save vector image
  file.path(OUTPUT_DIRECTORY, "recruitment_time_delta.pdf"),
  height = 2,
  width = 2/.17
)


SummaryDeltaTable