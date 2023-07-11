SampleRecruitmentTime <-
  # StatTable %>%
  NormStatTable %>%
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25
  ) %>%
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY),
    TIME_ADJUSTED = mean(TIME_ADJUSTED)
  ) %>% 
  group_by(
    COHORT
  ) %>%
  mutate(
    # REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    # QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>% 
  mutate(
    # REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    # QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  mutate(
    # REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    # QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  group_by(
    COHORT
  ) %>%
  mutate(
    # REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    # QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY),
    
    COHORT = factor(COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER)),
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = PROTEIN_ORDER),
    REFERENCE_PROTEIN = factor(REFERENCE_PROTEIN, levels = PROTEIN_ORDER)
  ) %>% 
  as.data.table()


ggplot()+
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  # geom_point(
  #   data = RecruitmentTimePoints,
  #   aes(
  #     x = REFERENCE_TIME_ADJUSTED,
  #     y = 0.5,
  #     color = REFERENCE_PROTEIN
  #   ),
  #   size = 3
  # ) +
  # geom_point(
  #   data = RecruitmentTimePoints,
  #   aes(
  #     x = QUERY_TIME_ADJUSTED,
  #     y = 0.5,
  #     color = QUERY_PROTEIN
  #   ),
  #   size = 3
  # ) +
  # scale_x_continuous(limits = c(0, max(SampleRecruitmentTime$TIME_ADJUSTED))) +
  # scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = ColorTable
  ) +
  labs(
    x = "Cluster Time (s)",
    y = "Int. (a.u.)",
    # y = "Avg. Scaled Int. ",
    color = "Protein"
  ) +
  theme_void(
    base_size = 6
  ) +
  theme(
    strip.background = element_blank()
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(
  "avg_time.pdf",
  height = .5,
  width = 1.5
)





SampleRecruitmentTime <-
  # StatTable %>% 
  NormStatTable %>%
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25
  ) %>%
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY),
    TIME_ADJUSTED = mean(TIME_ADJUSTED)
  ) %>% 
  group_by(
    COHORT
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  mutate(
    # REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    # QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  group_by(
    COHORT
  ) %>%
  mutate(
    # REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    # QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY),
    
    COHORT = factor(COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER)),
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = PROTEIN_ORDER),
    REFERENCE_PROTEIN = factor(REFERENCE_PROTEIN, levels = PROTEIN_ORDER)
  ) %>% 
  as.data.table()


ggplot()+
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  # geom_point(
  #   data = RecruitmentTimePoints,
  #   aes(
  #     x = REFERENCE_TIME_ADJUSTED,
  #     y = 0.5,
  #     color = REFERENCE_PROTEIN
  #   ),
  #   size = 3
  # ) +
  # geom_point(
  #   data = RecruitmentTimePoints,
#   aes(
#     x = QUERY_TIME_ADJUSTED,
#     y = 0.5,
#     color = QUERY_PROTEIN
#   ),
#   size = 3
# ) +
# scale_x_continuous(limits = c(0, max(SampleRecruitmentTime$TIME_ADJUSTED))) +
# scale_y_continuous(limits = c(0, 1)) +
scale_color_manual(
  values = cols[c(1, 4)],
  guide = "none"
) +
  labs(
    x = "Cluster Time (s)",
    y = "Int. (a.u.)",
    # y = "Avg. Scaled Int. ",
    color = "Protein"
  ) +
  theme_void(
    base_size = 6
  ) +
  theme(
    strip.background = element_blank()
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(
  "smooth.pdf",
  height = .5,
  width = 1.5
)


SampleRecruitmentTime <-
  # StatTable %>% 
  NormStatTable %>%
  filter(
    COHORT == "MyD88 TRAF6",
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
    LIGAND_DENSITY_CAT == 32,
    FPS == 0.25
  ) %>%
  group_by(
    COHORT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = median(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = median(QUERY_TOTAL_INTENSITY),
    TIME_ADJUSTED = mean(TIME_ADJUSTED)
  ) %>% 
  group_by(
    COHORT
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 7),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 7),
  ) %>% 
  mutate(
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 15),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 15)
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY - min(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY - min(QUERY_TOTAL_INTENSITY)
  ) %>%
  group_by(
    COHORT
  ) %>%
  mutate(
    REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY/max(REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY/max(QUERY_TOTAL_INTENSITY),
    
    COHORT = factor(COHORT, levels = paste(USE_REFERENCE_PROTEIN, PROTEIN_ORDER)),
    QUERY_PROTEIN = factor(QUERY_PROTEIN, levels = PROTEIN_ORDER),
    REFERENCE_PROTEIN = factor(REFERENCE_PROTEIN, levels = PROTEIN_ORDER)
  ) %>% 
  as.data.table()


ggplot()+
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
  # geom_point(
  #   data = RecruitmentTimePoints,
  #   aes(
  #     x = REFERENCE_TIME_ADJUSTED,
  #     y = 0.5,
  #     color = REFERENCE_PROTEIN
  #   ),
  #   size = 3
  # ) +
  # geom_point(
  #   data = RecruitmentTimePoints,
#   aes(
#     x = QUERY_TIME_ADJUSTED,
#     y = 0.5,
#     color = QUERY_PROTEIN
#   ),
#   size = 3
# ) +
# scale_x_continuous(limits = c(0, max(SampleRecruitmentTime$TIME_ADJUSTED))) +
# scale_y_continuous(limits = c(0, 1)) +
scale_color_manual(
  values = cols[c(1, 4)],
  guide = "none"
) +
  labs(
    x = "Cluster Time (s)",
    y = "Int. (a.u.)",
    # y = "Avg. Scaled Int. ",
    color = "Protein"
  ) +
  theme_void(
    base_size = 6
  ) +
  theme(
    strip.background = element_blank()
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(
  "scale.pdf",
  height = .5,
  width = 1.5
)




ggplot()+
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = REFERENCE_TOTAL_INTENSITY,
      color = REFERENCE_PROTEIN
    )
  ) +
  geom_path(
    data = SampleRecruitmentTime,
    aes(
      x = TIME_ADJUSTED,
      y = QUERY_TOTAL_INTENSITY,
      color = QUERY_PROTEIN
    )
  ) +
# geom_segment(
#   data = RecruitmentTimePoints %>% filter(COHORT=="MyD88 TRAF6") %>% as.data.table(),
#   aes(
#     x = REFERENCE_TIME_ADJUSTED,
#     xend = QUERY_TIME_ADJUSTED,
#     y = 0.5,
#     yend = 0.5
#   )
# ) +
geom_point(
  data = RecruitmentTimePoints %>% filter(COHORT=="MyD88 TRAF6") %>% as.data.table(),
  aes(
    x = REFERENCE_TIME_ADJUSTED,
    y = 0.5,
    color = REFERENCE_PROTEIN
  ),
  size = 3
) +
geom_point(
  data = RecruitmentTimePoints %>% filter(COHORT=="MyD88 TRAF6") %>% as.data.table(),
  aes(
    x = QUERY_TIME_ADJUSTED,
    y = 0.5,
    color = QUERY_PROTEIN
  ),
  size = 3
) +
scale_x_continuous(limits = c(0, max(SampleRecruitmentTime$TIME_ADJUSTED))) +
scale_y_continuous(limits = c(0, 1)) +
scale_color_manual(
  values = cols[c(1, 4)],
  guide = "none"
) +
  labs(
    x = "Cluster Time (s)",
    y = "Scaled Int. (a.u.)",
    # y = "Avg. Scaled Int. ",
    color = "Protein"
  ) +
  theme_void(
    base_size = 6
  ) +
  theme(
    strip.background = element_blank()
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(
  "assemble50.pdf",
  height = .5,
  width = 1.5
)



ggplot()+
  # geom_path(
  #   data = SampleRecruitmentTime,
  #   aes(
  #     x = TIME_ADJUSTED,
  #     y = REFERENCE_TOTAL_INTENSITY,
  #     color = REFERENCE_PROTEIN
  #   )
  # ) +
  # geom_path(
  #   data = SampleRecruitmentTime,
  #   aes(
#     x = TIME_ADJUSTED,
#     y = QUERY_TOTAL_INTENSITY,
#     color = QUERY_PROTEIN
#   )
# ) +
geom_segment(
  data = RecruitmentTimePoints %>% filter(COHORT=="MyD88 TRAF6") %>% as.data.table(),
  aes(
    x = REFERENCE_TIME_ADJUSTED,
    xend = QUERY_TIME_ADJUSTED,
    y = 0.5,
    yend = 0.5
  )
) +
  geom_point(
    data = RecruitmentTimePoints %>% filter(COHORT=="MyD88 TRAF6") %>% as.data.table(),
    aes(
      x = REFERENCE_TIME_ADJUSTED,
      y = 0.5,
      color = REFERENCE_PROTEIN
    ),
    size = 3
  ) +
  geom_point(
    data = RecruitmentTimePoints %>% filter(COHORT=="MyD88 TRAF6") %>% as.data.table(),
    aes(
      x = QUERY_TIME_ADJUSTED,
      y = 0.5,
      color = QUERY_PROTEIN
    ),
    size = 3
  ) +
  scale_x_continuous(limits = c(0, max(SampleRecruitmentTime$TIME_ADJUSTED))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = cols[c(1, 4)],
    guide = "none"
  ) +
  labs(
    x = "Cluster Time (s)",
    y = "Scaled Int. (a.u.)",
    # y = "Avg. Scaled Int. ",
    color = "Protein"
  ) +
  theme_void(
    base_size = 6
  ) +
  theme(
    strip.background = element_blank()
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(
  "delta50.pdf",
  height = .5,
  width = 1.5
)



