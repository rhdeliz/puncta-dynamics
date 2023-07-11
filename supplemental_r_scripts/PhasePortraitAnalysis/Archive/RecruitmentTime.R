setwd(OUTPUT_DIRECTORY)

RecruitmentTime <-
  AverageIntByTime %>% 
  filter(
    QUERY_TOTAL_INTENSITY >= .75
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN
  ) %>% 
  summarize(
    RECRUITMENT_TIME = min(FRAMES_ADJUSTED)
  )

ggplot(
  RecruitmentTime,
  aes(
    x = LIGAND_DENSITY_CAT,
    y = RECRUITMENT_TIME,
    group = QUERY_PROTEIN,
    color = QUERY_PROTEIN
  )
) +
  geom_line() +
  scale_x_continuous(
    # Log transform axis. The secondary axis reinforces that it's in log as it'll say 2^x
    trans = "log2",
    sec.axis = sec_axis(
      trans = ~.,
      breaks = trans_breaks("log2", function(x) 2^x),
      labels = trans_format("log2", math_format(2^.x))
    )
  ) +
  # scale_y_continuous(
  #   trans = "log2",
  #   sec.axis = sec_axis(
  #     trans = ~.,
  #     breaks = trans_breaks("log2", function(x) 2^x),
  #     labels = trans_format("log2", math_format(2^.x))
  #   )
  # ) +
  scale_color_brewer(
    palette = "Set1"
  ) +
  dark_theme_classic() +
  labs(
    y = "Recruitment Time (s)",
    x = "Ligand Density (mol. um^-2)",
    color = "Protein"
  ) +
  theme(
    legend.position = "right"
  )
  ggsave(
    # Save vector image
    "RecruitmentTime.pdf",
    height = 4.03,
    width = 5.64
  )