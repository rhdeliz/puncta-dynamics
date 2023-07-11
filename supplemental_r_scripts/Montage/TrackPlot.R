# Plot
ggplot(
  TrackTable
) +
  geom_vline(
    xintercept = unique(TrackTable$TIME)[SELECT_FRAMES_ADJSUTED],
    color = "white"
  ) +
  # geom_path(
  #   aes(
  #     x = TIME,
  #     y = NORMALIZED_CHANGEPOINT_INTENSITY,
  #     group = PROTEIN
  #   ),
  #   color = "darkgrey"
  # ) +
  geom_path(
    aes(
      x = TIME,
      y = NORMALIZED_INTENSITY,
      color = PROTEIN,
      group = PROTEIN
    )
  ) +
  # geom_path(
  #   aes(
  #     x = TIME,
  #     y = NORMALIZED_INTENSITY,
  #     group = PROTEIN
  #   ),
  #   color = "darkgrey"
  # ) +
  # geom_path(
  #   aes(
  #     x = TIME,
  #     y = NORMALIZED_CHANGEPOINT_INTENSITY,
  #     color = PROTEIN,
  #     group = PROTEIN
  #   )
  # ) +
  scale_color_manual(
    values = c("green", "magenta"),
    breaks = c(REFERENCE_PROTEIN, QUERY_PROTEIN)
  ) +
  facet_wrap(
    ~PROTEIN,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    x = "Time (s)",
    y = "Size",
    color = "Protein"
  ) +
  dark_theme_classic(
    base_size = 18
  ) +
  theme(
    legend.position = "bottom"
  )

setwd(OUTPUT_PATH)

# Save vector image
ggsave(
  "Track.pdf",
  height = 3*2,
  width = 4*2
)
