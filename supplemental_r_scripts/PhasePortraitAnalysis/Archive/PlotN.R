setwd(OUTPUT_DIRECTORY)

# N plot
NTable <-
  ppMTAll %>%
  group_by(
    LIGAND_DENSITY_CAT,
    CATEGORY,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>%
  mutate(
    N = N/max(N)*100
  ) %>%
  filter(
    CATEGORY == "All"
  )

# Plot
ggplot(
) +
  geom_point(
    data = NTable,
    aes(
      x = REFERENCE_TOTAL_INTENSITY,
      y = QUERY_TOTAL_INTENSITY,
      color = N,
      size = N
    )
  ) +
  # Color of magnitude. Goes from gray to white
  scale_color_distiller(
    palette = "RdPu",
    # trans = "log2"
  ) +
  scale_size(
    range = c(0.01, 1),
    trans = "log2"
  ) +
  # Use angle of arrows to color the triangle. Full spectrum
  scale_x_continuous(
    # Log transform axis. The secondary axis reinforces that it's in log as it'll say 2^x
    trans = "log2",
    sec.axis = sec_axis(
      trans = ~.,
      breaks = trans_breaks("log2", function(x) 2^x),
      labels = trans_format("log2", math_format(2^.x))
    )
  ) +
  scale_y_continuous(
    trans = "log2",
    sec.axis = sec_axis(
      trans = ~.,
      breaks = trans_breaks("log2", function(x) 2^x),
      labels = trans_format("log2", math_format(2^.x))
    )
  ) +
  labs(
    x = "MyD88 Size",
    y = "Query Size (Box has Name)",
    color = "N",
    size = "N"
  ) +
  facet_grid(
    # QUERY_PROTEIN~CATEGORY#~LIGAND_DENSITY_CAT
    CATEGORY~QUERY_PROTEIN+LIGAND_DENSITY_CAT
  ) +
  # Make dark to make the colors pop in a monitor
  # Use theme classic for printing
  dark_theme_classic(
    base_size = 18
  ) +
  theme(
    legend.position = "right"
  )
  ggsave(
    # Save vector image
    "N.pdf",
    height = 4.76,
    width = 11.5*2
    # height = 4.03,
    # width = 5.64
  )
