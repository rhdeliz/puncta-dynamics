# Plot
ppMT %>%
  filter(
    CATEGORY=="All",
    REFERENCE_TOTAL_INTENSITY <= 8,
    QUERY_TOTAL_INTENSITY <= 5.33
  ) %>% 
  ggplot(
  ) +
  # Puts arrows inside triangles so that we know which tip of the triangle to look at
  geom_quiver(
    aes(
      x = REFERENCE_TOTAL_INTENSITY,
      y = QUERY_TOTAL_INTENSITY,
      u = DELTA_REFERENCE_TOTAL_INTENSITY,
      v = DELTA_QUERY_TOTAL_INTENSITY,
      # Arrow color is the magnitude
      # Same as the triangle so that they can be visually paired
      color = MAGNITUDE,
      vecsize = 1/MAGNITUDE
    ),
    # Arrow size blowup
    # Center the arrow
    center = T
  ) +
  scale_color_distiller(
    palette = "RdPu",
    trans = "log2",
    breaks = trans_breaks("log2", function(x) 2^x),
    labels = trans_format("log2", math_format(2^.x))
  ) +
  labs(
    x = "MyD88 Size",
    y = "Query Size (Box has Name)",
    color = "Growth (proteins/s)"
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
    legend.position = "bottom",
    legend.key.width = unit(0.5,"in")
  ) +
  coord_fixed()