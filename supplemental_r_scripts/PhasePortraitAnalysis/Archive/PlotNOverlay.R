setwd(OUTPUT_DIRECTORY)

NTable <-
  statTable %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN,
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY
  ) %>% 
  summarize(
    N = n()
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    N = N/max(N)
  )

# Plot
ggplot(
) +
  geom_raster(
    data = NTable,
    aes(
      x = REFERENCE_TOTAL_INTENSITY,
      y = QUERY_TOTAL_INTENSITY,
      fill = N
    ),
    interpolate = T
  ) +
  scale_fill_gradient(
    low = "#303030",
    high = "white"
  ) +
  ggnewscale::new_scale_fill()+
  geom_vline(
    xintercept = c(1, 6)
  ) +
  geom_hline(
    yintercept = c(1, 4)
  ) +
  # # Make triangles of the arrows
  # geom_regon(
  #   # data = ppMT,
  #   data = ppMT %>% filter(CATEGORY=="Increase"),
  #   aes(
  #     x0 = REFERENCE_TOTAL_INTENSITY,
  #     y0 = QUERY_TOTAL_INTENSITY,
  #     # Size of triangle is the magnitude
  #     r = MAGNITUDE/2,
  #     # Angle of triangle is the derivative
  #     angle =  -(RAD_ANGLE+pi)+2*pi,
#     # 3 sides for a triangle
#     sides = 3,
#     # Angle of the derivatives is the color
#     fill = hsv(DEG_ANGLE),
#     # Redundancy so that the outline of the triangle also tells you the magnitude
#     # Also helps pair the arrow with the triangle
#     color = MAGNITUDE
#   ),
#   # Make transparent to see what's happening behind
#   alpha = 0.75
# ) +
# # Puts arrows inside triangles so that we know which tip of the triangle to look at
# geom_quiver(
#   # data = ppMT,
#   data = ppMT %>% filter(CATEGORY=="Increase"),
#   aes(
#     x = REFERENCE_TOTAL_INTENSITY,
#     y = QUERY_TOTAL_INTENSITY,
#     u = DELTA_REFERENCE_TOTAL_INTENSITY,
#     v = DELTA_QUERY_TOTAL_INTENSITY,
#     # Arrow color is the magnitude
#     # Same as the triangle so that they can be visually paired
#     color = MAGNITUDE
#   ),
#   # Arrow size blowup
#   vecsize = 2,
#   # Center the arrow
#   center = T
# ) +
# Color of magnitude. Goes from gray to white
scale_color_gradient(
  low = "#303030",
  high = "white"
) +
  # Use angle of arrows to color the triangle. Full spectrum
  scale_fill_identity(guide = "none") +
  scale_x_continuous(
    # Log transform axis. The secondary axis reinforces that it's in log as it'll say 2^x
    trans = "log2",
    sec.axis = sec_axis(
      trans = ~.,
      breaks = trans_breaks("log2", function(x) 2^x),
      labels = trans_format("log2", math_format(2^.x))
    ),
    # limits = c(2^-3, 2^6)
  ) +
  scale_y_continuous(
    trans = "log2",
    sec.axis = sec_axis(
      trans = ~.,
      breaks = trans_breaks("log2", function(x) 2^x),
      labels = trans_format("log2", math_format(2^.x))
    ),
    # limits = c(2^-4, 2^3.5)
  ) +
  labs(
    x = "MyD88 Size",
    y = "Query Size (Box has Name)"
  ) +
  facet_grid(
    # QUERY_PROTEIN~CATEGORY#~LIGAND_DENSITY_CAT
    QUERY_PROTEIN~LIGAND_DENSITY_CAT
  ) +
  # Make dark to make the colors pop in a monitor
  # Use theme classic for printing
  dark_theme_classic(
    base_size = 24
  ) +
  theme(
    legend.position = "none"
  ) 

ggsave(
  # Save vector image
  "N Overlay.pdf",
  height = 4.76*2,
  width = 11.5*2
  # height = 4.03,
  # width = 5.64
)
