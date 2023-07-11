PhasePortrait <-
  Table %>% 
  filter(
    FRAME %in% IMPORT_FRAMES
  )

PhasePortraitPlot <- NULL
PhasePortraitPlot$ref_int = PhasePortrait$NORMALIZED_INTENSITY[2]
PhasePortraitPlot$qry_int = PhasePortrait$COMPLEMENTARY_NORMALIZED_INTENSITY_1[2]

PhasePortraitPlot$d_ref_int <- (PhasePortrait$NORMALIZED_INTENSITY[3]-PhasePortrait$NORMALIZED_INTENSITY[1])/PhasePortrait$NORMALIZED_INTENSITY[2]
PhasePortraitPlot$d_qry_int <- (PhasePortrait$COMPLEMENTARY_NORMALIZED_INTENSITY_1[3]-PhasePortrait$COMPLEMENTARY_NORMALIZED_INTENSITY_1[1])/PhasePortrait$COMPLEMENTARY_NORMALIZED_INTENSITY_1[2]

PhasePortraitPlot <- as_tibble(PhasePortraitPlot)

PhasePortraitPlot <-
  PhasePortraitPlot %>% 
  mutate(
    # # Get absolute magnitude
    magnitude = abs(d_ref_int) - abs(d_qry_int),
    # Get angle to account for negative magnitudes
    rad_angle = atan2(-d_ref_int, -d_qry_int),
    # Make angle
    deg_angle = rad_angle*180/pi+180,
    deg_angle = deg_angle/360,
    N = 1
  )


ggplot(
  PhasePortraitPlot
) +
  # Make triangles of the arrows
  geom_regon(
    aes(
      x0 = ref_int,
      y0 = qry_int,
      # Size of triangle is the magnitude
      r = magnitude, # log(N/max(N)+1), # (MAGNITUDE)/10
      # Angle of triangle is the derivative
      angle = -(rad_angle+pi)+2*pi,
      # 3 sides for a triangle
      sides = 3,
      # Angle of the derivatives is the color
      fill = hsv(deg_angle),
      # Redundancy so that the outline of the triangle also tells you the magnitude
      # Also helps pair the arrow with the triangle
      color = N
    ),
    # Make transparent to see what's happening behind
    alpha = 0.75,
    size = 0.1
  ) +
  # Puts arrows inside triangles so that we know which tip of the triangle to look at
  geom_quiver(
    aes(
      x = ref_int,
      y = qry_int,
      u = d_ref_int,
      v = d_qry_int,
      # Arrow color is the magnitude
      # Same as the triangle so that they can be visually paired
      color = N,
      # vecsize = 1/magnitude
    ),
    # Arrow size blowup
    vecsize = 2,
    # Center the arrow
    center = T,
    size = 0.25
  ) +
  # Color of magnitude. Goes from gray to white
  scale_color_gradient(
    low = "white",
    high = "white",
    trans = "log2",
    breaks = trans_breaks("log2", function(x) 2^x),
    labels = trans_format("log2", math_format(2^.x))
  ) +
  # Use angle of arrows to color the triangle. Full spectrum
  scale_fill_identity(guide = "none") +
  labs(
    x = paste(REFERENCE_PROTEIN, "Size"),
    y = paste(QUERY_PROTEIN, "Size")
  ) +
  # Make dark to make the colors pop in a monitor
  # Use theme classic for printing
  dark_theme_classic(
    base_size = 18
  ) +
  theme(
    legend.key.width= unit(1.5, 'cm'),
    # legend.key.width= unit(1.5, 'cm'),
    legend.position = "none"
  ) +
  coord_fixed()

setwd(OUTPUT_DIRECTORY)
ggsave(
  "PhasePortrait.pdf",
  height = 3*2,
  width = 6
)