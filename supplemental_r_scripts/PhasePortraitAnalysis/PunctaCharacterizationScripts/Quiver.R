


TempTable <-
  LinearPhasePortrait %>% 
  group_by(
    QUERY_PROTEIN, COHORT
  ) %>% 
  mutate(
    TEST_COHORT = grepl(QUERY_PROTEIN, COHORT)
  ) %>% 
  filter(
    TEST_COHORT == TRUE
  ) %>% 
  select(-c(
    TEST_COHORT
  )) %>% 
  filter(
    FRAMES_SINCE_LANDING_CAT == 100,
    # Plot filters
    ROUNDED_REFERENCE_TOTAL_INTENSITY >= 1,
    # ROUNDED_REFERENCE_TOTAL_INTENSITY <= 12,
    PLOT_FACETS == PlotList$Facets[FacetX],
    # N_TEST == TRUE,
    
    # Category filters
    MAX_REFERENCE_TOTAL_INTENSITY_CAT == "All"
    # MAX_REFERENCE_TOTAL_INTENSITY_CAT == "≥4.5x MyD88"
  ) %>%
  group_by(
    PLOT_FACETS, IMAGENUMBER
  ) %>% 
  # mutate(
  #   # To select image with most spots
  #   IMAGE_N = n()
  # ) %>% 
  # ungroup() %>% 
  # filter(
  #   # Select image with most bins
  #   IMAGE_N == max(IMAGE_N)
  # ) %>% 
  # filter(
  #   IMAGENUMBER == min(IMAGENUMBER)
# ) %>% 
group_by(
  COHORT, LIGAND_DENSITY_CAT, FPS, MAX_REFERENCE_TOTAL_INTENSITY_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN, IMAGENUMBER, IMAGE,
  FRAMES_SINCE_LANDING_CAT,
  FACET, PLOT_FACETS
) %>% 
  mutate(
    FPS = paste(FPS, "Hz"),
    # DELTA_CAT = gsub("ease", ".", DELTA_CAT),
    # DELTA_CAT = factor(DELTA_CAT, levels = c("All", "Decr.", "Incr.")),
    # 95th percentile to make sure large triangles don't alter results
    # SD_MAGNITUDE = SD_MAGNITUDE/quantile(SD_MAGNITUDE, .95, na.rm = T),
    
    SD_ADJUSTED_MAGNITUDE = SD_ADJUSTED_MAGNITUDE/quantile(SD_ADJUSTED_MAGNITUDE, .95, na.rm = T),
    MAD_MAGNITUDE = MAD_MAGNITUDE/quantile(MAD_MAGNITUDE, .95, na.rm = T),
    MAD_ADJUSTED_MAGNITUDE = MAD_ADJUSTED_MAGNITUDE/quantile(MAD_ADJUSTED_MAGNITUDE, .95, na.rm = T),
    MAGNITUDE = MAGNITUDE/quantile(MAGNITUDE, .95, na.rm = T),
    ADJUSTED_MAGNITUDE = ADJUSTED_MAGNITUDE/quantile(ADJUSTED_MAGNITUDE, .95, na.rm = T)
  ) %>% 
  group_by(
    PLOT_FACETS, IMAGENUMBER
  ) %>% 
  mutate(
    # Define plot axis replacing variables with generic names
    ref_int = ROUNDED_REFERENCE_TOTAL_INTENSITY,
    qry_int = ROUNDED_QUERY_TOTAL_INTENSITY,
    
    d_ref_int = (!!as.name(TempVariables[1])),
    d_qry_int = (!!as.name(TempVariables[2])),
    magnitude = (!!as.name(TempVariables[3])),
    rad_angle = (!!as.name(TempVariables[4])),
    deg_angle = (!!as.name(TempVariables[5]))
  )


# Plot quiver
ggplot(
  TempTable
) +
  # Make triangles of the arrows
  geom_regon(
    aes(
      x0 = ref_int,
      y0 = qry_int,
      # Size of triangle is the magnitude
      r = (magnitude)/4, # log(N/max(N)+1), # (MAGNITUDE)/10
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
    title = PlotTitle,
    x = "MyD88 Size",
    y = paste(TempTable$QUERY_PROTEIN[1], "Size")
  ) +
  facet_grid(
    ~IMAGENUMBER+FRAMES_SINCE_LANDING_CAT
  ) +
  # Make dark to make the colors pop in a monitor
  # Use theme classic for printing
  dark_theme_classic(
    base_size = 20
  ) +
  theme(
    legend.key.width= unit(1.25, 'cm'),
    # legend.key.width= unit(1.5, 'cm'),
    legend.position = "none"
  ) +
  coord_fixed()


# File name
SaveName <-
  paste0(
    PlotList$Variables[FacetX],
    "PhasePortrait - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", round(1/STEP_SIZE, 2), " - ",
    ifelse(MOVING_AVERAGE_WINDOW > 0, paste("MvAvg", MOVING_AVERAGE_WINDOW, "- "), ""),
    ifelse(LOG_SCALE == T, "Log - ", ""),
    "Triangle.pdf"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)
# Save
ggsave(
  # Save vector image
  file.path(COHORT_FOLDER, SaveName),
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)

ggplot(
  TempTable
) +
  geom_arrow(
    aes(
      x = ref_int,
      y = qry_int,
      mag = .5,
      angle = atan2(d_qry_int, d_ref_int)*180/pi,
      # color = N
      color = hsv(deg_angle)
    ),
    # color = "white",
    size = 1,
    arrow.length = 1
  ) +
  scale_color_identity(guide = "none") +
  scale_mag(
    max = 1,
    guide = 'none'
  ) +
  # Color of magnitude. Goes from gray to white
  # scale_color_gradient(
  #   low = "black",
  #   high = "white",
  #   trans = "log2",
  #   breaks = trans_breaks("log2", function(x) 2^x),
  #   labels = trans_format("log2", math_format(2^.x))
  # ) +
  # scale_color_distiller(
  #   palette = "RdPu",
  #   trans = "log2",
#   breaks = trans_breaks("log2", function(x) 2^x),
#   labels = trans_format("log2", math_format(2^.x))
# ) +
# Use angle of arrows to color the triangle. Full spectrum
labs(
  title = PlotTitle,
  x = "MyD88 Size",
  y = paste(TempTable$QUERY_PROTEIN[1], "Size")
) +
  # facet_grid(
  #   ~DELTA_CAT
  # ) +
  # Make dark to make the colors pop in a monitor
  # Use theme classic for printing
  dark_theme_classic(
    base_size = 20
  ) +
  theme(
    legend.key.width= unit(1.25, 'cm'),
    legend.position = "none"
  ) +
  coord_fixed()


COHORT_FOLDER = gsub("\n", " - ", PlotTitle)
COHORT_FOLDER = file.path(OUTPUT_DIRECTORY, COHORT_FOLDER)

if(!file.exists(COHORT_FOLDER)){
  dir.create(COHORT_FOLDER)
}

# File name
SaveName <-
  paste0(
    PlotList$Variables[FacetX],
    "PhasePortrait - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", round(1/STEP_SIZE, 2), " - ",
    ifelse(MOVING_AVERAGE_WINDOW > 0, paste("MvAvg", MOVING_AVERAGE_WINDOW, "- "), ""),
    ifelse(LOG_SCALE == T, "Log - ", ""),
    "Arrow Direction.pdf"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)

ggsave(
  # Save vector image
  file.path(COHORT_FOLDER, SaveName),
  # height = 4.76*1.25,
  # width = 11.5*1.25
  height = 4.76,
  width = 11.5
)

