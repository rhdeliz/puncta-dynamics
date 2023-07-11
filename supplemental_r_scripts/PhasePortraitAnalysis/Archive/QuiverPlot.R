# CLEAN ENVIRONMENT----
remove(list = ls())
gc(reset = TRUE)
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Table location
TABLE_PATH = "/Users/u_deliz/PhaseCPResults.csv.gz"
# TABLE_PATH = "/Users/u_deliz/Downloads/Victoria_s preliminary trials/CombinedTables.csv.gz"

# Order of protein appearance
PROTEIN_ORDER = c("IRAK4", "IRAK1")

# TIME PARAMETERS
# Max number of frames included
MAX_FRAMES = 100 # s
# Lifetime minimum
LIFETIME_MIN = 6 # s
# Lead/Lag factor
LEAD_LAG = 1 #frame(s) before and __ frame(s) after

# SIGNAL INTENSITY PARAMETERS
# Max starting intensity allowed (Should starting intensity be low)
MAX_STARTING_INTENSITY = NA # molecules. NA for no limit
# Max intensity must be high
HIGH_MAX_INT = F #TRUE or FALSE
# Grid size for binning and plots
STEP_SIZE = .5 # molecules
# Log scale for plots
LOG_SCALE = T #FALSE for linear scale
# Limit intensity rage in plot
SCALE_LIMITS = F#TRUE or FALSE

# Load libraries
library(dplyr)
library(ggquiver)
library(ggplot2)
library(ggdark)
library(scales)
library(parallel)
library(ggforce)
library(data.table)

# Transform step size
STEP_SIZE = 1/STEP_SIZE

# Import compressed table
Table <- data.table::fread(TABLE_PATH)
Table <- Table %>% distinct()

# Filter table
FilteredTable <-
  Table %>% 
  mutate(
    DATE = substr(IMAGE, 0, 8)
  ) %>% 
  filter(
    # Keep only MyD88 as reference
    REFERENCE_PROTEIN == "MyD88",
    # Keep only ligand density at 32 mol./µm^-2
    # LIGAND_DENSITY_CAT == 10,
    # DATE != 20190718 | LIGAND_DENSITY_CAT != 0.32,
    DATE != 20200205 | LIGAND_DENSITY_CAT != 3.20,
    # DATE != 20190705 | LIGAND_DENSITY_CAT != 32.00,
    # QUERY_PROTEIN == "IRAK1"
  )

# Split table to run faster
SplitTable <-
  FilteredTable %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste(IMAGE, CELL, sep = "...")
  ) %>% 
  ungroup() %>% 
  group_split(
    UNIVERSAL_CELL_ID
  )

# Keep only tibbles
SplitTable <- unclass(SplitTable)
SplitTable <- SplitTable[c(1:NROW(SplitTable))]

# Small to large
SplitTable <- SplitTable[order(sapply(SplitTable,nrow))]
SplitTable <- SplitTable[c(1:NROW(SplitTable))]

# Loop in parallel to speed it up
TableFx <- function(TableX){
  tryCatch({
    library(tictoc)
    tic()
    # Temporary table
    TempTable <- as.data.table(TableX)
    TempTable <- TempTable[, LIFETIME := max(FRAMES_ADJUSTED), by = UNIVERSAL_TRACK_ID]
    
    TempTable <-
      TempTable %>% 
      filter(
        FRAMES_ADJUSTED <= MAX_FRAMES,
        LIFETIME >= LIFETIME_MIN
      ) #%>% 
    # mutate(
    #   # Get max intensity of track
    #   MAX_REFERENCE_INTENSITY = max(REFERENCE_TOTAL_INTENSITY),
    #   MAX_QUERY_INTENSITY = max(QUERY_TOTAL_INTENSITY),
    #   # Does it pass a threshold? (if applicable)
    #   MAX_REF_INT_TEST = ifelse(HIGH_MAX_INT, MAX_REFERENCE_INTENSITY >= 4, T),
    #   MAX_QRY_INT_TEST = ifelse(HIGH_MAX_INT, MAX_QUERY_INTENSITY >= 2, T),
    #   # Get starting intensity
    #   REF_START_INTENSITY = sum(case_when(FRAMES_ADJUSTED == 0 ~ ORIG_REFERENCE_TOTAL_INTENSITY), na.rm = T),
    #   QRY_START_INTENSITY = sum(case_when(FRAMES_ADJUSTED == 0 ~ ORIG_QUERY_TOTAL_INTENSITY), na.rm = T),
    #   # Is it small enough? Only if applicable
    #   REF_START_INTENSITY_TEST = ifelse(is.na(MAX_STARTING_INTENSITY), T, REF_START_INTENSITY <= MAX_STARTING_INTENSITY),
    #   QRY_START_INTENSITY_TEST = ifelse(is.na(MAX_STARTING_INTENSITY), T, QRY_START_INTENSITY <= MAX_STARTING_INTENSITY),
    # ) %>%
    # filter(
    #   MAX_REF_INT_TEST == T,
    #   MAX_QRY_INT_TEST == T,
    #   REF_START_INTENSITY_TEST == T,
    #   QRY_START_INTENSITY_TEST == T
    # )
    
    # Sort protein names
    TempTable$QUERY_PROTEIN = factor(TempTable$QUERY_PROTEIN, levels = PROTEIN_ORDER)
    
    TempTable <- as.data.table(TempTable)
    TempTable <- TempTable[, LAG_REFERENCE_TOTAL_INTENSITY := shift(.(REFERENCE_TOTAL_INTENSITY),LEAD_LAG, type = "lag"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, LEAD_REFERENCE_TOTAL_INTENSITY := shift(.(REFERENCE_TOTAL_INTENSITY),LEAD_LAG, type = "lead"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, LAG_QUERY_TOTAL_INTENSITY := shift(.(QUERY_TOTAL_INTENSITY),LEAD_LAG, type = "lag"), by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, LEAD_QUERY_TOTAL_INTENSITY := shift(.(QUERY_TOTAL_INTENSITY),LEAD_LAG, type = "lead"), by = UNIVERSAL_TRACK_ID]
    
    # Calculate change in intensity
    TempTable <- TempTable[, DELTA_REFERENCE_TOTAL_INTENSITY := LEAD_REFERENCE_TOTAL_INTENSITY-LAG_REFERENCE_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    TempTable <- TempTable[, DELTA_QUERY_TOTAL_INTENSITY := LEAD_QUERY_TOTAL_INTENSITY-LAG_QUERY_TOTAL_INTENSITY, by = UNIVERSAL_TRACK_ID]
    
    # Save for later
    TempTable$ORIG_REFERENCE_TOTAL_INTENSITY <- TempTable$REFERENCE_TOTAL_INTENSITY
    TempTable$ORIG_QUERY_TOTAL_INTENSITY <- TempTable$QUERY_TOTAL_INTENSITY
    
    # Get bin for reference
    TempTable <- TempTable[, REFERENCE_TOTAL_INTENSITY := fifelse(rep(LOG_SCALE, NROW(TempTable)), log(REFERENCE_TOTAL_INTENSITY, 2),REFERENCE_TOTAL_INTENSITY)]
    TempTable <- TempTable[, REFERENCE_TOTAL_INTENSITY := round(REFERENCE_TOTAL_INTENSITY*STEP_SIZE)/STEP_SIZE]
    TempTable <- TempTable[, REFERENCE_TOTAL_INTENSITY := fifelse(rep(LOG_SCALE, NROW(TempTable)), 2^REFERENCE_TOTAL_INTENSITY, REFERENCE_TOTAL_INTENSITY)]
    # Get bin for query
    TempTable <- TempTable[, QUERY_TOTAL_INTENSITY := fifelse(rep(LOG_SCALE, NROW(TempTable)), log(QUERY_TOTAL_INTENSITY, 2),QUERY_TOTAL_INTENSITY)]
    TempTable <- TempTable[, QUERY_TOTAL_INTENSITY := round(QUERY_TOTAL_INTENSITY*STEP_SIZE)/STEP_SIZE]
    TempTable <- TempTable[, QUERY_TOTAL_INTENSITY := fifelse(rep(LOG_SCALE, NROW(TempTable)), 2^QUERY_TOTAL_INTENSITY, QUERY_TOTAL_INTENSITY)]
    
    # Check if within limits (if applicable)
    TempTable <- TempTable[, REF_INT_WITHIN_RANGE :=
                             fifelse(
                               rep(SCALE_LIMITS, NROW(TempTable)),
                               REFERENCE_TOTAL_INTENSITY >= .5 & REFERENCE_TOTAL_INTENSITY <= 7.5,
                               T
                             )
    ]
    TempTable <- TempTable[, QRY_INT_WITHIN_RANGE :=
                             fifelse(
                               rep(SCALE_LIMITS, NROW(TempTable)),
                               QUERY_TOTAL_INTENSITY >= .5 & QUERY_TOTAL_INTENSITY <= 5,
                               T
                             )
    ]
    
    TempTable <-
      TempTable %>% 
      filter(
        REF_INT_WITHIN_RANGE == T,
        QRY_INT_WITHIN_RANGE == T
      ) %>%
      # Delete filtering variables
      select(
        LIGAND_DENSITY_CAT,
        
        REFERENCE_PROTEIN,
        QUERY_PROTEIN,
        
        FRAMES_ADJUSTED,
        
        ORIG_REFERENCE_TOTAL_INTENSITY,
        ORIG_QUERY_TOTAL_INTENSITY,
        
        REFERENCE_TOTAL_INTENSITY,
        QUERY_TOTAL_INTENSITY,
        
        DELTA_REFERENCE_TOTAL_INTENSITY,
        DELTA_QUERY_TOTAL_INTENSITY
      )
    
    toc()
    return(TempTable)
    # Error message
  }, error = function(e) {print("Error with CellFx")})}
# Run in parallel
statTable <- mclapply(SplitTable, TableFx, mc.cores = detectCores(logical = F))
# Combine results from all cells
statTable <- data.table::rbindlist(statTable)
# Remove NA values (those that don't have any lead/lag)
statTable <- statTable %>% tidyr::drop_na()

# Generate a table of what's the average reference and query intensities over time
AverageIntByTime <- 
  statTable %>% 
  filter(
    # Focus only on growth events
    # DELTA_REFERENCE_TOTAL_INTENSITY > 0,
      # DELTA_QUERY_TOTAL_INTENSITY >  0
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  summarize(
    REFERENCE_TOTAL_INTENSITY = mean(ORIG_REFERENCE_TOTAL_INTENSITY),
    QUERY_TOTAL_INTENSITY = mean(ORIG_QUERY_TOTAL_INTENSITY),
  ) %>% 
  arrange(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    FRAMES_ADJUSTED
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    # Smooth results
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 51),
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 21),
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 1, n = 3),
    REFERENCE_TOTAL_INTENSITY = signal::sgolayfilt(REFERENCE_TOTAL_INTENSITY, p = 3, n = 51),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 3, n = 21),
    QUERY_TOTAL_INTENSITY = signal::sgolayfilt(QUERY_TOTAL_INTENSITY, p = 1, n = 3)
  )

# Summarize results dealing with either reference OR query increasing
ppMTIncrease <-
  statTable %>% 
  filter(
    # Filter to have either reference OR query increasing
    DELTA_REFERENCE_TOTAL_INTENSITY > 0 |
    DELTA_QUERY_TOTAL_INTENSITY >  0
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY
  ) %>% 
  summarize(
    # Get number of spots in bin so that bins with few spots can be filtered out
    N = n(),
    # Get median of intensity changes. Mean could skew the data, thus not used
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T),
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
  ) %>% 
  mutate(
    # Keep only N's in the 50th percentile
    N_TEST = ifelse(N >= quantile(N, 0.50), T, F),
    # Get absolute magnitude so that the log can be taken later
    MAGNITUDE = abs(DELTA_REFERENCE_TOTAL_INTENSITY) + abs(DELTA_QUERY_TOTAL_INTENSITY),
    # Get angle to account for negative magnitudes
    RAD_ANGLE = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
    # Make angle
    DEG_ANGLE = RAD_ANGLE*180/pi+180,
    DEG_ANGLE = DEG_ANGLE/360,
    # Give facet name
    CATEGORY = "Increase"
  ) %>% 
  filter(
    # Filter out bins with small N
    N_TEST == T
  ) %>%
  mutate(
    # 95th percentile to make sure large triangles don't alter results
    MAGNITUDE = MAGNITUDE/quantile(MAGNITUDE, .95, na.rm = T)
  ) %>% 
  tidyr::drop_na()

# Summarize results dealing with either reference OR query decreasing
# Detailed description above. Difference highlighted below
ppMTDecrease <-
  statTable %>% 
  filter(
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
  ) %>% 
  filter(
    # Only this part changed as well as the label
    DELTA_REFERENCE_TOTAL_INTENSITY < 0 |
      DELTA_QUERY_TOTAL_INTENSITY <  0
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY
  ) %>% 
  summarize(
    N = n(),
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T),
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
  ) %>% 
  mutate(
    N_TEST = ifelse(N >= quantile(N, .5), T, F),
    MAGNITUDE = abs(DELTA_REFERENCE_TOTAL_INTENSITY) + abs(DELTA_QUERY_TOTAL_INTENSITY),
    RAD_ANGLE = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
    DEG_ANGLE = RAD_ANGLE*180/pi+180,
    DEG_ANGLE = DEG_ANGLE/360,
    # This label
    CATEGORY = "Decrease"
  ) %>% 
  filter(
    N_TEST == T
  ) %>%
  mutate(
    MAGNITUDE = MAGNITUDE/quantile(MAGNITUDE, .95, na.rm = T)
  ) %>% 
  tidyr::drop_na()

# Summarize results dealing with either reference OR query decreasing
# Detailed description above. Difference highlighted below
ppMTAll <-
  statTable %>% 
  filter(
    REFERENCE_TOTAL_INTENSITY > 0,
    QUERY_TOTAL_INTENSITY > 0,
  ) %>% 
  filter(
    # 0s will be counted if this area is commented out
    # However, if you want to combine decrease and increase, then run it
    DELTA_REFERENCE_TOTAL_INTENSITY != 0 |
    DELTA_QUERY_TOTAL_INTENSITY !=  0
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN,
    REFERENCE_TOTAL_INTENSITY,
    QUERY_TOTAL_INTENSITY
  ) %>% 
  summarize(
    N = n(),
    DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
    DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T),
  ) %>% 
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    N_TEST = ifelse(N >= quantile(N, .5), T, F),
    MAGNITUDE = abs(DELTA_REFERENCE_TOTAL_INTENSITY) + abs(DELTA_QUERY_TOTAL_INTENSITY),
    RAD_ANGLE = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
    DEG_ANGLE = RAD_ANGLE*180/pi+180,
    DEG_ANGLE = DEG_ANGLE/360,
    CATEGORY = "All"
  ) %>% 
  filter(
    N_TEST == T
  ) %>%
  group_by(
    LIGAND_DENSITY_CAT,
    REFERENCE_PROTEIN,
    QUERY_PROTEIN
  ) %>% 
  mutate(
    MAGNITUDE = MAGNITUDE/quantile(MAGNITUDE, .95, na.rm = T)
  ) %>% 
  tidyr::drop_na()

# Combine tables that separate deltas as negative, positive and all
ppMT <- bind_rows(ppMTAll, ppMTDecrease, ppMTIncrease)
# To ensure facets apear in the right order
ppMT$CATEGORY = factor(ppMT$CATEGORY, levels = c("Decrease", "All", "Increase"))

if(LOG_SCALE == T){
  
  # Plot
  ggplot(
  ) +
    # # Average intensity binned by time
    # geom_path(
    #   data = AverageIntByTime,
    #   aes(
    #     x = REFERENCE_TOTAL_INTENSITY,
    #     y = QUERY_TOTAL_INTENSITY,
    #     # To account for time, I colored the segments
    #     # The later in time, the lighter the color
    #     color = FRAMES_ADJUSTED/MAX_FRAMES
    #   )
    # ) +
    geom_vline(
      xintercept = c(3, 6)
    ) +
    # Make triangles of the arrows
    geom_regon(
      # data = ppMT,
      data = ppMT %>% filter(CATEGORY=="Increase"),
      aes(
        x0 = REFERENCE_TOTAL_INTENSITY,
        y0 = QUERY_TOTAL_INTENSITY,
        # Size of triangle is the magnitude
        r = MAGNITUDE/2,
        # Angle of triangle is the derivative
        angle =  -(RAD_ANGLE+pi)+2*pi,
        # 3 sides for a triangle
        sides = 3,
        # Angle of the derivatives is the color
        fill = hsv(DEG_ANGLE),
        # Redundancy so that the outline of the triangle also tells you the magnitude
        # Also helps pair the arrow with the triangle
        color = MAGNITUDE
      ),
      # Make transparent to see what's happening behind
      alpha = 0.75
    ) +
    # Puts arrows inside triangles so that we know which tip of the triangle to look at
    geom_quiver(
      # data = ppMT,
      data = ppMT %>% filter(CATEGORY=="Increase"),
      aes(
        x = REFERENCE_TOTAL_INTENSITY,
        y = QUERY_TOTAL_INTENSITY,
        u = DELTA_REFERENCE_TOTAL_INTENSITY,
        v = DELTA_QUERY_TOTAL_INTENSITY,
        # Arrow color is the magnitude
        # Same as the triangle so that they can be visually paired
        color = MAGNITUDE
      ),
      # Arrow size blowup
      vecsize = 2,
      # Center the arrow
      center = T
    ) +
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
      y = "Query Size (Box has Name)"
    ) +
    facet_grid(
      # QUERY_PROTEIN~CATEGORY#~LIGAND_DENSITY_CAT
      CATEGORY~QUERY_PROTEIN#+LIGAND_DENSITY_CAT
    ) +
    # Make dark to make the colors pop in a monitor
    # Use theme classic for printing
    dark_theme_classic(
      base_size = 18
    ) +
    theme(
      legend.position = "none"
    )
    ggsave(
      # Save vector image
      ifelse(SCALE_LIMITS, "PhasePortrait_Log Limited.pdf", "PhasePortrait_Log.pdf"),
      height = 4.76*2,
      width = 11.5*2
      # height = 4.03,
      # width = 5.64
    )
} else{
  
  # Plot
  ggplot(
  ) +
    # Average intensity binned by time
    # geom_path(
    #   data = AverageIntByTime,
    #   aes(
    #     x = REFERENCE_TOTAL_INTENSITY,
    #     y = QUERY_TOTAL_INTENSITY,
    #     # To account for time, I colored the segments
    #     # The later in time, the lighter the color
    #     color = FRAMES_ADJUSTED/MAX_FRAMES
    #   )
    # ) +
    # Make triangles of the arrows
    geom_regon(
      data = ppMT,
      # data = ppMT %>% filter(CATEGORY=="Increase"),
      aes(
        x0 = REFERENCE_TOTAL_INTENSITY,
        y0 = QUERY_TOTAL_INTENSITY,
        # Size of triangle is the magnitude
        r = MAGNITUDE/2,
        # Angle of triangle is the derivative
        angle = -(RAD_ANGLE+pi)+2*pi,
        # 3 sides for a triangle
        sides = 3,
        # Angle of the derivatives is the color
        fill = hsv(DEG_ANGLE),
        # Redundancy so that the outline of the triangle also tells you the magnitude
        # Also helps pair the arrow with the triangle
        color = MAGNITUDE
      ),
      # Make transparent to see what's happening behind
      alpha = 0.75
    ) +
    # Puts arrows inside triangles so that we know which tip of the triangle to look at
    geom_quiver(
      data = ppMT,
      # data = ppMT %>% filter(CATEGORY=="Increase"),
      aes(
        x = REFERENCE_TOTAL_INTENSITY,
        y = QUERY_TOTAL_INTENSITY,
        v = DELTA_REFERENCE_TOTAL_INTENSITY,
        u = DELTA_QUERY_TOTAL_INTENSITY,
        # Arrow color is the magnitude
        # Same as the triangle so that they can be visually paired
        color = MAGNITUDE
      ),
      # Arrow size blowup
      vecsize = 2,
      # Center the arrow
      center = T
    ) +
    # Color of magnitude. Goes from gray to white
    scale_color_gradient(
      low = "#303030",
      high = "white"
    ) +
    # Use angle of arrows to color the triangle. Full spectrum
    scale_fill_identity(guide = "none") +
    labs(
      x = "MyD88 Size",
      y = "Query Size (Box has Name)"
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
      legend.position = "none"
    )
    ggsave(
      # Save vector image
      ifelse(SCALE_LIMITS, "PhasePortrait_Linear Limited.pdf", "PhasePortrait_Linear.pdf"),
      height = 4.76*2,
      width = 11.5*2
      # height = 4.03,
      # width = 5.64
    )
  
}
# 
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

