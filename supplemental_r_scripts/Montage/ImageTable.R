# Get image
IMAGE_PATH <- file.path(CELL_PATH, paste0(REFERENCE_PROTEIN, "_intensity_ref.tif"))
img = ijtiff::read_tif(IMAGE_PATH, frames = IMPORT_FRAMES)
FrameFx <- function(FrameX){
  tryCatch({
    
    Z = img[,,,FrameX]
    X = NROW(Z)
    X = rep(1:X, NCOL(Z))
    
    Y = NCOL(Z)
    Y = rep(1:Y, each = NROW(Z))
    
    Z = as.vector(Z)
    
    plot_img = cbind(X, Y, Z)
    plot_img = as_tibble(plot_img)
    
    plot_img <-
      plot_img %>%
      mutate(
        y = X,
        x = Y,
        t = IMPORT_FRAMES[FrameX] - min(Table$FRAME) + 1,
        z = Z,
        # z = Z/Table$TOTAL_INTENSITY[1]/Table$NORMALIZED_INTENSITY[1]*9
      ) %>%
      select(-c(
        X,
        Y,
        Z
      ))
    
    return(plot_img)
  }, error = function(e){print(paste("     ERROR with FrameFx FrameX =", FrameX))})
}
plot_img <- mclapply(1:NROW(IMPORT_FRAMES), FrameFx)
plot_img <- data.table::rbindlist(plot_img)

# Get table of track coordinates
plot_spot <-
  Table %>%
  filter(
    UNIVERSAL_TRACK_ID == SELECT_UNIVERSAL_TRACK_ID,
    FRAME %in% IMPORT_FRAMES
  ) %>%
  mutate(
    t = FRAME - min(Table$FRAME) + 1
  ) %>%
  select(
    t,
    NORMALIZED_INTENSITY,
    POSITION_X,
    POSITION_Y
  )

# Combine
ref_plot_img_spot <- merge(plot_spot, plot_img, by = "t")

ref_plot_img_spot <-
  ref_plot_img_spot %>%
  group_by(
    t
  ) %>%
  mutate(
    x = x - POSITION_X,
    x = round(x),
    y = y - POSITION_Y,
    y = round(y),
    protein = REFERENCE_PROTEIN
  ) %>%
  filter(
    x >= -BOX_SIZE,
    y >= -BOX_SIZE,
    x <= BOX_SIZE,
    y <= BOX_SIZE
  ) %>% 
  ungroup() %>% 
  mutate(
    z = z - quantile(z, 0.50),
    z = ifelse(z < 1, 0, z),
    z = z/quantile(z, 0.9999),
    z = ifelse(z > 1, 1, z),
    color = rgb(z,z,z)
  )

inferno_list = viridis::inferno(2^10+1)
ref_plot_img_spot$color = round(ref_plot_img_spot$z*2^10)+1
ref_plot_img_spot$color = inferno_list[ref_plot_img_spot$color]


QUERY_PROTEIN = Table$COMPLEMENTARY_PROTEIN_1[1]
# Get image
IMAGE_PATH <- file.path(CELL_PATH, paste0(QUERY_PROTEIN, "_intensity_ref.tif"))
img = ijtiff::read_tif(IMAGE_PATH, frames = IMPORT_FRAMES)

FrameFx <- function(FrameX){
  tryCatch({
    
    Z = img[,,,FrameX]
    X = NROW(Z)
    X = rep(1:X, NCOL(Z))
    
    Y = NCOL(Z)
    Y = rep(1:Y, each = NROW(Z))
    
    Z = as.vector(Z)
    
    plot_img = cbind(X, Y, Z)
    plot_img = as_tibble(plot_img)
    
    plot_img <-
      plot_img %>%
      mutate(
        y = X,
        x = Y,
        t = IMPORT_FRAMES[FrameX] - min(Table$FRAME) + 1,
        z = Z,
        # z = Z/Table$TOTAL_INTENSITY[1]/Table$NORMALIZED_INTENSITY[1]*9
      ) %>%
      select(-c(
        X,
        Y,
        Z
      ))
    
    return(plot_img)
  }, error = function(e){print(paste("     ERROR with FrameFx FrameX =", FrameX))})
}
plot_img <- mclapply(1:NROW(IMPORT_FRAMES), FrameFx)
plot_img <- data.table::rbindlist(plot_img)

# Get table of track coordinates
plot_spot <-
  Table %>%
  filter(
    UNIVERSAL_TRACK_ID == SELECT_UNIVERSAL_TRACK_ID,
    FRAME %in% IMPORT_FRAMES
  ) %>%
  mutate(
    t = FRAME - min(Table$FRAME) + 1,
    NORMALIZED_INTENSITY = COMPLEMENTARY_NORMALIZED_INTENSITY_1
  ) %>%
  select(
    t,
    NORMALIZED_INTENSITY,
    POSITION_X,
    POSITION_Y
  )

# Combine
qry_plot_img_spot <- merge(plot_spot, plot_img, by = "t")

qry_plot_img_spot <-
  qry_plot_img_spot %>%
  group_by(
    t
  ) %>%
  mutate(
    x = x - POSITION_X,
    x = round(x),
    y = y - POSITION_Y,
    y = round(y),
    protein = QUERY_PROTEIN
  ) %>%
  filter(
    x >= -BOX_SIZE,
    y >= -BOX_SIZE,
    x <= BOX_SIZE,
    y <= BOX_SIZE
  ) %>% 
  ungroup() %>% 
  mutate(
    z = z - quantile(z, 0.50),
    z = ifelse(z < 1, 0, z),
    z = z/quantile(z, 0.9999),
    z = ifelse(z > 1, 1, z),
    # For white
    # color = rgb(z,z,z)
  )


inferno_list = viridis::inferno(2^10+1)
qry_plot_img_spot$color = round(qry_plot_img_spot$z*2^10)+1
qry_plot_img_spot$color = inferno_list[qry_plot_img_spot$color]
  
plot_img_spot <- rbindlist(list(ref_plot_img_spot, qry_plot_img_spot))

# Get merger color
merge_plot_img_spot <-
  plot_img_spot %>% 
  group_by(
    t,
    protein
  ) %>% 
  mutate(
    ref_z = ifelse(protein == USE_REFERENCE_PROTEIN, z, NA),
    qry_z = ifelse(protein == QUERY_PROTEIN, z, NA)
  ) %>% 
  group_by(
    t,
    x,
    y
  ) %>% 
  tidyr::fill(
    qry_z,
    .direction = "updown"
  ) %>% 
  tidyr::drop_na() %>% 
  group_by(
    t
  ) %>% 
  mutate(
    z = 0,
    color = rgb(qry_z, ref_z, qry_z),
    protein = "Merge"
  ) %>% 
  select(-c(
    ref_z,
    qry_z
  ))

# Merge all tables
plot_img_spot <- rbindlist(list(merge_plot_img_spot, plot_img_spot))

plot_img_spot <-
  plot_img_spot %>% 
  mutate(
    protein = factor(protein, levels = c(USE_REFERENCE_PROTEIN, QUERY_PROTEIN, "Merge"))
  )

RealTime <- NULL
RealTime$t = sort(unique(plot_img_spot$t))
if(is.null(SELECT_FRAMES_ADJSUTED)){
  RealTime$real_t <- sort((Table$TIME - min(Table$TIME)))
} else{
  RealTime$real_t <- sort((Table$TIME - min(Table$TIME))[SELECT_FRAMES_ADJSUTED])
}
RealTime$real_t <- round(RealTime$real_t, 0)
RealTime <- as.data.table(RealTime)

plot_img_spot <- merge(plot_img_spot, RealTime, by = "t")

temp_real_t <- unique(plot_img_spot$real_t)
temp_real_t <- sort(temp_real_t)
temp_real_t <- paste(temp_real_t, "s")

plot_img_spot$real_t <- paste(plot_img_spot$real_t, "s")
plot_img_spot$real_t <- factor(plot_img_spot$real_t, levels = temp_real_t)
