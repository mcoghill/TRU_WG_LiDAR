# Custom Functions

# Get clipping shapes for catalog
catalog_shapes <- function(ctg, clip_shp, id) {
  output <- do.call(rbind, lapply(engine_chunks(ctg), function(xx) {
    st_as_sf(st_as_sfc(st_bbox(xx))) |>
      cbind(t(as.matrix(st_bbox(xx)))) |>
      mutate(filename = paste0(paste(id, xmin, ymin, sep = "_"), ".las"),
             id = id) |>
      relocate(c(id, filename), 1)
  })) |> st_set_agr("constant") |>
    st_intersection(st_geometry(clip_shp))
}

# LAS cleaning and normalization - classification done in DJI Terra
ctg_clean <- function(las) {
  las <- classify_noise(las, ivf(res = 4, n = 15))
  las <- filter_poi(las, Classification != LASNOISE)
  las$Intensity <- as.integer(las$Intensity / 256L)
  return(las)
}

# Intensity normalization
# ctg_range <- function(las, sensor) {
#   las <- readLAS(las)
#   sensor <- sensor %>% 
#     filter(between(gpstime, floor(min(las$gpstime)), ceiling(max(las$gpstime))))
#   R <- get_range(las, sensor)
#   return(R)
# }

ctg_trees <- function(las) {
  las_s <- filter_nth(las, 1)
  las_s <- segment_snags(las_s, wing2015())
}

lmfx <- function (ws, hmin = 2, dist_2d = 3) {
  f = function(las) {
    context <- tryCatch({
      get("lidR.context", envir = parent.frame())
    }, error = function(e) {
      return(NULL)
    })
    lidR:::assert_is_valid_context(lidR:::LIDRCONTEXTITD, 
                                   name = "lmfx")
    # . <- X <- Y <- Z <- treeID <- NULL
    dist_2d = dist_2d^2
    if (is.numeric(ws) & length(ws) == 1) {
    } else if (is.function(ws)) {
      n = nrow(las@data)
      ws = ws(las@data$Z)
      if (!is.numeric(ws)) 
        stop("The function 'ws' did not return correct output. ", 
             call. = FALSE)
      if (any(ws <= 0)) 
        stop("The function 'ws' returned negative or nul values.", 
             call. = FALSE)
      if (anyNA(ws)) 
        stop("The function 'ws' returned NA values.", 
             call. = FALSE)
      if (length(ws) != n) 
        stop("The function 'ws' did not return correct output.", 
             call. = FALSE)
    }  else stop("'ws' must be a number or a function", call. = FALSE)
    # . <- X <- Y <- Z <- treeID <- NULL
    las = lidR::decimate_points(las, lidR::highest(1))
    is_maxima = lidR:::C_lmf(las, ws, hmin, TRUE, lidR:::getThread())
    LM = las@data[is_maxima, .(X, Y, Z)]
    data.table::setnames(LM, c("X", "Y", "Z"), c("x", "y", "z"))
    data.table::setorder(LM, -z)
    detected = logical(nrow(LM))
    detected[1] = TRUE
    for (i in 2:nrow(LM)) {
      distance2D = (LM$x[i] - LM$x[detected])^2 + (LM$y[i] - 
                                                     LM$y[detected])^2
      if (!any(distance2D < dist_2d)) {
        detected[i] = TRUE
      }
    }
    xoffset <- las[["X offset"]]
    yoffset <- las[["Y offset"]]
    zoffset <- las[["Z offset"]]
    
    xscale  <- las[["X scale factor"]]
    yscale  <- las[["Y scale factor"]]
    zscale  <- las[["Z scale factor"]]
    
    xscaled <- as.integer((LM[["x"]] - xoffset)/xscale)
    yscaled <- as.integer((LM[["y"]] - yoffset)/yscale)
    
    LM[, `:=` (Z = z, treeID = lidR:::bitmerge(xscaled, yscaled))]
    
    detected = LM[detected]
    # detected[, `:=`(treeID, 1:.N)]
    
    output = st_as_sf(detected, coords = c("x", "y", "z"), crs = sf::st_crs(las))
    return(output)
  }
  # class(f) <- c(lidR:::LIDRALGORITHMITD)
  f <- plugin_itd(f, omp = TRUE, raster_based = FALSE)
  return(f)
}

ws_fun <- function(x) {
  y <- 2.8 * (-(exp(-0.08*(x-5)) - 1)) + 2
  y[x < 5] <- 2
  y[x > 20] <- 4
  return(y)
}

# my_metrics <- function (x, y, z, i, rn, class, dz = 1, th = 2, zmin = 0, R, G, B) {
#   mts <- stdmetrics(x, y, z, i, rn, class, dz, th, zmin)
#   shp <- stdshapemetrics(x, y, z)
#   R <- as.integer(median(R, na.rm = TRUE))
#   G <- as.integer(median(G, na.rm = TRUE))
#   B <- as.integer(median(B, na.rm = TRUE))
#   metrics <- list(R = R, G = G, B = B)
#   metrics <- c(mts, shp, metrics)
#   return(metrics)
# }

# my_metrics <- function (x, y, z, R, G, B) {
#   mts <- stdtreemetrics(x, y, z)
#   width <- max(x) - min(x)
#   height <- max(y) - min(y)
#   R <- as.integer(median(R, na.rm = TRUE))
#   G <- as.integer(median(G, na.rm = TRUE))
#   B <- as.integer(median(B, na.rm = TRUE))
#   metrics <- list(width = width, height = height, R = R, G = G, B = B)
#   metrics <- c(mts, metrics)
#   return(metrics)
# }

my_metrics <- function (x, y, z, i, rn, class, dz = 1, th = 2, zmin = 0, R, G, B) {
  mts <- c(
    stdmetrics_z(z, dz, th, zmin),
    stdmetrics_i(i, z, class, rn)
  )
  npoints <- length(x)
  convhull_area <- round(lidR:::area_convex_hull(x, y), 3)
  crown_width_x <- max(x) - min(x)
  crown_width_y <- max(y) - min(y)
  R <- as.integer(median(R, na.rm = TRUE))
  G <- as.integer(median(G, na.rm = TRUE))
  B <- as.integer(median(B, na.rm = TRUE))
  metrics <- list(
    npoints = npoints,
    convhull_area = convhull_area,
    crown_width_x = crown_width_x, 
    crown_width_y = crown_width_y, 
    R = R, G = G, B = B
  )
  metrics <- c(mts, metrics)
  return(metrics)
}

# Get a list of tree ID's that are considered snags
snag_fun <- function(chunk, BBPRthrsh_mat) {
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  las <- segment_snags(las, wing2015(BBPRthrsh_mat = BBPRthrsh_mat))
  snag_ids <- las@data[
    , snagCls := fifelse(snagCls > 0, 1, 0)
  ][
    , .(ratio = sum(snagCls) / .N), by = treeID
  ][ratio > 0.5, treeID]
  return(snag_ids)
}

snag_id_fun <- function(las, snags) {
  suppressWarnings(las@data[, snagCls := treeID %in% snags])
  return(las)
}


# Set input directories and create output directories
mission <- "East"
shape_dir <- file.path("./00_Shapes")
dji_dir <- file.path("D:/DJI/PCGSPRO_1736208110/jchurch@tru.ca")
dji_proj <- file.path(dji_dir, paste0("WG_", mission, "_color"))
las_fname <- file.path(dji_proj, "lidars/terra_las/cloud_merged.las")
# sens_fname <- file.path(dji_proj, "lidars/flight_trajectory_recons.json")
traj_fname <- list.files(file.path(dji_proj, "lidars"), pattern = "L_sbet.txt$", full.names = TRUE)
tile_dir <- file.path("./01_tile")
cln_dir <- file.path("./02_clean")
norm_dir <- file.path("./02_clean_norm")
snag_dir <- file.path("./03_snags")
tree_dir <- file.path("./04_tree_seg")
