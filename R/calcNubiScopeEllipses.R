#' calcNubiScopeEllipses
#' 
#' Calculate ellipses/SPDF from observations
#' 
#' @param RawFile File location of ASCII coded temperature readings
#' @param NubiScopePosition SPDF representing position of NubiScope (1 point)
#' @param atmosphere
#' @param cs_radiances Clear sky radiances for NubiScope zenith viewing angles (output of Streamer)
#' @param cs_threshold Clear sky threshold as a function of housing temperature and zenith angle
#' @param horizon Boolean matrix expressing which NubiScope observations see the sky (TRUE)
#' @param default_cloud_height Default cloud height in case no clouds were detected
#' @param vertexes Number of vertexes of final ellipse polygons
#' 
#' @return SpatialPolygonDataFrame with cloud ellipses
#' 
#' @export
#' 
#' @import sp

calcNubiScopeEllipses <- function(RawFile,
                                  NubiScopePosition,
                                  atmosphere,
                                  cs_radiances, 
                                  cs_threshold, 
                                  horizon, 
                                  default_cloud_height, 
                                  vertexes = 50) {
  
  if (length(cs_radiances) != 30)
    stop('radiances vector does not meet NubiScope zenith angles')
  if (is.null(NubiScopePosition@proj4string))
    stop('projection of NubiScope position undefined')
  
  ### reading data
  raw_contents <- readLines(RawFile)
  # measured brightness temperatures [K]
  # columns = zenith, rows = azimuth
  Tbs <- as.matrix(read.table(text = raw_contents, header = FALSE, dec = ".", sep = "\t", skip = 7)[,4:33]) + 273.15
  # reference Temperature (NubiScope housing temperature) in line 4 of the file
  T_ref <- as.numeric(strsplit(raw_contents[4], ":")[[1]][3])
  
  ### constants + functions
  # FOV / 2
  phi <- deg2rad(3 / 2)
  # correction for turning direction:
  # - NubiScope starts N and turns clockwise
  # - mathematical positive starts E and turns anti-clockwise
  azimuth <- deg2rad((440 - seq(5, 355, 10)) %% 360)
  zenith <- deg2rad(seq(1.5, 88.5, 3))
  
  # physical constants
  planck <- 6.62606957e-34 #[J s]
  boltzmann <- 1.3806488e-23 #[J K-1]  
  speedOfLight <- 299792458.0 #[m s-1]
  # NubiScope frequency (wavelength = 8 to 14 µm)
  v_low <- speedOfLight / 14e-6
  v_high <- speedOfLight / 8e-6
  
  # temperature in K
  Plancks_law <- function(temperature, v) {
    (2 * planck * v ^ 3) / (speedOfLight ^ 2 * exp(planck * v / (boltzmann * temperature)) - 1)
  }
  
  # numerical approximation of temperature from Planck's law
  temp2radiance <- function(temperature){
    sapply(temperature, function(t) integral(function(v) Plancks_law(t, v), v_low, v_high))
  }
  
  ### actions
  
  # clear sky thresholds for all zenith angles [W m-2 sr-1]
  cs_thr <- cs_threshold(T_ref, seq_along(zenith))
  
  radiances <- apply(Tbs, 2, temp2radiance)
  is_clear <- t(apply(radiances, 1, function(r) r - cs_radiances < cs_thr))
  
  # substitute all values that shall not be calculated (clouds + horizon) with NA
  Tbs_gaps <- Tbs * horizon * apply(is_clear, 2, function(clear) ifelse(clear, NA, 1))
  
  # Tb -> height
  Tb2height <- function(Tb) 
    approx(x = atmosphere$temp, y = atmosphere$alt, xout = Tb)$y + coordinates(NubiScopePosition)[1, 3]
  # heights with horizon = NA, clear sky = NA
  heights <- apply(Tbs_gaps, 2, Tb2height)
  height_avg <- ifelse(all(is.na(heights)), default_cloud_height, mean(heights, na.rm = TRUE))
  # heigths with clear sky substituted by -average cloud height
  heights_w_clear <- apply(heights, 2, function(i) ifelse(is.na(i), -height_avg, i)) * horizon
  
  theights_w_clear <- t(heights_w_clear)
  zenithTransformed <- cos(phi) ^ 2 * cos(zenith) ^ 2 - sin(phi) ^ 2 * sin(zenith) ^ 2
  
  ellipse_a <- t(theights_w_clear * (sin(phi) * cos(phi) / zenithTransformed))
  ellipse_b <- t(theights_w_clear * (sin(phi) * sin(zenith) / sqrt(zenithTransformed)))
  # numerical problem with zenith 88.5° -> quick fix
  # ellipse_b[,30] <- heights_w_clear[,30] * (sin(phiRad) * sin(alpha[30]) / sqrt(1.1058862159352145E-17))
  ellipse_xp <- t(theights_w_clear * (tan(zenith) + sin(phi) * cos(phi) / zenithTransformed)) * cos(azimuth)
  ellipse_yp <- t(theights_w_clear * (tan(zenith) + sin(phi) * cos(phi) / zenithTransformed)) * sin(azimuth)
  
  
  # build one ellipse from parameters
  build_ellipse <- function(xp, yp, a, b, azimuth, height, radiance, Tb, ID) {
    # ellipse phase angles, reverse order (counterclockwise) for polygon writing
    # polygons have to end at the same point where they start, so add first angle again
    theta <- rev(c(1:vertexes, 1) * 2 * pi / vertexes)
    # calc ellipse points with ellipse equation in parameter form
    coords <- matrix(c(coordinates(NubiScopePosition)[1, 1] + xp + a * cos(theta) * cos(azimuth) - b * sin(theta) * sin(azimuth),
                       coordinates(NubiScopePosition)[1, 2] + yp + a * cos(theta) * sin(azimuth) + b * sin(theta) * cos(azimuth)),
                     ncol = 2, nrow = vertexes + 1)
    polys <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID)), proj4string = NubiScopePosition@proj4string)
    SpatialPolygonsDataFrame(polys, data.frame(height = height, radiance = radiance, Tb = Tb, row.names = row.names(polys)))     
  }
  
  # build up SpatialPolygonsDataFrame
  ellipses_list <- lapply(seq_along(Tbs), function(i) {
    tryCatch({
      if (!is.na(horizon[i]))
        build_ellipse(xp = ellipse_xp[i], 
                      yp = ellipse_yp[i], 
                      a = ellipse_a[i], 
                      b = ellipse_b[i],
                      azimuth = azimuth[(i - 1) %% length(azimuth) + 1],
                      height = heights_w_clear[i],
                      radiance = radiances[i],
                      Tb = Tbs[i],
                      ID = i)
      else
        NULL
    }, error = function(e) NULL)
  })
  
  #if (any(!sapply(ellipses_list, function(x) is.null(x)))) 
  # shapefile(do.call(rbind, Filter(Negate(is.null), ellipses_list)), "poly.shp", overwrite = TRUE)
  do.call(rbind, Filter(Negate(is.null), ellipses_list))
  #else
  #NULL
  # ellipse_b
}