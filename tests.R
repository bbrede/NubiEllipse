# tests
setwd('D:/Test/streamer/N110829/')

# rm(list=ls())

streamer_path <- 'streamer_win.exe'

atmos_tab <- read.table('N1108290.950_profiles.txt', sep = '\t', dec = '.', header = TRUE)
atmos <- makeAtmosphere(altitude = atmos_tab$height, temperature = atmos_tab$temperature, 
                        pressure = atmos_tab$pressure, watervapour = atmos_tab$moisture)
atmos$temp <- atmos$temp + 273.15

nubipos <- SpatialPointsDataFrame(matrix(c(439925, 5784598, 127), ncol = 3, byrow = TRUE), data.frame(ID = 1), 
                                  proj4string = CRS('+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))

stream <- runStreamer(overwrite = TRUE, strptime('N1108290950', 'N%y%m%d%H%M', 'UTC'), nubipos, 'N1108290.955.inp', 'N1108290.955.des', atmos, streamer_path)

cs_coeffs <- read.table('clearsky_coeffs.txt', header = TRUE, sep = '\t')
cs_threshold <- function(tref, i) cs_coeffs[i, 1] + tref * cs_coeffs[i, 2]
horizon <- matrix(1, nrow = 36, ncol = 30)
horizon[,30] <- NA

ell <- calcNubiScopeEllipses(raw_file = 'N1108290.955', 
                             NubiScope_Position = nubipos, 
                             atmosphere = atmos, 
                             cs_radiances = stream, 
                             cs_threshold = cs_threshold, 
                             horizon = horizon, 
                             default_cloud_height = 2000)



