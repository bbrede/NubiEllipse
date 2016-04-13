#' makeAtmosphere
#' 
#' Create named list with atmospheric profiles and/or column values
#' 
#' @return Named list
#' 
#' @export

makeAtmosphere <- function(altitude = NULL, 
                           pressure = NULL, 
                           temperature = NULL, 
                           watervapour = NULL,
                           ozone = NULL,
                           haze = NULL) {
  
  # make sure only same lengths for profiles exist
  #   ls <- sapply(list(altitude, pressure, temperature, watervapour, ozone, haze), length)
  #   possibles <- unique(ls)
  #   if ()
  
  list(alt = altitude, pres = pressure, temp = temperature, wv = watervapour, oz = ozone, hz = haze)
  
}