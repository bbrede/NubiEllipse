#' Run streamer
#' 
#' Wrapper for atmospheric radiative transfer model Streamer
#' 
#' @param dateTime POSIX with DateTime (any time zone possible)
#' @param NubiScopePosition SPDF representing position of NubiScope (1 point)
#' @param inpFile inpFile filename
#' @param desFile desFile filename
#' @param overwrite Overwrite inp and desFiles if they already exist?
#' @param atmosphere Named list with atmospheric profiles and/or column values (created with makeAtmosphere)
#' @param streamerPath File path to Streamer executable
#' 
#' @return Numeric vector with modeled clear sky radiances
#' 
#' @references \link{http://stratus.ssec.wisc.edu/streamer/streamer.html}
#' 
#' @export
#' 
#' @import rgdal
#' @import sp
#' @import pracma
#' @import lubridate


runStreamer <- function(dateTime, 
                        location, 
                        inpFile, 
                        desFile,
                        overwrite = FALSE, 
                        atmosphere,
                        streamerPath, 
                        cloudFile = NULL,
                        clouds = NULL,
                        nover = NULL,
                        NSTRLONG = 4, 
                        NCOEF = 24, 
                        STDPROF = 2, 
                        IZD = 2, 
                        ITD = 1, 
                        IWV = 1, 
                        IO3 = 1, 
                        ICTOP = 4, 
                        ICTHK = 2, 
                        theta = seq(1.5, 88.5, 3)) {
  
  library(rgdal)
  library(sp)
  library(pracma)
  library(lubridate)
  
  if (!file.exists(streamer_path))
    stop('streamer exe not found!')
  
  # prepare input file and run streamer in case description file is not available
  if (!file.exists(desFile) | overwrite) {
    if (!file.exists(inpFile) | overwrite) {
      
      # compose TRUE string as read by Streamer ('.TRUE.' or '.FALSE.')
      compBool <- function(x) paste0('.', x, '.')
      
      inpOptions <- c(
        'OPTIONS',
        compBool(FALSE), # L2 FLUXES
        compBool(FALSE), # L3 IR106
        compBool(FALSE), # L4 CLDFORCE
        paste('2', NSTRLONG, collapse = ' '), # L5 NSTRLONG
        NCOEF, # L6 NCOEF
        compBool(TRUE), # L7 GASABS
        compBool(FALSE), # L8 RAYLISHRT
        '0', # L9 ALBTYPE
        '4', # L10 EMISSTYPE
        paste(STDPROF, compBool(TRUE), collapse = ' '), # L11 STDPROF SPACE
        '2 1', # L12 AERMOD AERVERT
        paste(IZD, ITD, IWV, IO3, ICTOP, ICTHK, 3, collapse = ' '), # L13
        2, # L14 OUTLEVS
        compBool(TRUE), # L15 DESCRIP
        desFile, # L16
        compBool(FALSE), # L17
        '<user file name>', # L18
        '<bandweights file name>', # L19
        compBool(!is.null(cloudFile)),
        ifelse(!is.null(cloudFile), cloudFile, '<cloud file name>'),
        '<brdf file name>')
      
      # sort thetas with increasing order
      slice <- function(x, n) {
        N <- length(x)
        lapply(seq(1, N, n), function(i) x[i:min(i + n - 1, N)])
      }
      
      # maximum 20 thetas per CASE
      thetaBins <- slice(-cos(deg2rad(theta)), 20)
      
      dateTimeUTC <- with_tz(dateTime, 'UTC')
      pos <- spTransform(location, CRS('+proj=longlat +datum=WGS84 +ellps=WGS84'))
      
      # which form of values do the variables represent?
      profileCode <- Vectorize(function(atmosname) {
        x <- atmosphere[[atmosname]]
        len <- length(x)
        if (len > 1) 1 # profile
        else if (len == 0) 2 # case: NULL or NA, standard profile
        else if (is.na(x)) stop('This option is not implemented here!') # read from previous CASE not implemented
        else if (len == 1) 0 # total column value
      })
      
      atmosnames <- c('alt', 'pres', 'temp', 'wv', 'oz', 'haze')
      
      # print blanks if NULL
      null_print <- function(x) ifelse(is.null(x), '', x)
      
      inpCases <- lapply(thetaBins, function(tbin) {
        header <- c(
          'CASE',
          'casename',
          paste(year(dateTimeUTC), month(dateTimeUTC), day(dateTimeUTC), round(hour(dateTimeUTC) + minute(dateTimeUTC) / 60, 2), coordinates(pos)[1,2], coordinates(pos)[1,1], '-1.0', collapse = ' '),
          paste(length(tbin), paste(tbin, collapse = ' '), '1 0.0'),
          '36 62', # bstart, bend
          '1.0 1 1 1.0', # channel
          '-999.0 0', # tsurf, emis
          null_print(clouds), 
          null_print(nover), 
          paste(paste(sapply(atmosnames, function(name) profileCode(name)), collapse = ' '),
                max(sapply(atmosphere, length)), 
                paste(rep('1.0', 6), collapse = ' '), collapse = ' ')
        )
        
        atmosTable <- data.frame(do.call(cbind, lapply(atmosnames, function(name) if (profileCode(name) == 1) atmosphere[[name]])))
        names(atmosTable) <- atmosnames[profileCode(atmosnames) == 1]
        
        profiles <- paste(sapply(1:nrow(atmosTable), function(i) paste('\t', atmosTable[i,], collapse = '\t')), collapse = '\n')
        
        columnValues <- paste(
          ifelse(
            # case: no column values/only profiles
            all(profileCode(atmosnames) != 0), '',
            # case: at least one column value
            atmosphere[[atmosnames[profileCode(atmosnames) == 0]]]), 
          '\n', collapse = '')
        
        paste(paste(header, collapse = '\n'), profiles, columnValues, sep = '\n')
      }
      )
      
      inpCases$sep <- ''
      
      inpContents <- paste(c(inpOptions, do.call(paste, inpCases)), collapse = '\n')
      
      writeLines(inpContents, inpFile)
    }
    
    if (!file.exists(streamerPath))
      stop("streamer exe not found")
    
    # run Streamer
    system(paste('cmd.exe /c', streamerPath, inpFile), show.output.on.console = FALSE)
  }
  
  # read results
  
  # read.table(text = ...) => alternative
  #pattern <- '[:digit:]{1,4}[:print:]*-[:digit:]{1,2}[.][:digit:]{1,2}[:print:]*[.]00[:print:]*'
  radiance_pattern <- '[[:digit:]]*[.][[:digit:]]*E[+][[:digit:]]{2}'
  
  des_contents <- readLines(desFile)
  
  # return radiances (from description file)
  as.numeric(regmatches(des_contents, regexpr(radiance_pattern, des_contents)))
}
