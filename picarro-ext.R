args = commandArgs( trailingOnly = TRUE )

#INPUT_DIR <- args[ 1 ]
#EXP_FILE <- args[ 2 ]
INPUT_DIR <- "data"
EXP_FILE <- "exp.csv"
OUTPUT_DIR <- "out"

CHAMBER_RAD <- 9.9 # cm
Pa <- 101 #kPA
R <- 8.3145e-3 # m-3 kPa mol-1 K-1
Kelvin <- 273.15 # C to K conversion
PICARRO_V <- 482 # cm3 or sccm

DATE_FORMAT3 <- "%Y%m%d%H%M%S"
FILE_FORMAT <- "JFAADS2012-20[0-9]{6}-\\d+-DataLog_User"
SEP <- ifelse(Sys.info()['sysname'] != "Windows","/","\\")

alpha <- 2
gamma <- 3
first_measurement <- 1410742840

printlog <- function( msg="", ..., ts=TRUE, cr=TRUE, pre_cr=FALSE ) {
  
  if ( pre_cr ) cat( "\n" )
  if( ts ) cat( date(), " " )
  cat( msg, ... )
  if( cr ) cat( "\n")
  
} # printlog

# -----------------------------------------------------------------------------
# loadlibs
# Load a list of libraries
#     Arguments: liblist -- A character vector of library names.
# 
loadlibs <- function( liblist ) {
  
  #printlog( "Loading libraries..." )
  loadedlibs <- vector()
  not_installed <- vector()
  
  for( lib in liblist ) {
    
    #printlog( "Loading", lib )
    loadedlibs[ lib ] <- require( lib, character.only=TRUE, warn.conflicts=FALSE )
    if( !loadedlibs[ lib ] ){
      
      warning( "this package is not installed!" )
      not_installed <- c(not_installed,lib)
      
    }
  }
  
  if ( length(not_installed) != 0){
    
    #print("Installing packages.")
    chooseCRANmirror(graphics=FALSE,ind=85) # Choose USA (MD)
    install.packages(not_installed)
    
    for( lib in not_installed ) {
      
      #printlog( "Loading", lib )
      loadedlibs[ lib ] <- require( lib, character.only=TRUE, warn.conflicts=FALSE )
      
    }	
  }
  invisible( loadedlibs )
  
} # loadlibs
proper_file <- function( file ){ 
  
  file.split <- strsplit( file,"-" )
  datetime <- paste0( file.split[[1]][2],file.split[[1]][3] )
  
  str <- strptime( datetime, format=DATE_FORMAT3 )
  return ( as.POSIXct(str) )
  
}

get_filenames <- function(){   
    
    f <- list.files( INPUT_DIR,FILE_FORMAT )
    
    for ( i in seq( 1:length( f ) ) ){
      
      proper <- proper_file( f[i] ) 
      f[i] <- ifelse( !is.na( proper ),paste0( INPUT_DIR,SEP,f[i] ),NA )
      
    }

  return( f ) 
}

my_list <- function( l ){ vector("list",l) }

read_files <- function(filenames, alldata, exp){
  table_list <- my_list( length( filenames ) ) 
  alpha <- 300 # 300s : 5 minutes
  current_time <- 1410742840
  measurements <- data.table( EPOCH_TIME=numeric( 0 ),.I=numeric( 0 ) )
  alldata <- data.table()
  measurements_so_far <- 0
  pb <- txtProgressBar(min = 0, max = nrow(exp), style = 3)
  
  for ( i in 1:length( filenames ) ) {
    
    d <- as.data.table( read.table( filenames[i],header=T ) )
    d <- d[,c( "EPOCH_TIME","N2O_dry","N2O_dry30s","CO2_dry","CH4_dry","H2O","NH3" ),with=F]
    n_data <- nrow(d)
    
    setkey( d,EPOCH_TIME )
    d.experiment_indices <- d[J( seq.int( current_time,max( d$EPOCH_TIME ),by=alpha ) ),.I,roll="nearest"]
    n_experiments <- nrow(d.experiment_indices)
    
    current_time <- as.numeric(tail(d.experiment_indices$E,n=1))+alpha # first measurement of the next file
    for ( j in 1:n_experiments ){
      measurements_so_far = measurements_so_far + 1
      start <- d.experiment_indices$I[j]
      end <- ifelse(start + 200 > n_data,n_data,start+200)
      exp.temperature <- exp[measurements_so_far,Temperature]
      if ((end-start) <= 50){
        break
      }
      
      alldata <- rbindlist(list(alldata,calc(d[start:end,2:7,with=F],exp.temperature,end-start,d.experiment_indices$E[j])))
      setTxtProgressBar(pb, measurements_so_far)
    }
  
  }
  return(alldata)
  
}

get_r2 <- function(lm){
  return(round(summary(lm)$r.squared,2))
}

get_p <- function(lm){
  return(summary(lm)$coefficients[2,4])
}

get_slope <- function(lm){
  return(lm$coefficients[[2]])
}

flux <- function( resp, temp ){
  
  # We want to convert raw respiration (d[CO2]/dt) to a flux using
  # A = dC/dt * V/S * Pa/RT (e.g. Steduto et al. 2002), where
  # A is CO2 flux (umol/m2/s)
  # dC/dt is raw respiration as above (mole fraction/s)
  # V is total chamber volume (m3)
  # ...we are correcting for varying headspaces in the cores
  # S is ground surface area (m2)
  # ...but we're computing per kg of soil, so using dry mass instead
  # Pa is atmospheric pressure (kPa)
  # R is universal gas constant (8.3 x 10-3 m-3 kPa mol-1 K-1)
  # T is air temperature (K)
  
  S <- pi * CHAMBER_RAD ^ 2
  depth <- 10
  chamber_v <- depth * S
  V <- chamber_v + PICARRO_V
  Temp <- temp + Kelvin
  
  flux <- resp * ( V / S ) * Pa / ( R * Temp )
  
  return( flux )
  
}

calc <- function( d, temp, n_data, start_time){
  time <- seq(0,n_data)
  d.respirations <- d[,list(lms = lapply(.SD,function(gas)lm(gas ~ time)))]
  d.data <- d.respirations[,`:=` (lms = sapply(lms,get_slope),
                                  r2 = sapply(lms,get_r2),
                                  p = sapply(lms,get_p),
                                  gas = seq(1,6))]
  setkey(d.data,gas)
  
  d.data <- as.matrix(d.data[,lms := flux(lms,temp)])
  
  row <- data.table("Measurement Time"=as.POSIXct(start_time, origin="1970-01-01"))
  row <- row[,`:=`("Temp [C]"=temp,
                   "N2O_dry_flux [umol m-2 s-1]"=d.data[1,1],
                   "N2O_dry30s_flux [umol m-2 s-1]"=d.data[2,1],
                   "CO2_dry_flux [umol m-2 s-1]"=d.data[3,1],
                   "CH4_dry_flux [mol m-2 s-1]"=d.data[4,1]*1.0e3,
                   "H2O_flux [umol m-2 s-1]"=d.data[5,1],
                   "NH3_flux [umol m-2 s-1]"=d.data[6,1],
                   "N2O_dry_r2"=d.data[1,2],
                   "N2O_dry30s_r2"=d.data[2,2],
                   "CO2_dry_r2"=d.data[3,2],
                   "CH4_dry_r2"=d.data[4,2],
                   "H2O_r2"=d.data[5,2],
                   "NH3_r2"=d.data[6,2],
                   "N2O_dry_p"=d.data[1,3],
                   "N2O_dry30s_p"=d.data[2,3],
                   "CO2_dry_p"=d.data[3,3],
                   "CH4_dry_p"=d.data[4,3],
                   "H2O_p"=d.data[5,3],
                   "NH3_p"=d.data[6,3])]
  return(row)
}

savedata <- function( d, extension=".csv" ) {
  
  my.write <- function(x, file, header, f = write.csv, ...){
    datafile <- file(file, open='wt')
    on.exit(close(datafile))
    if (!missing(header)) writeLines(header,con=datafile)
    return( f(x, datafile, ...) )
  }
  
  stopifnot( file.exists( OUTPUT_DIR ) )
  fn <- paste0( OUTPUT_DIR, SEP, format( Sys.time(),"%d%B%Y_%H%M%S" ),"_fluxes",extension )
  printlog( "Saving", fn )
  my.write( d, fn, f = write.csv, row.names=FALSE )
  
} # savedata


# main ##########################################################################

stopifnot( file.exists( INPUT_DIR ) )
stopifnot( file.exists( EXP_FILE ) ) # make sure the folders and EXP_FILE exist 

if( !file.exists( OUTPUT_DIR ) ) { # create output dir if it does not exist
  
  printlog( "Creating", OUTPUT_DIR )
  dir.create( OUTPUT_DIR )
  
}

loadlibs( "data.table" )
options(datatable.old.bywithoutby=TRUE)

printlog("Reading metadata...")
exp <- as.data.table(read.csv(EXP_FILE,header=T))
printlog("Reading files...")
f <- get_filenames()
alldata <- data.table()
printlog("Calculating fluxes, quality control...")
alldata <- read_files( f,alldata,exp )

printlog("Saving...")
savedata(alldata)

#if (i == 1){
#  plot(d$EPOCH_TIME,d$CO2,type="l",xlim=c(1410742892,1410742892+85000),ylim=c(350,850))
#}
#else {
#  lines(d$EPOCH_TIME,d$CO2)
#}
#points(d$EPOCH_TIME[d.experiment_indices$I],d$CO2[d.experiment_indices$I],col="red")
#points(d$EPOCH_TIME[d.experiment_indices$I+200],d$CO2[d.experiment_indices$I+200],col="blue")


