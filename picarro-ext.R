args = commandArgs( trailingOnly = TRUE )

#INPUT_DIR <- args[ 1 ]
#EXP_FILE <- args[ 2 ]
INPUT_DIR <- "data"
EXP_FILE <- "exp.csv"

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
  
  printlog( "Loading libraries..." )
  loadedlibs <- vector()
  for( lib in liblist ) {
    
    printlog( "Loading", lib )
    loadedlibs[ lib ] <- require( lib, character.only=T, warn.conflicts=F )
    if( !loadedlibs[ lib ] )
      warning( "this package is not installed!" )
    
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
      start <- d.experiment_indices$.I[j]
      end <- ifelse(start + 200 > n_data,n_data,start+200)
      exp.temperature <- exp[measurements_so_far,Temperature]
      
      if ((end-start) <= 50){
        break
      }
      start_time <- d.experiment_indices$E[j]
      alldata <- rbindlist(list(alldata,
                                calc(apply(d[start:end,2:7,with=F],
                                                      2,
                                                      function(x){
                                                        summary(lm( x ~ seq(0,
                                                                            end-start)) ) 
                                                        }),
                                                exp.temperature)))
      
      ## TODO: 1. Split d into experiments. See LGR.
      ##       2. Find slope for each gas in this table.
      ##       3. Run QC on slopes.
      ##       4. Apply flux() with slope and temp.
      ##       5. Concat with timestamp into single table.
      ##       6. rbindlist() with current alldata.
      
    }
  
  }
  
}

quality_control <- function(resp,temperature){
  return(lapply(resp,function(x){
                                    c(unlist(round( x$r.squared, 2 )),
                                      unlist(coef(x)[2,4])) }))
}

calc <- function( resp,temperature ){
  # Perform QC for this measurement. 
  #print(resp)
  names_r2 <- c( "N2O_dry_r2","N2O_dry30s_r2","CO2_dry_r2","CH4_dry_r2","H2O_r2","NH3_r2" )
  names_p <- c( "N2O_dry_p","N2O_dry30s_p","CO2_dry_p","CH4_dry_p","H2O_p","NH3_p" )
  names <- c( names_r2,names_p )
  
  qc <- data.table(quality_control(resp))
  #setnames( qc,c( "Measurement_Time",paste0( "V",as.character(1:12) ) ),c( "Measurement_Time",names ) )
  print(colnames(qc))
  
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
  chamber_v <- 10 * S
  V <- chamber_v + PICARRO_V
  Pa <- 101 #kPA
  R <- 8.3145e-3 # m-3 kPa mol-1 K-1
  Kelvin <- 273.15 # C to K conversion
  Temp <- temp + Kelvin
  
  flux <- resp * ( V / S ) * Pa / ( R * Temp )
  
  return( flux )
  
}


# main ##########################################################################


loadlibs( "data.table" )
exp <- as.data.table(read.csv(EXP_FILE,header=T))
f <- get_filenames()
alldata <- data.table("Measurement_Time"=numeric(0),
                      "N2O_dry_flux [umol m-2 s-1]"=numeric(0),
                      "N2O_dry30s_flux [umol m-2 s-1]"=numeric(0),
                      "CO2_dry_flux [umol m-2 s-1]"=numeric(0),
                      "CH4_dry_flux [nmol m-2 s-1]"=numeric(0),
                      "H2O_flux [umol m-2 s-1]"=numeric(0),
                      "NH3_flux [umol m-2 s-1]"=numeric(0),
                      "Temp [C]"=numeric(0),
                      "N2O_dry_r2"=numeric(0),
                      "N2O_dry30s_r2"=numeric(0),
                      "CO2_dry_r2"=numeric(0),
                      "CH4_dry_r2"=numeric(0),
                      "H2O_r2"=numeric(0),
                      "NH3_r2"=numeric(0),
                      "N2O_dry_p"=numeric(0),
                      "N2O_dry30s_p"=numeric(0),
                      "CO2_dry_p"=numeric(0),
                      "CH4_dry_p"=numeric(0),
                      "H2O_p"=numeric(0),
                      "NH3_p"=numeric(0))

alldata <- read_files( f,alldata,exp )


#if (i == 1){
#  plot(d$EPOCH_TIME,d$CO2,type="l",xlim=c(1410742892,1410742892+85000),ylim=c(350,850))
#}
#else {
#  lines(d$EPOCH_TIME,d$CO2)
#}
#points(d$EPOCH_TIME[d.experiment_indices$.I],d$CO2[d.experiment_indices$.I],col="red")
#points(d$EPOCH_TIME[d.experiment_indices$.I+200],d$CO2[d.experiment_indices$.I+200],col="blue")


