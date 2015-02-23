args = commandArgs( trailingOnly = TRUE )

#INPUT_DIR <- args[ 1 ]
#EXP_FILE <- args[ 2 ]
INPUT_DIR <- "data"

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
      f[i] <- ifelse( !is.na( proper ),paste0(INPUT_DIR,SEP,f[i]),NA )
      
    }

  return( f ) 
}

my_list <- function( l ){ vector("list",l) }

read_files <- function(filenames, raw){
  table_list <- my_list( length( filenames ) ) 
  alpha <- 300 # 5 minutes
  current_time <- 1410742840
  measurements <- data.table(EPOCH_TIME=numeric(0),.I=numeric(0))
  
  for ( i in 1:length( filenames ) ) {
    d <- as.data.table( read.table( filenames[i],header=T ) )
    setkey(d,EPOCH_TIME)
    times <- seq.int(current_time,max(d$EPOCH_TIME),by=alpha)
    indices <- d[J(times),.I,roll="nearest"]
    measurements <- rbindlist(list(measurements,indices))
    print(measurements)
    current_time <- as.numeric(tail(measurements$E,n=1))+alpha # first measurement of the next file
    if (i == 1){
      plot(d$EPOCH_TIME,d$CO2,type="l",xlim=c(1410742892,1410742892+80000),ylim=c(100,750))
    }
    else {
      lines(d$EPOCH_TIME,d$CO2)
    }
    points(d$EPOCH_TIME[indices$.I],d$CO2[indices$.I],col="red")
    points(d$EPOCH_TIME[indices$.I+200],d$CO2[indices$.I+200],col="blue")
  }
  #return(measurements)
  
}
loadlibs("data.table")
f <- get_filenames()
read_files(f,data.table())

