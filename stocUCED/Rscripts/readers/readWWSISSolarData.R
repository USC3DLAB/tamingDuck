readWWSISSolarData <- function(state, stateCode, dataset = NULL, plotmap = FALSE) {
  require(ggmap)
  
  # Make sure that the dataset directory exists
  if ( !dir.exists("../dataset") )
    dir.create("../dataset")
  
  # Download the file, extract and delete the zip file
  inputDir = paste('../dataset/',stateCode,'-pv-2006/', sep = '')
  if ( !dir.exists(inputDir) ) {
    downloadLink = paste('https://www.nrel.gov/grid/assets/downloads/',stateCode,'-pv-2006.zip', sep = '');
    download.file(url = downloadLink, destfile = 'temp.zip', quiet = TRUE);
    unzip(zipfile = 'temp.zip', exdir = inputDir);
    file.remove('temp.zip');
    print("Downloaded and extracted the data into ../dataset folder.") }
  else
    print("Data already exists in ../dataset folder.")
  
  # Read the data from all the files in the folder.
  dataType <- NULL;
  if ( is.null(dataset) ) {
    fileNames = list.files(path = inputDir);
    dataset <- array(data = list(NULL)); i = 1;
    for (l in 1:length(fileNames)) {
      if ( fileNames[l] != "README" ) {
        temp <- unlist(strsplit(x = fileNames[l], split = "_"))
        tempName <- paste(temp[1],"_freq",sep = '')
        
        if ( is.na(as.double(temp[2])) ) 
          tempList <- list(latitude = temp[2])
        else 
          tempList <- list(latitude = as.double(temp[2]))
        
        if ( is.na(as.double(temp[3])) )
          tempList <- c(tempList, longitude = temp[3])
        else
          tempList <- c(tempList, longitude = as.double(temp[3]))
        
        if ( is.na(as.integer(temp[4])) )
          tempList <- c(tempList, year = temp[4])
        else ( is.na(as.integer(temp[4])) )
        tempList <- c(tempList, year = as.integer(temp[4]))
        
        tempList <- c(tempList, srcType = temp[5], capacity = as.numeric(gsub("\\D", "", temp[6])))
        
        if ( is.na(as.double(temp[7])) )
          tempFreq <- temp[7]
        else
          tempFreq <- 60/as.double(temp[7])
        
        # Identify the data type of the current file
        j = 1;
        while ( j <= length(dataType)) {
          if ( dataType[j] == temp[1] ){
            break;
          }
          j <- j+1;
        }
        if ( j > length(dataType) ) {
          dataType[j] = temp[1] 
        }
        
        # Check to see if there is an entry for the location 
        if ( j == 1 ) {
          dataset[[i]] <- tempList
          temp <- read.csv(file = paste(inputDir, fileNames[l], sep = ""), header = TRUE);
          if ( dim(temp)[1] == 0 )
            stop('No data availale for: ', fileNames[l])
          dataset[[i]]$timeStamp <- temp[,1]
          dataset[[i]]$ts[[dataType[[j]]]] <- temp[,2] 
          dataset[[i]]$freq[[dataType[j]]] <- tempFreq
        }
        else {
          k = 1;
          while ( k <= i ) {
            if ( dataset[[k]]$latitude == tempList$latitude && dataset[[k]]$longitude == tempList$longitude && 
                 dataset[[k]]$capacity == tempList$capacity) {
              break; 
            }
            k <- k+1 
          }
          if ( k < i ) {
            temp <- read.csv(file = paste(inputDir, fileNames[l], sep = ""), header = TRUE);
            if ( dim(temp)[1] == 0 )
              stop('No data availale for: ', fileNames[l])
            dataset[[k]]$ts[[dataType[[j]]]] <- temp[,2] 
            dataset[[k]]$freq[[dataType[j]]] <- tempFreq
          }
          else {
            stop("Failed to identify the location.")
          }
        }
        i = i+1 
      } 
    }
    print("Decomposed the data into useful components") 
  }
  
  if ( plotmap ) {
    # Plot the data on the maps
    df <- NULL;
    for ( l in 1:length(dataset)) {
      df = abind::abind(df, c(dataset[[l]]$longitude, dataset[[l]]$latitude), along = 2)
    }
    df <- data.frame(lon = df[1,], lat = df[2,])
    df <- unique(df)
    
    # Download the maps and plot the locations
    png(filename = paste(getwd(), "/", state,"_pv-2006", sep = ""))
    map.background <- get_googlemap(state, zoom = 6, key = "AIzaSyDi-HXKpj26zQonCSyn91geb59nj0vwHUU")
    ggmap(map.background, extent = 'panel') + geom_point(data = df, aes(x = lon, y = lat), color = "orange")
    dev.off()
  }
  
  return(dataset)
}