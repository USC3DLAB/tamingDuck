readNRELData <- function(dataSet = NULL) {
  # Reads the NREL dataset for wind, solar and load and puts it into a standard dataSet  
  
  # Make sure that the dataSet directory exists
  if ( !dir.exists("../dataSets") )
    dir.create("../dataSets")
  
  # Setup the full path for the files and read them.
  inputDir = sprintf('../dataSets/nrelData/');
  
  # Read the data from all the files in the folder.
  dataType <- NULL;
  if ( is.null(dataSet) ) {
    fileNames = list.files(path = inputDir);
    dataSet <- array(data = list(NULL)); i = 1;
    for (l in 1:length(fileNames)) {
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
        temp <- read.csv(file = paste(inputDir, fileNames[l], sep = ""), header = TRUE);
        dataSet[[i]] <- tempList
        if ( dim(temp)[1] == 0 )
          stop('No data availale for: ', fileNames[l])
        dataSet[[i]]$timeStamp <- temp[,1]
        dataSet[[i]]$ts[[dataType[[j]]]] <- temp[,2]
        dataSet[[i]]$freq[[dataType[j]]] <- tempFreq
      }
      else {
        temp <- read.csv(file = paste(inputDir, fileNames[l], sep = ""), header = TRUE);
        k = 1;
        while ( k <= i ) {
          if ( dataSet[[k]]$latitude == tempList$latitude && dataSet[[k]]$longitude == tempList$longitude && 
               dataSet[[k]]$capacity == tempList$capacity) {
            break;
          }
          k <- k+1 
        }
        if ( k < i ) {
          temp <- read.csv(file = paste(inputDir, fileNames[l], sep = ""), header = TRUE);
          if ( dim(temp)[1] == 0 )
            stop('No data availale for: ', fileNames[l])
          dataSet[[k]]$ts[[dataType[j]]] <- temp[,2]
          dataSet[[k]]$freq[[dataType[j]]] <- tempFreq
        }
        else {
          stop("Failed to identify the location.")
        }
      }
      i = i+1 
    }
    print("Decomposed the data into useful components") 
  }
  
  return(dataSet)
}