decomposeTSData <- function(inputTS, freq, decomposeBy = 'seasonal', identifier = 'summer') {
  # Input: 
  #   - inputTS     = A one-dimensioanl time series (assumed to be annual data)
  #   - freq        = Frequency of data
  #   - decomposeBy = 'seasonal', 'monthly', 'weekly', or 'daily'
  #   - identifier  = ('summer', 'autumn', 'winter', 'spring') or name of the month
  # Output:
  #   - Decomposed time series in the form of a matrix where each row corresponds to one days data. 
  
  Nhour <- freq
  
  if ( decomposeBy == 'seasonal' ) {
    # Session windows
    winterDays <- c(as.Date('06/12/21'), as.Date('06/03/20'))
    springDays <- c(as.Date('06/03/21'), as.Date('06/06/20'))
    summerDays <- c(as.Date('06/06/21'), as.Date('06/09/22'))
    autumnDays <- c(as.Date('06/09/23'), as.Date('06/12/20'))
    
    if ( identifier == 'winter') {
      iStart <- (as.numeric(strftime(x = as.POSIXlt(paste(winterDays[1], sep = '/')), format = '%j'))-1)*Nhour*24 + 1
      iEnd <- (as.numeric(strftime(x = as.POSIXlt(paste(winterDays[2], sep = '/')), format = '%j')))*Nhour*24
      idx = c(iStart:length(inputTS), 1:iEnd)
    }
    else if ( identifier == 'spring') {
      iStart <- (as.numeric(strftime(x = as.POSIXlt(paste(springDays[1], sep = '/')), format = '%j'))-1)*Nhour*24 + 1
      iEnd <- (as.numeric(strftime(x = as.POSIXlt(paste(springDays[2], sep = '/')), format = '%j')))*Nhour*24
      idx = iStart:iEnd
    }
    else if ( identifier == 'summer') {
      iStart <- (as.numeric(strftime(x = as.POSIXlt(paste(summerDays[1], sep = '/')), format = '%j'))-1)*Nhour*24 + 1
      iEnd <- (as.numeric(strftime(x = as.POSIXlt(paste(summerDays[2], sep = '/')), format = '%j')))*Nhour*24
      idx = iStart:iEnd
    }
    else if ( identifier == 'autumn') {
      iStart <- (as.numeric(strftime(x = as.POSIXlt(paste(autumnDays[1], sep = '/')), format = '%j'))-1)*Nhour*24 + 1
      iEnd <- (as.numeric(strftime(x = as.POSIXlt(paste(autumnDays[2], sep = '/')), format = '%j')))*Nhour*24
      idx = iStart:iEnd
    }
    else {
      stop('Unknown identifier,')
    }
    tsData <- inputTS[idx]
    tsData <- t(matrix(data = tsData, nrow = 24*Nhour))
  }
  else if ( decomposeBy == 'monthly' ) {
    # Monthly data
    if ( is.null(m <- grep(pattern = identifier, x = month.name)) ) 
      stop('Could not find the name of the month.')
    iStart <- (as.numeric(strftime(x = as.POSIXlt(paste('06',m,'01', sep = '/')), format = '%j'))-1)*Nhour*24+1
    if ( m != 12 ) {
      iEnd   <- (as.numeric(strftime(x = as.POSIXlt(paste('06',m+1,'01', sep = '/')), format = '%j'))-1)*Nhour*24
    }else {
      iEnd <- length(inputTS)
    }
    tsData <- inputTS[iStart:iEnd]
    tsData <- t(matrix(data = tsData, nrow = 24*Nhour))
  }
  else if ( decomposeBy == 'daily') {
    # Daily data
    tsData <- t(matrix(data = inputTS, nrow = Nhour*24))
  }
  else if ( decomposeBy == 'weekly' ) {
    # Weekly data
    tsData <- t(matrix(data = inputTS, nrow = Nhour*24*7))
  }
  else {
    stop('Unknown decomposeBy.')
  }
  
  return(tsData)
}