#' map.ROI
#' @description Function for mapping all ROI names in a series a studies
#' @param path A \code{character} value containing the path to browse for DICOM studies
#' @details DICOM RT structure set files often contain multiple \emph{ROI names} referencing to the same structure. \code{moddicom} is desinged to analyze only
#' one structure at a time. Using this function a \code{ROImap} class object is created, listing all ROI names for each case inside the borwsed directory.
#' Using the method \code{print.ROImap} a table showing all ROI names and their frequency in the explored path is shown.
#' @return An object of class \code{ROImap} that is a list containing all mapped ROIs in a series of studies
#' @export
#' @examples ## NOT RUN
#' ListOfROI <- map.ROI(path = '/my_path')
#' ## to get the table of ROIs in the explored path
#' print(ListOfROI)
map.ROI <- function(path) {
  patList <- list.dirs(path = path,recursive = FALSE)
  ROI.map <- list()
  n <- 0
  for(patID in patList) {
    n <- n + 1
    cat("\n ---------------------------------\n Computing Patient: ", patID, "\n")
    ogg.geoLet <- geoLet()
    ogg.geoLet$openDICOMFolder(pathToOpen = patID)  
    ROI.map[[n]] <- GLT.getROIList(obj.geoLet = ogg.geoLet)
    names(ROI.map)[[n]] <- patID
    rm(ogg.geoLet)
    gc()
  }
  attr(ROI.map, "class") <- "ROImap"
  return(ROI.map)
}

#' print.ROImap
#' @export print.ROImap
print.ROImap <- function(obj) {
  cat('\nNumber of patients:', length(x = obj))
  ROI.Table <- c()
  # browser()
  for (N in obj) {
    ROI.Table <- c(ROI.Table, N[2, ])
  }
  cat('\n')
  print(table(ROI.Table))
}

#' summary.ROImap
#' @export summary.ROImap
summary.ROImap <- function(obj) {
  cat('\nPatients and ROIs\n')
  cat('\n')
  y <- as.data.frame(x = table(names(unlist(unname(lapply(X = obj, FUN = function(x) return(table(length(x[2,])))))))))
  colnames(y) <- c('Number of ROIs', 'Number of Patients')
  print(y, row.names = FALSE)
  cat('\n\nPatients details:\n\n')
  print(lapply(X = obj, FUN = function(x) return(x[2,])))
  cat('\n')
}
