#' Setup geoLet derived objects for deep learning analysis
#' 
#' @param path The \code{path} of cases structured as \pkg{moddicom} usually requires (one subfolder with a DICOM study containing images and a RT structure DICOM file)
#' @param outcome A \code{data.frame} containing in one column the same IDs of patients listed in \code{path}, usually corresponding to the folder names, and a column with the given outcome
#' @param ratio a number between 0 and 1 giving the ratio between training set and verification set (default 2/3)
#' @param ROIname a \code{character} vector containing the ROI(s) name(s) used for extracting images from studies
#' @param ROINameNorm a \code{character} vector containing the ROI(s) name(s) used for used for normalizing the values, if \code{NULL} no normalization is applied
#' @param TrainTest \code{logical} value (default is \code{TRUE}) that gives the outcome already splitted for \emph{training} and \emph{test}. 
#' @param per_patient_sampling a \code{logical} value (default is \code{TRUE}) that sets the sampling for training and testing set according the patients IDs. If \code{FALSE} the sampling is obtained by images slices.
#' @param New_Pixel_Spacing a \code{numeric} vector, length must be 2 values, that sets the new pixel spacing for all the images in order to homogenize the spatial resolution among all patients.
#' @param ROIBitMask a \code{logical} value, default \code{TRUE}, return only the voxels within the ROI, if \code{FALSE} returns the whole VoxelCube.
#' @param TargetSize default value is \code{'auto'} that calculates the size in rows and columns of final array as mean of the single VoxelCubes sizes, otherwise user can set a numerical vector of length = 2, that sets number of rows and columns manually.
#' @param VirtualBiopsy \code{logical} value, default is \code{FALSE}. If \code{TRUE} returns for each patient returns a number maximum number of \code{VirtualBiopsyMaxNum} of small sampling images with size \code{VirtualBiopsyRay * 2 + 1}.
#' @param VirtualBiopsyRay a \code{numeric} vector of length 2, giving the width os extension of \emph{Virtual Biopsy} around the centroid. Final dimension of biopsy is \code{VirtualBiopsyRay * 2 + 1}.
#' @param VirtualBiopsyMaxNum a \code{integer} value setting the maqximum number of \emph{Virtual Biopsies} to be sampled.
#' @param backgroud the value of the pixel surrounding the ROI in the background: \code{'zero'} is 0; \code{'min'} is the minimum value in the ROI; \code{'mean'} is the mean value in the ROI; \code{'runif'} is a random uniform distribution between \code{min} and \code{max} values in the ROI; \code{'rnorm'} is a random normal distrbution with \code{mean} is the mean value of ROI, \code{sd} is the standard distribution of ROI values.
#' @param bckwidth a \code{numeric} value between 0 and 1: the width in the range of ROI values for \code{runif} background option or the multiplier of \code{sd} for \code{rnorm} background option.
#' @param threshold_filter The threshold number of pixels per slice \code{NA} below which images are deleted. Default is 150 pixels.
#' @param ROImap If available, the result of a \link{map.ROI} function, as an object of class ROImap given by \pkg{moddicom} package.
#' @description Provide the output needed for starting a deep learning model using \code{\link{keras}}:
#' @return 
#' \itemize{
#' \item{\code{ROIList}}{: The list of ROIs as achieved by \pkg{moddicom} \code{geoLet$getROIVoxels( )} function.}
#' \item{\code{cases_data.frame}}{: A \code{data.frame} containing the patients IDs (\code{pt_name}), each \code{ROI_name} mapped in each patient and the \code{pt_outcome} that reports the patient outcome.}
#' \item{\code{images_array}}{: The final \code{array} of selected images structured  \link{keras} requires.}
#' \item{\code{keras_outcome}}{: A \code{numeric} vector as required by \link{keras} for the classification of outcome}
#' \item{\code{df_outcome}}{: A \code{data.frame} containing the labels of the outcome as provided by user and the corresponding labels (numeric values) used for modeling in \link{keras}.}
#' \item{\code{cases_data.frame}}{: A \code{data.frame} containing the patients IDs (\code{pt_name}), each \code{ROI_name} mapped in each patient and the \code{pt_outcome} that reports the patient outcome.}
#' \item{\code{images_array_train}}{: The final \code{array} of selected images structured as \link{keras} requires. This is the \emph{training} dataset. Available if \code{TrainTest=TRUE}.}
#' \item{\code{images_array_test}}{: The final \code{array} of selected images structured as \link{keras} requires. This is the \emph{test} dataset. Available if \code{TrainTest=TRUE}.}
#' \item{\code{keras_outcome_train}}{: A \code{numeric} vector as required by \link{keras} for the classification of outcome. This is the outcome for the \emph{training} dataset. Available if \code{TrainTest=TRUE}.}
#' \item{\code{keras_outcome_test}}{: A \code{numeric} vector as required by \link{keras} for the classification of outcome. This is the outcome for the \emph{test} dataset. Available if \code{TrainTest=TRUE}.}
#' \item{\code{df_outcome}}{: A \code{data.frame} containing the labels of the outcome as provided by user and the corresponding labels (numeric values) used for modeling in \link{keras}.}
#' \item{\code{total_data.frame_train}}{: The \code{data.frame} with all cases and slices corresponding to the images in the \code{images_array_train} array. Available if \code{TrainTest=TRUE}.}
#' \item{\code{total_data.frame_test}}{: The \code{data.frame} with all cases and slices corresponding to the images in the \code{images_array_test} array. Available if \code{TrainTest=TRUE}.}
#' }
#' @export
#' @import keras
#' @import fields
Setup_For_Keras <- function(path, outcome = NULL, ratio = 2/3, ROIname, ROINameNorm = NULL, TrainTest = TRUE, per_patient_sampling = TRUE, New_Pixel_Spacing = c(),
                            ROIBitMask = TRUE, TargetSize = 'auto', VirtualBiopsy = FALSE, VirtualBiopsyRay = c(10, 10), VirtualBiopsyMaxNum = 25,
                            backgroud = c('zero', 'min', 'mean', 'runif', 'rnorm'), bckwidth = .5, threshold_filter = 150, ROImap = NULL){
  # check the ratio
  if ((ratio < 0) || (ratio > 1)) stop('\nratio must be > 0 AND < 1')
  # check the backgroud
  if ((bckwidth < 0) || (bckwidth > 1)) stop('\nbckwidth must be > 0 AND < 1')
  # match the choice
  backgroud <- match.arg(backgroud)
  # list of ROIs to be analyzed
  ROIList <- list()    # list of ROIs
  normList <- list()   # list of normalization ROIs
  pt_name <- c()       # output vector for data.frame of patients-structures: patient name
  ROI_name <- c()      # output vector for data.frame of structures-patients: ROI name
  pt_outcome <- c()    # output vector for data.frame of structures-patients: patient outcome
  final_pt_name <- c() # final vector containing patients name, one record for each slice
  final_ROI_name <- c()# final vector containing ROI names, one record for each slice
  final_outcome <- c() # final outcome repeating each pt_outcome for each image
  baseline <- c()      # baseline for zero level
  # counter of patients
  n <- 1
  # paths of patients images
  pathLIST <- list.dirs(path = path, recursive = F)
  if (is.null(ROImap)) ROIs<-map.ROI(path = path) else ROIs <- ROImap
  
  # function for normalizing the content of an array without having a reference normalization ROI
  normalize.biopsy <- function(inputROIList, normROIList = NULL) {
    outputROIList <- inputROIList
    if (is.null(normROIList)) {  # internal normalization from min to max with ROINameNorm = NULL
      minROIList <- min(unlist(lapply(X = inputROIList, FUN = function(x) return(min(x$cube)))))
      maxROIList <- max(unlist(lapply(X = inputROIList, FUN = function(x) return(max(x$cube)))))
      widthROIList <- maxROIList - minROIList
      for (z in 1:length(inputROIList)) outputROIList[[z]]$cube <- (inputROIList[[z]]$cube - minROIList) / widthROIList
    }
    else { # normalization to a ROINameNorm structure
      for (z in 1:length(inputROIList)) outputROIList[[z]]$cube <- inputROIList[[z]]$cube / mean(normROIList$masked.images$voxelCube, na.rm = TRUE)
      minROIList <- min(unlist(lapply(X = outputROIList, FUN = function(x) return(min(x$cube)))))
      maxROIList <- max(unlist(lapply(X = outputROIList, FUN = function(x) return(max(x$cube)))))
      widthROIList <- maxROIList - minROIList
      # normalize 0 to 1
      for (z in 1:length(outputROIList)) outputROIList[[z]]$cube <- (outputROIList[[z]]$cube - minROIList) / widthROIList
    }
    return(outputROIList)
  }
  
  ###########################################
  ###    function with virtual biopsies   ###
  ###########################################
  # loop inside directories for building image series
  # ROIname:  the name given by the user to extract voxel cube
  # ROIs:     the mapROI result, a list where each element is a patient
  # ROI:      the ROI being mapped now
  if (VirtualBiopsy == TRUE) {
    if (length(VirtualBiopsyRay) != 2) stop('\nVirtualBiopsyRay MUST be an integer vector of length 2.\n')
    for (patients_folder in pathLIST) { # cicle through patients
      if (any(ROIname %in% ROIs[[n]]))  # check if there is a ROIname chosen by user in the actual patient
      { objT <- geoLet()
        objT$openDICOMFolder(pathToOpen = patients_folder)
        for (ROI in ROIname) 
          if (ROI %in% ROIs[[n]]) 
          {
            cat('\n\nSampling biopsies of', ROI, 'in patient', patients_folder)
            m <- length(ROIList)
            # temporary virtual biopsy image
            ROIList[[m + 1]] <- objT$getVirtualBiopsy(ROIName = ROI, BiopsyCubeDim = c(VirtualBiopsyRay[1], VirtualBiopsyRay[2], 0), 
                                                      numOfSamples = VirtualBiopsyMaxNum, new.pixelSpacing = New_Pixel_Spacing)
            cat('\nEnd biopsies sampling')
            if (!is.null(ROINameNorm)) normList[[m + 1]] <- objT$getROIVoxels(Structure = ROINameNorm, new.pixelSpacing = New_Pixel_Spacing) # get the structure for normalization
            pt_name <- c(pt_name, patients_folder)
            ROI_name <- c(ROI_name, ROI)
            if (!is.null(outcome)) 
              pt_outcome <- c(pt_outcome, outcome[n])
          }
        rm(objT)
        gc()
      } else { 
        ROIList[[n]] <- NULL
        normList[[n]] <- NULL
      }
      n <- n + 1
    }
    # check for not using normalization in contours
    if (length(normList) == 0) normList <- NULL
    # code after cicling among virtual biopsies
    # set the final array size
    imgArr <- array(data = NA, dim = c(sum(sapply(X = ROIList, FUN = length)) , VirtualBiopsyRay[1] *2 + 1, VirtualBiopsyRay[2] * 2 + 1, 1))
    a <- 0
    cat('\n\nSetting image array\n')
    bp_n <- 0 # number of biopsies (and number of images in final array)
    for (n in 1:length(ROIList)) {
      temp_depth <- length(ROIList[[n]])
      bp_n <- bp_n + temp_depth 
      if (temp_depth > 0) { # for including a new series of biospies
        cat('\nAdding patient: ', pt_name[n], ';   ROI: ', ROI_name[n], '\n', sep = '')
        ROIList[[n]] <- normalize.biopsy(inputROIList = ROIList[[n]], normROIList = normList[[n]])
        if (!is.null(outcome)) 
          final_outcome <- c(final_outcome, rep.int(x = pt_outcome[n], times = temp_depth))
        final_pt_name   <- c(final_pt_name, rep.int(x = pt_name[n], times = temp_depth))
        final_ROI_name  <- c(final_ROI_name, rep.int(x = ROI_name[n], times = temp_depth))
        for (m in 1:temp_depth) imgArr[a+m,,,] <- ROIList[[n]][[m]]$cube # copy of the single biopsy    
        # start of copy of next patient
        a <- a + m
      } else cat('\nRemoving patient: ', pt_name[n], ';   ROI: ', ROI_name[n], '\n', sep = '')
    }
    
    if (!is.null(outcome)) 
      df_outcome <- as.data.frame(cbind('keras_outcome' = 1:length(levels(as.factor(final_outcome))), 'original_outcome' = levels(as.factor(final_outcome))))
    else {
      df_outcome <- NULL
      pt_outcome <- NULL
      final_outcome <- NULL
    }
    if (TrainTest == TRUE) {
      if (per_patient_sampling == TRUE) {
        # training set
        pt_selected <- sample(x = names(table(pt_name)), size = round(ratio * length(names(table(pt_name)))), replace = FALSE)   # sample the patients
        img_train <- which(final_pt_name %in% pt_selected)   # selecttion of the images index number
        final_pt_name_train <- final_pt_name[img_train]      # vector of patients names training set
        final_ROI_name_train <- final_ROI_name[img_train]    # vector of ROI names training set
        final_outcome_train <- final_outcome[img_train]      # vector of outcome training set
        img_arr_train <- imgArr[img_train,,,]                # final array of images training set     
        # testing set
        final_pt_name_test <- final_pt_name[-img_train]      # vector of patients names testing set
        final_ROI_name_test <- final_ROI_name[-img_train]    # vector of ROI names testing set
        final_outcome_test <- final_outcome[-img_train]      # vector of outcome testing set
        img_arr_test  <- imgArr[-img_train,,,]               # final array of images testing set
      } else {
        img_train <- sample(x = c(1:bp_n), size = round(ratio * bp_n), replace = FALSE) # sample the images
        final_pt_name_train <- final_pt_name[img_train]      # vector of patients names training set
        final_ROI_name_train <- final_ROI_name[img_train]    # vector of ROI names training set
        final_outcome_train <- final_outcome[img_train]      # vector of outcome training set
        img_arr_train <- imgArr[img_train,,,]                # final array of images training set     
        # testing set
        final_pt_name_test <- final_pt_name[-img_train]      # vector of patients names testing setA <- 
        final_ROI_name_test <- final_ROI_name[-img_train]    # vector of ROI names testing set
        final_outcome_test <- final_outcome[-img_train]      # vector of outcome testing set
        img_arr_test  <- imgArr[-img_train,,,]               # final array of images testing set
      }
      keras_outcome_train <- c()
      for(n in names(table(final_outcome_train))) keras_outcome_train <- cbind(keras_outcome_train, as.numeric(final_outcome_train == n))
      colnames(keras_outcome_train) <- names(table(final_outcome_train))
      keras_outcome_test <- c()
      for(n in names(table(final_outcome_test)))  keras_outcome_test  <- cbind(keras_outcome_test,  as.numeric(final_outcome_test  == n))
      colnames(keras_outcome_test)  <- names(table(final_outcome_test))
      if (!is.null(outcome)) {
        total_data.frame_train <- as.data.frame(cbind(final_pt_name_train, final_ROI_name_train, final_outcome_train))
        total_data.frame_test  <- as.data.frame(cbind(final_pt_name_test,  final_ROI_name_test,  final_outcome_test))
      }
      else {
        total_data.frame_train <- as.data.frame(cbind(final_pt_name_train, final_ROI_name_train))
        total_data.frame_test  <- as.data.frame(cbind(final_pt_name_test,  final_ROI_name_test ))
      }
      dim(img_arr_train) <- c(dim(img_arr_train) , 1)
      dim(img_arr_test)  <- c(dim(img_arr_test) , 1)
      return(list('ROIList' = ROIList, 'cases_data.frame' = as.data.frame(cbind(pt_name, ROI_name, pt_outcome)), 
                  'keras_images_array_train' = img_arr_train, 'keras_images_array_test' = img_arr_test, 
                  'keras_outcome_train' = keras_outcome_train, 'keras_outcome_test' = keras_outcome_test, 'df_outcome' = df_outcome, 
                  'total_data.frame_train' = total_data.frame_train, 'total_data.frame_test'  = total_data.frame_test))
    } else { # TrainTest = FALSE
      keras_outcome <- c()
      for(n in names(table(final_outcome))) keras_outcome <- cbind(keras_outcome, as.numeric(final_outcome == n))
      colnames(keras_outcome) <- names(table(final_outcome))
      if (!is.null(outcome)) 
        total_data.frame <- as.data.frame(cbind(final_pt_name, final_ROI_name, final_outcome))
      else 
        total_data.frame <- as.data.frame(cbind(final_pt_name, final_ROI_name))
      return(list('ROIList' = ROIList, 'cases_data.frame' = as.data.frame(cbind(pt_name, ROI_name, pt_outcome)), 
                  'keras_images_array' = imgArr, 'keras_outcome' = keras_outcome, 'total_data.frame' = total_data.frame))
    }
  }
  
  

  
  ###########################################
  ###  function without virtual biopsies  ###
  ###########################################
  # loop inside directories for building image series
  # ROIname:  the name given by the user to extract voxel cube
  # ROIs:     the mapROI result, a list where each element is a patient
  # ROI:      the ROI being mapped now
  for (patients_folder in pathLIST) { # cicle through patients
    if (any(ROIname %in% ROIs[[n]]))  # check if there is a ROIname chosen by user in the actual patient
    { 
      objT <- geoLet()
      objT$openDICOMFolder(pathToOpen = patients_folder)
      #browser()
      for (ROI in ROIname) 
        if (ROI %in% ROIs[[n]]) 
        {
          cat('\n\nMapping', ROI, 'in patient', patients_folder)
          m <- length(ROIList)
          if (ROIBitMask == TRUE)
            TEMPimg <- objT$getROIVoxels(Structure = ROI, ROIBitMask = TRUE, new.pixelSpacing = New_Pixel_Spacing)
          else
            TEMPimg <- objT$getROIVoxels(Structure = ROI, ROIBitMask = FALSE, new.pixelSpacing = New_Pixel_Spacing)
          
          if (min(TEMPimg$masked.images$voxelCube, na.rm = TRUE) < 0) baseline <- c(baseline, -2^15) else baseline <- c(baseline, 0) # baseline for preventing negative ROIs
          
          ROIList[[m + 1]] <- TEMPimg
          # normalize structure 
          if (!is.null(ROINameNorm)) { 
            normList[[m + 1]] <- objT$getROIVoxels(Structure = ROINameNorm, ROIBitMask = TRUE)$masked.images$voxelCube - baseline[m + 1]
            TEMPimg$masked.images$voxelCube <- TEMPimg$masked.images$voxelCube / mean(x = normList[[m + 1]], na.rm = TRUE) - baseline[m + 1]
          }
          pt_name <- c(pt_name, patients_folder)
          ROI_name <- c(ROI_name, ROI)
          if (!is.null(outcome)) 
            pt_outcome <- c(pt_outcome, outcome[n])
        }
      rm(objT)
      gc()
    } else ROIList[[n]] <- NULL
    n <- n + 1
  }
  
  # find mazimum dimensions of images and depth of final array
  if (ROIBitMask == TRUE)
  {
    # takes the maximum of voxel cubes
    rows_F  <- max(unlist(lapply(X = lapply(X = ROIList, FUN = function(x) return(dim(x$masked.images$voxelCube))), FUN = function(y) return(y[1]))))
    cols_F  <- max(unlist(lapply(X = lapply(X = ROIList, FUN = function(x) return(dim(x$masked.images$voxelCube))), FUN = function(y) return(y[2]))))
  }
  if ((TargetSize == 'auto') && (ROIBitMask == FALSE))
  { 
    # takes the mean of voxel cubes useful when you don't have masked voxel cubes and you want to prevent to take too noise from sourranding pixels
    rows_F  <- round(mean(unlist(lapply(X = lapply(X = ROIList, FUN = function(x) return(dim(x$masked.images$voxelCube))), FUN = function(y) return(y[1])))))
    cols_F  <- round(mean(unlist(lapply(X = lapply(X = ROIList, FUN = function(x) return(dim(x$masked.images$voxelCube))), FUN = function(y) return(y[2])))))
  }
  if ((TargetSize != 'auto') && (ROIBitMask == FALSE)) {
    if (!is.numeric(TargetSize) || (length(TargetSize) != 2)) stop('TargetSize MUST be either `auto` OR a numeric vector of length 2 for rows and columns number')
    rows_F <- TargetSize[1]
    cols_F <- TargetSize[2]
  }
  depth_F <- sum(unlist(lapply(X = lapply(X = ROIList, FUN = function(x) return(dim(x$masked.images$voxelCube))), FUN = function(y) return(y[3]))))
  
  # function for copying a smaller array in the center of a larger new one
  # input_array dim : number of rows, number of columns, number of images
  # output_array dim: number of images, number of rows, number of columns, number of channels (1: black and white, 3: RGB colors)
  copyARR <- function(input_array, output_array, Start, Stop) {
    # check the size of the array
    if ((dim(input_array)[1] > dim(output_array)[1]) || (dim(input_array)[2] > dim(output_array)[2]) || (dim(input_array)[3] > dim(output_array)[3]))
      stop('size of output_array must be bigger than in input_array')
    # calculate center of array with distance from border
    center <- round(c(dim(output_array)[1] - dim(input_array)[1], dim(output_array)[2] - dim(input_array)[2]) / 2)
    cat('\nCenter rows: ', center[1], ';  Center cols: ', center[2], ';  Start: ', Start, ';  Stop: ', Stop, sep = '')
    cat('\ndim input array:', dim(input_array), '\n')
    # copy the array in the center and between Start and Stop values
    output_array[(center[1]:(center[1] + dim(input_array)[1] - 1)), (center[2]:(center[2] + dim(input_array)[2] - 1)), Start:Stop] <- input_array
    return(output_array)
  }
  
  # function for copying and resizing matrices by interpolation
  rescale <- function(x, newrange=range(x)){
    xrange <- range(x)
    mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
    newrange[1]+(x-xrange[1])*mfac
  }
  
  ResizeMat <- function(mat, ndim = dim(mat)){
    if(!require(fields)) stop("`fields` package required.")
    
    # input object
    odim <- dim(mat)
    obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
    
    # output object
    ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
    ndim <- dim(ans)
    
    # rescaling
    ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
    loc <- ncord
    loc[,1] = rescale(ncord[,1], c(1,odim[1]))
    loc[,2] = rescale(ncord[,2], c(1,odim[2]))
    
    # interpolation
    ans[ncord] <- interp.surface(obj, loc)
    
    return(ans)
  }
  
  # function for resizing and copying an array into another
  resizeCopyARR <- function(input_array, output_array, Start, Stop) {
    cat('\ndim input array:', dim(input_array), '\n')
    cat('\nResize scale:', dim(output_array)[1] / dim(input_array)[1], 'x', dim(output_array)[2] / dim(input_array)[2], '\n')
    # copy the array in the center and between Start and Stop values
    output_array[,,Start:Stop] <- apply(X = input_array, MARGIN = 3, FUN = ResizeMat, ndim = c(rows_F, cols_F))
    return(output_array)
  }
  # empty array
  if (ROIBitMask == TRUE)
    img_arr <- array(data = NA, dim = c(rows_F + 2, cols_F + 2, depth_F))  # tradiational shape: rows, columns, depth
  else 
    img_arr <- array(data = NA, dim = c(rows_F, cols_F, depth_F))  # tradiational shape: rows, columns, depth
  
  a <- 0
  to_delete <- c() # vector with rows to be deleted in the output data.frame
  # three options available: zero: flatten to background, min: flatten to minimum value for each patient
  # background <- c('mean') 
  cat('\n\n\nAssembling new image array...')
  cat('\nInitial array size is:', dim(img_arr))
  for (n in 1:length(ROIList)) {
    if (!is.null(ROIList[[n]])) {
      # create the background
      tempARR <- ROIList[[n]]$masked.images$voxelCube - baseline[n]
      # filter the temporary array selecting images over the threshold
      selected_images <- which(apply(X = tempARR, MARGIN = 3, FUN = function(x) return(length(which(x > 0)))) > threshold_filter)
      # check length of images to keep, if no image kept delete the case
      if (length(selected_images) == 0) {
        cat('\n\n', n, ') Excluded patient: ', pt_name[n], '; Structure: ', ROI_name[n], sep = '')
        to_delete <- c(to_delete, n)
        rm(tempARR)
        gc()
        next
      }
      cat('\n\n', n, ') Included patient: ', pt_name[n], '; Structure: ', ROI_name[n], sep = '')
      cat('\nSelected slices:', selected_images)
      # takes only images over the threshold number of voxels
      tempARR <- tempARR[,,selected_images] 
      # check if the array is actually a single matrix (one image selected)
      if (length(dim(tempARR)) < 3) tempARR <- array(data = tempARR, dim = c(dim(tempARR), 1))
      # set the range of values between 0 and 1
      tempARR <- (tempARR - min(tempARR, na.rm = TRUE)) / (max(tempARR, na.rm = TRUE) - min(tempARR, na.rm = TRUE))
      # start of the copy
      a <- a + 1 
      # end of copy
      b <- a + dim(tempARR)[3] - 1
      # prepare values for background
      min_tempARR <- min(tempARR, na.rm = TRUE)
      mean_tempARR <- mean(x = tempARR, na.rm = TRUE)
      if (ROIBitMask == TRUE)
        # copy the array in the bigger one
        img_arr <- copyARR(input_array = tempARR, output_array = img_arr, Start = a, Stop = b)
      else 
        img_arr <- resizeCopyARR(input_array = tempARR, output_array = img_arr, Start = a, Stop = b)
      # define the outcome for tensorflow
      if (!is.null(outcome)) 
        final_outcome <- c(final_outcome, rep.int(x = pt_outcome[n], times = dim(tempARR)[3]))
      final_pt_name   <- c(final_pt_name, rep.int(x = pt_name[n], times = dim(tempARR)[3]))
      final_ROI_name  <- c(final_ROI_name, rep.int(x = ROI_name[n], times = dim(tempARR)[3]))
      # set the background
      if (backgroud == 'zero') img_arr[,,a:b][which(is.na(img_arr[,,a:b]))] <- 0
      if (backgroud == 'min')  img_arr[,,a:b][which(is.na(img_arr[,,a:b]))] <- min_tempARR
      if (backgroud == 'mean') img_arr[,,a:b][which(is.na(img_arr[,,a:b]))] <- mean_tempARR
      if (backgroud == 'runif') img_arr[,,a:b][which(is.na(img_arr[,,a:b]))] <- runif(n = length(img_arr[,,a:b][which(is.na(img_arr[,,a:b]))]), 
                                                                                      min = min_tempARR + bckwidth/2 * min_tempARR, max = max(tempARR, na.rm = TRUE) - bckwidth/2 * max(tempARR, na.rm = TRUE))
      if (backgroud == 'rnorm') img_arr[,,a:b][which(is.na(img_arr[,,a:b]))] <- rnorm(n = length(img_arr[,,a:b][which(is.na(img_arr[,,a:b]))]), 
                                                                                      mean = mean_tempARR, sd = sd(x = tempARR, na.rm = TRUE)*bckwidth)
      a <- b # set the new index for continuing copy
      rm(tempARR)
      gc()
    } else  # if NULL ROIList[[n]] remove elements from vectors of final data.frame
      to_delete <- c(to_delete, n)
  }
  
  # resize the final array
  final_dim <- c(dim(img_arr)[1:2], b)
  if (!is.null(outcome)) 
    df_outcome <- as.data.frame(cbind('keras_outcome' = 1:length(levels(as.factor(final_outcome))), 'original_outcome' = levels(as.factor(final_outcome))))
  else {
    df_outcome <- NULL
    pt_outcome <- NULL
    final_outcome <- NULL
  }
  img_arr <- (img_arr - min(img_arr, na.rm = TRUE)) / (max(img_arr, na.rm = TRUE) - min(img_arr, na.rm = TRUE)) # final scale array between 0 and 1
  img_arr <- array(data = img_arr[which(!is.na(img_arr))], dim = final_dim)
  cat('\n\nFinal number of images selected is:', dim(img_arr)[3], '\n\n')
  # componing final dataset
  if (!is.null(outcome)) 
    total_data.frame <- as.data.frame(cbind(final_pt_name, final_ROI_name, final_outcome))
  else 
    total_data.frame <- as.data.frame(cbind(final_pt_name, final_ROI_name))
  
  if (TrainTest == FALSE) # the case of TrainTest = FALSE 
  { 
    # array in the shape that keras requires
    final_array <- array(data = NA, dim = c(dim(img_arr)[3],dim(img_arr)[1],dim(img_arr)[2],1))
    for(i in 1:dim(final_array)[1]) final_array[i,,,] <- img_arr[,,i]
    
    return(list('ROIList' = ROIList, 'cases_data.frame' = as.data.frame(cbind(pt_name, ROI_name, pt_outcome)), 
                'keras_images_array' = final_array, 
                'keras_outcome' = as.numeric(as.factor(final_outcome)), 'df_outcome' = df_outcome, 'total_data.frame' = total_data.frame))
  }
  else { # TrainTest is TRUE so two datasets are prepared
    # create the two series of indeces for images and outcome data
    if (per_patient_sampling == FALSE) { # samples for each image
      train_seq <- sample(x = 1:b, size = round(ratio * b), replace = FALSE)
      test_seq  <- c(1:b)[-train_seq]
    } else { # sample for each patient for preventing overfitting
      train_pt <- sample(x = names(table(pt_name)), size = round(ratio * length(table(pt_name))), replace = FALSE)
      test_pt  <- names(table(pt_name))[-which(pt_name %in% train_pt)]
      train_seq <- c(1:b)[which(final_pt_name %in% train_pt)]
      test_seq <- c(1:b)[-train_seq]
    }
    # shaping the final arrays
    cat('\nShaping the final training array...\n')
    img_arr_train <- array(data = NA, dim = c(length(train_seq), dim(img_arr)[1], dim(img_arr)[2],1))
    total <- dim(img_arr_train)[1]
    for(i in 1:dim(img_arr_train)[1]) {
      # create progress bar
      pb <- txtProgressBar(min = 0, max = total, style = 3)
      # update progress bar
      setTxtProgressBar(pb, i)
      img_arr_train[i,,,] <- img_arr[,,train_seq[i]]
    }
    close(pb)
    
    cat('\nShaping the final testing array...\n')
    img_arr_test  <- array(data = NA, dim = c(length(test_seq),  dim(img_arr)[1], dim(img_arr)[2],1))
    total <- dim(img_arr_test)[1]
    for(i in 1:dim(img_arr_test)[1])  {
      # create progress bar
      pb <- txtProgressBar(min = 0, max = total, style = 3)
      # update progress bar
      setTxtProgressBar(pb, i)
      img_arr_test[i,,,]  <- img_arr[,,test_seq[i]]
    }
    total_data.frame_train <- total_data.frame[train_seq,]
    total_data.frame_test  <- total_data.frame[test_seq,]
    return(list('ROIList' = ROIList, 'cases_data.frame' = as.data.frame(cbind(pt_name, ROI_name, pt_outcome)), 
                'keras_images_array_train'  = img_arr_train, 'keras_images_array_test' =  img_arr_test, 
                'keras_outcome_train' = as.numeric(as.factor(final_outcome))[train_seq], 
                'keras_outcome_test'  = as.numeric(as.factor(final_outcome))[test_seq], 
                'df_outcome' = df_outcome, 'total_data.frame_train' = total_data.frame_train, 'total_data.frame_test' = total_data.frame_test))
  }
}
