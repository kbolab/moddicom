#' class for loading and presenting NIFTI data
#' 
#' @description  Instantiate an object of the class \code{geoNIFTI}.This represents just the classname, 
#'               methods are exposed with the technique of 'closure'.
#'               In order to see manuals for the single mathods, consider the vignette or use the 
#'               available for the following wrapping functions:
#'               \itemize{
#'               \item \code{GLT.getImageVoxelCube( );} : to get the ImageVoxelCube stored into a geoNIFTI object
#'               \item \code{GLT.getROIVoxels( );} : to get the IMAGE Voxels geometrically located into a ROI, for a given geoNIFTI object
#'               }
#'               The original methods for the class geoLet can also be invocked using the same name without the previx 'GTL.', i.e.:
#' @examples \dontrun{
#' # first of all create an object obj
#' obj.geoNIFTI<-geoNIFTI()
#' tmp<-obj.geoNIFTI$getROIVoxels("/Users/davpas/Desktop/bis/BOOGLO","*metastatic*")
#' tmp2<- obj.geoNIFTI$getImageVoxelCube("/Users/davpas/Desktop/bis/BOOGLO")
#' matrixImageROI<-tmp$masked.images$voxelCube
#' 
#' }
#' @export 
#' @import oro.nifti
#' @useDynLib moddicom
#' 



geoNIFTI <- function(){
  
  # rotate makes a rotation of images 
  rotate <- function(x) t(apply(x, 2, rev))
  
  readNIfTI_frompath_rotate<- function(matrixNIFTI,matrixROI){
    ##########################################
    # readNIfTI_frompath_rotate
    # this function takes in input a n-dimensional matrix of data (dimension NxMxZ row, columns and slices)
    # return in output:
    # 1) matrixRegionOfInterest: a NxMxZ matrix given by matrixNIFTI*matrixROI
    # 2) matrixROI: the rotated matrix of ROI with the same dimension (NxMxZ)
    # It uses the function rotate because original image and matrixROI are rotated.
    #########################################
    print("Computing matrixNIFTI*matrixROI")
    nslices= dim(matrixROI)
    for (i in 1:nslices[3]){
      matrixROI[,,i]<-apply(rotate(rotate(matrixROI[,,i])),2,rev)
      
    }
    
    
    matrixRegionOfInterest<-matrixNIFTI*matrixROI
    
    return(list("matrixRegionOfInterest"=matrixRegionOfInterest,"matrixROI_rotated"=matrixROI))
  }
  
  
  
  ImageVoxelSelection<-function(matrixRegionOfInterest){
    #######################################################
    # ImageVoxelSelection takes in input matrixRegionOfInterest (NxMxZ)
    # returns in output c(start_roi_slice,end_roi_slice), a list of two number:
    # 1) start_roi_slice: the number of first slice with a ROI inserted
    # 2) end_roi_slice: the number of last slice with a ROI inserted
    #
    #######################################################
    print("Estimating Voxel selection")
    #browser()
    nslices= dim(matrixRegionOfInterest)
    #print(nslices)
    start_roi_slice<-0
    end_roi_slice<-0
    for (i in 1:nslices[3]){
      if (start_roi_slice==0){
        if(sum(matrixRegionOfInterest[,,i])!=0){
          start_roi_slice<- i
        }
      }
    }
    
    for (i in 1:nslices[3]){
      if (end_roi_slice==0){
        if(sum(matrixRegionOfInterest[,,nslices[3]-i])!=0){
          end_roi_slice<- nslices[3]-i
        }
      }
    }
    
    
    return(c(start_roi_slice,end_roi_slice))
    
  }
  
  
  resizeROIVoxel<-function(matrixRegionOfInterest){
    ################################################
    # resizeROIVoxel takes in input matrixRegionOfInterest (NxMxZ)
    # returns as output c(min_row,max_row, min_column,max_column)
    # this function computes the min and max row and columns of the Voxel cube 
    ##################################################
    
    print("Estimating ROI Voxel cube")
    nslices= dim(matrixRegionOfInterest)
    min_row<-10000
    max_row<- -100
    min_column<-10000
    max_column<- -100
    for (i in 1:nslices[3]){
      list_row_column_pixel_roi <-which(matrixRegionOfInterest[,,i]!=0, arr.ind = TRUE)
      tmp_min_row<-min(list_row_column_pixel_roi[,1])
      tmp_max_row<- max(list_row_column_pixel_roi[,1])
      tmp_min_column<- min(list_row_column_pixel_roi[,2])
      tmp_max_column<- max(list_row_column_pixel_roi[,2])
      if (tmp_min_row<min_row){
        min_row<- tmp_min_row 
      }
      if (tmp_max_row>max_row){
        max_row<- tmp_max_row 
      }
      
      if (tmp_min_column<min_column){
        min_column<- tmp_min_column
      }
      
      if (tmp_max_column>max_column){
        max_column<- tmp_max_column
      }
      
    }
    return(c(min_row,max_row, min_column,max_column)) 
  }
  
  
  getImageVoxelCube<- function(path){
    ##############################################
    # getImageVoxelCube, given a path that is associated to a patient,
    # returns as output a NxMxZ n-dimensional matrix with the original 
    # NIFTI data.
    ###############################################
    aa <- getwd()
    setwd(path)
    list_of_files_i<-list.files()
    list_lenght_names= nchar(list_of_files_i)
    index_minimum_filename <- which(list_lenght_names==min(list_lenght_names))
    filenameNIFTI= list_of_files_i[index_minimum_filename[2]]
    print("Opening NIFTI image")
    print(filenameNIFTI)
    imgNIFTI <- readNIfTI(filenameNIFTI)
    list_of_values<- imgNIFTI@pixdim
    matrixNIFTI<-imgNIFTI@.Data[,,]
    
    #ritorno al path iniziale
    setwd(aa)
    
    return(matrixNIFTI)
    
  }
  
  
  getROIVoxels<-function(pathNIFTI, regulexp_roi){
    
    ##############################################
    # getROIVoxels given a path that is associated to a patient 
    # and a regular expression used in order to specify the ROI
    # It returns as output the following list:
    # list("geometricalInformationOfImages"=geometricalInformationOfImages, "masked.images"=masked.images)) 
    # 1) geometricalInformationOfImages contains: 
    # pixelSpacing
    # SliceThickness
    # supposedNumberOfSlices
    # 2) masked.images cointaining: 
    # - a voxel cube,  a UxWxT n-dimensional matrix that is a subset with the ROI of the original image
    # - coordinates of the VoxelCube: min.x, max.x, min.y, max.y, min.z, max.z
    ###############################################

    aa <- getwd()
    setwd(pathNIFTI)
    list_of_files_i<-list.files()
    #print(list_of_files_i)
    
    #find shortest filename
    list_lenght_names= nchar(list_of_files_i)
    index_minimum_filename <- which(list_lenght_names==min(list_lenght_names))
    filenameNIFTI= list_of_files_i[index_minimum_filename[2]]
    tmpfilenameROI=list.files(pattern=regulexp_roi)
    
    if (length(tmpfilenameROI)==0) {
      cat(paste("No ROI:",regulexp_roi))
      cat("\n")
      setwd(aa)
      return(list("geometricalInformationOfImages"=NULL, "masked.images"=NULL))  
    }
    if (length(tmpfilenameROI)==1) {
      filenameROI=tmpfilenameROI
    }
    if (length(tmpfilenameROI)==2) {
      filenameROI=tmpfilenameROI[1]
    }
    if (length(tmpfilenameROI)>2) {
      stop("Error: regular expression is associated to more than 2 files (.img and .hdr)")
    }
    
  
    print("reading NIFTI")
    print(filenameNIFTI)
    imgNIFTI <- readNIfTI(filenameNIFTI)
    print("reading ROI NIFTI")
    print(filenameROI)
    imgROI <-readNIfTI(filenameROI)
    

    if(length(dim(imgNIFTI@.Data)) == 4){
      image_voxels <- imgNIFTI@scl_slope * imgNIFTI@.Data[,,,1] + imgNIFTI@scl_inter
    }
    else{
      image_voxels <- imgNIFTI@scl_slope * imgNIFTI@.Data + imgNIFTI@scl_inter
    }
    
    matrixNIFTI <- image_voxels 
    #matrixNIFTI<- image_voxels + abs(min(image_voxels,na.rm = T))
    
    
    list_of_values<- imgNIFTI@pixdim
    #image_parameters<- imgNIFTI@pixdim
    matrixROI<-imgROI@.Data[,,]
    
    
    matrixImageRoi_and_matrixROI<-readNIfTI_frompath_rotate(matrixNIFTI,matrixROI)
    matrixImageRoi<-matrixImageRoi_and_matrixROI$matrixRegionOfInterest
    
    start_stop_slice <- ImageVoxelSelection(matrixImageRoi)
    
    resized_index = resizeROIVoxel(matrixImageRoi[,,start_stop_slice[1]:start_stop_slice[2]])
    
    pixelSpacing<<-list_of_values[2:3]
    SliceThickness<<-toString(list_of_values[4])
    supposedNumberOfSlices=start_stop_slice[2]-start_stop_slice[1]+1
    #geometricalinfo<-list()
    geometricalInformationOfImages<<-list("pixelSpacing"=pixelSpacing,"SliceThickness"=SliceThickness, "supposedNumberOfSlices"=supposedNumberOfSlices)
    
    list_pixel_dati<<-list_of_values
    matrixROI<-matrixImageRoi_and_matrixROI$matrixROI
    matrixImageRoi[matrixROI== 0] <- NA
    matrice_dati_ndimensionale<-matrixImageRoi[resized_index[1]:resized_index[2],resized_index[3]:resized_index[4],start_stop_slice[1]:start_stop_slice[2]]
    location<-list("min.x"=resized_index[1],"max.x"=resized_index[2],"min.y"=resized_index[3], "max.y"=resized_index[4], "min.z"=start_stop_slice[1],"max.z"=start_stop_slice[2])
    masked.images<<-list("voxelCube"=matrice_dati_ndimensionale,"location"=location)
    
    #ritorno al path iniziale
    setwd(aa)
    
    return(list("geometricalInformationOfImages"=geometricalInformationOfImages, "masked.images"=masked.images))    
    
  }
  
  return(list(
    "getImageVoxelCube" = getImageVoxelCube,
    "getROIVoxels"= getROIVoxels
  ))
}
