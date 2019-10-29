#'  Filter to select a particular subpopulation in the structure
#' 
#' @description  Instatiated an object of the class geoLet, Select works on a single patient containing the extracted ROI
#' N.B: pay attention, Select requires a geoLet object in which the ROI is already extracted. See the examples for more details.
#' Select applies two image preprocessing steps on each study you want to analyze:
#' 1) Normalize the values inside ROI, considering the histogram of gray levels
#'    and setting as extremes fstPerc and lst Perc
#'    NormValue = (OriginalValue-fstPerc)/(lstPerc-fstPerc)
#' 2) Binarize the images inside ROI by applying a Threshold, and setting to 
#'    normalised values, the pixels with intensities between ThresholdDown and ThresholdUp
#'    NaN all values with intensities outside the threshold values
#' @param fstPerc is the lowest Percentile of the histogram that you use to normalize the pixel values
#' @param lstPerc is the highest Percentile of the histogram that you use to normalize the pixel values
#' @param ThDown is the lowest threshold value  
#' @param ThUp is the highest threshold value  
#' @param ThStep is the step threshold value
#' @param without.na is the option to decide to if analyze the image with or without na. Set this parameter to TRUE is suggested.
#' @return  #' Select returns a mmButo oject with the values selected
#' 
#'@examples \dontrun{
#'#First of all create an object obj of the class geoLet()
#'obj <- geoLet()
#'
#'#now load the DICOM serie of the patient that you want analyze
#'#using the method openDICOMFolder and indicating the path
#'# obj$openDICOMFolder(Path = "/home/kboaria/Desktop/Davide/Retto/Validation/Patient1/")
#'
#'#Before running Select you have to extract the ROI that you want analyze by using the method getROIVoxels of the geoLet object
#'geoLetVoxelList<-obj$getROIVoxels(Structure = "GTV")
#'# Now you can run Select on GTV object;#'
#'
#'}
#' @export
Select<-function(geoLetVoxelList, fstPerc=0, lstPerc=100, ThDown=0, ThUp=100, ThStep=10, without.na=FALSE) {
  
  a <- Normer(geoLetVoxelList = geoLetVoxelList, fstPerc=fstPerc, lstPerc=lstPerc)
  b <- SThresholder(geoLetVoxelList = a, ThDown=ThDown, ThUp=ThUp, ThStep=ThStep, without.na=without.na) 
  return(b)
}

#' 
#' @description  Thresholder prende in ingresso un geoLetVoxelList, cioè un oggetto con ROI 
#' già estratta
#' Questa funzione restituisce la stessa geoLetVoxelList, ma con dentro i VoxelCube sogliati
SThresholder<-function(geoLetVoxelList, ThDown=ThDown, ThUp=ThUp, ThStep = ThStep, without.na=FALSE) {

  
  objS<-services();
  
  ROIVoxelData<-geoLetVoxelList
  
    # prendi i voxelData ed estendili
    voxelData.ready<-ROIVoxelData$masked.images$voxelCube
    
    # salvo le varie soglie rispetto al valore minore, valore maggiore e lo step
    ListaSoglie <- seq(from = ThDown, to =  ThUp, by = ThStep)
    ListaRange <- t(combn(ListaSoglie, 2))
    FilteredROI <- list()
    
    for(RangeTh in 1:nrow(ListaRange)){
      
      for (q in 1:dim(voxelData.ready)[3])  {
        ####################Binarizzazione da ThDown a ThUp#################################
        arrayRM2_3d<-voxelData.ready[,,q]
        arrayRM2_3d[which(arrayRM2_3d<ListaRange[RangeTh,1])]<-NA
        #-im Cambio i minori uguali 07-12-16DC#
        arrayRM2_3d[which(arrayRM2_3d>ListaRange[RangeTh,2])]<-NA
        #-fm Cambio i minori uguali 07-12-16DC#
        
        if (without.na ==T){
          #Sostituisco i Na con 0
          arrayRM2_3d[is.na(arrayRM2_3d)] <- 0
        }
        
        ROIVoxelData$masked.images$voxelCube[,,q]<-arrayRM2_3d
      }
      FilteredROI[[paste0(ListaRange[RangeTh,], collapse = ',')]]$masked.images$voxelCube <- ROIVoxelData$masked.images$voxelCube
      FilteredROI[[paste0(ListaRange[RangeTh,], collapse = ',')]]$masked.images$location <- ROIVoxelData$masked.images$location
    }
    
  return(FilteredROI) 
}  

#'  badali'3
#' 
#' @description  function that apply the normalization
Normer<-function(geoLetVoxelList, fstPerc=fstPerc, lstPerc=lstPerc) 
  #Normer prende in ingresso un mButoVoxel List, cioè un oggetto mButo che ha già caricato
  #gli studi e ha già estratto la ROI 
  #Questa funzione restituisce la stessa geoLetVoxelList, ma con dentro i VoxelCube le immagini normalizzate
{
  
  objS<-services();
  
  ROIVoxelData<-geoLetVoxelList
  

    # prendi i voxelData 
    voxelData.ready<-ROIVoxelData$masked.images$voxelCube
    newStudio<-voxelData.ready+100000
    q<-c(fstPerc,lstPerc)/100
    Quantils<-quantile(newStudio, probs = q, na.rm = T)
    MaxStudio<-Quantils[2]
    MinStudio<-Quantils[1]
    newStudio[which(newStudio>MaxStudio)]<-NaN
    newStudio[which(newStudio<MinStudio)]<-NaN
    MaxStudio<-max(newStudio, na.rm = T)
    MinStudio<-min(newStudio, na.rm = T)
    NormStudio<-100*((newStudio-MinStudio)/(MaxStudio-MinStudio))
    ROIVoxelData$masked.images$voxelCube<-NormStudio

  return(ROIVoxelData) 
}  