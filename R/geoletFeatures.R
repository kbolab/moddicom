#' Function extractor for multiple filter parameter
#'
#' @description  It extracts the feature values by applying the filter at multiple LoG sigma.
#' @param path The path to be browsed for folders containing DICOM images and Structure Set of patients
#' @param ROIName A \code{character} object containing the name(s) (as vector) of ROIs to be analyzed. If there is more than one ROI (multiple references) or no ROI is listed the analysis is stopped
#' @param feature.family The list of features to be extracted, by default all. The possible families are: statistical, grey level co-occurrence matrix, run-length matrix, size-zone matrix and fractal
#' @param interpolate A \code{logical} value, by default \code{FALSE}, for interpolating sigma values in filters using sigma for computation as LoG
#' @param px The pixel spacing in \emph{x} direction as string for interpolating sigma, by default empty.
#' @param py The pixel spacing in \emph{y} direction as string for interpolating sigma, by default empty.
#' @param threshold A \code{numeric} value setting the threshold for applying the interpolation, by default .001
#' @param discretize ROI discretization for feature extractions, by default \code{none}. The possible values are \code{fixed.bin.size} and \code{fixed.bin.number}.
#' @param bin.size Size of bin for ROI discretization for feature extractions, by default 25.
#' @param bin.number Number of bin for ROI discretization for feature extractions, by default 25.
#' @param fileName The Rdata file name saved locally during the running, by default 'tmp.f.extractor.sing.par.RData'.
#' @param forceRecalculus A \code{logical} value, by defaul \code{TRUE}, for forcing recalculus by deleting all the local Rdata previously created.
#' @param from.sigma The first sigma to be analyzed
#' @param to.sigma The last sigma to be analyzed
#' @param def.by The sigma step to be analyzed
#' @param strategy The strategy to be applied for the radiomics analysis, by default "fixed". For now no other strategy have been implemented.
#' @export 
#' @import radiomics data.table plyr
#' @examples \dontrun{
#' f.extractor.pluri.LoG.par(path = './RadiomicsFolder', ROIName = 'GTV', feature.family = c('stat', 'glcm'), forceRecalculus = T, from.sigma = 0.3, to.sigma = 2, def.by = 0.1)
#' }
f.extractor.pluri.LoG.par<- function(path, ROIName,
                                 feature.family=c("stat","morph","glcm","rlm","szm"),
                                 interpolate = FALSE, px="", py="", threshold = .1, discretize="", bin.size = 25, bin.number=25,
                                 fileName = "tmp.f.extractor.pluri.par.RData" ,
                                 forceRecalculus = TRUE, from.sigma=from.sigma, to.sigma=to.sigma, def.by = def.by, strategy="fixed") {
  
  cat(path)
  cat("\n")
  sigma.array <-c()
  big.matrice<-list()
  if(strategy=="fixed") sigma.array <- seq( from.sigma, to.sigma, by = def.by)
  
  if(forceRecalculus==FALSE) {
    if(file.exists(fileName)==TRUE) {
      load(file = fileName)  
    }
  } else {

    if(file.exists(fileName)==TRUE) {
      file.remove(fileName)
    }
  }  
  
  for( iterazione in sigma.array) {
    
    if( !(iterazione %in% names(big.matrice))) {
    
      filterPipeline<- list()
      filterPipeline[[1]]<-list("kernel.type"="LoG", "sigma"=iterazione)
  
      nomeFile <- paste(c( "f.extractor.single.par_",iterazione,".RData" ),collapse='')
      
      tmp.matrice <- f.extractor.sing.par(path = path,ROIName = ROIName,feature.family = feature.family,filterPipeline = filterPipeline,interpolate = interpolate,px = px,py = py, 
                                          fileName = nomeFile, forceRecalculus = forceRecalculus, discretize = discretize, bin.size = bin.size, bin.number = bin.number)
      big.matrice[[as.character(iterazione)]]<-tmp.matrice
      save(big.matrice, file = fileName)
      file.remove(nomeFile)
      
    }
  }
  return(big.matrice);
}

#' Function extractor for single filter parameter
#'
#' @description  It extracts the feature values by applying the filter at a single LoG sigma. 
#' @param path The path to be browsed for folders containing DICOM images and Structure Set of patients
#' @param ROIName A \code{character} object containing the name(s) (as vector) of ROIs to be analyzed. If there is more than one ROI (multiple references) or no ROI is 
#' listed the analysis is stopped
#' @param feature.family The list of features to be extracted, by default all. The possible families are: statistical, grey level co-occurrence matrix, run-length matrix, size-zone matrix and fractal
#' @param filterPipeline A \code{list} contains the filter pipeline to be applied, by default empty, to images before feature extraction. The first list argument is "kernel.type" as character and the second argument is "sigma" as numeric.
#' @param interpolate A \code{logical} value, by default \code{FALSE}, for interpolating sigma values in filters using sigma for computation as LoG
#' @param px The pixel spacing in \emph{x} direction as string for interpolating sigma, by default empty.
#' @param py The pixel spacing in \emph{y} direction as string for interpolating sigma, by default empty.
#' @param threshold A \code{numeric} value setting the threshold for applying the interpolation, by default .001
#' @param discretize ROI discretization for feature extractions, by default \code{none}. The possible values are \code{fixed.bin.size} and \code{fixed.bin.number}.
#' @param fileName The Rdata file name saved locally during the running, by default 'tmp.f.extractor.sing.par.RData'.
#' @param bin.size Size of bin for ROI discretization for feature extractions, by default 25.
#' @param bin.number Number of bin for ROI discretization for feature extractions, by default 25.
#' @param forceRecalculus A \code{logical} value, by defaul \code{TRUE}, for forcing recalculus by deleting all the local Rdata previously created.
#' @details Using an object of class \code{ROImap} it is possible to list all the ROIs contained in an explored directory for the analysis. It is possible to analyze
#' different ROIs names referencing to the same structure in the different studies by putting al their names a a vector of \code{character} objects in the argument \code{ROIName}
#' @export
#' @examples \dontrun{
#' # Create a filter pipeline list
#' filterPip <- list()
#' filterPip[[1]] <- list('kernel.type'='LoG', 'sigma' = 0.3)
#' f.extractor.sing.par(path = './RadiomicsFolder', ROIName = 'GTV', feature.family = c('stat', 'glcm'), filterPipeline = filterPip, forceRecalculus = T)
#' }
f.extractor.sing.par<- function(path, ROIName ,
                          feature.family=c("stat","morph","glcm","rlm","szm"),
                          filterPipeline=list(), interpolate = FALSE, px="", py="", threshold = .001, discretize="", bin.size = 25, bin.number=25,
                          fileName = "tmp.f.extractor.sing.par.RData" ,
                          forceRecalculus = TRUE) {
  
  if((px == "" | py=="") & !(px == "" & py == "" ))  {
    stop("ERRORE: px e py devono essere specificati entrambi!")
  }

  patList <- list.dirs(path = path,recursive = FALSE)
  matrice <- c()
  # -im Carlotta 23112017
  rigaFeatures <- c()
  patientsPath <- c()
  # -fm Carlotta 23112017
  
  if(forceRecalculus==FALSE) {
    if(file.exists(fileName)==TRUE) {
      load(file = fileName)  
    }   
  } else {
    if(file.exists(fileName)==TRUE) {
      file.remove(fileName)
    }    
  }
  
  for( patID in patList) {
    
    if(!(patID %in% matrice[,1])) {
      cat(patID)
      cat("\n")

      temp.ROIName <- ROIName
      ogg.geoLet <- geoLet()
      ogg.geoLet$openDICOMFolder(pathToOpen = patID)
      #browser();
      cat("\n ---------------------------------\n Computing Patient: ",patID,"\n")
      
      if (length(which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)) > 1) {  # check for multiple ROI mapped
        str <- ogg.geoLet$getROIList()[2, ][which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)] # multiple ROI names mapped
        stop(paste('ERROR: more than one ROI mapped in a single structure set:', paste(str, collapse = ', ')))
      }
      if (length(which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)) == 0) 
        stop(paste('ROI', paste(ROIName, collapse = ', '), 'not mapped into structure set')) # check for not mapped ROI
      ROIName <- ROIName[which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)]
      riga <- computeFeatures.geoLet( obj.geoLet = ogg.geoLet, ROIName = ROIName ,
                                      feature.family=feature.family, filterPipeline=filterPipeline,
                                      px = px, py = py, discretize=discretize, bin.size=bin.size,bin.number=bin.number) 
      # aggiungi il nome del paziente
      
      # -im Carlotta 23112017
      # riga <- c( patID , riga)
      # matrice <- rbind(matrice,riga)
      rigaFeatures <- rbind( rigaFeatures , unlist(riga))
      patientsPath <- rbind(patientsPath, as.character(patID))
      matrice <- data.frame("patID"=patientsPath, rigaFeatures)
      # -fm Carlotta 23112017
      
      # salva la matrice
      save(matrice,file = fileName)
      ROIName <- temp.ROIName
    }
    
  }
  # -im Carlotta 23112017
  matrice$patID <- as.character(matrice$patID)
  # -fm Carlotta 23112017
  return(matrice)
}

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

#' computeFeatures.geoLet lavora su piu
#'
#' @description  Intiate an object of the class \code{geoLet}.This represents just the classname (f.extractor.sing.par)
#' @export
computeFeatures.geoLet<- function( obj.geoLet, ROIName , feature.family=c("stat","morph","glcm","rlm","szm"), filterPipeline=c(),
                                   px = "", py ="", threshold = .1, discretize="", bin.size = "", bin.number="", n_grey = 100) {

  objS <- services()
  if( (px == "" | py=="") & !(px == "" & py == "" ))  {
    stop("ERRORE: px e py devono essere specificati entrambi!")
  }
  if(px == "" & py == "" )
    ROI <- obj.geoLet$getROIVoxels(Structure = ROIName)
  else{
    pixelSpacing <- obj.geoLet$getPixelSpacing()
    if( (abs(pixelSpacing[1]-as.numeric(px)))>threshold | (abs(pixelSpacing[2]-as.numeric(py))>threshold) ){
      ROI <- obj.geoLet$getROIVoxels(Structure = ROIName, new.pixelSpacing=c(px,py))
    }
    else{
      ROI <- obj.geoLet$getROIVoxels(Structure = ROIName)
    }
  }

  # Ragionamento sul pixelspaxing legato alle esigenze del filtro di convoluzione
  old.px <- ROI$geometricalInformationOfImages$pixelSpacing[1]
  old.py <- ROI$geometricalInformationOfImages$pixelSpacing[2]
  old.pz <- as.numeric(ROI$geometricalInformationOfImages$SliceThickness)

  if( px == "") px <- ROI$geometricalInformationOfImages$pixelSpacing[1]
  if( py == "") py <- ROI$geometricalInformationOfImages$pixelSpacing[2]
  pz <-  as.numeric(ROI$geometricalInformationOfImages$SliceThickness)


  #filtro se richiesto
  if(length(filterPipeline)!=0) {
    res <- FIL.2D.conv.Filter.geoLet(obj.geoLet = obj.geoLet,ROIName = ROIName,
                                     filter.pipeline = filterPipeline,
                                     scaleFactor = 'space', px = px, py = py)

  } else {
    res <- ROI$masked.images
  }

  naVALUE <- as.integer(min(voxelCube = res$voxelCube,na.rm = T)- 10)
  res.noNA <- res$voxelCube
  res.noNA[ which( is.na(res.noNA),arr.ind = T  )  ] <- naVALUE
  

  #discretizzo la ROI in base ai parametri in entrata
  original <- res$voxelCube
  
  if(discretize=="fixed.bin.size"){
    cat("\n")
    cat("discretizing at fix bin size...")
    cat("\n")
    #### DISCRETIZE FIXED BIN SIZE
    for(j in 1:dim(res$voxelCube)[3]){
        discr <- res$voxelCube[,,j]
        wb <- bin.size
        Nv <- length(discr[which(!is.na(discr))])
        Xgl_min <- min(discr,na.rm = T)
        for(i in 1:Nv){
          if(discr[which(!is.na(discr))][i] == Xgl_min) discr[which(!is.na(discr))][i] <- 1
          else if(discr[which(!is.na(discr))][i] > Xgl_min) {
            discr[which(!is.na(discr))][i] <- 1 + floor((discr[which(!is.na(discr))][i] - Xgl_min) / wb)
            #discr[which(!is.na(discr))][i] <- 1 + (discr[which(!is.na(discr))][i] - Xgl_min) / wb
            #discr[which(!is.na(discr))][i] <- round(discr[which(!is.na(discr))][i],digits = 0)
          }
        }
        res$voxelCube[,,j] <- discr
    }
  }
  
  if(discretize=="fixed.bin.number"){
    cat("\n")
    cat("discretizing at fix bin number...")
    cat("\n")
    
    #### DISCRETIZE FIXED BIN NUMBER

    for(j in 1:dim(res$voxelCube)[3]){
      discr <- res$voxelCube[,,j]
      Ng <- bin.number
      Nv <- length(discr[which(!is.na(discr))])
      Xgl_min <- min(discr,na.rm = T)
      Xgl_max <- max(discr,na.rm = T)
      for(i in 1:Nv){
        if(discr[which(!is.na(discr))][i] == Xgl_min) discr[which(!is.na(discr))][i] <- 1
        else if(discr[which(!is.na(discr))][i] > Xgl_min) {
          discr[which(!is.na(discr))][i] <- 1 + (Ng*(discr[which(!is.na(discr))][i] - Xgl_min) / (Xgl_max - Xgl_min))
          discr[which(!is.na(discr))][i] <- round(discr[which(!is.na(discr))][i],digits = 0)
        }
      }
      res$voxelCube[,,j] <- discr
    }
  }

  # stat.df<-c()
  def <- c()
  if("stat" %in% feature.family) {
    cat("computing stat features...\n")
    stat.df <- statisticalFeatures(res$voxelCube)
    def <- c(def,stat.df)
  }

  if("morph" %in% feature.family) {
    cat("computing morph features...\n")
    morph.df <- morphologicalFeatures(original,px=px,py=py,pz=old.pz)
    def <- c(def,morph.df)
  }

  if("glcm" %in% feature.family){
    cat("computing glcm features...\n")
    F_cm <-  glcmTexturalFeatures(res$voxelCube, n_grey=n_grey)
    F_cm <- do.call(data.frame,lapply(F_cm, function(x) replace(x, is.infinite(x),NA)))
    cm.df <- colMeans(F_cm,na.rm = T)
    def <- c(def,cm.df)
  }
  
  if("glcm" %in% feature.family){
    cat("computing glcm merged features...\n")
    F_cm_merged <-  glcmTexturalFeaturesMerged(res$voxelCube, n_grey=n_grey)
    F_cm_merged <- do.call(data.frame,lapply(F_cm_merged, function(x) replace(x, is.infinite(x),NA)))
    cm_merged.df <- colMeans(F_cm_merged,na.rm = T)
    def <- c(def,cm_merged.df)
  }
  
  if("glcm" %in% feature.family){
    cat("computing glcm 2.5D features...\n")
    F_cm_25D <-  glcmTexturalFeatures25D(res$voxelCube, n_grey=n_grey)
    F_cm_25D <- do.call(data.frame,lapply(F_cm_25D, function(x) replace(x, is.infinite(x),NA)))
    cm_25d.df <- colMeans(F_cm_25D,na.rm = T)
    def <- c(def,cm_25d.df)
  }
  
  if("glcm" %in% feature.family){
    cat("computing glcm 2.5D merged features...\n")
    F_cm_25D_merged <-  glcmTexturalFeatures25Dmerged(res$voxelCube, n_grey=n_grey)
    F_cm_25D_merged <- do.call(data.frame,lapply(F_cm_25D_merged, function(x) replace(x, is.infinite(x),NA)))
    def <- c(def,F_cm_25D_merged)
  }

  if("rlm" %in% feature.family) {
    cat("computing glrm features...\n")
    F_rlm <-glrlmTexturalFeatures(res$voxelCube, n_grey=n_grey)
    F_rlm <- do.call(data.frame,lapply(F_rlm, function(x) replace(x, is.infinite(x),NA)))
    rlm.df <- colMeans(F_rlm,na.rm = T)
    def <- c(def,rlm.df)
  }
  
  if("rlm" %in% feature.family) {
    cat("computing glrm merged features...\n")
    F_rlm_merged <-glrlmTexturalFeaturesMerged(res$voxelCube, n_grey=n_grey)
    F_rlm_merged <- do.call(data.frame,lapply(F_rlm_merged, function(x) replace(x, is.infinite(x),NA)))
    rlm_merged.df <- colMeans(F_rlm_merged,na.rm = T)
    def <- c(def,rlm_merged.df)
  }
  
  if("rlm" %in% feature.family) {
    cat("computing glrm 2.5D features...\n")
    F_rlm_25D <-glrlmTexturalFeatures25D(res$voxelCube, n_grey=n_grey)
    F_rlm_25D <- do.call(data.frame,lapply(F_rlm_25D, function(x) replace(x, is.infinite(x),NA)))
    rlm_25D.df <- colMeans(F_rlm_25D,na.rm = T)
    def <- c(def,rlm_25D.df)
  }
  
  if("rlm" %in% feature.family) {
    cat("computing glrm 2.5D merged features...\n")
    F_rlm_25D_merged <-glrlmTexturalFeatures25Dmerged(res$voxelCube, n_grey=n_grey)
    F_rlm_25D_merged <- do.call(data.frame,lapply(F_rlm_25D_merged, function(x) replace(x, is.infinite(x),NA)))
    def <- c(def,F_rlm_25D_merged)
  }

  if("szm" %in% feature.family) {
    cat("computing szm features...\n")
    F_szm <- glszmTexturalFeatures(res$voxelCube, n_grey=n_grey)
    F_szm <- do.call(data.frame,lapply(F_szm, function(x) replace(x, is.infinite(x),NA)))
    szm.df <- colMeans(F_szm,na.rm = T)
    def <- c(def,szm.df)
  }
  
  if("szm" %in% feature.family) {
    cat("computing szm 2.5D features...\n")
    F_szm_25D <- glszmTexturalFeatures25D(res$voxelCube, n_grey=n_grey)
    F_szm_25D <- do.call(data.frame,lapply(F_szm_25D, function(x) replace(x, is.infinite(x),NA)))
    szm_25D.df <- colMeans(F_szm_25D,na.rm = T)
    def <- c(def,szm_25D.df)
  }

  if("fractal" %in% feature.family) {
    cat("computing fractal features...\n")
    logObj<-logHandler()
    logObj$handle(type="error", msg=c("Fractal features has to be extracted using fractalFeatures function. Please see the help for further information"))
    # F_fractal <- fractalFeatures(original)
    # fractal.df <- colMeans(F_fractal)
    # def <- c(def,fractal.df )
  }


  return(def)
}
