#' @export

f.extractor.sing.par.nifti<- function(path, regulexp_roi ,
                                feature.family=c("stat","morph","glcm","rlm","szm"),
                                filterPipeline=list(), interpolate = FALSE, threshold = .001, discretize="", bin.size = 25, bin.number=25,
                                fileName = "tmp.f.extractor.sing.par.RData" ,
                                forceRecalculus = TRUE) {


  patList <- list.dirs(path = path,recursive = FALSE)
  matrice <- c()
  rigaFeatures <- c()
  patientsPath <- c()
  
  if(forceRecalculus==FALSE) {
    if(file.exists(fileName)==TRUE) {
      load(file = fileName)  
    }   
  } else {
    if(file.exists(fileName)==TRUE) {
      file.remove(fileName)
    }    
  }
  
  patListNew <- c()
  for(pat in patList){
    ll <- list.files(path = pat, pattern = regulexp_roi)
    if(length(ll) == 1 | length(ll) == 2){
      patListNew <- c(patListNew,pat)
    }
  }
  

  for( patID in patListNew) {

    if(!(patID %in% matrice[,1])) {
      cat(patID)
      cat("\n")
      
      riga <- computeFeatures.geoNifti( path.nifti = patID, ROIName = regulexp_roi ,
                                      feature.family=feature.family, filterPipeline=filterPipeline,
                                      discretize=discretize, bin.size=bin.size,bin.number=bin.number) 
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
      #ROIName <- temp.ROIName
    }
    
  }
  
  matrice$patID <- as.character(matrice$patID)
  return(matrice)
}

# qui entrano path che produrranno voxelcube non vuoti (hanno la roi richiesta)
computeFeatures.geoNifti<- function( path.nifti, ROIName , feature.family=c("stat","morph","glcm","rlm","szm"), filterPipeline=c(),
                                     threshold = .1, discretize="", bin.size = "", bin.number="", n_grey = 100) {

  ogg.geoNIFTI <- geoNIFTI()
  ROI <- ogg.geoNIFTI$getROIVoxels(path.nifti,"*met*")
  
  px <- ROI$geometricalInformationOfImages$pixelSpacing[1]
  py <- ROI$geometricalInformationOfImages$pixelSpacing[2]
  pz <- as.numeric(ROI$geometricalInformationOfImages$SliceThickness)
  
  
  
  
  #filtro se richiesto (non implementato per nifti)
  if(length(filterPipeline)!=0) {
    res <- FIL.2D.conv.Filter.geoLet(obj.geoLet = obj.geoLet,ROIName = ROIName,
                                     filter.pipeline = filterPipeline,
                                     scaleFactor = 'space', px = px, py = py)

  } else {
    res <- ROI$masked.images
  }
  
  naVALUE <- as.integer(min(voxelCube = res$voxelCube,na.rm = T) - 10)
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
    morph.df <- morphologicalFeatures(original,px=px,py=py,pz=pz)
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
