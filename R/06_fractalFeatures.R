#' Fractal features
#'
#' @description  This function extracts fractal features using the box counting technique.
#' @param imgObj The imgObj get as input the GTV object extracted using getROIVoxels method.
#' @param fstPerc A \code{numeric} value for normalized the image. It is the first percentile.
#' @param lstPerc A \code{numeric} value for normalized the image. It is the last percentile.
#' @param ThDown  A \code{numeric} value for masked image filter. It is the minimum value of the threshold
#' @param ThUp A \code{numeric} value for masked image filter. It is the maximum value of the threshold
#' @export
#' @examples \dontrun{
#' # Create an instance of geoLet object and load a case
#' obj <- geoLet()
#' obj$openDICOMFolder(pathToOpen = './patient1')
#' 
#' # get the ROI
#' GTV <- obj$getROIVoxels(Structure = 'GTV')
#' 
#' # extract the fractal feature
#'  Fractal <- fractalFeatures(imgObj = GTV, fstPerc = 1, lstPerc = 99, ThDown = 0, ThUp = 100)
#' 
#' }
fractalFeatures <- function(imgObj,fstPerc=0, lstPerc=100, ThDown=0, ThUp=100, onepxwide=T) {
  
  if( length(dim(imgObj)) == 2 ) {
    n.imgObj <- array(0, dim = c(dim(imgObj),3) )
    n.imgObj[,,2] <- imgObj
    imgObj <- n.imgObj
  }
  
  a <- NormerObj(imgObj = imgObj, fstPerc=fstPerc, lstPerc=lstPerc)
  b <- ThresholderObj(imgObj = a, ThDown=ThDown, ThUp=ThUp)
  ciuccia <- DFCalculator(imgObj = b , onepxwide=onepxwide)
  report <-c()
  res.mean <- mean(ciuccia$FractalDimension$dimF.patient,na.rm = TRUE)
  res.median <- median(ciuccia$FractalDimension$dimF.patient,na.rm = TRUE)
  res.min <- min(ciuccia$FractalDimension$dimF.patient,na.rm = TRUE)
  res.max <- max(ciuccia$FractalDimension$dimF.patient,na.rm = TRUE)
  res.sd <- sd(ciuccia$FractalDimension$dimF.patient,na.rm = TRUE)
  report <- rbind(report, c(res.mean,res.median,res.min, res.max,res.sd))
  colnames(report)<-c("meanFD","medianFD","minFD","maxFD","sdFD")
  risultato <- list("report"=report,"tutte_slice"=ciuccia$FractalDimension, "cubi_tutte_slice"=ciuccia$cubi_boxCounting, "sottpop"=a)
  return(risultato)
}

#' @useDynLib moddicom
DFCalculator<-function(imgObj, onepxwide=T){
  #DFCalculator prende in ingresso un oggetto GeoLet. L'immagine è gia stata passata da Thresholder
  #Questa funzione restituisce lo stesso Geolet, ma con dentro i VoxelCube le immagini normalizzate
  
  #Caricamento delle library e del file C che serve per il boxcounting
  # library(moddicomV2)
  # dyn.load("/home/davide/Scrivania/FractalRadiomics/fractal2D.so")
  
  #Inizializzazione Variabili
  ROIVoxelData<-imgObj
  risultati.tot<-list(); j<-1
  ID <- c(); dimF <- c(); slice <- c()

  #Interlock che valuta se lo studio è binarizzato. Se non lo è non consente il calcolo di DF
  # if (summary(ROIVoxelData$voxelCube)[6]!=1)
  if (summary(ROIVoxelData)[6]!=1)
  {
    Attention<-"DF Calculation not possible: Study not binarized "
    return( list("result"=NA,"error"=TRUE, "nonBinarized"=patient)   )
  }
  
  objS<-services();
  
  # voxelData.ready<-ROIVoxelData$voxelCube
  voxelData.ready<-ROIVoxelData
  dimF.patient<-c()
  slice.Patient<-c()

  for (q in 1:(dim(voxelData.ready)[3])){
    risultati<-list()
    i<-1
    
    # Estrazione del contorno dell'immagine se è richiesta la condizione onepxwide
    
    if (onepxwide==T) {
      
      arrayRM2_3d<-voxelData.ready[,,q]
      
      #########Matrice per erosion######
      arrayRM_erosion<-voxelData.ready[,,q]
      arrayRM_erosion[which(arrayRM_erosion==0)]<-NA
      
      obj<-services()
      numeroPixelNonZero <- sum(!is.na(arrayRM_erosion))
      
      if(numeroPixelNonZero>0)
      {
        Mask_NormFetta<-obj$applyErosion(ROIVoxelData = arrayRM_erosion,margin.x = 1,margin.y = 1)
      }
      else Mask_NormFetta <- arrayRM_erosion
      
      Mask_NormFetta[is.na(Mask_NormFetta)] <- 0
      Bordo<-arrayRM_erosion-Mask_NormFetta
      Bordo[Bordo==0]<-NA
      #Sostituisco i Na con 0
      arrayRM2_3d<-Bordo
    }
    
    if (onepxwide==F) {
      arrayRM2_3d<-voxelData.ready[,,q]
      arrayRM2_3d[which(arrayRM2_3d==0)]<-NA
    }
    # arrayRM2_3d è l'immagine sulla quale si calcola il box-counting: può essere piena o solo il contorno
    #calcolo degli estremi su cui deve lavorare il box-counting
    dim1<-dim(arrayRM2_3d)[1]
    dim2<-dim(arrayRM2_3d)[2]
    if(sum(!is.na(arrayRM2_3d))!=0) {
      estremi<-as.data.frame(which(arrayRM2_3d!=0,useNames = dim1,dim2))
      min.row<-min(estremi$row);      max.row<-max(estremi$row)
      min.col<-min(estremi$col);      max.col<-max(estremi$col)
      dim1.r<-(max.row-min.row)+1;    dim2.r<-(max.col-min.col)+1
      arrayRM2_3d.r<-arrayRM2_3d[min.row:max.row,min.col:max.col]
      
      #Determinazione delle dimensioni iniziali del Cubo locale che spazzola la slice
      dimCubo<-min(dim1.r,dim2.r)
      NumGrigi<-1
      iter_ux<-seq(from=1, to=((dimCubo)-1))
      
      ########For in cui calcola il numero di Grigi al variare della dim del cubo scritto in C #########
      
      for (ux in iter_ux)
      {
        # Metto il contatore dei Grigi a zero, lo devo fissare ogni volta che entro nel for
        N<-0
        
        #Inizializzo il Quadrato che spazzola la matrice con degli zeri, poi lo riassegno localmente nel programma in C
        QuadratoFrattale<-rep(0,((ux*2)+1)*((ux*2)+1))
        t<-length(QuadratoFrattale)
        
        #Sistemo la matrice arrayRMborders così da non perdermi mai i bordi
        arrayRMborders<-matrix(data=0,nrow = dim1.r+(ux*2)+1, ncol = dim2.r+(ux*2)+1)
        dim1.b<-dim(arrayRMborders)[1]; dim2.b<-dim(arrayRMborders)[2]
        arrayRMborders[1:dim1.r,1:dim2.r]<-arrayRM2_3d.r
        arrayRMborders[is.na(arrayRMborders)] <- 0
        
        #Invoco la funzione in C e rimetto i risultati di sta funzione nel vettore C_Results
        C_Results<- .C("fractal2D", as.double (arrayRMborders), as.integer (dim1.b), as.integer (dim2.b),
                       as.double(QuadratoFrattale), as.integer(ux), as.integer(ux), as.integer(N))
        
        #Popolo la lista che si chiama risultati, comandata dal contatore i
        #names(risultati)<-c("Threshold","LatoCubo","NumGrigi")
        risultati[[i]]<-c(2*(C_Results[[6]])+1,C_Results[[7]])
        i=i+1
        #Popolo una lista globale
        risultati.tot[[j]]<-c(q,2*(C_Results[[6]])+1,C_Results[[7]])
        j=j+1
      }
      
      
      #Trasferisco i risultati della lista risultati in output dove ho il report dei grigi per una soglia al variare della dimensione cubo
      output<-matrix(data = unlist(risultati),nrow = length(risultati),ncol = 2,byrow = TRUE)
      output.tot<-matrix(data = unlist(risultati.tot),nrow = length(risultati.tot),ncol = 3,byrow = TRUE)
      colnames(output.tot)<-c("Fetta","DimCubo","NGrigi")
      #CALCOLO DIMENSIONE FRATTALE
      #preprocessing per levare gli uni e i valori ripetuti che fanno schifio
      output[which(output==0)]<-1
      vect<-duplicated(x = output[,2])
      index<-which(vect==FALSE)
      output2<-output[index,1:2]
      
      # #la condizione chiede di avere almeno 4 punti per calcolare la dimF
      # if (length(output2)>9) {
        # fit<-(lm(formula = log(output2[,2]) ~ log(output2[,1])))
        # y<-log(output2[,2])
        # x<-log(output2[,1])
        
        # commentato il 12/11/2018 per inserire vincoli sul plot
        # # con i ripetuti
        # fit<-(lm(formula = log(output[,2]) ~ log(output[,1])))
        # y<-log(output[,2])
        # x<-log(output[,1])
        
        # calcolo la dimensione frattale con i vincoli: 4 < ln(N) < 100
        ind_lower <- which(output[,2] > 4)
        ind_upper <- which(output[,2] < 100)
        ind <- ind_lower[which(ind_lower %in% ind_upper)]
        output_new <- output[ind,]
        # la condizione chiede di avere almeno 4 punti per calcolare la dimF
        if (length(output_new)>9 && length(unique(output_new[,2])) >1 ) {
        fit<-(lm(formula = log(output_new[,2]) ~ log(output_new[,1])))
        ss <- summary(fit)
        y<-log(output_new[,2])
        x<-log(output_new[,1])
        FractalDim<-(-fit$coefficients[[2]])
        if (FractalDim <1)  {FractalDim<-NaN}
        #if (ss$r.squared < .99) {FractalDim<-NaN}
      }
      else  {FractalDim<-NaN}
    } else  {FractalDim<-NaN}
    
    #ID <- c(ID, Paziente)
    #dimF <- c(dimF, FractalDim)
    #slice <- c(slice, q)
    dimF.patient<-c(dimF.patient,FractalDim)
    slice.Patient<-c(slice.Patient,q)
  }
  ReportData.Paziente<-as.data.frame(cbind(slice.Patient,dimF.patient))
  ROIVoxelData$FractalDimension<-ReportData.Paziente
  ind_lower <- which(output.tot[,3] > 4)
  ind_upper <- which(output.tot[,3] < 100)
  ind <- ind_lower[which(ind_lower %in% ind_upper)]
  output.tot_new <- output.tot[ind,]
  ROIVoxelData$cubi_boxCounting <- output.tot_new
  return(ROIVoxelData)
}

# L'Output è una lista dove sono inseriti anche i valori di Dimensione Frattale
#'  Funzione 2
#'
#' @description  ThresholderObj prende in ingresso un geoLet che ha già caricato l'immagine ed  estratto la ROI
#' Questa funzione restituisce lo stesso Geolet, ma con dentro i VoxelCube le immagini binarizzate
ThresholderObj<-function(imgObj, ThDown=0, ThUp=100){
  
  #objS<-services();
  ROIVoxelData<-imgObj
  # voxelData.ready<-ROIVoxelData$masked.images$voxelCube
  voxelData.ready<-ROIVoxelData
  
  for (q in 1:dim(voxelData.ready)[3])
  {
    ####################Binarizzazione da ThDown a ThUp#################################
    arrayRM2_3d<-voxelData.ready[,,q]
    arrayRM2_3d[which(arrayRM2_3d<ThDown)]<-NA
    arrayRM2_3d[which(arrayRM2_3d>ThUp)]<-NA
    arrayRM2_3d[which(arrayRM2_3d<=ThUp)]<-1
    
    #Sostituisco i NaN con zero
    arrayRM2_3d[is.na(arrayRM2_3d)] <- 0
    
    # ROIVoxelData$masked.images$voxelCube[,,q]<-arrayRM2_3d
    ROIVoxelData[,,q] <- arrayRM2_3d
  }
  return(ROIVoxelData)
}

#'  Funzione 3
#'
#' @description  NormerObj prende in ingresso un geoLet che ha già caricato l'immagine ed estratto la ROI
#' Questa funzione restituisce la stessa mButoVoxelList, ma con dentro i VoxelCube le immagini normalizzate
NormerObj<-function(imgObj, fstPerc=0, lstPerc=100){
  #objS<-services();
  
  ROIVoxelData<-imgObj
  # prendi i voxelData
  # voxelData.ready<-imgObj$masked.images$voxelCube
  voxelData.ready<-imgObj
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
  # ROIVoxelData$masked.images$voxelCube<-NormStudio
  ROIVoxelData<-NormStudio
  return(ROIVoxelData)
}


#' Fractal features extractor from a folder based on threshold filter
#'
#' @description  This function extracts fractal features from a folder using the threshold filter.
#' @param path The path to be browsed for folders containing DICOM images and Structure Set of patients
#' @param ROIName A \code{character} object containing the name(s) (as vector) of ROIs to be analyzed. If there is more than one ROI (multiple references) or no ROI is 
#' listed the analysis is stopped
#' @param fstPerc A \code{numeric} value for normalized the image. It is the first percentile.
#' @param lstPerc A \code{numeric} value for normalized the image. It is the last percentile.
#' @param ThDown  A \code{numeric} value for masked image filter. It is the minimum value of the threshold
#' @param ThUp A \code{numeric} value for masked image filter. It is the maximum value of the threshold
#' @param filename The Rdata file name saved locally during the running, by default 'Fractals.Rdata'.
#' @export
#' @examples \dontrun{
#' prova <- fractal.Extr.threshold(path = './Images', ROIName = 'GTV', fstPerc = 1, lstPerc = 99, 
#' ThDown = 0, ThUp = 100, ThStep = 50, filename = 'Fractals_Lung.Rdata')
#' }
fractal.Extr.threshold <- function(path, ROIName, fstPerc, lstPerc, ThDown,
                              ThUp, ThStep, filename){
  Report<-list()
  ReportGlobale<-data.frame()
  row<-c()
  patList <- list.dirs(path = path,recursive = FALSE)
  for (i in patList){
    
    obj<-geoLet()
    obj$openDICOMFolder(pathToOpen = i)
    geoLetVoxelList<-obj$getROIVoxels(Structure = ROIName)
    
    Filtrate<-Select(geoLetVoxelList = geoLetVoxelList, fstPerc = fstPerc, 
                     lstPerc = lstPerc, ThDown = ThDown, ThUp = ThUp, 
                     ThStep = ThStep )
    
    vect<-c()
    nomi<-c()
    row<-c(row,i)
    for (j in names(Filtrate)){
      a<-fractalFeatures(imgObj = Filtrate[[j]]$masked.images$voxelCube)
      
      vect<-c(vect,a$report)
      nomi<-c(nomi,paste(colnames(a$report),j))
    }
    Report[[i]]=vect
    #ReportFrame=as.data.frame(t(Report))
    ReportGlobale=rbind(ReportGlobale, Report[[i]])
  }
  colnames(ReportGlobale)<-nomi
  
  # rownames(ReportGlobale)<-row
  save(ReportGlobale,file = filename)
  ReportGlobale <- cbind(row, ReportGlobale)
  names(ReportGlobale)[1] <- "patID"
  ReportGlobale$patID <- as.character(ReportGlobale$patID) 
  
  return(ReportGlobale);
}
