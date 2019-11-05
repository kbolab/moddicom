#' class for loading and presenting DICOM data
#'
#' @description  Create a report using as input a path containing directories with patient,
#'               \itemize{
#'               \item \code{setPath(), getPath(), pixelSpacing;} : to load a DICOM series into an geoLet object
#'               }
#'      
#' @examples \dontrun{
#' # first of all create an object aa 
#' aa <- RAD.scoutNIFTI(path = '/Users/davpas/Desktop/bis', '*metastatic*')
#' aa$setPath("/Users/davpas/Desktop/moddicom/moddicom/R/")
#' aa$pixelSpacing()
#'
#' }
#' @export
#' @useDynLib moddicom
#' @import rmarkdown
RAD.scoutNIFTI<- function (path,regexp_roi){
  
  # RAD.scoutNIFTI is a class that takes in input:
  # a path of a folder with all the patient directory 
  # a regular expression given to specify the ROI (metastatic, primary, ecc..)
  
  attr.G.def.RMark.dir<-''
  attr.G.def.output.report<-''
  attr.G.def.DICOM.pattern.ext<-''
  attr.G.def.report.fileName.scout<-''
  attr.G.def.report.fileName.sigma<-''
  attr.G.path.final<-''
  attr.G.fileName<-''
  #attr.G.regexp_roi<-regexp_roi
  attr.G.regexp_roi<-''
  


  # --------------------------------------------------
  # setPath
  # --------------------------------------------------
  setPath<-function(def.RMark.dir = NA, def.output.report = NA,
                    def.DICOM.pattern.ext = NA, def.report.fileName.scout=NA, def.report.fileName.sigma=NA,
                    fileName=NA) {

    if(!is.na(def.RMark.dir)) attr.G.def.RMark.dir<<-def.RMark.dir
    if(!is.na(def.output.report)) attr.G.def.output.report<<-def.output.report
    if(!is.na(def.DICOM.pattern.ext)) attr.G.def.DICOM.pattern.ext<<-def.DICOM.pattern.ext
    if(!is.na(def.report.fileName.scout)) attr.G.def.report.fileName.scout<<-def.report.fileName.scout
    if(!is.na(def.report.fileName.sigma)) attr.G.def.report.fileName.sigma<<-def.report.fileName.sigma
    if(!is.na(fileName)) attr.G.fileName<<-fileName

  }
  # --------------------------------------------------
  # getPath
  # --------------------------------------------------
  getPath<-function() {
    return(list("def.RMark.dir"=attr.G.def.RMark.dir,"def.output.report"=attr.G.def.output.report,"def.DICOM.pattern.ext"=attr.G.def.DICOM.pattern.ext,
                "def.report.fileName.scout"=attr.G.def.report.fileName.scout,
                "def.report.fileName.sigma"=attr.G.def.report.fileName.sigma, "fileName"=attr.G.fileName))
  }
  # --------------------------------------------------
  # pixelSpacing
  # --------------------------------------------------
  # This method generates the report using .Rmd file
  
  pixelSpacing <- function(){

    tmp_list_of_patient <- list()
    patList <- list.dirs(path = attr.G.path.final,recursive = FALSE)
    print(path)
    for(pat in list.files(attr.G.path.final)){
      tmp_check_files= list.files(path = paste(path,pat,sep = "/"), pattern = regexp_roi)
      print(regexp_roi)
      #browser()
      if (length(tmp_check_files)==1 || length(tmp_check_files)==2){
        obj.geoNIFTI<-geoNIFTI()
        tmp_list_of_patient[[pat]] <-obj.geoNIFTI$getROIVoxels(paste(path,pat,sep = "/"),regexp_roi)
        gc()
        rm(obj.geoNIFTI)
      }
    }
    pixelX <- c()
    pixelY <- c()
    SpessoreFetta<-c()
    lista_nomi<-names(tmp_list_of_patient)

    #browser()
    print(lista_nomi)
    for (i in tmp_list_of_patient){
      print(names(i))
      print(i$geometricalInformationOfImages$pixelSpacing)
      pixelX<-c(pixelX,i$geometricalInformationOfImages$pixelSpacing[1])
      pixelY<-c(pixelY,i$geometricalInformationOfImages$pixelSpacing[2])
      SpessoreFetta<-c(SpessoreFetta,i$geometricalInformationOfImages$SliceThickness)
    }
    
    # modifing patList including only directories that have files with regexp given in input
   patListNew<-c()
   #browser()
   for (i in names(tmp_list_of_patient)){
     print(i)
     folder_to_add_patList<-grep(i,patList, value=TRUE)
     if (nchar(folder_to_add_patList)>1){
       patListNew<-c(patListNew,folder_to_add_patList)
     }
   }
   patList<-patListNew
   #browser()
      
      rm(tmp_list_of_patient)
      gc()
    #}
    browser()
    dir.x.RMARKDOWN <- paste( c(attr.G.def.RMark.dir,"scoutPixelSpacingNIFTI.Rmd")  , collapse='')
    render(dir.x.RMARKDOWN, "pdf_document", output_file=attr.G.def.report.fileName.scout, output_dir=attr.G.def.output.report)


  }


  constructor<-function(path, regexp_roi) {
    attr.G.def.RMark.dir<<-'./'
    attr.G.def.output.report<<-'./'
    attr.G.def.DICOM.pattern.ext <<-'\\.dcm$'
    attr.G.def.report.fileName.scout<<-'scoutPixelSpacingNIFTI.pdf'
    #attr.G.def.report.fileName.sigma<<-'scoutSigmaNIF.pdf'
    attr.G.path.final<<-path
    attr.G.fileName<<-"tmp.f.extractor.pluri.par.RData"
    attr.G.regexp_roi<-regexp_roi
  }
  constructor(path = path, regexp_roi=regexp_roi)

  #return(list("pixelSpacing"=pixelSpacing,"getPath"=getPath,"setPath"=setPath,"scoutSigma"=scoutSigma))
  return(list("pixelSpacing"=pixelSpacing,"getPath"=getPath,"setPath"=setPath)) #,"scoutSigma"=scoutSigma))
  
}
