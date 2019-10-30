# Moddicom
Radiomics toolbox for R

[WORK IN PROGRESS]

## Table of contents

* [Installation (Windows)](#installation-windows "Goto Installation(Windows)")
* [Installation (Ubuntu)](#installation-ubuntu "Goto Installation(Ubuntu)")
* [Installation (Mac)](#installation-mac "Goto Installation(Mac)")
* [Usage](#usage "Goto Usage")
  * [Use case 1: Extract a ROI for a single patient DICOM series](#use-case-1-extract-a-roi-for-a-single-patient-dicom-series "")
  * [Use case 2: Extract features for a single patient DICOM series](#use-case-2-extract-features-for-a-single-patient-dicom-series "")
  * [Use case 3: Extract features for multiple patients DICOM series](#use-case-3-extract-features-for-multiple-patients-dicom-series "")
  * [Use case 4: Extract features from multiple patients DICOM series with LoG filter](# "")
  * [Use case 5: working with Analyize/NIFTI standard](# "")
 * [Contributing](#contributing "")
 * [Credits](#credits "")
 * [License](#license "")
  
  
 

## Installation (Windows)

1) Download and install R: https://cran.r-project.org/bin/windows/base/
2) Download and install Rstudio: https://www.rstudio.com/
3) Download and install RTools: https://cran.r-project.org/bin/windows/Rtools/ 
4) Download and unzip DCMTK libraries: https://dicom.offis.de/dcmtk.php.en (look for “DCMTK 3.6.4 - executable binaries” and choose you OS architecture)
5) Add the .\bin folder of the unzipped DCMTK libraries to the PATH system environment variable
6) Launch Rstudio and run the script to install all required R packages (install_packages.R)
7) Install moddicom from github

```library(devtools)```

```install_github("kbolab/moddicom")```

## Installation (Ubuntu)

1) Download and install R: https://cran.r-project.org/bin/windows/base/
2) Download and install Rstudio: https://www.rstudio.com/
3) Download and unzip DCMTK libraries: https://dicom.offis.de/dcmtk.php.en (look for “DCMTK 3.6.4 - executable binaries” and choose you OS architecture)
4) Launch Rstudio and run the script to install all required R packages (install_packages.R)
5) Install moddicom from github

```library(devtools)```

```install_github("kbolab/moddicom")```

## Installation (Mac)

1) Download and install R: https://cran.r-project.org/bin/windows/base/
2) Download and install Rstudio: https://www.rstudio.com/
3) Download and unzip DCMTK libraries: https://dicom.offis.de/dcmtk.php.en (look for “DCMTK 3.6.4 - executable binaries” and choose you OS architecture)
4) Launch Rstudio and run the script to install all required R packages (install_packages.R)
5) Install moddicom from github

```library(devtools)```

```install_github("kbolab/moddicom")```


## Usage

Requisites for a single patient folder:
The class geoLet allows to open a folder of a single patient. The folder must not have subfolders and it must cointain only DICOM files with extension .dcm
In the folder there must be only one series of images (CT or MRI) and one RTStruct; optionally there can be one RTDose and/or one RTPlan. The examples below can be run on the test images at https://github.com/kbolab/testImagesModdicom

### Use case 1: Extract a ROI from a single patient DICOM series

Launch Rstudio.

Load Moddicom library

```library(moddicom)```

Create a geoLet class instance.

```obj.geolet <- geoLet()```

Load the patient image series:

```obj.geolet$openDICOMFolder(pathToOpen= ".\testImages\pat1")```

moddicom will load one by one the images and file with the structures, and will align the geometries.

Typing:

```obj.geolet$getROIList()```

you will see the ROI names for that patient.

TO extract a particular ROI, for instance “GTV-1”, you just type

```roiVoxels <- obj.geolet$getROIVoxels(Structure = "GTV-1")```

This will return a “list” data type. Since for our puroposes we’re just interested in the ROI VoxelCube, we can extract that as follows: 

```GTV_voxelcube <- roiVoxels$masked.images$voxelCube```

Now is a 3D matrix with coordinates [x,y,z]
To visualize a particular slice (the fifth in this example):

```image(GTV_voxelcube[,,5], col = grey.colors(256))```

### Use case 2: Extract features from a single patient DICOM series

```library(moddicom)```

```folder <- ".\testImages\pat1"```

```roi <- "GTV-1"```

```features <- f.extractor.sing.par(path = path, ROIName = roi, feature.family=c("stat","morph","glcm","rlm","szm"))```

### Use case 3: Extract features from multiple patients DICOM series

Note: all patients must have a ROI with the same name (“GTV” in the example below)

```library(moddicom)```

```folder <- ".\testImages\"```

```roi <- "GTV-1"```

```Features <- f.extractor.sing.par(path = folder, ROIName = roi, feature.family=c("stat","morph","glcm","rlm","szm"),fileName = paste0("features",".Rdata"), forceRecalculus = FALSE)```

### Use case 4: Extract features from multiple patients DICOM series with LoG filter

### Use case 5: working with Analyize/NIFTI standard

## Contributing

## Credits

## License
