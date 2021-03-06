% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_fractalFeatures.R
\name{fractal.Extr.threshold}
\alias{fractal.Extr.threshold}
\title{Fractal features extractor from a folder based on threshold filter}
\usage{
fractal.Extr.threshold(path, ROIName, fstPerc, lstPerc, ThDown, ThUp,
  ThStep, filename)
}
\arguments{
\item{path}{The path to be browsed for folders containing DICOM images and Structure Set of patients}

\item{ROIName}{A \code{character} object containing the name(s) (as vector) of ROIs to be analyzed. If there is more than one ROI (multiple references) or no ROI is 
listed the analysis is stopped}

\item{fstPerc}{A \code{numeric} value for normalized the image. It is the first percentile.}

\item{lstPerc}{A \code{numeric} value for normalized the image. It is the last percentile.}

\item{ThDown}{A \code{numeric} value for masked image filter. It is the minimum value of the threshold}

\item{ThUp}{A \code{numeric} value for masked image filter. It is the maximum value of the threshold}

\item{filename}{The Rdata file name saved locally during the running, by default 'Fractals.Rdata'.}
}
\description{
This function extracts fractal features from a folder using the threshold filter.
}
\examples{
\dontrun{
prova <- fractal.Extr.threshold(path = './Images', ROIName = 'GTV', fstPerc = 1, lstPerc = 99, 
ThDown = 0, ThUp = 100, ThStep = 50, filename = 'Fractals_Lung.Rdata')
}
}
