% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_EWData.R
\name{load_EWData}
\alias{load_EWData}
\title{Load EW Mortality Data}
\usage{
load_EWData(ages, years, series = c("Male", "Female", "Total"))
}
\arguments{
\item{ages}{A vector specifying the age range to subset. Allowed range is 0 to 110.}

\item{years}{A vector specifying the year range to subset. Allowed range is 1841 to 2021.}

\item{series}{A string specifying which series to select. Allowed values are "Female", "Male", "Total".}
}
\value{
A matrix containing the subsetted data for the specified series.
}
\description{
This function loads the Mx_1x1_EW data from the package, subsets it based on age and year range, and selects the specified series (Female, Male, Total).
}
