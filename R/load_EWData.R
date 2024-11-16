#' Load EW Mortality Data
#'
#' This function loads the Mx_1x1_EW data from the package, subsets it based on age and year range, and selects the specified series (Female, Male, Total).
#'
#' @param ages A vector specifying the age range to subset. Allowed range is 0 to 110.
#' @param years A vector specifying the year range to subset. Allowed range is 1841 to 2021.
#' @param series A string specifying which series to select. Allowed values are "Female", "Male", "Total".
#' @return A matrix containing the subsetted data for the specified series.
#' 
#' @examples
#' # Load EW male mortality data with ages between 60 and 89, and years between 1960 and 2009.
#' Data <- load_EWData(ages = 60:89, years = 1960:2009, series = "Male")
#' @export
load_EWData <- function(ages, years, series = c("Male", "Female", "Total")) {
  # Check ages
  if (any(ages < 0 | ages > 110)) {
    stop("Error: Ages must be in the range 0 to 110.")
  }

  # Check years
  if (any(years < 1933 | years > 2021)) {
    stop("Error: Years must be in the range 1933 to 2021.")
  }

  # Check series
  if (!series %in% c("Female", "Male", "Total")) {
    stop('Error: Series must be one of "Female", "Male", "Total".')
  }

  # Load the data file
  data_file <- system.file("data", "Mx_1x1_EW.txt", package = "RHals")
  data <- read.table(data_file, skip = 2, header = TRUE)

  # Convert Age to numeric and filter out "110+" rows
  data$Age <- as.numeric(gsub("\\+", "", data$Age))
  data <- data[data$Age <= 110, ]

  # Subset the data based on ages and years
  data <- data[data$Age %in% ages & data$Year %in% years, ]
  data <- na.omit(data)

  # Convert columns 3 to 5 to numeric
  for (i in 3:5) {
    data[, i] <- as.numeric(data[, i])
  }

  # Select the specified series and return as a matrix
  if (series == "Female") {
    return(matrix(data[, 3], nrow = length(unique(data$Age))))
  } else if (series == "Male") {
    return(matrix(data[, 4], nrow = length(unique(data$Age))))
  } else if (series == "Total") {
    return(matrix(data[, 5], nrow = length(unique(data$Age))))
  }
}
