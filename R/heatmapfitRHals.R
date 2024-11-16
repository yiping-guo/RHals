#' Define a generic method for drawing a heat-map of Residuals for RH Models
#' @export

heatmap <- function(object, ...) {
  UseMethod("heatmap")
}

#' Draw a Heat-map of Residuals for RH Models
#'
#' This function generates a heat-map (color map) of the residuals 
#' from the fitted model object with class \code{fitRHals} returned by \code{RHfit}. 
#' The residuals are defined as:
#' \deqn{d_{x,t} = \log m_{x,t} - \log \hat{m}_{x,t}}
#' where \eqn{m_{x,t}} are the observed mortality rates, 
#' and \eqn{\hat{m}_{x,t}} are the fitted values from the model.
#' Positive residuals are shown in blue, while negative residuals are shown in red.
#'
#' @param object An \code{fitRHals} object returned by \code{RHfit()}.
#'
#' @return This function draws a heat-map (color map) of the residuals using the 
#' \code{plot()} method from the \code{StMoMo} package.
#'
#' @examples
#' # Load example data
#' Data <- load_EWData(ages = 60:89, years = 1960:2009, series = "Male")
#'
#' # Fit the RH model using the RHfit function
#' fit_rh <- RHfit(Data, ages = 60:89, years = 1960:2009)
#'
#' # Generate the heatmap of residuals
#' RHheatmap(fit_rh)
#'
#' @export
heatmap.fitRHals <- function(object) {
  # Compute residuals: log observed - log fitted mortality rates
  res <- fitted(object) - log(object$mxt)
  
  # Assign row and column names for better visualization
  colnames(res) <- object$years
  rownames(res) <- object$ages
  
  # Create a residual structure compatible with StMoMo plot function
  X <- structure(
    list(residuals = res, ages = object$ages, years = object$years),
    class = "resStMoMo"
  )
  
  # Plot the residuals as a heatmap (color map)
  plot(X, type = "colourmap")
}
?heatmap
