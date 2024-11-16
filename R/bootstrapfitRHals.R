#' Bootstrap a fitted Generalized Renshaw-Haberman Model
#' 
#' A generic \code{S3} method that uses residual bootstrap to produce bootstrap  
#' parameters described in Guo and Li (2024) to account for parameter uncertainty.
#' @param object An \code{fitRHals} object returned by \code{RHalsfitting()}.
#' @param nBoot number of bootstrap samples to produce.
#' @return A list with class \code{"bootStMoMo"} with components:
#' 
#' \item{bootParameters}{ a list of of length \code{nBoot} with the fitted 
#' parameters for each bootstrap replication.} 
#' \item{model}{ the model fit that has been bootstrapped.}
#' 
#' @examples
#' # Load example data
#' Data <- load_EWData(ages = 60:89, years = 1960:2009, series = "Male")
#'
#' # Fit the RH model using the RHfit function
#' fit_rh <- RHfit(Data, ages = 60:89, years = 1960:2009)
#' 
#' # Perform 1000 residual bootstraps
#' boot_rh <- bootstrap(fit_rh, nBoot = 1000)
#' 
#' # Extract the standard error of bx
#' seBoot(boot_rh)$bx
#' 
#' # Draw the fan charts
#' plot(boot_rh, nCol = 3)
#' @export

bootstrap.fitRHals<- function(object, nBoot = 1,...) {
  p <- length(object$ages)
  n <- length(object$years)
  
  # Calculate residual of the log rates
  res <- log(object$mxt) - fitted(object)
  bootSamples <- lapply(1:nBoot, function(i) {
    exp(matrix(sample(res, size = p * n, replace = TRUE), nrow = p, ncol = n) + fitted(object))
  })
  
  #Fit the model to each of bootstrapped sample
  refit <- function(bootdata) {
    RHalsfitting(bootdata, ages = object$ages, years = object$years, lc = object$lc,
             m = object$m, const.b0x = object$const.b0x, approxConst = object$approxConst)
  }
  bootmodels <- lapply(bootSamples, FUN = refit)  
  
  boot <- create_bootRHals()
  boot$bootParameters <- bootmodels
  boot$model <- object
  
  return(boot)
}


#' Define a generic method for drawing a calculating the standard error
#' of bootstrapped parameters
#' @export
seBoot <- function(object, ...) {
  UseMethod("seBoot")
}

#' Calculate Standard Errors of Bootstrapped Parameters
#'
#' This function calculates the standard error (SE) of the bootstrapped parameters 
#' from an \code{bootRHals} object.
#'
#' @param object An \code{bootRHals} object.
#' 
#' @return A list containing the calculated standard error of the parameters:
#' \code{ax}, \code{bx}, \code{kt}, \code{b0x}, and \code{gc}, each retaining their original structure.
#' @export

seBoot.bootRHals <- function(object, ...) {
  nBoot <- length(object$bootParameters)  # Number of bootstraps
  
  # Helper function to extract parameters from each bootstrap
  extract_parameter <- function(param_name) {
    lapply(object$bootParameters, function(model) model[[param_name]])
  }
  
  # Helper function to calculate the standard error for each element
  calculate_se <- function(param_values) {
    # Stack all bootstrapped versions of the parameter into a 3D array
    param_array <- simplify2array(param_values)  # Creates an array
    
    # Calculate the standard error: sqrt(var / nBoot)
    apply(param_array, MARGIN = c(1, 2), function(x) {
      sqrt(var(x) / nBoot)
    })
  }
  
  # Calculate the SE for each parameter
  se_ax <- calculate_se(extract_parameter("ax"))
  se_bx <- calculate_se(extract_parameter("bx"))
  se_kt <- calculate_se(extract_parameter("kt"))
  se_b0x <- calculate_se(extract_parameter("b0x"))
  se_gc <- calculate_se(extract_parameter("gc"))
  
  # Return a list of SE results, retaining the original structure
  return(list(
    ax = se_ax,
    bx = se_bx,
    kt = se_kt,
    b0x = se_b0x,
    gc = se_gc
  ))
}