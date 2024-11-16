#' Fit a Generalized Renshaw-Haberman Model using Least Squares
#'
#' This function fits a RH mortality model using the least squares method.
#' For the RH model, the fitting is via alternating least squares.
#' Updating the cohort effects at each iteration uses an iterative SVD.
#' The model is defined as:
#' \deqn{
#' \log m_{x,t} = a_x + \sum_{i=1}^{m}b^{(i)}_x k^{(i)}_t + b_x^{(0)}\gamma_{t-x}
#' }
#'
#' @param data A matrix of central mortality rates \eqn{m_{x,t}}, with rows representing ages and columns representing years.
#' @param ages A vector of ages corresponding to the rows of the mortality matrix.
#' @param years A vector of years corresponding to the columns of the mortality matrix.
#' @param m A positive integer corresponding to the number of age-period terms.
#' @param lc A logical value controlling whether the cohort effects are not included. If \code{TRUE}, a Lee-Carter model with \eqn{m} age-period terms will be fitted using the order-\eqn{m} singular value decomposition (SVD).
#' @param const.b0x A logical value controlling whether \eqn{b_x^{(0)}=1}. If \code{TRUE}, an H1 model will be fitted. Only effective when \code{lc = FALSE}.
#' @param approxConst A logical value defining whether the approximate identifiability constraint \eqn{\sum (s-\bar{s})\gamma_s=0} is applied. If \code{TRUE}, a Lagrange multiplier method proposed by Guo and Li (2024) will be implemented. Only effective when \code{lc = FALSE}.
#' @param maxits A positive integer corresponding to the maximum number of the main iterations. Only effective when \code{lc = FALSE}.
#' @param tol A small positive number corresponding to the tolerance level of the main iterations. Only effective when \code{lc = FALSE}.
#' @return An \code{fitRHals} object containing the following components:
#' \item{\code{model}}{An \code{StMoMo} object containing the fitted model.}
#' \item{\code{ax}}{A vector with the fitted values of the static age function \eqn{a_x}.}
#' \item{\code{bx}}{A matrix with the fitted values of the age-period function \eqn{b_x}.}
#' \item{\code{kt}}{A matrix with the fitted values of the period indexes \eqn{k_t}.}
#' \item{\code{b0x}}{A vector with the fitted values of the age-cohort function \eqn{b_x^{0}}.}
#' \item{\code{gc}}{A vector with the fitted values of the cohort indexes \eqn{\gamma_{t-x}}.}
#' \item{\code{ssr}}{The sum of squares error of the fitted model.}
#' \item{\code{npar}}{A positive integer representing the effective number of parameters in the model.}
#' \item{\code{loglik}}{The Gaussian log-likelihood of the fitted model.}
#' \item{\code{AIC}}{The Akaike's Information Criterion (AIC) of the fitted model.}
#' \item{\code{BIC}}{The Bayesian Information Criterion (BIC) of the fitted model.}
#' \item{\code{ages}}{A vector of ages of the data set.}
#' \item{\code{years}}{A vector of ages of the data set.}
#' \item{\code{cohorts}}{A vector of cohorts of the data set.}
#' \item{\code{mxt}}{The input mortality matrix.}
#' \item{\code{m}}{The input parameter \eqn{m}.}
#' \item{\code{lc}}{The input logical variable \code{lc}.}
#' \item{\code{const.b0x}}{The input logical variable \code{const.b0x}.}
#' \item{\code{approxConst}}{The input logical variable \code{approxConst}.}
#' 
#' @examples
#' # Example data
#' Data <- load_EWData(ages = 60:89, years = 1960:2009, series = "Male")
#' 
#' # Fit the RH model using the least squares method
#' fit_rh <- RHfit(Data, 60:89, 1960:2009)
#' 
#' # Extract the estimated gc
#' fit_rh$gc
#' 
#' # Plot the parameters estimates
#' plot(fit_rh$model)
#' 
#' # Other models
#' 
#' # Fit a Lee-Carter model
#' fit_lc <- RHfit(Data, 60:89, 1960:2009, lc = TRUE)
#' # Fit a RH model with two age-period terms
#' fit_rh2 <- RHfit(Data, 60:89, 1960:2009, m = 2)
#' # Fit an H1 model
#' fit_h1 <- RHfit(Data, 60:89, 1960:2009,const.b0x == TRUE)
#' # Fit an H1 model with the additional identification constraint
#' fit_h1_appro <- RHfit(Data, 60:89, 1960:2009,const.b0x == TRUE,approxConst == TRUE)
#' 
#' # Extract the effective number of parameters, the log-likelihood, AIC and BIC.
#' fit_rh$npar
#' fit_rh$loglik
#' fit_rh$AIC
#' fit_rh$BIC
#' 
#' # Forecasting -- fitted RH model as the example
#' # Forecasting window = 10 years
#' # Multivaraite drift random walk to model gc, and ARIMA(1,1,0) to model gc.
#' for_rh <- forecast(fit_rh, h = 10, gc.order = c(1, 1, 0)) 
#' # Fan chart for project kt
#' plot(for_rh, only.kt = TRUE)
#' @export
RHfit <- function(data, ages, years, m = 1, lc = FALSE, const.b0x = FALSE, 
                     approxConst = FALSE, maxits=10000, tol=1e-6) {
  # Check if mortality is specified and is a numeric matrix
  if (missing(data) || is.null(data)) {
    stop("Error: `data` is not specified.")
  } else if (!is.matrix(data) || !is.numeric(data)) {
    stop("Error: `data` must be a numeric matrix.")
  }

  # Check if ages are specified and are a numerical vector
  if (missing(ages) || is.null(ages)) {
    stop("Error: Ages are not specified.")
  } else if (!is.numeric(ages)) {
    stop("Error: Ages must be a numerical vector.")
  }

  # Check if years are specified and are a numerical vector
  if (missing(years) || is.null(years)) {
    stop("Error: Years are not specified.")
  } else if (!is.numeric(years)) {
    stop("Error: Years must be a numerical vector.")
  }

  # Check compatibility of ages and years with the dimensions of the mortality matrix.
  if (length(ages) != nrow(data)) {
    stop("Error: Length of ages vector does not match the number of rows in data.")
  }

  if (length(years) != ncol(data)) {
    stop("Error: Length of years vector does not match the number of columns in data.")
  }
  
  if (approxConst && !const.b0x) {
    stop("Error: `approxConst` not yet available for the non-parametric age-cohort terms.")
  }

  # Re-define ages and years for simplicity
  # Also use p and x to denote the corresponding dimensions
  x <- ages
  p <- length(x)
  t <- years
  n <- length(t)
  # Define and calculate the length for cohorts
  c <- (t[1] - x[p]):(t[n] - x[1])
  nc <- length(c)

  # Define and initialize bx and kt
  bx <- array(0, dim = c(p, m), dimnames = list(x, 1:m))
  kt <- array(0, dim = c(m, n), dimnames = list(1:m, t))

  # Define and initialize b^0_x
  if (lc == TRUE) {
    b0x <- array(0, dim = p, dimnames = list(x))
  } else {
    b0x <- array(1, dim = p, dimnames = list(x))
  }
  

  # Define and initialize gamma_c
  gc <- array(0, dim = nc, dimnames = list(c))

  # Log-transform the mortality rates
  Y <- log(data)

  # Define and initialize ax
  ax<- rowMeans(Y)
  
  # Define and Initial Ya_hat and Yc_hat 
  Ya_hat <- bx %*% kt
  Yc_hat <- Convert_AC_AT(b0x %*% t(gc))
  
  ### Fit the Lee-Carter model (when lc = TRUE)
  if (lc == TRUE){
    svd_Y <- svd(Y-ax,m,m)
    bx <- sweep(svd_Y$u, 2, colSums(svd_Y$u), FUN = "/")
    if (m == 1){
      kt <- sweep(svd_Y$d[1] * t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")
    } else if (m > 1){
      kt <- sweep(diag(svd_Y$d[1:m]) %*% t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")
    }
    Ya_hat <- bx %*% kt
  } else {
    ### Fit the Renshaw-Haberman model (when lc = FALSE)
    
    # Initialize estimating errors
    rel_err <- 1 # Relative error
    iter <- 0 # Index of the iterations
    sse_old <- sum((Y-ax-Ya_hat-Yc_hat)^2)
    
    # Iterations
    while((rel_err > tol)& (iter <= maxits)) {
      ### Step 1 (update ax)
      ax<- rowMeans(Y-Ya_hat-Yc_hat)
      
      ### Step 2 (update bx and kt)
      svd_Y <- svd(Y-ax-Yc_hat,m,m)
      bx <-sweep(svd_Y$u, 2, colSums(svd_Y$u), FUN = "/")
      if (m == 1){
        kt <- sweep(svd_Y$d[1] * t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")
      } else if (m > 1){
        kt <- sweep(diag(svd_Y$d[1:m]) %*% t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")
      }
      
      Ya_hat<- bx%*%kt
      Z<- Convert_AT_AC(Y-ax-Ya_hat)
      
      # Step 3 (update b0x and gc)
      # For general Renshaw-Haberman model (non-parametric b0x)
      if (const.b0x == F){
        start.Z_rec <- b0x%*%t(gc)
        Second <- iterative_SVD(Z,start.Z_rec = start.Z_rec)
        b0x <- Second$b0x
        gc <- Second$gc
      } else {
        # For H1 model (b0x=1) without additional constraint
        if (approxConst == F){
          for (i in 1:nc) {
            ismiss <- is.na(Z[,i])
            z_obs<- Z[!ismiss,i]
            gc[i]<- mean(z_obs)
          }
        }
        else {
          # For H1 model (b0x=1) with additional approximate constraint \sum c*(gc-mean(gc))=1
          coeff_mu<- rep(0,nc)
          rhs<- rep(0,nc)
          for (i in 1:nc) {
            s<- i-(n+p)/2
            ismiss <- is.na(Z[,i])
            z_obs<- Z[!ismiss,i]
            coeff_mu[i]<- s^2/sum(b0x[!ismiss]^2)
            rhs[i]<- s*sum(b0x[!ismiss]*z_obs)/sum(b0x[!ismiss]^2)
          }
          mu<- sum(rhs)/sum(coeff_mu)
          for (i in 1:nc) {
            s<- i-(n+p)/2
            ismiss <- is.na(Z[,i])
            z_obs<- Z[!ismiss,i]
            gc[i]<- (sum(b0x[!ismiss]*z_obs)-s*mu)/sum(b0x[!ismiss]^2)
          }
        }
      }
      
      ##### Force sum_gc=0 ######
      d<- mean(gc)
      gc<- gc-d # normalize gc
      ax<- ax+as.vector(d*b0x)  # adjust ax
      
      Yc_hat<- Convert_AC_AT(b0x%*%t(gc))
      
      sse <- sum((Y-ax-Ya_hat-Yc_hat)^2)
      rel_err <- (sse_old-sse)/sse_old
      sse_old <- sse
      iter <- iter + 1
    }
  }
  
  # Get the fitted values of m_{x,t}
  mhat <- exp(Y-ax-Ya_hat-Yc_hat)
  
  ### Calculate the AIC and BIC ###
  # Get the log-likelihood (Gaussian)
  
  
  ssr <- sum((Y-ax-Ya_hat-Yc_hat)^2)
  sigma2 <- ssr/(n*p)
  loglik <- -ssr/(2*sigma2) - n*p*log(sqrt(2*pi*sigma2))
  
  # Get the effective number of parameters
  if (lc == TRUE){
    npar <- p+m*(p+n-2)
  } else {
    if (const.b0x == F){
      npar <- n+3*p-3+m*(n+p-2)
    } else {
      if (approxConst == F){
        npar <- n+2*p-2+m*(n+p-2)
      } else {
        npar <- n+2*p-3+m*(n+p-2)
      }
    }
  }
  
  # Get the AIC and BIC
  AIC <- 2*npar-2*loglik 
  BIC <- log(n*p)*npar-2*loglik 

  # Create an fitRHals object
  fit <- create_fitRHals()
  
  # Fill in the necessary components
  if (lc == TRUE){
    fit$model <- StMoMo::StMoMo(link = "log-Gaussian", staticAgeFun = TRUE, periodAgeFun = rep("NP",m))
  } else {
    if (const.b0x == F){
      fit$model <- StMoMo::StMoMo(link = "log-Gaussian", staticAgeFun = TRUE, periodAgeFun = rep("NP",m),
                                  cohortAgeFun = "NP")
    } else {
      fit$model <- StMoMo::StMoMo(link = "log-Gaussian", staticAgeFun = TRUE, periodAgeFun = rep("NP",m),
                                  cohortAgeFun = "1")
    }
  }
  
  fit$ax <- ax
  fit$bx <- bx
  fit$kt <- kt
  fit$b0x <- b0x
  fit$gc <- gc
  fit$ssr <- ssr
  fit$npar <- npar
  fit$loglik <- loglik
  fit$AIC <- AIC
  fit$BIC <- BIC
  fit$ages <- x
  fit$years <- t
  fit$cohorts <- c
  fit$mxt <- data
  fit$m <- m
  fit$lc <- lc
  fit$const.b0x <- const.b0x
  fit$approxConst <- approxConst
  fit$oxt <- matrix(0, nrow = p, ncol = n)
  fit$wxt <- matrix(1, nrow = p, ncol = n)

  return(fit)
}

