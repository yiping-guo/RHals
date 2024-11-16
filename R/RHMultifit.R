#' Fit a Generalized Common Age Effect Model using Alternating Least Squares
#'
#' This function fits a generalized common age effect (GCAE) model using the least squares method. It accommodates multiple settings for sharing parameters across populations and allows for different model specifications.
#'
#' The general model structure is defined as:
#' \deqn{
#' \log m_{x,t} = a_x + \sum_{i=1}^{m}b^{(i)}_x k^{(i)}_t + b_x^{(0)}\gamma_{t-x}
#' }
#'
#' @param data A list of matrices, where each matrix represents the central mortality rates \eqn{m_{x,t}} for a population, with rows corresponding to ages and columns to years.
#' @param ages A vector of ages corresponding to the rows of each mortality matrix.
#' @param years A vector of years corresponding to the columns of each mortality matrix.
#' @param m A positive integer indicating the number of age-period terms.
#' @param lc A logical value controlling whether the cohort effects are excluded. If \code{TRUE}, a Lee-Carter model with \eqn{m} age-period terms will be fitted using the order-\eqn{m} singular value decomposition (SVD).
#' @param const.b0x A logical value specifying whether \eqn{b_x^{(0)}=1}. If \code{TRUE}, an H1 model will be fitted. Only effective when \code{lc = FALSE}.
#' @param common.bx A logical value indicating whether the age-period parameters \eqn{b_x^{(i)}} should be shared across populations.
#' @param common.b0x A logical value indicating whether the age-cohort parameter \eqn{b_x^{(0)}} should be shared across populations. Only effective when \code{lc = FALSE} and \code{const.b0x = FALSE}.
#' @param maxits A positive integer representing the maximum number of main iterations. Only effective when \code{lc = FALSE}.
#' @param tol A small positive number representing the tolerance level for convergence in the main iterations. Only effective when \code{lc = FALSE}.
#' @return A list containing the following components:
#' \item{\code{fits}}{A list of \code{fitRHals} objects, one for each population. Each \code{fitRHals} object contains:}
#'   \itemize{
#'     \item{\code{model}}{An \code{StMoMo} object containing the fitted model.}
#'     \item{\code{ax}}{A vector with the fitted values of the static age function \eqn{a_x}.}
#'     \item{\code{bx}}{A matrix or array with the fitted values of the age-period function \eqn{b_x}.}
#'     \item{\code{kt}}{A matrix with the fitted values of the period indexes \eqn{k_t}.}
#'     \item{\code{b0x}}{A vector or array with the fitted values of the age-cohort function \eqn{b_x^{(0)}}.}
#'     \item{\code{gc}}{A vector with the fitted values of the cohort indexes \eqn{\gamma_{t-x}}.}
#'     \item{\code{ages}}{A vector of ages of the data set.}
#'     \item{\code{years}}{A vector of years of the data set.}
#'     \item{\code{cohorts}}{A vector of cohorts of the data set.}
#'     \item{\code{mxt}}{The input mortality matrix for the population.}
#'     \item{\code{m}}{The input parameter \eqn{m}.}
#'     \item{\code{lc}}{The input logical variable \code{lc}.}
#'     \item{\code{const.b0x}}{The input logical variable \code{const.b0x}.}
#'   }
#' \item{\code{ssr}}{The total sum of squares error across all populations.}
#' \item{\code{loglik}}{The cumulative Gaussian log-likelihood across all populations.}
#' \item{\code{npar}}{The effective number of parameters in the multi-population model, adjusted for shared parameters.}
#' \item{\code{AIC}}{The Akaike's Information Criterion (AIC) for the overall multi-population model.}
#' \item{\code{BIC}}{The Bayesian Information Criterion (BIC) for the overall multi-population model.}
#' 
#' @examples
#' # Example two-population data (EW male and female)
#' Data1 <- load_EWData(ages = 60:89, years = 1960:2009, series = "Male")
#' Data2 <- load_EWData(ages = 60:89, years = 1960:2009, series = "Female")
#' MData <- list(Data1, Data2)
#'
#' # Fit a RH model for each of the two sub-populations (not common age effects)
#' multifit_RH <- RHMultifit(MData, ages = 60:89, years = 1960:2009) 
#' # Fit a RH model with common age effects (both bx and b0x)
#' multifit_GCAE <- RHMultifit(MData, ages = 60:89, years = 1960:2009,common.bx = TRUE, common.b0x = TRUE)
#' 
#' # Extract the estimated gc from the GCAE model for the female data (indexed the second in Mdata)
#' multifit_GCAE$fits[[2]]$gc
#' @export

RHMultifit <- function(data, ages, years, m = 1, lc = FALSE, const.b0x = FALSE, 
                               common.bx = FALSE, common.b0x = FALSE, 
                               maxits = 10000, tol = 1e-6) {
  
  # Check if data is specified and is a list
  if (missing(data) || is.null(data)) {
    stop("Error: `data` is not specified.")
  } else if (!is.list(data) || !all(sapply(data, is.matrix))) {
    stop("Error: `data` must be a list of numeric matrices.")
  }
  
  # Check if all matrices in data list have the same dimensions
  data_dims <- sapply(data, dim)
  if (any(apply(data_dims, 1, function(x) length(unique(x)) > 1))) {
    stop("Error: All matrices in `data` must have the same dimensions.")
  }
  
  # Check if all matrices in data list are numeric
  if (!all(sapply(data, is.numeric))) {
    stop("Error: All matrices in `data` must be numeric.")
  }
  
  # Check if ages are specified and are a numeric vector
  if (missing(ages) || is.null(ages)) {
    stop("Error: `ages` is not specified.")
  } else if (!is.numeric(ages)) {
    stop("Error: `ages` must be a numeric vector.")
  }
  
  # Check if years are specified and are a numeric vector
  if (missing(years) || is.null(years)) {
    stop("Error: `years` is not specified.")
  } else if (!is.numeric(years)) {
    stop("Error: `years` must be a numeric vector.")
  }
  
  # Check compatibility of ages and years with the dimensions of the mortality matrices in the list
  if (length(ages) != nrow(data[[1]])) {
    stop("Error: Length of `ages` vector does not match the number of rows in each matrix in `data`.")
  }
  
  if (length(years) != ncol(data[[1]])) {
    stop("Error: Length of `years` vector does not match the number of columns in each matrix in `data`.")
  }
  
  # Check that common.b0x can only be set TRUE when lc = FALSE and const.b0x = FALSE
  if (common.b0x && (lc || const.b0x)) {
    stop("Error: `common.b0x` can only be set to TRUE when `lc` is FALSE and `const.b0x` is FALSE.")
  }
  
  
  # Define variables based on dimensions
  I <- length(data)  # Number of populations
  p <- length(ages)  # Number of ages
  n <- length(years) # Number of years
  x <- ages
  t <- years
  c <- (t[1] - x[p]):(t[n] - x[1])
  nc <- length(c)
  
  # Log-transform the mortality data
  Y <- lapply(data, log)
  
  # Initialize lists for each population
  ax <- lapply(1:I, function(i) array(0, dim = p, dimnames = list(x)))
  kt <- lapply(1:I, function(i) array(0, dim = c(m, n), dimnames = list(1:m, t)))
  gc <- lapply(1:I, function(i) array(0, dim = nc, dimnames = list(c)))
  Ya_hat <- lapply(1:I, function(i) array(0, dim = c(p, n))) # Initialize Ya_hat for each population
  Yc_hat <- lapply(1:I, function(i) array(0, dim = c(p, n))) # Initialize Yc_hat for each population
  
  # Initialize bx
  bx <- if (common.bx) {
    array(0, dim = c(p, m), dimnames = list(x, 1:m))  # shared array
  } else {
    lapply(1:I, function(i) array(0, dim = c(p, m), dimnames = list(x, 1:m)))  # list of arrays
  }
  
  # Initialize b0x
  b0x <- if (lc == TRUE){
    array(0, dim = p, dimnames = list(x))
  } else {
    if (const.b0x == TRUE){
      b0x <- array(1, dim = p, dimnames = list(x))
    } else {
      if (common.b0x == TRUE){
        array(1, dim = p, dimnames = list(x))
      } else {
        lapply(1:I, function(i) array(1, dim = p, dimnames = list(x)))
      }
    }
  }
  
  # Initial values for ax, Yc_hat, and Ya_hat
  for (i in 1:I) {
    ax[[i]] <- rowMeans(Y[[i]])
    
    if (common.bx) {
      Ya_hat[[i]] <- bx %*% kt[[i]]
    } else {
      Ya_hat[[i]] <- bx[[i]] %*% kt[[i]]
    }
    
    if (lc == FALSE && const.b0x == FALSE && common.b0x == FALSE) {
      Yc_hat[[i]] <- Convert_AC_AT(b0x[[i]] %*% t(gc[[i]]))
    } else {
      Yc_hat[[i]] <- Convert_AC_AT(b0x %*% t(gc[[i]]))
    }
  }
  
  
  if (lc) {  
    # Lee-Carter Model with closed-form solution via SVD
    if (common.bx) {
      combined_Y <- do.call(cbind, lapply(1:I, function(i) Y[[i]] - ax[[i]] - Yc_hat[[i]]))
      svd_Y <- svd(combined_Y, nu = m, nv = m)
      bx <- sweep(svd_Y$u, 2, colSums(svd_Y$u), FUN = "/")
      
      for (i in 1:I) {
        if (m == 1) {
          kt[[i]] <- sweep(svd_Y$d[1] * t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")[, (1 + (i-1) * n):(i * n)]
        } else if (m > 1) {
          kt[[i]] <- sweep(diag(svd_Y$d[1:m]) %*% t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")[, (1 + (i-1) * n):(i * n)]
        }
        Ya_hat[[i]] <- bx %*% kt[[i]]
      }
    } else {
      for (i in 1:I) {
        svd_Yi <- svd(Y[[i]] - ax[[i]] - Yc_hat[[i]], nu = m, nv = m)
        bx[[i]] <- sweep(svd_Yi$u, 2, colSums(svd_Yi$u), FUN = "/")
        
        if (m == 1) {
          kt[[i]] <- sweep(svd_Yi$d[1] * t(svd_Yi$v), 1, colSums(svd_Yi$u), FUN = "*")
        } else if (m > 1) {
          kt[[i]] <- sweep(diag(svd_Yi$d[1:m]) %*% t(svd_Yi$v), 1, colSums(svd_Yi$u), FUN = "*")
        }
        Ya_hat[[i]] <- bx[[i]] %*% kt[[i]]
      }
    }
  } else {  
    
    # Iterative approach for the RH model
    # Begin iterations for model fitting
    rel_err <- 1
    iter <- 0
    sse_old <- sum(sapply(1:I, function(i) sum((Y[[i]] - ax[[i]] - Ya_hat[[i]] - Yc_hat[[i]])^2)))
    
    while ((rel_err > tol) & (iter <= maxits)) {
      # Step 1: Update ax for each population
      for (i in 1:I) {
        ax[[i]] <- rowMeans(Y[[i]] - Ya_hat[[i]] - Yc_hat[[i]])
      }
      
      # Step 2: Update bx and kt using SVD for each model type
      if (common.bx) {
        combined_Y <- do.call(cbind, lapply(1:I, function(i) Y[[i]] - ax[[i]] - Yc_hat[[i]]))
        svd_Y <- svd(combined_Y, nu = m, nv = m)
        bx <- sweep(svd_Y$u, 2, colSums(svd_Y$u), FUN = "/")
        
        for (i in 1:I) {
          if (m == 1) {
            kt[[i]] <- sweep(svd_Y$d[1] * t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")[, (1 + (i-1) * n):(i * n)]
          } else if (m > 1) {
            kt[[i]] <- sweep(diag(svd_Y$d[1:m]) %*% t(svd_Y$v), 1, colSums(svd_Y$u), FUN = "*")[, (1 + (i-1) * n):(i * n)]
          }
          # Update Ya_hat after updating bx and kt
          Ya_hat[[i]] <- bx %*% kt[[i]]
        }
      } else {
        for (i in 1:I) {
          svd_Yi <- svd(Y[[i]] - ax[[i]] - Yc_hat[[i]], nu = m, nv = m)
          bx[[i]] <- sweep(svd_Yi$u, 2, colSums(svd_Yi$u), FUN = "/")
          
          if (m == 1) {
            kt[[i]] <- sweep(svd_Yi$d[1] * t(svd_Yi$v), 1, colSums(svd_Yi$u), FUN = "*")
          } else if (m > 1) {
            kt[[i]] <- sweep(diag(svd_Yi$d[1:m]) %*% t(svd_Yi$v), 1, colSums(svd_Yi$u), FUN = "*")
          }
          # Update Ya_hat after updating bx and kt
          Ya_hat[[i]] <- bx[[i]] %*% kt[[i]]
        }
      }
      
      
      # Step 3: Update b0x and gc using iterative SVD for each model type
      if (const.b0x == FALSE) {
        # General RH model (non-parametric b0x)
        if (common.b0x){
          # Shared b0x across populations
          combined_Z <- do.call(cbind, lapply(1:I, function(i) Convert_AT_AC(Y[[i]] - ax[[i]] - Ya_hat[[i]])))
          start.Z_rec <- do.call(cbind, lapply(1:I, function(i) b0x %*% t(gc[[i]])))
          Second <- iterative_SVD(combined_Z, start.Z_rec = start.Z_rec)
          b0x <- Second$b0x
          for (i in 1:I) {
            gc[[i]] <- Second$gc[(1 + (i-1) * nc):(i * nc)]
          }
        } else {
          # General RH model (non-parametric b0x) with different b0x for each population
          for (i in 1:I) {
            Z <- Convert_AT_AC(Y[[i]] - ax[[i]] - Ya_hat[[i]])
            start.Z_rec <- b0x[[i]] %*% t(gc[[i]])
            Second <- iterative_SVD(Z, start.Z_rec = start.Z_rec)
            b0x[[i]] <- Second$b0x
            gc[[i]] <- Second$gc
          }
        }
      } else {
        # General H1 model (const b0x)
        for (i in 1:I) {
          Z <- Convert_AT_AC(Y[[i]] - ax[[i]] - Ya_hat[[i]])
          for (j in 1:nc) {
            ismiss <- is.na(Z[,j])
            z_obs <- Z[!ismiss, j]
            gc[[i]][j] <- mean(z_obs)
          }
        }
      }
      
      # Normalize gc to enforce sum(gc) = 0 for each population
      # And update Yc_hat
      for (i in 1:I) {
        d <- mean(gc[[i]])
        gc[[i]] <- gc[[i]] - d
        if (lc == FALSE && common.b0x == FALSE && const.b0x == FALSE){
          ax[[i]] <- ax[[i]] + as.vector(d * b0x[[i]])
          Yc_hat[[i]] <- Convert_AC_AT(b0x[[i]] %*% t(gc[[i]]))
        } else {
          ax[[i]] <- ax[[i]] + as.vector(d * b0x)
          Yc_hat[[i]] <- Convert_AC_AT(b0x %*% t(gc[[i]]))
        }
      }
      
      # Calculate the residuals and SSE for convergence check
      sse <- sum(sapply(1:I, function(i) sum((Y[[i]] - ax[[i]] - Ya_hat[[i]] - Yc_hat[[i]])^2)))
      rel_err <- (sse_old - sse) / sse_old
      sse_old <- sse
      iter <- iter + 1
    }
  }
  
  # Multi-population Output Processing
  
  # Initialize empty list to store fitRHals objects for each population
  fit_list <- vector("list", I)
  total_ssr <- 0  # Total SSR
  total_loglik <- 0  # Total log-likelihood
  
  for (i in 1:I) {
    # Get the fitted values for each population
    mhat <- exp(Y[[i]] - ax[[i]] - Ya_hat[[i]] - Yc_hat[[i]])
    
    # Calculate SSR for the population
    ssr <- sum((Y[[i]] - ax[[i]] - Ya_hat[[i]] - Yc_hat[[i]])^2)
    total_ssr <- total_ssr + ssr  # Accumulate total SSR
    
    # Calculate log-likelihood for the population
    sigma2 <- ssr / (n * p)
    loglik <- -ssr / (2 * sigma2) - n * p * log(sqrt(2 * pi * sigma2))
    total_loglik <- total_loglik + loglik  # Accumulate total log-likelihood
    
    # Create an individual fitRHals object for the population
    fit <- create_fitRHals()
    
    # Fill in components for the individual fitRHals object
    if (lc == TRUE) {
      fit$model <- StMoMo::StMoMo(link = "log-Gaussian", staticAgeFun = TRUE, periodAgeFun = rep("NP", m))
    } else {
      if (const.b0x == FALSE) {
        fit$model <- StMoMo::StMoMo(link = "log-Gaussian", staticAgeFun = TRUE, periodAgeFun = rep("NP", m),
                                    cohortAgeFun = "NP")
      } else {
        fit$model <- StMoMo::StMoMo(link = "log-Gaussian", staticAgeFun = TRUE, periodAgeFun = rep("NP", m),
                                    cohortAgeFun = "1")
      }
    }
    
    # Fill in remaining components
    fit$ax <- ax[[i]]
    
    if (common.bx) {
      fit$bx <- bx  # Shared bx across populations
    } else {
      fit$bx <- bx[[i]]  # Population-specific bx
    }
    
    fit$kt <- kt[[i]]
    
    if (lc == FALSE && const.b0x == FALSE && common.b0x == FALSE) {
      fit$b0x <- b0x[[i]]  # Population-specific b0x (only suitted for RH model with shared b0x)
    } else {
      fit$b0x <- b0x  # Shared b0x across populations
    }
    
    fit$gc <- gc[[i]]
    fit$ages <- x
    fit$years <- t
    fit$cohorts <- c
    fit$mxt <- data[[i]]
    fit$m <- m
    fit$lc <- lc
    fit$const.b0x <- const.b0x
    fit$oxt <- matrix(0, nrow = p, ncol = n)
    fit$wxt <- matrix(1, nrow = p, ncol = n)
    
    # Add fitRHals object for this population to the list
    fit_list[[i]] <- fit
  }
  
  
  # Calculating overall metrics for the multi-population model
  # Total SSR and log-likelihood have been accumulated
  
  # Calculate the effective number of parameters (npar) based on the model specifications
  npar <- NA  
  
  if (lc == TRUE) { 
    # LC Model
    if (common.bx) {
      npar <- I * p + m * p + I * m * n - m * (I + 1)
    } else {
      npar <- I * (p + m * (p + n - 2))
    }
  } else {
    # H1 Model
    if (const.b0x == TRUE){
      if (common.bx) {
        npar <- I * p + m * p + I * m * n + I * (n + p - 1) - (m + m * I + I)
      } else {
        npar <- I * (n + 2 * p - 2 + m * (n + p - 2))
      }
    } else {
      # RH Model
      if (common.bx && common.b0x) {
        npar <- I * p + m * p + I * m * n + p + I * (n + p - 1) - (m + m * I + 1 + I)
      } else if (common.bx && !common.b0x) {
        npar <- I * p + m * p + I * m * n + I * p + I * (n + p - 1) - (m + m * I + 2 * I)
      } else if (!common.bx && common.b0x) {
        npar <- I * p + I * m * p + I * m * n + p + I * (n + p - 1) - (2 * m * I + 1 + I)
      } else {
        npar <- I * (p + m * p + m * n + p + (n + p - 1) - 2 * (m + 1))
      }
    }
  }
  
  
  # Calculate AIC and BIC using the accumulated log-likelihood and total npar
  AIC <- 2 * npar - 2 * total_loglik
  BIC <- log(n * p * I) * npar - 2 * total_loglik
  
  # Combine all information into a summary list
  multi_fit <- list(
    fits = fit_list,
    ssr = total_ssr,
    loglik = total_loglik,
    npar = npar,
    AIC = AIC,
    BIC = BIC
  )
  
  # Return the final list containing all individual fits and overall metrics
  return(multi_fit)
}
