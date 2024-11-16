#' Convert the Age-Time matrix to the Age-Cohort matrix
#' @keywords internal
Convert_AT_AC <- function(X) {
  p <- nrow(X)
  n <- ncol(X)
  Y <- matrix(rep(NA, p * (n + p - 1)), nrow = p)
  for (i in 1:p) {
    Index_fill <- (p - i + 1):(n + p - i)
    Y[i, Index_fill] <- X[i, ]
  }
  return(Y)
}

#' Convert the Age-Cohort matrix to the Age-Time matrix
#' @keywords internal
Convert_AC_AT <- function(Z) {
  p <- nrow(Z)
  n <- ncol(Z) - p + 1
  Y <- matrix(rep(NA, p * n), nrow = p)
  for (i in 1:p) {
    Index_fill <- (p - i + 1):(n + p - i)
    Y[i, ] <- Z[i, Index_fill]
  }
  return(Y)
}

#' Create als_fitStMoMo Object
#' @keywords internal
create_fitRHals <- function() {
  fit <- list(
    model = NULL,
    ax = NULL,
    bx = NULL,
    kt = NULL,
    b0x = NULL,
    gc = NULL,
    ssr = NULL,
    npar = NULL,
    loglik = NULL,
    AIC = NULL,
    BIC = NULL,
    ages = NULL,
    years = NULL,
    cohorts = NULL,
    mxt = NULL,
    m = NULL,
    lc = NULL,
    const.b0x = NULL,
    approxConst = NULL,
    oxt = NULL,
    wxt = NULL
  )
  class(fit) <- c("fitRHals", "fitStMoMo")
  return(fit)
}

#' Print Method for fitRHals Objects
#' @export
print.fitRHals <- function(x, ...) {
  cat("Generalized Renshaw-Haberman Model fit")
  cat("\n\n")
  print(x$model)
  cat(paste("\nYears in fit:", min(x$years), "-", max(x$years)))
  cat(paste("\nAges in fit:", min(x$ages), "-", max(x$ages), "\n"))
}

#' Create bootRHals Object
#' @keywords internal
create_bootRHals <- function() {
  boot <- list(
    bootParameters = NULL,
    model = NULL
  )
  class(boot) <- c("bootRHals", "bootStMoMo")
  return(boot)
}


#' Print Method for bootRHals Objects
#' @export
print.bootRHals <- function(x, ...) {
  cat("Bootstrapped Generalized Renshaw-Haberman Model")
  cat("\n\n")
  cat("Residual bootstrap based on\n")
  print(x$model$model)
  cat(paste("\n\nNumber of bootstrap samples:", length(x$bootParameters)))
  cat(paste("\nYears in fit:", min(x$model$years), "-", max(x$model$years)))
  cat(paste("\nAges in fit:", min(x$model$ages), "-", max(x$model$ages), "\n"))
}




#' Function for iterative SVD (used for estimating the Renshaw-Haberman model)
#' @keywords internal
iterative_SVD<- function(Z,start.Z_rec=NULL,maxits=1000,tol=1e-6){
  p<- nrow(Z) # The number of age bands
  n<- ncol(Z)-p+1 # The number of time bands
  
  ### Step 1: Initialization
  rel_err <- 1 # Relative error
  iter <- 0 # Index of the iterations
  # Set initial values for the imputation matrix (Z_rec means Z recov)
  if (is.null(start.Z_rec)){
    Z_rec<- matrix(rep(rowMeans(Z,na.rm=T),n+p-1),nrow = p)
  } else {
    Z_rec<- start.Z_rec
  }
  
  ismiss <- is.na(Z)
  sse_old <- sum(((Z-Z_rec)[!ismiss])^2)
  
  ### Step 2: Iterations
  while((rel_err > tol)& (iter <= maxits)) {
    Z_imp<- Z_rec
    Z_imp[!ismiss]<- Z[!ismiss]
    SVD_Z_imp<- svd(Z_imp,1,1)
    Z_rec<- SVD_Z_imp$u%*%t(SVD_Z_imp$v)*SVD_Z_imp$d[1]
    
    sse <- sum(((Z_rec-Z)[!ismiss])^2)
    rel_err <- (sse_old-sse)/sse_old
    sse_old <- sse
    iter <- iter + 1
  }
  b0x=SVD_Z_imp$u/sum(SVD_Z_imp$u)
  gc<- SVD_Z_imp$v*SVD_Z_imp$d[1]*sum(SVD_Z_imp$u)
  return(list("b0x"=b0x,"gc"=gc,"Z_rec"=Z_rec))
}
