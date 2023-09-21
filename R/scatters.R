# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


# TCOV scatter matrix -----


#' Pairwise one-step M-estimate of scatter for ICS
#' 
#' Wrapper function for the pairwise one-step M-estimator of scatter with 
#' weights based on pairwise Mahalanobis distances, as computed by 
#' \code{\link{tcov}()}. Note that this estimator is based on pairwise 
#' differences and therefore no location estimate is returned.
#' 
#' @param x  a numeric matrix or data frame.
#' @param beta  a positive numeric value specifying the tuning parameter of the 
#' pairwise one-step M-estimator (default to 2), see \code{\link{tcov}()}.
#'
#' @return An object of class \code{"ICS_scatter"} with the following
#' components:
#' \item{location}{this is \code{NULL} as the estimator is based on pairwise 
#' differences and does not use a location estimate.}
#' \item{scatter}{a numeric matrix giving the estimate of the scatter matrix.}
#' \item{label}{a character string providing a label for the scatter matrix.}
#' 
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link[ICS]{ICS}()}
#'
#' \code{\link{tcov}()}, \code{\link{ucov}()}, \code{\link{ICS_ucov}()}
#' 
#' 
#' @export

ICS_tcov <- function(x, beta = 2) {
  # compute scatter estimate
  out <- list(location = NULL, scatter = tcov(x, beta = beta), label = "TCOV")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


#' Pairwise one-step M-estimate of scatter
#' 
#' Computes a pairwise one-step M-estimate of scatter with weights based on 
#' pairwise Mahalanobis distances. Note that it is based on pairwise 
#' differences and therefore does not require a location estimate.
#' 
#' For a sample \eqn{\boldsymbol{X}_{n} = (\mathbf{x}_{1}, \dots, 
#' \mathbf{x}_n)^{\top}}, a positive and decreasing weight function \eqn{w}, 
#' and a tuning parameter \eqn{\beta > 0}, the pairwise one-step M-estimator 
#' of scatter is defined as
#' \deqn{\mathrm{TCOV}_{\beta}(\boldsymbol{X}_{n}) =
#' \frac{\sum_{i=1}^{n-1} \sum_{j=i+1}^{n} 
#' w(\beta \, r^{2}(\mathbf{x}_{i}, \mathbf{x}_{j})) 
#' (\mathbf{x}_{i} - \mathbf{x}_{j}) 
#' (\mathbf{x}_{i} - \mathbf{x}_{j})^{\top}}{\sum_{i=1}^{n-1} \sum_{j=i+1}^{n} 
#' w(\beta \, r^{2}(\mathbf{x}_{i}, \mathbf{x}_{j}))},}
#' where 
#' \deqn{r^{2}(\mathbf{x}_{i}, \mathbf{x}_{j}) = 
#' (\mathbf{x}_{i} - \mathbf{x}_{j})^{\top}
#' \mathrm{COV}(\boldsymbol{X}_n)^{-1} 
#' (\mathbf{x}_{i} - \mathbf{x}_{j})}
#' denotes the squared pairwise Mahalanobis distance between observations 
#' \eqn{\mathbf{x}_{i}} and \eqn{\mathbf{x}_{j}} based on the sample 
#' covariance matrix \eqn{\mathrm{COV}(\boldsymbol{X}_n)}. Here, the weight 
#' function \eqn{w(x) = \exp(-x/2)} is used.
#' 
#' @param x  a numeric matrix or data frame.
#' @param beta  a positive numeric value specifying the tuning parameter of the 
#' pairwise one-step M-estimator (defaults to 2), see \sQuote{Details}.
#' 
#' @return A numeric matrix giving the pairwise one-step M-estimate of scatter.
#' 
#' @author Andreas Alfons and Aurore Archimbaud
#' 
#' @references 
#' Caussinus, H. and Ruiz-Gazen, A. (1993) Projection Pursuit and Generalized 
#' Principal Component Analysis. In Morgenthaler, S., Ronchetti, E., Stahel, 
#' W.A. (eds.) *New Directions in Statistical Data Analysis and Robustness*, 
#' 35-46. Monte Verita, Proceedings of the Centro Stefano Franciscini Ascona 
#' Series. Springer-Verlag.
#' 
#' Caussinus, H. and Ruiz-Gazen, A. (1995) Metrics for Finding Typical 
#' Structures by Means of Principal Component Analysis. In *Data Science and 
#' its Applications*, 177-192. Academic Press.
#'
#' @seealso
#' \code{\link{ICS_tcov}()}, \code{\link{ucov}()}, \code{\link{ICS_ucov}()}
#' 
#' 
#' @useDynLib ICSClust, .registration = TRUE
#' @export

tcov <- function(x, beta = 2) {
  # initializations
  x <- as.matrix(x)
  cn <- colnames(x)
  # call internal function (which in turn calls C++ function)
  V <- tcov_cpp(x, beta = beta)
  # set row and column names and return scatter matrix
  dimnames(V) <- list(cn, cn)
  V
}


# # reference implementation using package 'amap'
# #' @importFrom amap W
# #' @importFrom stats var
# tcov_amap <- function(x, beta = 2) {
#   # initializations
#   x <- as.matrix(x)
#   cn <- colnames(x)
#   # compute inverse of maximum likelihood estimate of covariance matrix
#   # n <- nrow(x)
#   # S_inv <- solve(var(x) * ((n-1)/n))
#   S_inv <- solve(var(x))
#   V <- amap::W(x, h = 1/sqrt(beta), D = S_inv, kernel = "gaussien")
#   # set row and column names and return scatter matrix
#   dimnames(V) <- list(cn, cn)
#   V
# }


# SCOV scatter matrix -----
#' @name ICS_ucov
#' @export
ICS_scov <- function(x, location = TRUE, beta = 0.2) {
  # initializations
  x <- as.matrix(x)
  n <- nrow(x)
  location <- isTRUE(location)
  # compute location and scatter estimates
  m <- colMeans(x)
  S_inv <- solve(var(x))
  scatter <- .scov(x, m = m, S_inv = S_inv, beta = beta)
  out <- list(location = if (location) m, scatter = scatter, label = "SCOV")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}

#' @name ucov
#' @export
scov <- function(x, beta = 0.2) {
  # initializations
  x <- as.matrix(x)
  n <- nrow(x)
  # compute sample means and inverse of covariance matrix
  m <- colMeans(x)
  S_inv <- solve(var(x))
  # call internal function
  .scov(x, m = m, S_inv = S_inv, beta = beta)
}

## internal function to avoid recomputation of sample mean
#' @useDynLib ICSClust, .registration = TRUE
.scov <- function(x, m, S_inv, beta = 0.2) {
  # initializations
  cn <- colnames(x)
  # call internal function (which in turn calls C++ function)
  V <- scov_cpp(x, m = m, S_inv = S_inv, beta = beta)
  # set row and column names and return scatter matrix
  dimnames(V) <- list(cn, cn)
  V
}


# UCOV scatter matrix -----


#' Simple robust estimates of scatter for ICS
#' 
#' Wrapper functions for the one-step M-estimator of scatter with weights based 
#' on Mahalanobis distances as computed by \code{\link{scov}()}, or the simple 
#' related estimator that is based on a transformation as computed by 
#' \code{\link{ucov}()}.
#' 
#' @param x  a numeric matrix or data frame.
#' @param location  a logical indicating whether to include the sample
#' mean as location estimate (defaults to \code{TRUE}).
#' @param beta  a positive numeric value specifying the tuning parameter of the 
#' estimator (default to 0.2), see \code{\link{ucov}()}.
#'
#' @return An object of class \code{"ICS_scatter"} with the following
#' components:
#' \item{location}{if requested, a numeric vector giving the location estimate.}
#' \item{scatter}{a numeric matrix giving the estimate of the scatter matrix.}
#' \item{label}{a character string providing a label for the scatter matrix.}
#' 
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link[ICS]{ICS}()}
#'
#' \code{\link{tcov}()}, \code{\link{ICS_tcov}()}, \code{\link{ucov}()}
#' 
#' 
#' @export

ICS_ucov <- function(x, location = TRUE, beta = 0.2) {
  # initializations
  x <- as.matrix(x)
  location <- isTRUE(location)
  # compute location and scatter estimates
  m <- colMeans(x)
  scatter <- .ucov(x, m = m, beta = beta)
  out <- list(location = if (location) m, scatter = scatter, label = "UCOV")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


#' Simple robust estimates of scatter
#' 
#' Compute a one-step M-estimator of scatter with weights based on Mahalanobis 
#' distances, or a simple related estimator that is based on a transformation.
#' 
#' For a sample \eqn{\boldsymbol{X}_{n} = (\mathbf{x}_{1}, \dots, 
#' \mathbf{x}_n)^{\top}}, a positive and decreasing weight function \eqn{w}, 
#' and a tuning parameter \eqn{\beta > 0}, the one-step M-estimator 
#' of scatter is defined as
#' \deqn{\mathrm{SCOV}_{\beta}(\boldsymbol{X}_{n}) =
#' \frac{\sum_{i=1}^{n} 
#' w(\beta \, r^{2}(\mathbf{x}_{i})) 
#' (\mathbf{x}_{i} - \mathbf{\bar{x}}_{n}) 
#' (\mathbf{x}_{i} -  \mathbf{\bar{x}}_{n})^{\top}}{\sum_{i=1}^{n} 
#' w(\beta \, r^{2}(\mathbf{x}_{i}))},}
#' where 
#' \deqn{r^{2}(\mathbf{x}_{i}) = 
#' (\mathbf{x}_{i} -  \mathbf{\bar{x}}_{n})^{\top}
#' \mathrm{COV}(\boldsymbol{X}_n)^{-1} 
#' (\mathbf{x}_{i} -  \mathbf{\bar{x}}_{n})}
#' denotes the squared Mahalanobis distance of observation \eqn{\mathbf{x}_{i}} 
#' from the sample mean \eqn{\mathbf{\bar{x}}_{n}} based on the sample 
#' covariance matrix \eqn{\mathrm{COV}(\boldsymbol{X}_n)}. Here, the weight 
#' function \eqn{w(x) = \exp(-x/2)} is used.
#' 
#' A simple robust estimator that is consistent under normality is obtained via 
#' the transformation
#' \deqn{\mathrm{UCOV}_{\beta}(\boldsymbol{X}_{n}) = 
#' (\mathrm{SCOV}_{\beta}(\boldsymbol{X}_{n})^{-1} - 
#' \beta \, \mathrm{COV}(\boldsymbol{X}_{n})^{-1})^{-1}.}
#' 
#' @param x  a numeric matrix or data frame.
#' @param beta  a positive numeric value specifying the tuning parameter of the 
#' estimator (defaults to 0.2), see \sQuote{Details}.
#' 
#' @return A numeric matrix giving the estimate of the scatter matrix.
#' 
#' @author Andreas Alfons and Aurore Archimbaud
#' 
#' @references 
#' Caussinus, H. and Ruiz-Gazen, A. (1993) Projection Pursuit and Generalized 
#' Principal Component Analysis. In Morgenthaler, S., Ronchetti, E., Stahel, 
#' W.A. (eds.) *New Directions in Statistical Data Analysis and Robustness*, 
#' 35-46. Monte Verita, Proceedings of the Centro Stefano Franciscini Ascona 
#' Series. Springer-Verlag.
#' 
#' Caussinus, H. and Ruiz-Gazen, A. (1995) Metrics for Finding Typical 
#' Structures by Means of Principal Component Analysis. In *Data Science and 
#' its Applications*, 177-192. Academic Press.
#' 
#' Ruiz-Gazen, A. (1996) A Very Simple Robust Estimator of a Dispersion Matrix. 
#' *Computational Statistics & Data Analysis*, **21**(2), 149-162. 
#' \doi{10.1016/0167-9473(95)00009-7}.
#'
#' @seealso
#' \code{\link{ICS_ucov}()}, \code{\link{tcov}()}, \code{\link{ICS_tcov}()}
#' 
#' 
#' @export

ucov <- function(x, beta = 0.2) {
  # initializations
  x <- as.matrix(x)
  # compute sample means
  m <- colMeans(x)
  # call internal function
  .ucov(x, m = m, beta = beta)
}


# internal function to avoid recomputation of sample mean
.ucov <- function(x, m, beta = 0.2) {
  # initializations
  n <- nrow(x)
  # compute inverse of sample covariance matrix and inverse of SCOV
  S_inv <- solve(var(x))
  scov_inv <- solve(.scov(x, m = m, S_inv = S_inv, beta = beta))
  # compute UCOV scatter matrix
  solve(scov_inv - beta * S_inv)
}

# ## reference implementation using package 'amap'
# #' @importFrom amap W
# #' @importFrom stats var
# ucov_amap <- function(x, beta = 0.2) {
#   # initializations
#   x <- as.matrix(x)
#   cn <- colnames(x)
#   # compute inverse of maximum likelihood estimate of covariance matrix
#   # n <- nrow(x)
#   # S_inv <- solve(var(x) * ((n-1) / n))
#   S_inv <- solve(var(x))
#   # compute inverse of SCOV (error in 'amap' package: should be 1/h^2)
#   h <- 1/sqrt(beta)
#   scov_inv <- solve(amap::varrob(x, h = h, D = S_inv, kernel = "gaussien"))
#   V <- solve(scov_inv + (h-1)/(h^2) * S_inv)
#   # set row and column names and return scatter matrix
#   dimnames(V) <- list(cn, cn)
#   V
# }


# MLC scatter matrix -----

#' Cauchy location and Scatter Estimates for ICS
#' 
#' It is a wrapper for the Cauchy estimator of location and scatter
#' for a multivariate t-distribution, as computed by [ICS::tM()].
#' 
#' @param x a numeric matrix or data frame.
#' @param location  a logical indicating whether to include the M-estimate of 
#' location (defaults to `FALSE`).
#' @param ... potential further arguments passed to [ICS::ICS_tM()].
#'
#' @return An object of class \code{"ICS_scatter"} with the following
#' components:
#' \item{location}{if requested, a numeric vector giving the location
#' estimate.}
#' \item{scatter}{a numeric matrix giving the estimate of the scatter matrix.}
#' \item{label}{a character string providing a label for the scatter matrix.}
#' @author Andreas Alfons and Aurore Archimbaud
#' @seealso [ICS::tM()], [ICS::ICS_tM()]
#' @export
#' @importFrom ICS tM
ICS_mlc <- function(x, location = FALSE, ...) {
  # compute scatter estimates
  mlc_out <- ICS::tM(x, df = 1, ...) # we fix the df to have only cauchy estimate
  location <- isTRUE(location)
  location <- if (location) mlc_out$mu
  # compute scatter estimate
  out <- list(location = location, scatter = mlc_out$V, label = "MLC")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}

# LCOV scatter matrix -----

#' Local Shape Scatter Estimates for ICS
#' 
#' It is a wrapper for the local shape estimator of scatter
#' as computed by [fpc::localshape()].
#' 
#' @param x a numeric matrix or data frame.
#' @param mscatter `"mcd"` or `"cov"` (default); specified minimum covariance 
#' determinant or classical covariance matrix to be used for Mahalanobis
#'  distance computation.
#' @param proportion proportion of points to be considered as neighbourhood.
#' @param ... potential further arguments passed to [fpc::localshape()].
#'
#' @return An object of class \code{"ICS_scatter"} with the following
#' components:
#' \item{location}{this is NULL as the estimator does not use a
#' location estimate.}
#' \item{scatter}{a numeric matrix giving the estimate of the scatter matrix.}
#' \item{label}{a character string providing a label for the scatter matrix.}
#' @author Andreas Alfons and Aurore Archimbaud
#' @seealso [fpc::localshape()]
#' 
#' @export
#' @importFrom fpc localshape
ICS_lcov <- function(x, mscatter = "cov", proportion = 0.1, ...) {
  # compute scatter estimate
  out <- list(location = NULL, 
              scatter = fpc::localshape(xdata = x, mscatter = mscatter,
                                        proportion = proportion, ...), 
              label = "LCOV")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}



# MCD scatter matrix -----

#' MCD location and Scatter Estimates for ICS
#' 
#' It is a wrapper for the (reweighted) MCD estimators of location and scatter
#' as computed by [rrcov::CovMcd()].
#' 
#' @param x a numeric matrix or data frame.
#' @param location  a logical indicating whether to include the MCD-estimate of
#' location (defaults to `FALSE`). 
#' @param nsamp number of subsets used for initial estimates or `"best"`, 
#' `"exact"` or `"deterministic"` (default).
#' @param alpha numeric parameter controlling the size of the subsets over 
#' which the determinant is minimized as in [rrcov::CovMcd()].
#' @param ... potential further arguments passed to [rrcov::CovMcd()].
#' 
#' @details
#' - [ICS_mcd_raw()]: computes the raw MCD estimates.
#' - [ICS_mcd_rwt()]: computes the reweighted MCD estimates.
#' 
#' 
#' 
#' @return An object of class \code{"ICS_scatter"} with the following
#' components:
#' \item{location}{if requested, a numeric vector giving the location
#' estimate.}
#' \item{scatter}{a numeric matrix giving the estimate of the scatter matrix.}
#' \item{label}{a character string providing a label for the scatter matrix.}
#' 
#' @seealso [rrcov::CovMcd()]
#' @author Andreas Alfons and Aurore Archimbaud
#' @name ICS_mcd
#' @export
#' @importFrom rrcov CovMcd
ICS_mcd_raw <- function(x, location = FALSE,
                    nsamp = "deterministic",  alpha = 0.5, ...) {
  # compute scatter estimates
  mcd_out <- rrcov::CovMcd(x, raw.only = TRUE, alpha = alpha, nsamp = nsamp, 
                           ...)
  location <- isTRUE(location)
  location <- if (location) mcd_out@center
  # compute scatter estimate
  out <- list(location = location, scatter = mcd_out@cov, label = "MCD")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}


#' @name ICS_mcd
#' @export
#' @importFrom rrcov CovMcd
ICS_mcd_rwt <- function(x, location = FALSE, nsamp = "deterministic",
                     alpha = 0.5, ...) {
  # compute scatter estimates
  rmcd_out <- rrcov::CovMcd(x, raw.only = FALSE, alpha = alpha, nsamp = nsamp,  
                            ...)
  location <- isTRUE(location)
  location <- if (location) rmcd_out@center
  # compute scatter estimate
  out <- list(location = location, scatter = rmcd_out@cov, 
              label = "RMCD")
  # add class and return object
  class(out) <- "ICS_scatter"
  out
}