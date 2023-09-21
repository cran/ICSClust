#' Simulation of a mixture of Gaussian distributions
#' 
#' Simulation of a \eqn{n \times p} data frame according to a mixture of \eqn{q}
#' Gaussian distributions with \eqn{q < p}, different location parameters
#' \eqn{\mu_1, \dots, \mu_q}, and the identity matrix as the covariance matrix.
#'
#' @param pct_clusters a vector of marginal probabilities for each group, i.e
#' mixture weights. 
#' Default is two 
#' balanced clusters.
#' @param n integer. The number of observations.
#' @param p integer. The number of variables.
#' @param delta integer. The location shift.
#' 
#' @details
#' Let \eqn{X} be a \eqn{p}-variate real random vector distributed according to
#' a mixture of \eqn{q} Gaussian distributions with \eqn{q < p}, 
#' different location parameters \eqn{\mu_1, \dots, \mu_q}, and the same positive 
#' definite covariance matrix \eqn{I_p}:
#'  \deqn{X \sim \sum_{h=1}^{q} \epsilon_h \, {\cal N}(\mu_h,I_p),}
#'  where \eqn{\epsilon_{1}, \dots, \epsilon_{q}} are mixture weights with 
#'  \eqn{\epsilon_1 + \cdots + \epsilon_q = 1},  \eqn{\mu_1 = 0_p},
#'   and  \eqn{\mu_{h+1} = \delta e_h} with \eqn{h = 1, \dots, q-1}.
#'  
#' 
#'
#' @return A dataframe of *n* observations and *p+1* variables with the first
#' variable
#' indicating the cluster assignment using a character string.
#' @export
#' 
#' @importFrom mvtnorm rmvnorm
#' @author Aurore Archimbaud
#' 
#' @references
#' Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2022). 
#' Tandem clustering with invariant coordinate selection. 
#'  \emph{arXiv preprint arXiv:2212.06108}..
#' 
#' @examples
#' X <- mixture_sim()
#' summary(X)
mixture_sim = function(pct_clusters = c(0.5,0.5) , n = 500, p = 10, delta = 10){
  
  # Checks and initialization of inputs
  if(sum(pct_clusters)!=1){
    stop("the sum of the groups is not equal to 1!") 
  }
  n_groups = floor(n*pct_clusters)
  if(sum(n_groups)!= n){
    n_groups[length(pct_clusters)] <- 
      rev(n-cumsum(n_groups)[1:(length(pct_clusters)-1)])[1]
    
  }
  if(sum(n_groups)!= n){
    warning(paste("the total number of observation is not equal to n",
                  paste(round(pct_clusters,2), collapse = " - ")))
  }
 
  # We simulate normal gaussian for each cluster
  clusters_means <- rep(0,p)
  X_list <- lapply(1:length(pct_clusters), function(i){
    n <-  n_groups[i]
    # we introduce the shift location 
    if(i>1){
      clusters_means[i-1] = delta
    }
    if(n>0){
      data.frame(mvtnorm::rmvnorm(n = n, mean = clusters_means, 
                                  sigma = diag(1,p)))
    }
  })
  
  data.frame(cluster = rep(paste0("Group", 1:length(pct_clusters)), n_groups),
                    do.call(rbind, (X_list)))
}

#' Uniform distribution outside a given range
#' 
#' Draw from a multivariate uniform distribution outside a given range.
#' Intuitively speaking, the observations are drawn from a multivariate 
#' uniform distribution on a hyperrectangle with a hole in the middle (in the 
#' shape of a smaller hyperrectangle). This is useful, e.g., for adding random 
#' noise to a data set such that the noise consists of large values that do not 
#' overlap the initial data.
#'
#' @param n  an integer giving the number of observations to generate.
#' @param min  a numeric vector giving the minimum of each variable of the 
#' initial data set (outside of which to generate random noise).
#' @param max  a numeric vector giving the maximum of each variable of the 
#' initial data set (outside of which to generate random noise).
#' @param mult  multiplication factor (larger than 1) to expand the 
#' hyperrectangle around the initial data (which is given by `min` and `max`). 
#' For instance, the default value 2 gives a hyperrectangle for which each 
#' side is twice as long as the range of the initial data.  The data are then 
#' drawn from a uniform distribution on the expanded hyperrectangle from which 
#' the smaller hyperrectangle around the data is cut out. See the examples for 
#' an illustration.
#' 
#' @return A matrix of generated points.
#' 
#' @author Andreas Alfons
#'
#' @examples
#' ## illustrations for argument 'mult'
#' 
#' # draw observations with argument 'mult = 2'
#' xy2 <- runif_outside_range(1000, min = rep(-1, 2), max = rep(1, 2), 
#'                            mult = 2)
#' # each side of the larger hyperrectangle is twice as long as 
#' # the corresponding side of the smaller rectanglar cut-out
#' df2 <- data.frame(x = xy2[, 1], y = xy2[, 2])
#' ggplot(data = df2, mapping = aes(x = x, y = y)) + 
#'   geom_point()
#' 
#' # draw observations with argument 'mult = 4'
#' xy4 <- runif_outside_range(1000, min = rep(-1, 2), max = rep(1, 2), 
#'                            mult = 4)
#' # each side of the larger hyperrectangle is four times as long 
#' # as the corresponding side of the smaller rectanglar cut-out
#' df4 <- data.frame(x = xy4[, 1], y = xy4[, 2])
#' ggplot(data = df4, mapping = aes(x = x, y = y)) + 
#'   geom_point()
#' 
#' @importFrom stats runif
#' @export
#' 
#' @references
#' #' Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2022). 
#' Tandem clustering with invariant coordinate selection. 
#' \emph{arXiv preprint arXiv:2212.06108}.

runif_outside_range <- function(n, min = 0, max = 1, mult = 2) {
  # check number of observations
  n <- rep(n, length.out = 1L)
  if(!is.numeric(n) || !is.finite(n) || n < 0L) {
    stop("number of observations 'n' must be a positive integer")
  } else n <- as.integer(n)
  # check supplied range
  min <- as.numeric(min)
  p <- length(min)
  max <- as.numeric(max)
  if (length(max) != p) stop("'min' and 'max' must have the same length")
  if (!all(is.finite(min)) || !all(is.finite(max))) {
    stop("'min' and 'max' may not contain missing or infinite values")
  }
  if (any(min >= max)) {
    stop("values of 'min' must be smaller than the corresponding values of ", 
         "'max'")
  }
  # check multiplication factor
  mult <- rep(mult, length.out = 1L)
  if(!is.numeric(mult) || !is.finite(mult) || mult <= 1L) {
    stop("multiplication factor 'mult' must be larger than 1")
  }
  # center and dimension of hyperrectangles
  m <- (min + max) / 2
  # length of sides of hyperrectangle around data
  a <- max - min
  # halved length of sides of expanded hyperrectangle
  h <- mult * a / 2
  # lower and upper bounds of expanded hyperrectangle
  lower <- m - h
  upper <- m + h
  # draw from uniform distribution on expanded hyperrectangle and reject the 
  # observation if it lies in the smaller hyperrectangle around the data
  n_ok <- 0L
  X <- matrix(NA_real_, nrow = n, ncol = p)
  while (n_ok < n) {
    # draw new point from uniform distribution on expanded hyperrectangle
    x <- runif(p, min = lower, max = upper)
    # reject if it lies inside of the smaller hyperrectangle around the data
    reject <- all(x >= min & x <= max)
    # add to matrix
    if (!reject) {
      n_ok <- n_ok + 1L
      X[n_ok, ] <- x
    }
  }
  # return matrix of generated points
  X
}
