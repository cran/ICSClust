#' *k*-means clustering
#' 
#' Wrapper for performing k-means clustering from [stats::kmeans()].
#'
#' @param X a numeric matrix or data frame of the data. It corresponds to the 
#' argument \code{x}.
#' @param k the number of clusters searched for. It corresponds to the argument
#'  \code{centers}.
  #' @param clusters_only boolean. If \code{TRUE} only the partition of the data
  #'  is returned as a vector. If \code{FALSE} the usual output of 
  #'  the \link[stats]{kmeans} function is returned.
#' @param iter.max the maximum number of iterations allowed.
#' @param nstart if `centers` is a number, how many random sets should be 
#' chosen.
#' @param ... other arguments to pass to the [stats::kmeans()] function.
#'
#' @return If \code{clusters_only} is \code{TRUE} a vector of the new partition
#'  of the data is returned, i.e a vector of integers (from \code{1:k}) 
#'  indicating the cluster to which each observation is allocated. 
#' 
#' 
#' Otherwise a list is returned with the following components:
#' \item{clust_method}{the name of the clustering method, i.e. "kmeans".}
#' \item{clusters}{the vector of the new partition of the data, i.e. a vector of 
#' integers (from \code{1:k}) indicating the cluster to which each observation 
#' is allocated.}
#' \item{...}{an object of class \code{"kmeans"}}.
#' 
#' @seealso [stats::kmeans()] 
#' @export
#' 
#' @importFrom stats kmeans
#' 
#' @author Aurore Archimbaud
#'
#' @examples
#' kmeans_clust(iris[,1:4], k = 3, clusters_only = TRUE)
#' 
kmeans_clust <- function(X, k, clusters_only = FALSE, iter.max = 100,
                         nstart = 20, ...){
  # Perform kmeans
  clust <- kmeans(x = X, centers = k, iter.max = iter.max, nstart = nstart, 
                  ...)
  # Get the new partition
  out <-  clust$cluster
  
  # Return the results
  if (!clusters_only) out <- append(list(clust_method = "kmeans", 
                                         clusters = out),
                                    clust)
  out
}


#' Trimmed k-means clustering
#' 
#' Wrapper for performing trimmed k-means clustering from 
#' [tclust::tkmeans()].
#'
#'
#'
#' @param X a numeric matrix or data frame of the data. It corresponds to the 
#' argument \code{x}.
#' @param k the number of clusters searched for. It corresponds to the argument
#'  \code{k}.
#' @param clusters_only boolean. If \code{TRUE} only the partition of the data 
#' is returned as a vector. If \code{FALSE} the usual output of the 
#' \link[tclust]{tkmeans} function is returned.
#' @param alpha the proportion of observations to be trimmed.
#' @param ... other arguments to pass to the [tclust::tkmeans()]
#'
#' @return If \code{clusters_only} is \code{TRUE} a vector of the new partition 
#' of the data is returned, i.e a vector of integers (from \code{1:k}) 
#' indicating the cluster to which each observation is allocated. 
#' 0 indicates trimmed observations.
#' 
#' 
#' Otherwise a list is returned with the following components:
#' \item{clust_method}{the name of the clustering method, i.e. "tkmeans".}
#' \item{clusters}{the vector of the new partition of the data, i.e. a vector of 
#' integers (from \code{1:k}) indicating the cluster to which each observation
#'  is allocated. 0 indicates trimmed observations.}
#' \item{...}{an object of class \code{"tkmeans"}}.
#' 
#' @seealso [tclust::tkmeans()]
#' @export
#' 
#' @author Aurore Archimbaud
#' @importFrom tclust tkmeans
#'
#' @examples
#' tkmeans_clust(iris[,1:4], k = 3, alpha = 0.1, clusters_only = TRUE)
tkmeans_clust <- function(X, k, clusters_only = FALSE, alpha = 0.05, ... ){
  # Perform tkmeans
  clust <- tclust::tkmeans(x = X, k = k, alpha = alpha, ...)
  
  # Get the new partition: the additional cluster of outliers are coded by 0
  out <-  clust$cluster
  
  # Return the results
  if (!clusters_only) out <- append(list(clust_method = "tkmeans",
                                         clusters = out), 
                                    clust)
  out
}

#' Partitioning Around Medoids clustering
#' 
#' Wrapper for performing Partitioning Around Medoids clustering from 
#' [cluster::pam()].
#'
#'
#'
#' @param X a numeric matrix or data frame of the data. It corresponds to the 
#' argument \code{x}.
#' @param k the number of clusters searched for. It corresponds to the argument
#'  \code{k}.
#' @param clusters_only boolean. If \code{TRUE} only the partition of the data 
#' is returned as a vector. If \code{FALSE} the usual output of the 
#' [cluster::pam()] function is returned.
#' @param ... other arguments to pass to the [cluster::pam()].
#'
#' @return If \code{clusters_only} is \code{TRUE} a vector of the new partition
#'  of the data is returned, i.e a vector of integers (from \code{1:k}) 
#'  indicating the cluster to which each observation is allocated.
#'   0 indicates trimmed observations.
#' 
#' 
#' Otherwise a list is returned with the following components:
#' \item{clust_method}{the name of the clustering method, i.e "clara_pam".}
#' \item{clusters}{the vector of the new partition of the data, i.e a vector of 
#' integers (from \code{1:k}) indicating the cluster to which each observation 
#' is allocated.
#'  0 indicates outlying observations.}
#' \item{...}{an object of class \code{"pam"}}.
#' 
#' @seealso [cluster::pam()]
#'
#' @export
#' @author Aurore Archimbaud
#' @importFrom cluster pam
#'
#' @examples
#' pam_clust(iris[,1:4], k = 3, clusters_only = TRUE)
#'  
pam_clust <- function(X, k, clusters_only = FALSE, ...){
  # Perform PAM
  clust <- cluster::pam(x = X, k = k, diss = FALSE, ...)
  
  # Get the new partition: the additional cluster of outliers are coded by 0
  out <- clust$clustering
  
  # Return the results
  if (!clusters_only) out <- append(list(clust_method = "clara_pam", 
                                         clusters = out), 
                                    clust)
  out
}




#' Robust Improper Maximum Likelihood Clustering
#' 
#' Wrapper for performing Robust Improper Maximum Likelihood Clustering 
#' clustering from [otrimle::rimle()].
#'
#' @param X a numeric matrix or data frame of the data. It corresponds to the 
#' argument \code{data}.
#' @param k the number of clusters searched for. It corresponds to the argument
#'  \code{G}.
#' @param clusters_only boolean. If \code{TRUE} only the partition of the data 
#' is returned as a vector. If \code{FALSE} the usual output of the 
#' [otrimle::rimle()] function is returned.
#' @param ... other arguments to pass to [otrimle::rimle()].
#'
#' @return If \code{clusters_only} is \code{TRUE} a vector of the new partition
#'  of the data is returned, i.e a vector of integers (from \code{1:k}) 
#'  indicating the cluster to which each observation is allocated.
#'   0 indicates trimmed observations.
#' 
#' 
#' Otherwise a list is returned with the following components:
#' \item{clust_method}{the name of the clustering method, i.e, "rimle".}
#' \item{clusters}{the vector of the new partition of the data, i.e. a vector of 
#' integers (from \code{1:k}) indicating the cluster to which each observation 
#' is allocated.
#'  0 indicates outlying observations.}
#' \item{...}{an object of class \code{"rimle"}}
#' 
#' @seealso [otrimle::rimle()]
#'
#' @export
#' @author Aurore Archimbaud
#' @importFrom otrimle rimle
#'
#' @examples
#' rimle_clust(iris[,1:4], k = 3, clusters_only = TRUE)
rimle_clust <- function(X, k, clusters_only = FALSE, ...){
  # Perform rimle
  clust <- otrimle::rimle(data = X, G = k,  ...)
  
  # Get the new partition: the additional cluster of outliers are coded by 0
  out <- clust$cluster
  
  # Return the results
  if (!clusters_only) out <- append(list(clust_method = "rimle", 
                                         clusters = out), 
                                    clust)
  out
}

#' Model-Based Clustering
#' 
#' Wrapper for performing Model-Based Clustering from [mclust::Mclust()]
#' allowing noise or not.
#'
#' @param X a numeric matrix or data frame of the data. It corresponds to the 
#' argument \code{data}.
#' @param k the number of clusters searched for. It corresponds to the argument
#'  \code{G} of function [mclust::Mclust()].
#' @param clusters_only boolean. If \code{TRUE} only the partition of the data 
#' is returned as a vector. If \code{FALSE} the usual output of the 
#' [mclust::Mclust()] function is returned.
#' @param ... other arguments to pass to [mclust::Mclust()].
#' 
#' @details
#' - [mclust_clust()]: does not allow noise
#' - [rmclust_clust()]: allows noise
#' 
#'
#' @return If \code{clusters_only} is \code{TRUE} a vector of the new partition
#'  of the data is returned, i.e a vector of integers (from \code{1:k}) 
#'  indicating the cluster to which each observation is allocated.
#'   0 indicates trimmed observations.
#' 
#' 
#' Otherwise a list is returned with the following components:
#' \item{clust_method}{the name of the clustering method, i.e "rimle".}
#' \item{clusters}{the vector of the new partition of the data, i.e a vector of 
#' integers (from \code{1:k}) indicating the cluster to which each observation 
#' is allocated.
#'  0 indicates outlying observations for [rmclust_clust()] only.}
#' \item{...}{an object of class "\code{mclust}"}
#' 
#' @seealso [mclust::Mclust()]
#'
#' @export
#' @author Aurore Archimbaud
#' @importFrom mclust Mclust mclustBIC
#' @rdname mclust_clust
#' @examples
#' mclust_clust(iris[,1:4], k = 3, clusters_only = TRUE)
mclust_clust <- function(X, k, clusters_only = FALSE, ...){
  # Perform rimle
  clust <- mclust::Mclust(X, G = k, ...)
  # Get the new partition: the additional cluster of outliers are coded by 0
  out <- clust$classification
  
  # Return the results
  if (!clusters_only) out <- append(list(clust_method = "mclust", 
                                         clusters = out), 
                                    clust)
  out
}
#' @rdname mclust_clust
#' @export
rmclust_clust <- function(X, k, clusters_only = FALSE, ...){
  # Perform rimle
  clust <- mclust::Mclust(X, G = k, initialization = list(noise = TRUE), ...)
  # Get the new partition: the additional cluster of outliers are coded by 0
  out <- clust$classification
  
  # Return the results
  if (!clusters_only) out <- append(list(clust_method = "rmclust", 
                                         clusters = out), 
                                    clust)
  out
}

