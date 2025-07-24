#' Tandem clustering with ICS
#' 
#' Sequential clustering approach: (i) dimension reduction through the Invariant 
#' Coordinate Selection method using the \code{\link[ICS]{ICS-S3}} function and (ii)
#' clustering of the transformed data. 
#'
#' @param X a numeric matrix or data frame containing the data.
#' @param nb_select the number of components to select. 
#' It is used only in case \code{criterion} is either `"med_crit"`, `"var_crit"` 
#' or `"discriminatory_crit"`.  By default it is set to \code{NULL}, i.e the number 
#' of components to select is the number of clusters minus one.
#' @param nb_clusters the number of clusters searched for. 
#' @param ICS_args list of \code{\link[ICS]{ICS-S3}} arguments. Otherwise, default 
#' values of \code{\link[ICS]{ICS-S3}} are used.
#' @param criterion criterion to automatically decide which invariant components
#'  to keep. Possible values are `"med_crit"`, `"normal_crit"`, `"var_crit"` and 
#'  `"discriminatory_crit"`. The default value is `"med_crit"`. 
#'  See [med_crit()], [normal_crit()], [var_crit()] or 
#' [discriminatory_crit()] for more details.
#' @param ICS_crit_args list of arguments passed to [med_crit()], [normal_crit()],
#' [var_crit()] or \cr [discriminatory_crit()] for choosing the components to keep.
#' @param method clustering method to perform. Currently implemented wrapper
#'  functions are `"kmeans_clust"`, `"tkmeans_clust"`, `"pam_clust"`,
#'  `"mclust_clust"`, `"rmclust_clust"` or `"rimle_clust"`.
#'  The default value is `"kmeans_clust"`. 
#' @param clustering_args list of [kmeans_clust()], 
#' [tkmeans_clust()],  [pam_clust()], [rimle_clust()], [mclust_clust()] or
#' [rmclust_clust()]
#' arguments for performing cluster analysis.
#' @param clusters a vector indicating the true clusters of the data. By default,
#' it is \code{NULL} but it is required to choose the components based on the 
#' discriminatory criterion \code{\link{discriminatory_crit}}.
#'
#' @details
#' Tandem clustering with ICS is a sequential method:
#' 
#' - \code{\link[ICS]{ICS-S3}} is performed.
#' 
#' - only a subset of the first and/or the last few components are
#'  selected based on a criterion.
#'  
#' - the clustering method is performed only on the subspace
#'  of the selected components.
#'
#' - wrapper for several different clustering methods are provided. Users can however
#'   also write wrappers for other clustering methods.
#'
#' @return  
#' An object of class `"ICSClust"` with the following components:
#' - `ICS_out`: An object of class \code{"ICS"}. 
#' See \code{\link[ICS]{ICS-S3}}.
#' - `select`: a vector of the names of the selected invariant
#'  coordinates.
#' - `clusters`: a vector of the new partition of the data, i.e a vector
#'  of integers (from \code{1:k}) indicating the cluster to which each
#'   observation is allocated. 0 indicates outlying observations.
#'   
#' [summary()][summary.ICSClust()] and [plot()][plot.ICSClust()] methods are available.
#'  
#' @references
#' Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2024). 
#' Tandem clustering with invariant coordinate selection. 
#' Econometrics and Statistics.
#'  \doi{10.1016/j.ecosta.2024.03.002}.
#' 
#' @export
#' 
#' @author Aurore Archimbaud
#' 
#' @seealso [med_crit()], [normal_crit()], 
#' [var_crit()], \code{\link[ICS]{ICS-S3}}, 
#' [discriminatory_crit()], [kmeans_clust()],
#' [tkmeans_clust()], [pam_clust()],
#' [rimle_clust()], [mclust_clust()]
#' [summary()] and [plot()] methods
#' 
#'
#' @examples
#' X <- iris[,1:4]
#' 
#' # indicating the number of components to retain for the dimension reduction
#' # step as well as the number of clusters searched for.
#' out <- ICSClust(X, nb_select = 2, nb_clusters = 3)
#' summary(out)
#' plot(out)
#' 
#' # changing the scatter pair to consider in ICS
#' out <- ICSClust(X, nb_select = 1, nb_clusters = 3,
#' ICS_args = list(S1 = ICS_mcd_raw, S2 = ICS_cov,S1_args = list(alpha = 0.5)))
#' summary(out)
#' plot(out)
#'  
#' # changing the criterion for choosing the invariant coordinates
#' out <- ICSClust(X, nb_clusters = 3, criterion = "normal_crit",
#' ICS_crit_args = list(level = 0.1, test = "anscombe.test", max_select = NULL))
#' summary(out)
#' plot(out)
#' 
#' # changing the clustering method
#' out <- ICSClust(X, nb_clusters = 3, method  = "tkmeans_clust", 
#' clustering_args = list(alpha = 0.1))
#' summary(out)
#' plot(out)
ICSClust <- function(X, nb_select = NULL, nb_clusters = NULL, ICS_args = list(),
                     criterion = c("med_crit", "normal_crit", "var_crit",
                                   "discriminatory_crit"), 
                     ICS_crit_args = list(),
                     method = c("kmeans_clust", "tkmeans_clust", "pam_clust",
                                "mclust_clust",  "rmclust_clust", "rimle_clust"),
                     clustering_args = list(),
                     clusters = NULL
){
  
  # Initialization
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  
  if(is.null(nb_clusters)){
    stop("You should specify the `nb_clusters` argument.")
  }
  
  # ICS ----
  ICS_out <-  do.call(ICS::ICS, append(list(X = X), ICS_args))
  
  # Choice of components ----
  if(criterion %in% c("med_crit", "var_crit", "discriminatory_crit")){
    # if nb_select is NULL we put the nb_clusters-1
    nb_select <- ifelse(is.null(nb_select), nb_clusters-1, nb_select)
    ICS_crit_args <- append(ICS_crit_args, c(nb_select = nb_select))
  }
  if (criterion %in% c("discriminatory_crit")){
    ICS_crit_args <- append(ICS_crit_args, list(clusters = clusters))
  }
 
  select <- do.call(criterion, append(list(object = ICS_out,
                                           select_only = TRUE),
                                      ICS_crit_args))
  
  # Clustering ----
  reduced_df <- ICS::components(ICS_out, select = select)
  if(ncol(reduced_df) == 0){
    clusters <- NULL
    warning("No component has been selected.")
  }else{
    selected <- paste(colnames(reduced_df), collapse = ",")
    clusters <- do.call(method,
                        append(list(X = reduced_df, k = nb_clusters,
                                    clusters_only = TRUE),
                               clustering_args))
    
  }
  
  
  out <- list(ICS_out = ICS_out, select = select, clusters = clusters)
  class(out) <- "ICSClust"
  out
}


#' Summary of an `ICSClust` object
#' 
#' Summarizes an `ICSClust` object in an informative way.
#'
#' @param object object of class `"ICSClust"`.
#' @param ... additional arguments passed to [summary()]
#'
#' @return An object of class `"ICSClust_summary"` with the following components:
#' - `ICS_out`: `ICS_out` object
#' - `nb_comp`: number of selected components
#' - `select`: vector of names of selected components
#' - `nb_clusters`: number of clusters
#' - `table_clusters`: frequency table of clusters
#' 
#' 
#' @export
#' 
#' @name summary.ICSClust
#' @method summary ICSClust
#' @author Aurore Archimbaud

summary.ICSClust <- function(object, ...) {
  out <- list(ICS_out = object$ICS_out,
              nb_comp =  length(object$select),
              select = object$select,
              nb_clusters = length(unique(object$clusters)),
              table_clusters = table(object$clusters))
  
  class(out) <- "ICSClust_summary"
  out
}

#' Print of an `ICSClust_summary` object
#' 
#' Prints an `ICSClust_summary` object in an informative way.
#'
#' @param x object of class `"ICSClust_summary"`.
#' @param info logical, either TRUE or FALSE. If TRUE, prints additional
#' information on arguments used for computing scatter matrices
#' (only named arguments that contain numeric, character, or logical scalars)
#' and information on the parameters of the algorithm.
#' Default is FALSE.
#' @param digits  number of digits for the numeric output.
#' @param ... additional arguments are ignored.
#' @return The supplied object of class `"ICSClust_summary"` is returned invisibly.
#' @export
#' 
#' @name print.ICSClust_summary
#' @method print ICSClust_summary
#' @author Aurore Archimbaud
print.ICSClust_summary <-  function(x, info = FALSE, digits = 4L, ...) {
  # print information on ICS
  print(x$ICS_out, info = info, digits = digits)
  
  # print information on selected components
  cat("\n", x$nb_comp, "components are selected:", x$select)
  
  # print information on clusters
  cat("\n\n", x$nb_clusters, "clusters are identified:\n")
  print(x$table_clusters)
  
  # return object invisibly
  invisible(x)
}



#' Scatterplot Matrix with densities on the diagonal
#' 
#' Wrapper for [component_plot()].
#' 
#' @param x an object of class `"ICSClust"`. 
#' @param \dots additional arguments to be passed down to [component_plot()]
#' 
#' @return An object of class [`"ggmatrix"`][GGally::ggmatrix()] (see 
#' [GGally::ggpairs()]).
#' 
#' @name plot.ICSClust
#' @method plot ICSClust
#' @author Aurore Archimbaud
#' @export
plot.ICSClust <- function(x, ...) {
  if(length(x$select) == 0){
    warning("No component has been selected.")
  }else{
    component_plot(x$ICS_out, select = x$select, 
                   clusters = factor(x$clusters), ...)
  }
}
