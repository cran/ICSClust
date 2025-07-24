#' Selection of Non-normal Invariant Components Using Marginal 
#' Normality Tests
#' 
#' Identifies invariant coordinates that are non normal using univariate 
#' normality tests as in the `comp.norm.test` function from the 
#' `ICSOutlier` package, with the difference that both the 
#' first and last few components are investigated.
#'
#' @param object object of class `"ICS"` or a data frame or matrix.
#' @param level the initial level used to make a decision based on the test
#' p-values. See details. Default is 0.05.
#' @param test name of the normality test to be used. Possibilities are 
#' `"jarque.test"`, `"anscombe.test"`, `"bonett.test"`, `"agostino.test"`,
#'  `"shapiro.test"`.
#' Default is `"agostino.test"`.
#' @param max_select the maximal number of components to select.
#' @param select_only boolean. If `TRUE` only the vector names of the selected 
#' invariant components is returned. If `FALSE` additional details are returned.
#' @param gen_kurtosis vector of generalized kurtosis values.
#' @param \dots  additional arguments are currently ignored.
#'  
#' @details
#' The procedure sequentially tests the first and the last components until 
#' finding no additional components as non-normal. The quantile levels are
#' adjusted for multiple testing by taking the level as `level`/*j* for the 
#' *j*th component.
#' 
#' @return If `select_only` is  `TRUE` a vector of the names of the invariant
#'  components or variables to select. If `FALSE` an object of class `"ICS_crit"`
#'  is returned with the following objects: 
#'  - `crit`: the name of the criterion "normal".
#'  - `level`: the level of the test.
#'  - `max_select`: the maximal number of components to select.
#'  - `test`: name of the normality test to be used.
#'  - `pvalues`: the p-values of the tests.
#'  - `adjusted_levels`: the adjusted levels.
#'  - `select`: the names of the invariant components or variables to select.
#'  - `gen_kurtosis`: the vector of generalized kurtosis values in case of 
#'  `ICS` object.
#'  
#' @export
#' 
#' 
#' @references
#' Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2024). 
#' Tandem clustering with invariant coordinate selection. 
#' Econometrics and Statistics.
#'  \doi{10.1016/j.ecosta.2024.03.002}.
#' 
#' Archimbaud, A., Nordhausen, K., and Ruiz-Gazen, A. (2018). 
#' ICSOutlier: Unsupervised Outlier Detection for Low-Dimensional Contamination Structure, 
#' The RJournal, Vol. 10(1):234–250. \doi{10.32614/RJ-2018-034}
#' 
#' Archimbaud, A., Nordhausen, K., and Ruiz-Gazen, A. (2016). 
#' ICSOutlier: Outlier Detection Using Invariant Coordinate Selection. 
#' R package version 0.3-0
#' 
#' @seealso [med_crit()], [var_crit()], [discriminatory_crit()],
#'  [moments::jarque.test()], [moments::anscombe.test()], 
#' [moments::bonett.test()], [moments::agostino.test()], [stats::shapiro.test()].
#'
#' @author Andreas Alfons, Aurore Archimbaud, Klaus Nordhausen and Anne Ruiz-Gazen
#' 
#' @rdname normal_crit
#' @examples
#' X <- iris[,-5]
#' out <- ICS(X)
#' normal_crit(out, level = 0.1, select_only = FALSE)
#' 
normal_crit <- function(object, ...) UseMethod("normal_crit")

#' @rdname normal_crit
#' @import ICS
#' @method normal_crit ICS
#' @export
normal_crit.ICS <- function(object, level = 0.05,  
                            test = c("agostino.test", "jarque.test", 
                                     "anscombe.test", "bonett.test", 
                                     "shapiro.test"), 
                            max_select = NULL, select_only = FALSE, ...){
  # get the generalized values
  gen_kurtosis <- ICS::gen_kurtosis(object, scale = FALSE)
  
  # apply the normal criterion
  normal_crit(ICS::components(object), gen_kurtosis = gen_kurtosis, 
              level = level, test = test, max_select = max_select,
              select_only = select_only)
}


#' @export
#' @rdname normal_crit
#' @method normal_crit default
#' @importFrom moments jarque.test anscombe.test bonett.test agostino.test
#' @importFrom stats shapiro.test
normal_crit.default <- function(object, level = 0.05,  
                                test = c("agostino.test", "jarque.test", 
                                         "anscombe.test", "bonett.test", 
                                         "shapiro.test"), 
                                max_select = NULL, select_only = FALSE,
                                gen_kurtosis = NULL, ...){
  
  # Initialization
  test <- match.arg(test)
  comp_select <- colnames(object)
  max_select <- ifelse(is.null(max_select), ncol(object)-1, max_select)
  
  
  # Apply marginal normality tests to all components and keep only the ones lower
  test_pvals <- apply(object, 2, test)
  test_pvals <- unlist(lapply(test_pvals, function(x) x$p.value))
  out <- list(crit = "normal", level = level,  max_select =  max_select, 
              test = test, pvalues = test_pvals)
  comp_signif <- (test_pvals <= level)
  
  if(sum(comp_signif) == 0){
    select <- vector()
    adjusted_levels <- level
  }else{
    # We select the components on the extreme left and right while 
    # they are not gaussian
    temp <- 1
    select <- vector()
    adjusted_levels <- level
    while (temp <= max_select) {
      
      # for each iteration we compare the first and last p value
      
      left_pval <- test_pvals[1]
      right_pval <- rev(test_pvals)[1]
      
      if (left_pval < level & left_pval < right_pval) {
        
        # we select the first component if its p value is the smallest and 
        # significant
        select <- c(select, comp_select[1])
        comp_select <- comp_select[-1]
        test_pvals <- test_pvals[-1]
        
      } else if (right_pval < level)  {
        
        # we retrieve the selected component to only look at the next one
        select <- c(select, comp_select[length(comp_select)])
        comp_select <- comp_select[-length(comp_select)]
        test_pvals <- test_pvals[-length(test_pvals)]
        
      } else {
        
        break
        
      }
      
      temp <- temp + 1
      # we adjust the alpha level by dividing by the old weight and multiplying 
      # by the new one.
      test_pvals <- test_pvals * temp/(temp-1)
      adjusted_levels <- c(adjusted_levels, level/temp)
      
    }
  }
  if(select_only){
    out <- select
  }else{
    out <- append(out, list(adjusted_levels =  adjusted_levels, select = select,
                            gen_kurtosis = gen_kurtosis))
    class(out) <- "ICS_crit"
  } 
  
  out
  
}



#' Selection of Invariant components using the med criterion
#'
#' Identifies as interesting invariant coordinates whose generalized eigenvalues are the furthermost away 
#' from the median of all generalized eigenvalues.
#'
#' @param object object of class `"ICS"`.
#' @param nb_select the exact number of components to select. By default it is set to
#' `NULL`, i.e the number of components to select is the number of variables minus one.
#' @param select_only boolean. If `TRUE` only the vector names of the selected 
#' invariant components is returned. If `FALSE` additional details are returned. 
#' @param \dots  additional arguments are currently ignored.
#' 
#' 
#' @details
#' If more than half of the components are "uninteresting" and have the same generalized eigenvalue then the median 
#' of all generalized eigenvalues corresponds 
#' to the uninteresting component generalized eigenvalue. The components of interest are the ones whose generalized eigenvalues differ 
#' the most from the median. The motivation of this criterion depends therefore on the assumption that at least half of the 
#' components have equal generalized eigenvalues.
#'                              
#' @return If `select_only` is  `TRUE` a vector of the names of the invariant
#'  components or variables to select. If `FALSE` an object of class `"ICS_crit"`
#'  is returned with the following objects: 
#'  - `crit`: the name of the criterion "med".
#'  - `nb_select`: the number of components to select.
#'  - `gen_kurtosis`: the vector of generalized kurtosis values.
#'  - `med_gen_kurtosis`: the median of the generalized kurtosis values.
#'  - `gen_kurtosis_diff_med`: the absolute differences between the generalized
#'   kurtosis values and the median.
#'  - `select`: the names of the invariant components or variables to select.
#' 
#'  
#' @export
#' 
#' @rdname med_crit
#' 
#' @references
#' Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2024). 
#' Tandem clustering with invariant coordinate selection. 
#' Econometrics and Statistics.
#'  \doi{10.1016/j.ecosta.2024.03.002}.
#' 
#' @seealso [normal_crit()], [var_crit()], [discriminatory_crit()].
#'
#' @author Andreas Alfons, Aurore Archimbaud and Klaus Nordhausen
#' @import ICS
#' @examples
#' X <- iris[,-5]
#' out <- ICS(X)
#' med_crit(out, nb_select = 2, select_only = FALSE)

med_crit <- function(object, ...) UseMethod("med_crit")

#' @import ICS
#' @method med_crit ICS
#' @rdname med_crit
#' @export
med_crit.ICS <- function(object, nb_select = NULL, select_only = FALSE, ...){
  med_crit(ICS::gen_kurtosis(object, scale = FALSE), nb_select = nb_select,
           select_only = select_only)
}

#' @method med_crit default
#' @rdname med_crit
#' @export
med_crit.default <- function(object, nb_select = NULL, select_only = FALSE, ...){
  # Initialization
  nb_select <- ifelse(is.null(nb_select), length(object)-1, nb_select)
  
  # we take the components associated to the furthest eigenvalues from the median
  med_gen_kurtosis <- median(object)
  gen_kurtosis_diff <- abs(object - med_gen_kurtosis)
  out <- names(sort(gen_kurtosis_diff, decreasing = TRUE))[seq(0, nb_select)]
  
  if (!select_only){
    out <- append(list(crit = "med", nb_select =  nb_select, 
                       gen_kurtosis = object, 
                       med_gen_kurtosis = med_gen_kurtosis,
                       gen_kurtosis_diff_med = gen_kurtosis_diff),
                  list(select = out))
    class(out) <- "ICS_crit"
  } 
  out
  
}


#' Selection of Invariant components using the var criterion
#' 
#' Identifies the interesting invariant coordinates based on the rolling
#' variance criterion as used in the `ICSboot` function of the `ICtest`
#' package. It computes rolling variances on the generalized eigenvalues
#' obtained through [ICS::ICS()].
#'   
#' @param object object of class `"ICS"`.
#' @param nb_select the exact number of components to select. By default it is set to
#' `NULL`, i.e the number of components to select is the number of variables minus one.
#' @param select_only boolean. If `TRUE` only the vector names of the selected 
#' invariant components is returned. If `FALSE` additional details are returned.
#' @param \dots  additional arguments are currently ignored. 
#' 
#' @details
#'  Assuming that the generalized eigenvalues of the uninformative components are all the same
#'  means that the variance of these generalized eigenvalues must be minimal. 
#'  Therefore when \code{nb_select} components should be selected, the method identifies 
#'  the \code{p - nb_select} neighboring generalized eigenvalues with minimal variance, 
#'  where \code{p} is the total number of components. The number of interesting components should be at
#'  most \code{p-2} as at least two uninteresting components are needed to compute a variance.
#' 
#'
#' @return If `select_only` is `TRUE` a vector of the names of the invariant
#'  components or variables to select. If `FALSE` an object of class `"ICS_crit"`
#'  is returned with the following objects: 
#'  - `crit`: the name of the criterion "var".
#'  - `nb_select`: the number of components to select.
#'  - `gen_kurtosis`: the vector of generalized kurtosis values.
#'  - `select`: the names of the invariant components or variables to select.
#'  - `RollVarX`: the rolling variances of order d-`nb_select`.
#'  - `Order`: indexes of the ordered invariant components such that the
#'  ones associated to the smallest variances of the eigenvalues are at the
#'  end.
#' 
#' @export
#' @references
#' Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2024). 
#' Tandem clustering with invariant coordinate selection. 
#' Econometrics and Statistics.
#'  \doi{10.1016/j.ecosta.2024.03.002}.
#' 
#' Radojicic, U., & Nordhausen, K. (2019). 
#' Non-gaussian component analysis: Testing the dimension of the signal subspace.
#'  In Workshop on Analytical Methods in Statistics (pp. 101–123). Springer. 
#'  \doi{10.1007/978-3-030-48814-7_6}.
#' 
#' @seealso [normal_crit()], [med_crit()], [discriminatory_crit()].
#' 
#' @rdname var_crit
#'
#' @author Andreas Alfons, Aurore Archimbaud and Klaus Nordhausen
#' @import ICS
#' @examples
#' X <- iris[,-5]
#' out <- ICS(X)
#' var_crit(out, nb_select = 2, select_only = FALSE)
#' 
var_crit <- function(object, ...) UseMethod("var_crit")

#' @import ICS
#' @method var_crit ICS
#' @rdname var_crit
#' @export
var_crit.ICS <- function(object,  nb_select = NULL, select_only = FALSE, ...){
  var_crit(ICS::gen_kurtosis(object, scale = FALSE), 
           nb_select = nb_select,
           select_only = select_only)
}

#' @method var_crit default
#' @rdname var_crit
#' @export
var_crit.default <- function(object, nb_select = NULL, select_only = FALSE, ...){
  # Initialization
  d <- length(object)
  nb_select <- ifelse(is.null(nb_select), d-1, nb_select)
  
  # if the number of non-spherical components is equal or higher to p-1,
  # it makes no sense to compute the rolling variance of one component
  if (nb_select >= (d-1)){
    stop("The nb_select is higher or equal to the number of variables
            minus one, so it makes no sense to select some components based on
            the rolling variance of only one invariant component.")
    out <- vector()
  }else{
    orderD <- fixOrder(object, d-nb_select)
    out <- as.vector(names(object)[orderD$Order[seq(0,nb_select)]])
    if (!select_only){
      out <- append(list(crit = "var", nb_select = nb_select,
                         gen_kurtosis = object, select = out), 
                    orderD)
      class(out) <- "ICS_crit"
    }
  }
  out
  
  
}

#' Fix Order based on rolling variances
#' 
#' Simplified version of the non exported function `fixOrder` from `ICtest`
#' package.
#'
#' @param x a vector of generalized kurtosis values.
#' @param nb_spherical the number of non spherical components to consider. 
#'
#'
#' @return A list with:
#'  - `RollVarX`: the rolling variances.
#'  - `Order`: indexes of the ordered invariant components.
#' 
#' @importFrom RcppRoll roll_var 
#' @noRd
fixOrder <- function (x, nb_spherical) 
{
  P <- length(x)
  Index <- seq(1, P)
  RollVarX <- RcppRoll::roll_var(x, nb_spherical)
  Start <- which.min(RollVarX)
  End <- Start + nb_spherical - 1
  Smallest <- seq(Start, End)
  Order <- c(Index[-Smallest], Index[Smallest])
  Ox <- x[Order]
  out <- list(RollVarX = RollVarX, Order = Order)
  out
}


#' Selection of ICS components based on discriminatory power
#' 
#' Identifies invariant coordinates associated to the highest discriminatory 
#' power. Currently, the implemented measure is "eta2" as quantified by the
#'  Wilks' partial eta-squared, computed using the [heplots::etasq()]
#'   function.
#'
#' @param object dataframe or object of class `"ICS"`.
#' @param clusters a vector of the same length as the number of
#'  observations, indicating the true clusters. It is used to compute
#' the discriminatory power based on it.
#' @param method the name of the discriminatory power. 
#' Only `"eta2"` is implemented.
#' @param nb_select the exact number of components to select. 
#' By default it is set to \code{NULL}, i.e the number 
#' of components to select is the number of clusters minus one.
#' @param select_only boolean. If `TRUE` only the vector names of the selected 
#' invariant components are returned. If `FALSE` additional details are returned. 
#' @param gen_kurtosis vector of generalized kurtosis values.
#' @param \dots  additional arguments are currently ignored.
#'
#' @details
#' The discriminatory power is evaluated for each combination of the
#'  first and/or last combinations of `nb_select` components. The combination
#'  achieving the highest discriminatory power is selected.
#'  
#' More specifically, we compute \eqn{\eta^{2} = 1 - \Lambda^{1/s}}, where \eqn{\Lambda}  
#' denotes Wilks' lambda: 
#' \deqn{
#' \Lambda = \frac{\det(E)}{\det(T)},
#' }
#' where \eqn{E} is the within-group sum of squares and cross-products matrix,
#' \eqn{H} is the between-group sum of squares and cross-products matrix and 
#'  \eqn{T} is the total sum of squares and cross-products matrix, with 
#' \eqn{T = H + E},  \eqn{s=min(p, df_h)} with \eqn{p} being the number of
#'  latent roots of \eqn{HE^{-1}}. See [heplots::etasq()] for more details.
#'
#' @return If `select_only` is `TRUE` a vector of the names of the invariant
#'  components or variables to select. 
#'  If `FALSE` an object of class `"ICS_crit"`
#'  is returned with the following objects: 
#'  - `crit`: the name of the criterion "discriminatory".
#'  - `method`: the name of the discriminatory power.
#'  - `nb_select`: the number of components to select.
#'  - `select`: the names of the invariant components or variables to select.
#'  - `power_combinations`: the discriminatory values for each of the considered
#'  combinations of `nb_select` components.
#'  - `gen_kurtosis`: the vector of generalized kurtosis values in case of 
#'  `ICS` object.
#' 
#'  
#' @export
#' 
#' @rdname discriminatory_crit
#' 
#' @references
#' Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2024). 
#' Tandem clustering with invariant coordinate selection. 
#' Econometrics and Statistics.
#'  \doi{10.1016/j.ecosta.2024.03.002}.
#'  
#' Muller, K. E. and Peterson, B. L. (1984). Practical methods for
#' computing power in testing the Multivariate General Linear Hypothesis
#' \emph{Computational Statistics and Data Analysis}, \bold{2}, 143-158.
#' 
#' @seealso [normal_crit()], [med_crit()], [var_crit()], [heplots::etasq()].
#'
#' @author Aurore Archimbaud and Anne Ruiz-Gazen
#' @import ICS
#' @examples
#' X <- iris[,-5]
#' out <- ICS(X)
#' discriminatory_crit(out, clusters = iris[,5], select_only = FALSE)

#' @export
discriminatory_crit <- function(object, ...) UseMethod("discriminatory_crit")

#' @import ICS
#' @method discriminatory_crit ICS
#' @rdname discriminatory_crit
#' @export
discriminatory_crit.ICS <- function(object, clusters, method = "eta2", 
                                    nb_select = NULL, select_only = FALSE, ...){
  gen_kurtosis <- ICS::gen_kurtosis(object, scale = FALSE)
  
  if (is.null(clusters)){
    stop("The 'clusters' argument is mandatory to compute the discriminatory 
            power of the reduced data frame.")
    out <- vector()
  }else{
    # if nb_select is NULL we put the nb_clusters-1
    nb_select <- ifelse(is.null(nb_select), length(unique(clusters))-1, nb_select)
    discriminatory_crit(ICS::components(object), 
                        clusters = clusters,
                        method = method, nb_select = nb_select, 
                        select_only = select_only, gen_kurtosis = gen_kurtosis)
  }
}


#' @export
#' @method discriminatory_crit default
#' @rdname discriminatory_crit
#' @importFrom heplots etasq
discriminatory_crit.default <- function(object, clusters, method = "eta2", 
                                        nb_select = NULL, select_only = FALSE,
                                        gen_kurtosis = NULL, ...){
  # Initialization
  method <- match.arg(method)
  nb_select <- ifelse(is.null(nb_select), ncol(object)-1, nb_select)
  # First we construct all potential combinations of first and last components
  # to analyze
  d <- ncol(object)
  IC_last <- rev(sapply(seq(0, nb_select), function(i){
    if(i == nb_select){
      IC_last <- 0
    }else{
      IC_last <- d-(seq((nb_select-(i+1)),0))
    }
  }, simplify = FALSE))
  
  IC_first <- sapply(seq(0, nb_select), function(i){seq(0,(nb_select-i))})
  
  all_comb <- sapply(seq(1, (nb_select+1)), function(i){
    IC_all = sort(c(IC_first[[i]], IC_last[[i]]))
    IC_all[IC_all != 0]
  })
  
  if(nb_select == 1){
    all_comb <- data.frame(all_comb[1], all_comb[2])
  }
  
  
  # we compute the discriminatory power for each combination to identify which
  # one has the highest power.
  if (method == "eta2"){
    all_power <- sapply(seq(1, ncol(all_comb)), function(i){
      eta2_power(object, clusters, select = all_comb[,i])
    })
    names(all_power) <- apply(all_comb, 2, function(x) 
      paste(paste("IC", x, sep = "."), collapse = ","))
    
    # We select the combination with the highest power
    ind <- which.max(all_power)
    select <- all_comb[,ind]
    power <- all_power[ind]
  }else{
    warning("No other method than 'eta2' has been implemented yet.")
    select <- vector()
    power <- 0
  }
  
  out <- colnames(object)[select]
  if (!select_only){
    out <- list(crit = "discriminatory",
                method = method, nb_select = nb_select, 
                select = out, power = power,
                power_combinations = all_power,
                gen_kurtosis = gen_kurtosis )
    class(out) <- "ICS_crit"
  }
  out
  
}




#' @importFrom stats as.formula lm manova median var
#' @noRd
eta2_power <- function(object, clusters, select){
  if(is.null(clusters)){
    warning("The 'clusters' argument is mandatory to compute the discriminatory 
            power of the reduced data frame.")
  }else{
    df <- data.frame(clusters = as.factor(clusters), object[, select])
    
    # Univariate case: ANOVA
    if(length(select) == 1){
      
      ICS_mod <- lm(as.formula(paste("cbind(", colnames(df)[-1], ") ~ clusters")),
                    data = df)
      etasq(ICS_mod)[1,1]
      
    }else{
      # Multivariate case: MANOVO with Wilks test
      
      ICS_mod <- manova(as.formula(paste("cbind(", 
                                         paste(colnames(df)[-1], collapse = ","), 
                                         ") ~ clusters")), data = df)
      etasq(ICS_mod, test = "Wilks")[1,1]
    }
  }
}
