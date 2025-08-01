% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_crit.R
\name{var_crit}
\alias{var_crit}
\alias{var_crit.ICS}
\alias{var_crit.default}
\title{Selection of Invariant components using the var criterion}
\usage{
var_crit(object, ...)

\method{var_crit}{ICS}(object, nb_select = NULL, select_only = FALSE, ...)

\method{var_crit}{default}(object, nb_select = NULL, select_only = FALSE, ...)
}
\arguments{
\item{object}{object of class \code{"ICS"}.}

\item{\dots}{additional arguments are currently ignored.}

\item{nb_select}{the exact number of components to select. By default it is set to
\code{NULL}, i.e the number of components to select is the number of variables minus one.}

\item{select_only}{boolean. If \code{TRUE} only the vector names of the selected
invariant components is returned. If \code{FALSE} additional details are returned.}
}
\value{
If \code{select_only} is \code{TRUE} a vector of the names of the invariant
components or variables to select. If \code{FALSE} an object of class \code{"ICS_crit"}
is returned with the following objects:
\itemize{
\item \code{crit}: the name of the criterion "var".
\item \code{nb_select}: the number of components to select.
\item \code{gen_kurtosis}: the vector of generalized kurtosis values.
\item \code{select}: the names of the invariant components or variables to select.
\item \code{RollVarX}: the rolling variances of order d-\code{nb_select}.
\item \code{Order}: indexes of the ordered invariant components such that the
ones associated to the smallest variances of the eigenvalues are at the
end.
}
}
\description{
Identifies the interesting invariant coordinates based on the rolling
variance criterion as used in the \code{ICSboot} function of the \code{ICtest}
package. It computes rolling variances on the generalized eigenvalues
obtained through \code{\link[ICS:ICS-S3]{ICS::ICS()}}.
}
\details{
Assuming that the generalized eigenvalues of the uninformative components are all the same
means that the variance of these generalized eigenvalues must be minimal.
Therefore when \code{nb_select} components should be selected, the method identifies
the \code{p - nb_select} neighboring generalized eigenvalues with minimal variance,
where \code{p} is the total number of components. The number of interesting components should be at
most \code{p-2} as at least two uninteresting components are needed to compute a variance.
}
\examples{
X <- iris[,-5]
out <- ICS(X)
var_crit(out, nb_select = 2, select_only = FALSE)

}
\references{
Alfons, A., Archimbaud, A., Nordhausen, K., & Ruiz-Gazen, A. (2024).
Tandem clustering with invariant coordinate selection.
Econometrics and Statistics.
\doi{10.1016/j.ecosta.2024.03.002}.

Radojicic, U., & Nordhausen, K. (2019).
Non-gaussian component analysis: Testing the dimension of the signal subspace.
In Workshop on Analytical Methods in Statistics (pp. 101–123). Springer.
\doi{10.1007/978-3-030-48814-7_6}.
}
\seealso{
\code{\link[=normal_crit]{normal_crit()}}, \code{\link[=med_crit]{med_crit()}}, \code{\link[=discriminatory_crit]{discriminatory_crit()}}.
}
\author{
Andreas Alfons, Aurore Archimbaud and Klaus Nordhausen
}
