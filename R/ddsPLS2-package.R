#' @details
#' The only function you're likely to need from ddsPLS2 is \code{\link{ddsPLS}}, a function which performs:
#' \itemize{
#'   \item{Bootstrap operations to optimally select the regularization coefficients per component.}
#'   \item{Automatically select the optimal number of components to reduce the over-fitting of the model.}
#'   \item{Build the optimal sparse PLS model, denoted as data-driven sparse PLS (`ddsPLS`) model.}
#' }
#' Apart from this, you would certainly need to use the three S3 methods which are:
#' \describe{
#'  \item{\code{\link{summary.ddsPLS}}}{To get principal variance and prediction descriptors over the model built.}
#'  \item{\code{\link{predict.ddsPLS}}}{To predict the `y` values for a new `X` data-set.}
#'  \item{\code{\link{plot.ddsPLS}}}{To get different visualizations corresponding to your model.}
#' }
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"
