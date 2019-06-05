# TidyComb
# Functions for fitting single drug dose-response curve
# Copyright Shuyu Zheng
#
# Functions on this page:
# Exported fuctions:
# * FitDoseResponse: Fitting single drug dose-response model
# Internal functions:
# * FindModelType: Get the model type which was used to fit dose-response curve.
# * CalculateIC50: Caluculate absolute IC50 from dose-response curve fitted
#                  coefficients.
# * PredictResponse: Get the predict values from fitted dose-response model.

#' Fitting single drug dose-response model
#'
#' Function \code{FitDoseResponse} fits dose-response model by using
#' \code{\link[drc]{drm}} function.
#'
#' Pre-fitting process:
#' 1. Change the 0 value in concentration into 10^-10 to avoide raising error
#' when taking log.
#' 2. If the variance of "response" values equal to 0, add 10^-10 to the last
#' "response" value.
#'
#' Model choice:
#' First use "L.4" model to fit the raw data. If error or waring occurs, use
#' "LL.4" model to fit \code{log(raw data)}.
#'
#' @param data A data frame. It contains two columns:
#' \itemize{
#'   \item \strong{conc} The concentration of drugs added in experiment.
#'   \item \strong{response} The response of cell lines to drug with different
#'   concentrations.
#' }
#'
#' @param Emin A numeric or \code{NA}. the minimal effect of the drug used in
#'    the 4-parameter log-logistic function to fit the dose-response curve. If
#'    it is not NA, it is fixed the value assigned by the user. Default setting
#'    is \code{NA}.
#'
#' @param Emax A numeric or \code{NA}. the maximal effect of the drug used in
#'    the 4-parameter log-logistic function to fit the dose-response curve. If
#'    it is not NA, it is fixed the value assigned by the user. Default setting
#'    is \code{NA}.
#'
#' @return An object of class 'drc'. It contains imformation of fitted model.
#'
#' @author \itemize{
#'    \item{Liye He \email{liye.he@helsinki.fi}}
#'    \item{Shuyu Zheng \email{shuyu.zheng@helsinki.fi}}
#' }
#'
#' @references Seber, G. A. F. and Wild, C. J (1989)
#' href{https://onlinelibrary.wiley.com/doi/book/10.1002/0471725315}{Nonlinear
#' Regression, New York}: Wiley \& Sons (p. 330).
#'
#' @export
FitDoseResponse <- function (data, Emin = NA, Emax = NA) {

  if (!all(c("dose", "response") %in% colnames(data))) {
    stop('The input must contain columns: "dose", "respone".')
  }

  # nonzero concentrations to take the log
  data$dose[which(data$dose == 0)] <- 10^-10

  if (nrow(data) != 1 & stats::var(data$response) == 0) {
    data$response[nrow(data)] <- data$response[nrow(data)] + 10 ^ -10
  }

  drug.model <- tryCatch({
    drc::drm(response ~ dose, data = data,
             fct = drc::LL.4(fixed = c(NA, Emin = Emin, Emax = Emax, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE,
                                 otrace = TRUE))
  }, warning = function(w) {
    drc::drm(response ~ log(dose), data = data,
             fct = drc::L.4(fixed = c(NA, Emin, Emax, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE,
                                 otrace = FALSE))
  }, error = function(e) {
    drc::drm(response ~ log(dose), data = data,
             fct = drc::L.4(fixed = c(NA, Emin, Emax, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE,
                                 otrace = FALSE))
  })
  return(drug.model)
}

FindModelType <- function(model) {
  type <- model$call$fct[[1]][[3]]
  return(type)
}

CalculateIC50 <- function(coef, type, max.conc){
  if (type == "LL.4") {
    ic50 <- coef[["e:(Intercept)"]]
  } else if (type == "L.4") {
    ic50 <- exp(coef[["e:(Intercept)"]])
  }

  if (ic50 > max.conc) {
    ic50 = max.conc
  }

  return (ic50)

}


#' Predict response value at certain drug dose
#'
#' \code{PredictResponse} uses \code{\link[drc]{drm}} function to fit the dose
#' response model and generate the predict response value at the given dose.
#'
#' \strong{Note}: Random number generator used in \code{AddNoise} with
#' \code{method = "random"}. If the analysis requires for reproductiblity,
#' plesase set the random seed before calling this function.
#'
#' @param df A data frame. It contains two variable:
#' \itemize{
#'   \item \strong{dose} a serial of concentration of drug;
#'   \item \strong{response} the cell line response to each concentration of
#'   drug. It should be the inhibition rate according to negative control.
#' }
#'
#' @param dose A numeric value. It specifies the dose at which user want to
#' predict the response of cell line to the drug.
#'
#' @return A numeric value. It is the response value of cell line to the drug at
#' inputted dose.
#'
#' @author Shuyu Zheng{shuyu.zheng@helsinki.fi}
#'
#' @export
PredictResponse <- function(df, dose) {
  if (stats::var(df$response, na.rm = TRUE) == 0) {
    pred <- df$response[1]
  } else {
    model <- FitDoseResponse(df)

    if (model$call$fct[[1]][[3]] == "LL.4") {
      pred <- stats::predict(model, data.frame(dose = dose))
    } else {
      pred <- stats::predict(model, data.frame(dose = log(dose)))# NB! use log
    }

    if (pred > 100) {
      pred <- 100 + stats::runif(1, -0.01, 0)
    }
  }
  return(pred)
}