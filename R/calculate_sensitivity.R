################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, November 2020
################################################################################

# TidyComb
# Functions for calculating cell line's sensitivity to drugs or drug combination
#
# Functions on this page:
#
# Exported:
#
# CalculateSens: Calculate sensitivity score (relative inhibition)
# ImputeIC50: Impute missing value at IC50 concentration of drug
# PredictResponse: Predict response value at certain drug dose
#
# Internal:
# scoreCurve/scoreCurve.L4: facility functions for CalculateSens
# own_log/own_log2: facility functions for CalculateSens
# CalculateIC50: Transform IC50 from coefficients from fitted dose-response model

#' Calculate relative inhibition (RI) for dose-response curve
#'
#' Function \code{CalculateSens} calculates cell line sensitivity to a drug or a
#' combination of drugs from dose response curve.
#'
#' This function measures the sensitivity by calculating the Area Under Curve
#' (AUC) according to the dose response curve. The lower bouder is chosen as
#' lowest non-zero concentration in the dose response data.
#'
#' @param df A data frame. It contains two variables:
#' \itemize{
#'   \item \strong{dose} the concentrations of drugs.
#'   \item \strong{response} the response of cell lines at crresponding doses.
#'   We use inhibition rate of cell line growth to measure the response.
#' }
#' @param pred A logical value. If it is \code{TRUE}, the function will
#' return one more table in the result. It contains the predicted response value
#' at input doses (according to fitted dose-response model) and corresponding
#' standard deviation. This table could be used in \code{TRUE}.
#'
#' \strong{Note}: The input data frame must be sorted by "dose" with ascending
#' order.
#'
#' @return If \code{pred} is \code{FALSE}, only the RI value will be return. If
#' \code{pred} is set to be \code{TRUE}, one more data frame which contains
#' predicted resposne values and corresponding standard deviations will be
#' return. It could be used to \code{\link{RIConfidenceInterval}} for confidence
#' interval calculation.
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' # LL.4
#' df <- data.frame(dose = c(0, 0.1954, 0.7812, 3.125, 12.5, 50),
#'                  response = c(2.95, 3.76, 18.13, 28.69, 46.66, 58.82))
#' RI <- CalculateSens(df)
#'
#' RI_with_pred <- CalculateSens(df, pred = TRUE)

CalculateSens <- function(df, pred = FALSE) {
  #options(show.error.messages = FALSE)
  df <- df[which(df$dose != 0),]
  if (nrow(df) == 1) {
    score <- df$response[1]
    if (pred) {
      warning("Standard deviation of relative inhibition (RI) score can not be",
      "calculated with data which contains only one dose.")
    }
    res <- score
  } else {
    tryCatch({
      model <- drc::drm(response ~ dose, data = df, fct = drc::LL.4(),
                           control = drc::drmc(errorm = FALSE, noMessage = TRUE,
                                               otrace = TRUE))
      fitcoefs <- model$coefficients
      names(fitcoefs) <- NULL
      # Calculate DSS
      score <- round(scoreCurve(d = fitcoefs[3] / 100,
                                c = fitcoefs[2] / 100,
                                b = fitcoefs[1],
                                m = log10(fitcoefs[4]),
                                c1 = log10(min(df$dose)),
                                c2 = log10(max(df$dose)),
                                t = 0), 3)
    }, error = function(e) {
      # Skip zero conc, log, drc::L.4()
      # message(e)
      model <<- drc::drm(response ~ log10(dose), data = df, fct = drc::L.4(),
                           control = drc::drmc(errorm = FALSE, noMessage = TRUE,
                                               otrace = FALSE))
      fitcoefs <- model$coefficients
      names(fitcoefs) <- NULL
      score <<- round(scoreCurve.L4(d = fitcoefs[3] / 100,
                                   c = fitcoefs[2] / 100,
                                   b = fitcoefs[1],
                                   e = fitcoefs[4],
                                   c1 = log10(min(df$dose)),
                                   c2 = log10(max(df$dose)),
                                   t = 0), 3)
    })

    if (pred) {
      # Calculate SD for DSS
      pred <- suppressWarnings(predict(model, model$data, interval="prediction"))
      pred <- cbind(model$data[,1:2], as.data.frame(pred))
      pred$sd <- (pred[,"Upper"] - pred[,"Prediction"])/(2*1.96)
      pred$sd[is.na(pred$sd)] <- sqrt(100 - pred$Prediction[is.na(pred$sd)]) # problematic?
      res <- list(RI = score, pred = pred[, c(1:3, 6)])
    } else {
      res <- score
    }
  }
  return(res)
  #Clean up
  gc()
}

#' CSS facilitate function - scoreCurve for curves fitted by LL.4
#'
#' New function used to score sensitivities given either a single-agent or a
#' fixed conc (combination) columns. The function calculates the AUC of the
#' log10-scaled dose-response curve. \strong{IMPORTANT:} note that with
#' \code{\link[drc]{LL.4()}} calls, this value is already logged since the
#' input concentrations are logged.
#'
#' @param b numeric, fitted parameter b from \code{\link[drc]{LL.4()}} model
#' @param c numeric, fitted parameter c from \code{\link[drc]{LL.4()}} model
#' @param d numeric, fitted parameter d from \code{\link[drc]{LL.4()}} model
#' @param m numeric, relative IC50 for the curve. log10(e), where e is the
#'   fitted parameter e from \code{\link[drc]{LL.4()}} model
#' @param c1 numeric, log10(min conc) (this is the minimal nonzero concentration)
#' @param c2 numeric, log10(max conc) (this is the maximal concentration)
#' @param t numeric, threshold (usually set to zero)
#'
#' @return numeric, RI or CSS scores
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
scoreCurve <- function(b, c, d, m, c1, c2, t) {
  # y <- c + (d - c) / (1 + (e / x) ^ (-b)) # LL.4
  # b <- coef[1]
  # c <- coef[2]
  # d <- coef[3]
  # e <- coef[4]
  # m <- log10(e)
  int_y <- (((((d - c) * own_log(-b, c2, m)) / ((-b) * log(10))) + c * c2) -
              ((((d - c) * own_log(-b, c1, m)) / ((-b) * log(10))) + c * c1))
  # int_y <- (((((a - d) * own_log(b, c, x2)) / (b * log(10))) + a * x2) -
  #            ((((a - d) * own_log(b, c, x1)) / (b * log(10))) + a * x1)) -
  #              (t * (x2 - x1))

  ratio <- int_y / ((1 - t) * (c2 - c1))
  sens <- ratio * 100 # scale by 100
  return(sens)
}

#' CSS facilitate function - log calculation (nature based) LL.4 version
#'
#' #' This function calculates ln(1+10^(b*(c-x))) to be used in \code{scoreCurve}
#' function
#'
#' @param b fitted parameter b from \code{\link[drc]{L.4()}} model
#' @param c fitted parameter c from \code{\link[drc]{L.4()}} model
#' @param x relative IC50 for the curve. log10(e), where e is the
#'   fitted parameter e from \code{\link[drc]{L.4()}} model
#'
#' @return ln(1+10^(b*(c-x)))
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
own_log = function(b, c, x)
{
  arg = 1 + 10^(b*(c-x))
  if(is.infinite(arg)==T) res = b*(c-x)*log(10) else res = log(arg)
  return(res)
}

#' CSS facilitate function - scoreCurve for curves fitted by L.4
#'
#' This function is used to score sensitivities given either a single-agent or a
#' fixed conc (combination) columns. The function calculates the AUC of the
#' log10-scaled dose-response curve.
#'
#' @param b numeric, fitted parameter b from \code{\link[drc]{L.4()}} model
#' @param c numeric, fitted parameter c from \code{\link[drc]{L.4()}} model
#' @param d numeric, fitted parameter d from \code{\link[drc]{L.4()}} model
#' @param e numeric, fitted parameter e from \code{\link[drc]{L.4()}} model
#' @param c1 numeric, log10(min conc) (this is the minimal nonzero concentration)
#' @param c2 numeric, log10(max conc) (this is the maximal concentration)
#' @param t numeric, threshold (usually set to zero)
#'
#' @return numeric, RI or CSS scores
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
scoreCurve.L4 <- function(b, c, d, e, c1, c2, t) {
  # y <- c + (d - c) / (1 + exp(b * (x - e))) # L4
  # b <- coef[1]
  # c <- coef[2]
  # d <- coef[3]
  # e <- coef[4]
  # m <- log10(e)

  int_y <- d * (c2 - c1) + ((c - d) / b) *
    (own_log2(b * (c2 - e)) - own_log2(b * (c1 - e)))
  ratio <- int_y / ((1 - t) * (c2 - c1))
  sens <- ratio * 100 # scale by 100
  return(sens)
}

#' CSS facilitate function - log (nature based) calculation L.4 version
#'
#' This function calculates ln(1+exp(x)) to be used in \link{scoreCurve.L4}
#' function
#'
#' @param x relative IC50 for the curve. The fitted parameter e from
#'  \code{\link[drc]{L.4()}} model
#'
#' @return ln(1+exp(x))
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
own_log2 <- function(x){
  arg = 1 + exp(x)
  if(is.infinite(arg)==T) res = x else res = log(arg)
  return(res)
}

#' Calculate Confidence Interval of RI
#'
#' Function \code{RIConfidenceInterval} is used to calculate The confidence
#' interval of RI (Relitave Inhibition) score for cell sensitivity to certain
#' drug.
#'
#' @details \code{RIConfidenceInterval} takes the prediction results from
#' \code{\link{CalculateSens}} to generate simmulated dose respons data. The
#' number of iteration is controlled by parameter \code{iter}. The RIs will
#' be computated for all the simulated data. Finally, the function will return:
#' Simulated RI (Median), and two boundaries of RI's confidence interval.
#'
#' @param pred A data frame. It must contain:
#'   \itemize{
#'     \item \strong{dose} The concentration of drugs used for model fitting
#'     and prediction.
#'     \item \strong{prediction} The predicted response value at certain
#'     doses.
#'     \item \strong{sd} The standard deviation of predicted response value.
#'   }
#' It can be generated by calling \code{\link{CalculateSens}(df, pred = TRUE)}.
#' @param iter A numeric value. It indicates the number of iterations for
#' simulation.
#'
#' @return A named numeric vector. It contains simulated RI, two boundaries of
#' the RI's confidence interval
#'
#' @author
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#' Jing Tang \email{jing.tang@helsinki.fi}
#'
#' @export
#'
#' @examples
#' df <- data.frame(dose = c(0, 0.1954, 0.7812, 3.125, 12.5, 50),
#'                  response = c(2.95, 3.76, 18.13, 28.69, 46.66, 58.82))
#'
#' RI_with_pred <- CalculateSens(df, pred = TRUE)
#'
RIConfidenceInterval <- function(pred, iter = 100){
  dose <- sort(pred$dose)
  n <- length(dose)
  score.simulated <- c()

  # Generate simmulated responses at certain doses.
  response_sim <- sapply(1:iter, function(x){
    set.seed(x)
    rnorm(n, mean = pred$Prediction, sd = pred$sd)
  })

  # Calculate RI for simmulated dose-response data frames.
  score_sim <- apply(response_sim, 2, function(x){
    df <- data.frame(dose = dose, response = x)
    tryCatch(
      score <- suppressWarnings(CalculateSens(df, pred = FALSE)),
      error = function(e){
        print(e)
        score <<- NA
      },
      finally = return(score)
    )

  })

  # confidence interval
  res <- c(simulated_RI = median(score_sim, na.rm = TRUE),
             lower = quantile(score_sim, probs = 0.025,
                              names = FALSE, na.rm = TRUE),
             upper = quantile(score_sim, probs = 0.975,
                              names = FALSE, na.rm = TRUE))

  return (res)
}
#' Impute missing value at IC50 concentration of drug
#'
#' \code{ImputeIC50} uses the particular experiment's values to predict the
#' missing values at the desired IC50 concentration of the drug.
#
#' This function is only called when trying to fix a drug at its selected IC50
#' concentration where the response values have not been tested in experiment.
#'
#' \code{ImputeIC50} fits dose-response models (with \code{\link[drc]{drm}}
#' function) by fixing the concentrations of the
#' \strong{other} drug successively, and uses each fit to predict the missing
#' value at the combination (missing IC50, fixed conc).
#'
#' @param response.mat A matrix. It contains response value of a block of drug
#' combination.
#'
#' @param row.ic50 A numeric. The IC50 value of drug added to rows.
#'
#' @param col.ic50 A numeric. The IC50 value of drug added to columns.
#'
#' @return a data frame contains all response value at the IC50 concentration
#' of certein drug. It could be directly passed to function
#' \code{CalculateSens} for scoring.
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
ImputeIC50 <- function(response.mat, col.ic50, row.ic50) {

  colconc <- as.numeric(colnames(response.mat))
  rowconc <- as.numeric(rownames(response.mat))
  n_col <- length(colconc)
  n_row <- length(rowconc)

  if (n_row == 2) {
    tempcf_c <- data.frame(dose = colconc, response = response.mat[2, ])
  } else {
    response <- apply(response.mat, 2, function(x){
          df <- data.frame(dose = rowconc, response = x)
          pred <- PredictResponse(df, row.ic50)
          return(pred)
        }
      )
    tempcf_c <- data.frame(dose = colconc, response = response)
  }

  if (n_col == 2) {
    tempcf_r <- data.frame(dose = rowconc, response = response.mat[, 2])
  } else {
    response <- apply(response.mat, 1, function(x){
        df <- data.frame(dose = colconc, response = x)
        pred <- PredictResponse(df, col.ic50)
        return(pred)
      }
    )
    tempcf_r <- data.frame(dose = rowconc, response = response)
  }

  tempres <- list(tempcf_c = tempcf_c, tempcf_r = tempcf_r)
  return(tempres)

  # Clean up
  gc()
}

CalculateIC50 <- function(coef, type, max.conc){
  if (type == "LL.4") {
    ic50 <- coef[["e_EC50"]]
  } else if (type == "L.4") {
    ic50 <- exp(coef[["e_EC50"]])
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
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
PredictResponse <- function(df, dose) {
  if (stats::var(df$response, na.rm = TRUE) == 0) {
    pred <- df$response[1]
  } else {
    model <- FitDoseResponse(df)
    if (model$call$fct[[1]][[3]] == "LL.4") {
      pred <- stats::predict(model, data.frame(dose = dose))
    } else if(model$call$fct[[1]][[3]] == "L.4") {
      pred <- stats::predict(model, data.frame(dose = log(dose)))# NB! use log
    } else {
      stop("Fitted model should be LL.4 or L.4.")
    }

    if (pred > 100) {
      pred <- 100 + stats::runif(1, -0.01, 0)
    }
  }
  return(pred)
}