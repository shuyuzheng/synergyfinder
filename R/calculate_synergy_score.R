# TidyComb
# Functions for calculating drug synergy scores.
#
# Functions in this page:
#
# CalculateZIP/Bliss/HSA/Loewe: Four functions to calculate synergy scores.
# eq.LL4/L4.LL4/L4: Four functions to calculate loewe score in CalculateLoewe.
# fun: Function used in CalculateLoewe

#' Calculate Delta synergy score based on ZIP model
#'
#' \code{CalculateZIP} calculates the \eqn{\Delta} score matrix for a block of
#' drug combination by using Zero Interaction Potency (ZIP) method.
#'
#' zero interaction potency (ZIP) is a reference model for evaluating the
#' interaction between two drugs. It captures the drug interaction relationships
#' by comparing the change in the potency of the dose-response curves between
#' individual drugs and their combinations. More details about this model could
#' be found in original publication.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @param drug.col.model (optional) a character. It indicates the model used for
#' fitting dose-response curve for drug added to columns.
#'
#' @param drug.row.model (optional) a character. It indicates the model type
#' used for fitting dose-response curve for drug added to rows.
#'
#' @param ... Other arguments from nested functions.
#' @author \itemize{
#'    \item{Liye He \email{liye.he@helsinki.fi}}
#'    \item{Jing Tang \email{jing.tang@helsinki.fi}}
#'    \item{Shuyu Zheng \email{shuyu.zheng@helsinki.fi}}
#' }
#' @references \itemize{
#'    \item{Yadav B, Wennerberg K, Aittokallio T, Tang J. (2015).
#'    \href{https://doi.org/10.1016/j.csbj.2015.09.001}{Searching for Drug
#'    Synergy in Complex Dose-Response Landscape Using an Interaction Potency
#'    Model.} Comput Struct Biotechnol J, 13:504– 513.}
#' }
#' @return A matrix of \eqn{\Delta} score calculated via Zero Interaction
#' Potency (ZIP) method.
#'
#' @export
CalculateZIP <- function(response.mat, drug.row.model = NULL, drug.col.model = NULL, ...) {
  if (is.null(drug.row.model)) {
    drug.row <- ExtractSingleDrug(response.mat, dim = "row")
    drug.row.model <- FitDoseResponse(drug.row)
  }
  drug.row.fit <- suppressWarnings(stats::fitted(drug.row.model))

  if (is.null(drug.col.model)) {
    drug.col <- ExtractSingleDrug(response.mat, dim = "col")
    drug.col.model <- FitDoseResponse(drug.col)
  }
  drug.col.fit <- suppressWarnings(stats::fitted(drug.col.model))

  n.row <- nrow(response.mat)
  n.col <- ncol(response.mat)

  # generate drug_row fitting matrix
  tmp <- data.frame(dose = as.numeric(rownames(response.mat))[-1])
  updated.col.mat <- matrix(nrow = n.row - 1, ncol = n.col - 1)

  for (i in 2:n.col) {
    # nonzero concentrations to take the log
    tmp$response <- response.mat[-1, i]
    if(nrow(tmp) == 1) {
      # # no fitting
      fitted.response <- tmp$response - 10 ^ -10
    } else {
      tmp.min <- drug.col.fit[i]
      tmp.model <- FitDoseResponse(data = tmp, Emin = tmp.min, Emax = 100)
      if (!is.null(tmp.model$convergence)){
        fitted.response <- tmp$response
      } else {
        fitted.response <- suppressWarnings(stats::fitted(tmp.model))
      }
    }

    # if (fitted.inhibition[length(fitted.inhibition)] < 0)
    #  fitted.inhibition[length(fitted.inhibition)] <- tmp.min
    updated.col.mat[, i - 1] <- fitted.response
  }

  # generate drug_col fitting matrix
  tmp <- data.frame(dose = as.numeric(colnames(response.mat))[-1])
  updated.row.mat <- matrix(nrow = n.row - 1, ncol = n.col - 1)

  for (i in 2:n.row) {
    # nonzero concentrations to take the log
    tmp$response <- response.mat[i, -1]
    if(nrow(tmp) == 1) {
      # # no fitting
      fitted.response <- tmp$response - 10 ^ -10
    } else {
      tmp.min <- drug.row.fit[i]
      tmp.model <- FitDoseResponse(data = tmp, Emin = tmp.min, Emax = 100)
      if (!is.null(tmp.model$convergence)){
        fitted.response <- tmp$response
      } else {
        fitted.response <- suppressWarnings(stats::fitted(tmp.model))
      }
    }
    # if (fitted.inhibition[length(fitted.inhibition)] < 0)
    #  fitted.inhibition[length(fitted.inhibition)] <- tmp.min
    updated.row.mat[i - 1 , ] <- fitted.response
  }

  fitted.mat <- (updated.col.mat + updated.row.mat) / 2

  zip.mat <- matrix(nrow = n.row -1, ncol = n.col - 1)
  for (i in 1:(n.row - 1)) {
    for (j in 1:(n.col - 1)) {
      zip.mat[i, j] <- drug.row.fit[i + 1] + drug.col.fit[j + 1] -
        drug.row.fit[i + 1] * drug.col.fit[j + 1] / 100
    }
  }

  delta.mat <- fitted.mat - zip.mat

  # add 0 to first column and first row
  delta.mat <- rbind(rep(0, times = n.row -1), delta.mat)
  delta.mat <- cbind(rep(0, times = n.col), delta.mat)

  # add colname and rowname for delta.mat
  colnames(delta.mat) <- colnames(response.mat)
  rownames(delta.mat) <- rownames(response.mat)

  return(delta.mat)

  # clean up
  gc()
}

#' Calculate Bliss synergy score
#'
#' \code{CalculateBliss} calculates the synergy score matrix for a block of
#' drug combination by using a druginteraction reference model introduced by
#' C. I. Bliss in 1939.
#'
#' This model is a reference model for evaluating the interaction between two
#' drugs. The basic assumption of this model is "The expected effect of two
#' drugs acting independently". More details about this model could be found in
#' original publication:
#' .
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @return A matrix for synergy score calculated via reference model introduced
#' by C. I. Bliss.
#'
#' @author \itemize{
#'    \item{Liye He \email{liye.he@helsinki.fi}}
#'    \item{Shuyu Zheng \email{shuyu.zheng@helsinki.fi}}
#' }
#'
#' @references \itemize{
#'    \item{Yadav B, Wennerberg K, Aittokallio T, Tang J. (2015).
#'    \href{https://doi.org/10.1016/j.csbj.2015.09.001}{Searching for Drug
#'    Synergy in Complex Dose-Response Landscape Using an Interaction Potency
#'    Model.} Comput Struct Biotechnol J, 13:504– 513.}
#'    \item{Bliss, C. I. (1939).
#'    \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1744-7348.1939.tb06990.x}{The
#'    toxicity of poisons applied jointly.} Annals of Applied Biology,
#'    26(3):585–615.}
#' }
#'
#' @export
CalculateBliss <- function (response.mat) {
  drug.row <- response.mat[, 1]
  drug.col <- response.mat[1, ]
  reference.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      reference.mat[i, j] <- drug.row[i] + drug.col[j] -
        drug.row[i] * drug.col[j]/100
    }
  }
  synergy.mat <- response.mat - reference.mat

  return(synergy.mat)

  # clean up
  gc()
}

#' Calculate HSA synergy score
#'
#' \code{CalculateHSA} calculates the synergy score matrix for a block of
#' drug combination by using Highest Single Agent (HSA) reference model.
#'
#' This model is a reference model for evaluating the interaction between two
#' drugs. The basic assumption of this model is "The reference effect of drug
#' combination is the maximal single drug effect".
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @return A matrix for synergy score calculated via Highest Single Agent (HSA).
#'
#' @author \itemize{
#'    \item{Liye He \email{liye.he@helsinki.fi}}
#'    \item{Shuyu Zheng \email{shuyu.zheng@helsinki.fi}}
#' }
#'
#' @references \itemize{
#'    \item{Yadav B, Wennerberg K, Aittokallio T, Tang J.(2015).
#'    \href{https://doi.org/10.1016/j.csbj.2015.09.001}{Searching for Drug
#'    Synergy in Complex Dose-Response Landscape Using an Interaction Potency
#'    Model.}Comput Struct Biotechnol J, 13:504– 513.}
#'    \item{Berenbaum MC. (1989).
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/2692037}{What is synergy?}
#'    Pharmacol Rev 1990 Sep;41(3):422.
#'    }
#' }
#'
#' @export
CalculateHSA <- function(response.mat) {
  drug.row <- response.mat[, 1]
  drug.col <- response.mat[1, ]
  reference.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      reference.mat[i, j] <- max(drug.row[i], drug.col[j])
    }
  }
  synergy.mat <- response.mat - reference.mat

  return(synergy.mat)

  #clean up
  gc()
}

# Four functions to calculate loewe
eq.LL4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / (drug.col.par[4] * (((x - drug.col.par[3]) /
                            (drug.col.par[2] - x)) ^ (1/drug.col.par[1]))) +
    x2 / (drug.row.par[4] * (((x - drug.row.par[3]) /
                              (drug.row.par[2] - x)) ^ (1/drug.row.par[1]))) - 1
}# Eq.8 in the ZIP paper

eq.L4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / exp((drug.col.par[4] + log((drug.col.par[3] - x) /
                                    (x - drug.col.par[2])) / drug.col.par[1])) +
    x2 / exp((drug.row.par[4] + log((drug.row.par[3] - x) /
                                  (x - drug.row.par[2])) / drug.row.par[1])) -1
}# x1, x2 to be log scaled

eq.LL4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / (drug.col.par[4] * (((x - drug.col.par[3]) /
                              (drug.col.par[2] - x)) ^ (1 / drug.col.par[1]))) +
    x2 / exp((drug.row.par[4] + log((drug.row.par[3] - x) /
                                  (x - drug.row.par[2])) / drug.row.par[1])) -1
}# x2 to be log-scaled

eq.L4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / exp((drug.col.par[4] + log((drug.col.par[3] - x) /
                                    (x - drug.col.par[2])) / drug.col.par[1])) +
    x2 / (drug.row.par[4] * (((x - drug.row.par[3]) /
                            (drug.row.par[2] - x)) ^ (1 / drug.row.par[1]))) - 1
}# x1 to be log-scaled

# function used to calculate loewe if termination code from 'nleqslv' function
# is -10, 1, or 2, which mean:
# * -10 User supplied Jacobian is most likely incorrect.
# * 1 Function criterion is near zero. Convergence of function values has been
# achieved.
# * 2 x-values within tolerance. This means that the relative distance between
# two consecutive x-values is smaller than xtol but that the function value
# criterion is still larger than ftol. Function values may not be near zero;
# therefore the user must check if function values are acceptably small.
#
fun <- function(col_conc, row_conc, drug.par, model) {
  # LL.4, conc must be raw
  if(model == "LL.4") {
    conc = col_conc + row_conc
    (drug.par[3] + drug.par[2] *
        (conc / drug.par[4]) ^ drug.par[1]) /
      (1 + (conc / drug.par[4]) ^ drug.par[1])
  } else if (model == "L.4"){# L.4, conc must be logscaled, ie. log(conc)
    conc = log(col_conc+row_conc)
    (drug.par[2] + (drug.par[3] - drug.par[2]) /
        (1 + exp(drug.par[1] * (conc - drug.par[4]))))
  } else {
    stop("Model type is incorrect. Available values are 'LL.4' or 'L.4' ")
  }
}

#' Calculate Loewe synergy score
#'
#' \code{CalculateLoewe} calculates the synergy score matrix for a block of
#' drug combination by using a druginteraction reference model introduced by
#' Loewe in 1953.
#'
#' This model is a reference model for evaluating the interaction between two
#' drugs. The basic assumption of this model is "The referece effect of drug
#' combination is the expected effect of a drug combined with itself".
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @param quiet A logical value. If it is \code{TRUE} then the warning message
#' will not show during calculation.
#'
#' @param drug.col.type (optional) a character. It indicates the model type used
#' for fitting dose-response curve for drug added to columns.
#'
#' @param drug.row.type (optional) a character. It indicates the model type used
#' for fitting dose-response curve for drug added to rows.
#'
#' @param drug.col.par (optional) a named vector. It contains the coeficients of
#' fitted dose-response model for drug added to columns.
#'
#' @param drug.row.par (optional) a named vector. It contains the coeficients of
#' fitted dose-response model for drug added to rows.
#'
#' @param ... Other arguments from nested functions.
#'
#' @return A matrix for Synergy score calculated via reference model introduced
#' by Loewe, S.
#'
#' @author \itemize{
#'    \item{Liye He \email{liye.he@helsinki.fi}}
#'    \item{Jing Tang \email{jing.tang@helsinki.fi}}
#'    \item{Shuyu Zheng \email{shuyu.zheng@helsinki.fi}}
#' }
#'
#' @references \itemize{
#'    \item{Yadav B, Wennerberg K, Aittokallio T, Tang J.(2015).
#'    \href{https://doi.org/10.1016/j.csbj.2015.09.001}{Searching for Drug
#'    Synergy in Complex Dose-Response Landscape Using an Interaction Potency
#'    Model.}Comput Struct Biotechnol J, 13:504– 513.}
#'    \item{[Loewe, 1953] Loewe, S. (1953).
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/13081480}{The problem of
#'    synergism and antagonism of combined drugs.} Arzneimittelforschung,
#'    3(6):285–290.
#'    }
#' }
#'
#' @export
CalculateLoewe <- function (response.mat, quiet = TRUE, drug.col.type = NULL,
                            drug.row.type = NULL, drug.col.par = NULL,
                            drug.row.par = NULL, ...) {
  if (quiet) {
    options(warn = -1)
  }

  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.col <- ExtractSingleDrug(response.mat, dim = "col")

  con <- sapply(list(drug.col.type, drug.col.par, drug.row.type, drug.row.par),
                is.null)

  if (!all(!con)) {

    drug.row.model <- FitDoseResponse(drug.row)
    drug.row.par <- stats::coef(drug.row.model)
    drug.row.type <- FindModelType(drug.row.model)

    drug.col.model <- FitDoseResponse(drug.col)
    drug.col.par <- stats::coef(drug.col.model)
    drug.col.type <- FindModelType(drug.col.model)
  }

  drug.row$dose[drug.row$dose == 0] = 10^-10 # avoid log(0)
  drug.col$dose[drug.col$dose == 0] = 10^-10 # avoid log(0)

  loewe.mat <- response.mat
  eq <- switch (paste(drug.col.type, drug.row.type),
                "LL.4 LL.4" = eq.LL4.LL4,
                "L.4 L.4"   = eq.L4.L4,
                "LL.4 L.4"  = eq.LL4.L4,
                "L.4 LL.4"  = eq.L4.LL4)

  x <- max(drug.col.par[2], drug.row.par[2]) + 1

  for (i in 1:(nrow(drug.col) - 1)) {
    for (j in 1:(nrow(drug.row) - 1)) {
      x1 <- drug.col$dose[i + 1]
      x2 <- drug.row$dose[j + 1]

      options(warn = -1)
      slv <- tryCatch({
        slv <- nleqslv::nleqslv(x, eq, method = "Newton", x1=x1, x2=x2,
                                drug.col.par = drug.col.par,
                                drug.row.par = drug.row.par)
        }, error = function(e){
          slv <- list(termcd = 999)
        }
      )

      if (slv$termcd < 3) {
        y.loewe <- slv$x
      } else {
        y.loewe1 <- fun(x1, x2, drug.par = drug.col.par,
                        model = drug.col.type) # x1 col, x2 row
        y.loewe2 <- fun(x1, x2, drug.par = drug.row.par,
                        model = drug.row.type) # x1 col, x2 row
        y.loewe <- max(y.loewe1, y.loewe2)
      }

      loewe.mat[j + 1, i + 1] <- ifelse(y.loewe > 100, 100, y.loewe)
    }
  }

  synergy.mat <- response.mat - loewe.mat

  # Output results
  return(synergy.mat)

  options(warn = 0)
  # clean up
  gc()
}