# SynergyFinder
# Functions for processing drug response matrix.
#
# Functions in this page:
#
# AdjustColumnName: Adjust column names of input data table
# ReshapeData: pre-process the response data for further calculation and plot.
# ImputeNear: Impute missing value with nearest values
# AddNoise: Add noise to response value
# ExtractSingleDrug: Extract single drug response from matrix
# CorrectBaseLine: Do base line correction to dose-response matrix.


#' Adjust column names of input data table
#'
#' This function changes the column names in other format into the style:
#' block_id, drug1, drug2, conc1, conc2, response, conc_unit1, conc_unit2.
#'
#' @param data
#'
#' @return The data frame with the changed column names.
#'
#' @export
AdjustColumnName <- function(data) {
  colnames <- colnames(data)
  colnames <- tolower(gsub("([a-z])([A-Z])", "\\1_\\L\\2",
    colnames,
    perl = TRUE
  ))
  colnames <- gsub("^conc_col$", "conc2", colnames, perl = TRUE)
  colnames <- gsub("^conc_row$", "conc1", colnames, perl = TRUE)
  colnames <- gsub("^pair_index$", "block_id", colnames, perl = TRUE)

  colnames <- gsub("^drug_row$", "drug1", colnames, perl = TRUE)
  colnames <- gsub("^drug_col$", "drug2", colnames, perl = TRUE)
  colnames <- gsub("^conc_r$", "conc1", colnames, perl = TRUE)
  colnames <- gsub("^conc_c$", "conc2", colnames, perl = TRUE)
  colnames <- gsub("^conc_unit$", "conc_unit1", colnames, perl = TRUE)
  colnames <- gsub("^conc_r_unit$", "conc_unit1", colnames, perl = TRUE)
  colnames <- gsub("^conc_c_unit$", "conc_unit2", colnames, perl = TRUE)


  colnames(data) <- colnames
  return(data)
}

#' Pre-process the response data for further calculation and plot
#'
#' A function to transform the response data from data frame format to
#' dose-response matrices. Several processes could be chose to add noise, impute
#' missing values or correct base line to the dose-response matrix.
#'
#' @details The input data must contain the following columns: (block_id/BlockId/PairIndex),
#' (drug_row/DrugRow/Drug1), (drug_col/DrugCol/Drug2), (response/Response),
#' (conc_r/ConcRow/Conc1), (conc_c/ConcCol/Conc2), and (ConcUnit/conc_r_unit,
#' conc_c_unit/ConcUnit1, ConcUnit2, ConcUnit3)
#'
#' @param data drug combination response data in a data frame format
#' @param impute a logical value. If it is \code{TRUE}, the \code{NA} values
#'   will be imputed by \code{\link{ImputeNA}}. Default is \code{TRUE}.
#' @param noise a logical value. It indicates whether or not adding noise to
#'   to the "response" values in the matrix. Default is \code{TRUE}.
#' @param seed a single value, interpreted as an integer, or NULL. It is the
#'   random seed for calculating the noise. Default setting is \code{NULL}
#' @param correction a character. This argument is extended from the argument
#'   \code{method} of \code{\link{CorrectBaseLine}} function. There are three
#'   available valuse: \code{non}, \code{part}, \code{all}.
#'   The default setting is \code{non}.
#' @param data.type a parameter to specify the response data type which can be
#'   either "viability" or "inhibition".
#'
#' @return a list of the following components:
#'   \itemize{
#'     \item \strong{dose.response.mats} a list of the dose-response matrices
#'       with \%inhibition as the response data. Row names and column names are
#'       drug concentrations.
#'   \item \strong{adjusted.response.mats} The dose response matrix adjusted.
#'     The processes are chosen by arguments \code{impute}, \code{noise}, and
#'     \code{correction}. If no process was chosen, the final result will not
#'     contain this result.
#'   \item \strong{drug.pairs} a data frame contains the name of the row drug,
#'     the name of the column drug, concentration unit and block IDs.
#'   }
#'
#' @author
#'   \itemize{
#'     \item Liye He \email{liye.he@helsinki.fi}
#'     \item Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'  }
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' data("mathews_screening_data")
#' # set a random number seed for generating the noises
#' set.seed(1)
#' data <- ReshapeData(mathews_screening_data)
ReshapeData <- function(data, impute = TRUE, noise = TRUE, seed = NULL,
                        correction = "non", data.type = "viability") {
  data <- AdjustColumnName(data)
  # 1.1 check column names

  if (!all(c(
    "block_id", "drug1", "drug2", "response", "conc1", "conc2",
    "conc_unit1"
  ) %in%
    colnames(data))) {
    stop(
      "The input data must contain the following columns: (block_id/BlockId/PairIndex), ",
      "(drug_row/DrugRow/Drug1/drug1), (drug_col/DrugCol/Drug2/drug2), (response/Response),",
      "(conc_r/ConcRow/Conc1/conc1), (conc_c/ConcCol/Conc2/conc2), 
         (ConcUnit/conc_r_unit/ConcUnit1/conc_unit1)"
    )
  }

  drugs <- grep("drug\\d+", colnames(data), value = TRUE)
  conc_units <- gsub("conc_unit", "drug", drugs)

  for (i in 1:length(conc_units)) {
    if (!conc_units[i] %in% colnames(data)) {
      data[conc_units[i]] <- data$conc_unit1
    }
  }

  # 1.2 Check missing values
  if (!impute & sum(is.na(data$response))) {
    stop(
      "There are missing values in input data. Please fill it up or run 'ReshapeData' ",
      "with 'impute=TRUE'."
    )
  }

  # 2. Split data
  data <- list(
    drug.pairs = unique(data[, !grepl("(response|conc\\d+)",
      colnames(data),
      perl = TRUE
    )]),
    response.df = unique(data[, grepl("(block_id|response|conc\\d+)",
      colnames(data),
      perl = TRUE
    )])
  )

  # 3. Check the data type
  if (data.type == "viability") {
    data$response.df$response <- 100 - data$response.df$response
  } else if (data.type == "inhibition") {
    data$response.df$response <- data$response.df$response
  } else {
    stop(
      "Please tell me the data type of response valuse: 'viability' or ",
      "'inhibition'."
    )
  }

  # 4. Dealing with replicates
  data$replicate <- sum(duplicated(dplyr::select(data$response.df, -response))) > 0
  if (data$replicate) {
    data <- repResponse(data)
  } else {
    data$drug.pairs$replicate <- FALSE
  }

  # 5. Mark multi-drug data
  multidrug <- sum(grepl("^drug\\d$", colnames(data$drug.pair), perl = TRUE)) > 2
  data$multidrug <- multidrug

  # 6. Adjust dose response values
  if (multidrug) { # more than two drugs combination
    if (noise) {
      if (data$replicate) {
        data$replicate.response$response_adj <- data$replicate.response$response +
          stats::rnorm(nrow(data$replicate.response), 0, 0.001)
      } else {
        data$response.df$response_adj <- data$response.df$response +
          stats::rnorm(nrow(data$response.df), 0, 0.001)
      }
    }
  } else { # 2 drugs combination
    if (impute | noise | correction != "non") {
      tmp.df <- NULL
      if (data$replicate) {
        response.df <- data$replicate.response
      } else {
        response.df <- data$response.df
      }
      blocks <- unique(data$drug.pairs$block_id)
      for (b in data$drug.pairs$block_id) {
        tmp.mat <- response.df %>%
          dplyr::filter(block_id == b)
        tmp.mat <- reshape2::acast(conc1 ~ conc2,
          data = tmp.mat,
          value.var = "response"
        )
        # process data according to setting of arguments

        if (impute) {
          tmp.mat <- ImputeNA(tmp.mat)
        }
        if (noise) {
          set.seed(seed)
          tmp.mat <- AddNoise(tmp.mat)
        }
        tmp.mat <- CorrectBaseLine(tmp.mat, method = correction)
        # change it back to data frame
        tmp.mat <- reshape2::melt(tmp.mat)
        colnames(tmp.mat) <- c("conc1", "conc2", "response_adj")
        tmp.mat$block_id <- b
        tmp.df <- rbind.data.frame(tmp.df, tmp.mat)
      }
      if (data$replicate) {
        data$replicate.response <- data$replicate.response %>%
          dplyr::left_join(tmp.df, by = c("block_id", "conc1", "conc2"))
      }
      data$response.df <- data$response.df %>%
        dplyr::left_join(tmp.df, by = c("block_id", "conc1", "conc2"))
    }
  }
  return(data)
}

#' Impute missing value with nearest values
#'
#' Function \code{ImputeNA} does missing value imputation by assigning the
#' average of values in nearest 4 cells (top, bottom, left, right) to the NA
#' cell. This pocess will be done repeadly until there is no missing values in
#' the matrix.
#'
#' @param response.mat A matrix which has missing value.
#'
#' @return A matrix which is same as input matrix except the missing values are
#' imputed.
#'
#' @author Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
#' response.mat <- data$dose.response.mats[[1]]
#' # introduce some NA values into matrix
#' response.mat[3:4, 3:5] <- NA
#' adjusted.mat <- ImputeNA(response.mat)
ImputeNA <- function(response.mat) {
  while (sum(is.na(response.mat))) {
    x <- array(c(
      rbind(response.mat[-1, ], NA),
      rbind(NA, response.mat[-nrow(response.mat), ]),
      cbind(response.mat[, -1], NA),
      cbind(NA, response.mat[, -ncol(response.mat)])
    ),
    dim = c(nrow(response.mat), ncol(response.mat), 4)
    )
    x.imp <- apply(x, c(1, 2), mean, na.rm = TRUE)
    index.na <- is.na(response.mat)
    response.mat[index.na] <- x.imp[index.na]
  }
  return(response.mat)
}

#' Add noise to response value
#'
#' Function \code{AddNoise} calculates and add a noise to values in response
#' matrix. The noises obey normal distribution ~N(0, 0.001) wich are generated
#' by fucntion \code{rnorm}.
#'
#' \strong{Note}: If the analysis requires for reproductiblity, plesase set the
#' random seed before calling this function.
#'
#' @param response.mat A matrix. It contains the response data for one drug
#' combination.
#'
#' @return A matrix. It contains the response value added with noises.
#'
#' @author Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
#' response.mat <- data$dose.response.mats[[1]]
#' set.seed(1)
#' adjusted.mat <- AddNoise(response.mat)
AddNoise <- function(response.mat) {
  noise <- matrix(stats::rnorm(
    nrow(response.mat) * ncol(response.mat),
    0, 0.001
  ),
  nrow = nrow(response.mat),
  ncol = ncol(response.mat)
  )
  response.mat <- response.mat + noise
  return(response.mat)
}

#' Extract single drug response from matrix
#'
#' \code{ExtractSingleDrug} extracts the dose-response values of single drug
#'  from a drug combination dose-response matrix.
#'
#' @param response A data frame. It must contain the columns: "conc1", "conc2",
#' ..., for the concentration of the combined drugs and "response" for the
#' observed %inhibition at certain combination.
#'
#' @return A list contains several data frames each of which contains two columns:
#'   \itemize{
#'     \item \strong{dose} The concertration of drug.
#'     \item \strong{response} The cell's response (inhibation rate) to
#'       corresponding drug concertration.
#' }
#'
#' @author Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
#' response <- data$response[data$response$block_id == 1, c("conc1", "conc2", "response")]
#' single <- ExtractSingleDrug(response.df)
ExtractSingleDrug <- function(response) {
  concs <- grep("(conc)", colnames(response), perl = TRUE, value = TRUE)
  single.drug <- vector("list", length(concs))
  names(single.drug) <- concs
  conc_sum <- rowSums(response[, concs])
  for (conc in concs) {
    index <- which(response[, conc] == conc_sum)
    single.drug[[conc]] <- data.frame(
      dose = response[index, conc],
      response = response[index, "response"]
    )
  }
  return(single.drug)
}

#' Base line correction
#'
#' \code{CorrectBaseLine} adjusts the base line of drug combination
#' dose-response matrix to make it closer to 0.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#'   and row name are representing the concerntrations of drug added to column
#'   and row, respectively. The values in matrix indicate the inhibition rate to
#'   cell growth.
#' @param method A character value to indicate using which method to do
#'   baseline correction. Available values ate:
#'   \itemize{
#'     \item \strong{non} means no baseline corection.
#'     \item \strong{part} means only adjust the negative values in the matrix.
#'     \item \strong{all} means adjust all values in the matrix.
#'   }
#'
#' @return A matrix which base line have been adjusted.
#'
#' @author \itemize{
#'    \item{Liye He \email{liye.he@helsinki.fi}}
#'    \item{Shuyu Zheng \email{shuyu.zheng@helsinki.fi}}
#' }
#'
#' @export
#'
#' @examples
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
#' response.mat <- data$dose.response.mats[[1]]
#' adjusted.mat <- CorrectBaseLine(response.mat, method = "part")
CorrectBaseLine <- function(response.mat, method = c("non", "part", "all")) {
  method <- match.arg(method)

  if (method == "non") {
    return(response.mat)
  } else if (method == "part") {
    negative.ind <- which(response.mat < 0, arr.ind = TRUE)
    if (length(negative.ind) == 0) {
      return(response.mat)
    }
    drug.row <- ExtractSingleDrug(response.mat, dim = "row")
    drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row)))

    drug.col <- ExtractSingleDrug(response.mat, dim = "col")
    drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col)))
    drug.3 <- ExtractSingleDrug(response.mat, dim = "col")
    drug.3.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col)))
    baseline <- min(c(
      min(as.numeric(drug.row.fit)),
      min(as.numeric(drug.col.fit))
    ))
    response.mat[negative.ind] <- vapply(
      response.mat[negative.ind],
      function(x) {
        x - ((100 - x) / 100 * baseline)
      }, numeric(1)
    )
    return(response.mat)
  } else if (method == "all") {
    drug.row <- ExtractSingleDrug(response.mat, dim = "row")
    drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row)))

    drug.col <- ExtractSingleDrug(response.mat, dim = "col")
    drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col)))

    baseline <- min(c(
      min(as.numeric(drug.row.fit)),
      min(as.numeric(drug.col.fit))
    ))
    response.mat <- response.mat - ((100 - response.mat) / 100 * baseline)
    return(response.mat)
  }
}
