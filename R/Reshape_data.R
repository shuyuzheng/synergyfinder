# SynergyFinder
# Functions for processing drug response matrix.
#
# Functions in this page:
#
# ReshapeData: pre-process the response data for further calculation and plot.
# ImputeNear: Impute missing value with nearest values
# AddNoise: Add noise to response value
# ExtractSingleDrug: Extract single drug response from matrix
# CorrectBaseLine: Do base line correction to dose-response matrix.

#' Pre-processe the response data for furhter calculation and plot
#'
#' A function to transform the response data from data frame format to 
#' dose-response matrices. Several processes could be chose to add noise, impute
#' missing values or correct base line to the dose-response matrix.
#'   
#' @details The input data must contain the following columns: block_id, 
#'   drug_row, drug_col, response, conc_r, conc_c, conc_r_unit, conc_c_unit.
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
#' @export
#' 
#' @examples
#' data("mathews_screening_data")
#' # set a random number seed for generating the noises
#' set.seed(1) 
#' data <- ReshapeData(mathews_screening_data)
ReshapeData <- function(data, impute=TRUE, noise=TRUE, seed = NULL, 
                        correction = "non", data.type = "viability") {
  # 1. Check the input data
  # Adjust column names
  colnames(data) <- tolower(gsub("([a-z])([A-Z])", "\\1_\\L\\2", 
                                 colnames(data), perl = TRUE))
  colnames(data) <- gsub("conc_col", "conc_c", colnames(data), perl = TRUE)
  colnames(data) <- gsub("conc_row", "conc_r", colnames(data), perl = TRUE)
  # 1.1 check column names
  if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r", "conc_c",
             "conc_r_unit", "conc_c_unit") %in%
           colnames(data)))
    stop("The input data must contain the following columns: block_id, ",
         "drug_row, drug_col, response, conc_r, conc_c, conc_r_unit,",
         "conc_c_unit")
  # 1.2 Check the data type
  if (data.type == "viability") {
    data$response <- 100 - data$response
  } else if (data.type == "inhibition") {
    data <- data
  } else {
    stop("Please tell me the data type of response valuse: 'viability' or ",
         "'inhibition'.")
  }
  # 1.3 Check missing values
  if (!impute & sum(is.na(data$response))) {
    stop("There are missing values in input data. Please run 'ReshapeData' ", 
         "with 'impute=TRUE'.")
  }

  # obtain block IDs
  blocks <- unique(data$block_id)
  
  # 2. Create containers
  # 2.1 List dose.response.mats for storing all the dose-response matrices.
  #     Setting the block_id as the name of each element.
  dose.response.mats <- vector(mode="list", length=length(blocks))
  names(dose.response.mats) <- blocks
  
  # 2.2 List adjusted.response.mats for storing all the adjusted dose-response 
  #     matrices. Setting the block_id as the name of each element.
  if (impute | noise | correction != "non") {
    adjusted.response.mats <- vector(mode="list", length=length(blocks))
    names(adjusted.response.mats) <- blocks
  }
  
  # 2.3 Data frame drug.pairs for storing all the drug name, concentration unit.
  drug.pairs <- unique(data[, colnames(data) %in% c("block_id", "drug_row", 
                                                    "drug_col", "conc_r_unit", 
                                                    "conc_c_unit")])
  
  # 3. Reshape the data
  for (block in blocks) {
    tmp.mat <- data[data$block_id == block, ]
    block <- as.character(block)

    # response matrix for one drug combination
    response.mat <- reshape2::acast(conc_r ~ conc_c, data = tmp.mat, 
                          value.var = "response")
    # save dose-response matrix
    dose.response.mats[[block]] <- response.mat
    # process data according to setting of arguments
    if (impute | noise | correction != "non") {
      if (impute) {
        response.mat <- ImputeNA(response.mat)
      }
      if (noise){
        set.seed(seed)
        response.mat <- AddNoise(response.mat)
      } 
        response.mat <- CorrectBaseLine(response.mat, method = correction)
      adjusted.response.mats[[block]] <- response.mat
    }
  }

  if (impute | noise | correction != "non") {
    return(list(dose.response.mats = dose.response.mats, 
                adjusted.response.mats = adjusted.response.mats,
                drug.pairs = drug.pairs))
  } else {
    return(list(dose.response.mats = dose.response.mats, 
                drug.pairs = drug.pairs))
  }
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
  x <- array(c(rbind(response.mat[-1,], NA),
               rbind(NA, response.mat[-nrow(response.mat), ]),
               cbind(response.mat[,-1], NA),
               cbind(NA, response.mat[, -ncol(response.mat)])),
             dim=c(nrow(response.mat), ncol(response.mat), 4))
  x.imp <- apply(x, c(1,2), mean, na.rm = TRUE)
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
    noise <- matrix(stats::rnorm(nrow(response.mat) * ncol(response.mat),
                                 0, 0.001),
                    nrow = nrow(response.mat),
                    ncol = ncol(response.mat))
  response.mat <- response.mat + noise
  return(response.mat)
}

#' Extract single drug response from matrix
#'
#' \code{ExtractSingleDrug} extracts the dose-response values of single drug (
#' drug added in column or row) from a drug combination dose-response matrix.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#'   and row name are representing the concerntrations of drug added to column 
#'   and row, respectively. The values in matrix indicate the inhibition rate to
#'   cell growth.
#' @param dim A character. It should be either "col" or "row" to indicate which
#'   drug's dose-response value will be extracted.
#'
#' @return A data frame. It contains two variables:
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
#' response.mat <- data$dose.response.mats[[1]]
#' drug.row <- ExtractSingleDrug(response.mat, dim = "row")
ExtractSingleDrug <- function(response.mat, dim = "row") {
  dose_col <- as.numeric(colnames(response.mat))
  dose_row <- as.numeric(rownames(response.mat))
  if (dim == "row") {
    single.drug <- data.frame(response = response.mat[, dose_col == 0],
                              dose = dose_row)
  } else if (dim == "col") {
    single.drug <- data.frame(response = response.mat[dose_row == 0, ],
                              dose = dose_col)
  } else {
    stop("Values for 'dim' should be eighther 'row' or 'col'!")
  }
  rownames(single.drug) <- NULL
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
CorrectBaseLine <- function(response.mat, method = c("non", "part", "all")){

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

    baseline <- mean(c(min(as.numeric(drug.row.fit)),
                       min(as.numeric(drug.col.fit))))
    response.mat[negative.ind] <- vapply(response.mat[negative.ind],
                                         function(x) {
                                           x - ((100 - x) / 100 * baseline)
                                         }, numeric(1))
    return(response.mat)
  } else if (method == "all"){
    drug.row <- ExtractSingleDrug(response.mat, dim = "row")
    drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row)))

    drug.col <- ExtractSingleDrug(response.mat, dim = "col")
    drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col)))

    baseline <- mean(c(min(as.numeric(drug.row.fit)),
                       min(as.numeric(drug.col.fit))))
    response.mat <- response.mat - ((100 - response.mat) / 100 * baseline)
    return(response.mat)
  }
}

