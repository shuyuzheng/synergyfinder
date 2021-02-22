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
  colnames <- tolower(
    gsub("([a-z])([A-Z])", "\\1_\\L\\2", colnames, perl = TRUE)
  )
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
#' @details The input data must contain the following columns:
#' (block_id/BlockId/PairIndex), (drug_row/DrugRow/Drug1),
#' (drug_col/DrugCol/Drug2), (response/Response), (conc_r/ConcRow/Conc1),
#' (conc_c/ConcCol/Conc2), and (ConcUnit/conc_r_unit, conc_c_unit/ConcUnit1,
#' ConcUnit2, ConcUnit3)
#'
#' @param data drug combination response data in a data frame format
#' @param impute a logical value. If it is \code{TRUE}, the \code{NA} values
#'   will be imputed by \code{\link[mice]{mice}}. Default is \code{TRUE}.
#' @param impute_method a single string. It sets the \code{method} parameter
#'   in function \code{\link[mice]{mice}} to specify the imputation method. 
#'   Please check the documentation of \code{\link[mice]{mice}} to find the
#'   available methods.
#' @param noise a logical value. It indicates whether or not adding noise to
#'   to the "response" values in the matrix. Default is \code{TRUE}.
#' @param seed a single value, interpreted as an integer, or NULL. It is the
#'   random seed for calculating the noise. Default setting is \code{NULL}
#' @param correction a character. This argument is extended from the argument
#'   \code{method} of \code{\link{CorrectBaseLine}} function. There are three
#'   available valuse: \code{non}, \code{part}, \code{all}.
#'   The default setting is \code{non}.
#' @param data_type a parameter to specify the response data type which can be
#'  either "viability" or "inhibition".
#'
#' @return a list of the following components:
#'   \itemize{
#'     \item \strong{drug_pairs} A data frame contains the name of all the
#'     tested drugs, concentration unit, block IDs and a logical column 
#'     "replicate" to indicate whether there are replicates in the corresponding
#'     block.
#'     \item \strong{response} A data frame contains the columns: "concX" 
#'     concentrations for drugs from input data; "response_origin" response
#'     values from input data; "response" % inhibition value for downstream
#'     analysis.
#'   \item \strong{response_statistics} A data frame. It will be output if
#'     there is block have replicated response values. It contains the
#'     block ID, the concentrations for all the tested drugs, and statistics for
#'     % inhibition values across replicates (including mean, standard
#'     deviation, standard error of mean and 95% confidence interval).
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
ReshapeData <- function(data, impute = TRUE, impute_method = NULL,
                        noise = TRUE, seed = NULL, data_type = "viability") {
  data <- AdjustColumnName(data)
  data <- dplyr::as_tibble(data)
  # 1 check column names

  if (!all(c(
    "block_id", "drug1", "drug2", "response", "conc1", "conc2",
    "conc_unit1"
  ) %in%
    colnames(data))) {
    stop(
      "The input data must contain the following columns:",
      "(block_id/BlockId/PairIndex), ",
      "(drug_row/DrugRow/Drug1/drug1), ",
      "(drug_col/DrugCol/Drug2/drug2), ",
      "(response/Response),",
      "(conc_r/ConcRow/Conc1/conc1), ",
      "(conc_c/ConcCol/Conc2/conc2), ",
      "(ConcUnit/conc_r_unit/ConcUnit1/conc_unit1)"
    )
  }
  
  # Complete the conc_unit columns
  drugs <- grep("drug\\d", colnames(data), value = TRUE)
  conc_unit <- sub("drug", "conc_unit", drugs)
  conc_unit <- setdiff(conc_unit, grep("conc_unit\\d", colnames(data), value = TRUE))
  if (length(conc_unit) > 0) {
    for (i in conc_unit) {
      data[[i]] <- data$conc_unit1
    }
  }
  
  if (!impute & sum(is.na(data$response))) {
    stop(
      "There are missing values in input data. Please fill it up or run ", 
      "'ReshapeData' with 'impute=TRUE'."
    )
  }
  
  drugs <- grep("drug\\d+", colnames(data), value = TRUE)
  conc_units <- gsub("conc_unit", "drug", drugs)

  for (i in 1:length(conc_units)) {
    if (!conc_units[i] %in% colnames(data)) {
      data[conc_units[i]] <- data$conc_unit1
    }
  }
  
  # 2. Split data
  drug_pairs <- data %>% 
    dplyr::select(block_id, dplyr::starts_with(c("drug", "conc_unit"))) %>% 
    dplyr::arrange(block_id) %>% 
    unique()
  
  response <- data %>% 
    dplyr::select(
      block_id,
      dplyr::matches("conc\\d", perl = TRUE),
      response
    ) %>% 
    dplyr::arrange(block_id) %>% 
    dplyr::mutate(response_origin = response) %>% 
    unique()
  
  # 3. make sure the response values are % inhibition
  if (data_type == "viability") {
    response$response <- 100 - response$response
  } else if (data_type == "inhibition") {
    response$response <- response$response
  } else {
    stop(
      "Please tell me the data type of response valuse: 'viability' or ",
      "'inhibition'."
    )
  }
  drug_pairs$input_type <- data_type
  
  # 4. Add random noise
  if (noise) {
    set.seed(seed)
    response$response <- response$response + 
      stats::rnorm(nrow(response), 0, 0.001)
  }
  
  # 5.Impute missing values
  # Whether all blocks are full matrices
  concs <- grep("conc\\d+", colnames(response), value = TRUE)
  combs <- response %>% 
    dplyr::select(-response, -response_origin) %>% 
    dplyr::group_by(block_id) %>% 
    tidyr::nest(data = all_of(concs)) %>% 
    dplyr::mutate(
      data = purrr::map(
        data, 
        function(d) {
          expand.grid(lapply(d, function(x) unique(x)))
        }
      )
    ) %>% 
    tidyr::unnest(cols = c(data)) %>% 
    dplyr::left_join(response, by = c("block_id", concs))
  
  non_complete_block <- combs$block_id[is.na(combs$response)]
  if (length(non_complete_block) > 0) {
    if (!imput){
      stop(
        "The blocks: ", paste(non_complete_block, sep = ", "),
        " are not full combination matrices.",
        "Please complete them manually or run 'ReshapeData' with 'impute=TRUE'."
      )
    } else {
      imp <- mice(combs, method = impute_method, printFlag = FALSE)
      combs <- complete(imp) %>% 
        dplyr::rename(response_adj = response) %>% 
        dplyr::left_join(response, by = c("block_id", concs))
    }
  }
  
  # 6. Dealing with replicates
  replicate_response <- response %>%
    dplyr::group_by(dplyr::across(c(-response, -response_origin)))%>%
    dplyr::summarise(
      response_sd = sd(response),
      response_mean = mean(response),
      response_origin_sd = sd(response_origin),
      response_origin_mean = mean(response_origin),
      n = dplyr::n(), .groups = "keep"
    ) %>%
    dplyr::filter(n > 1) %>% 
    dplyr::mutate(
      response_sem = response_sd / sqrt(n),
      response_CI95 = qt(0.975, df = n - 1) * response_sem,
      response_origin_sem = response_origin_sd / sqrt(n),
      response_origin_CI95 = qt(0.975, df = n - 1) * response_origin_sem,
    )
  dup_blocks <- replicate_response$block_id
  drug_pairs$replicate <- drug_pairs$block_id %in% dup_blocks

  # # 7. Extract single drug
  # single_drug_data <- response %>% 
  #   dplyr::group_by(block_id) %>% 
  #   dplyr::group_map(~ ExtractSingleDrug(.x))
  # drug_pairs$single_drug_data <- single_drug_data
  # 
  # # 8. Fit single drug dose response curve
  # single_drug_model <- lapply(
  #   single_drug_data, 
  #   function(x) lapply(x, FitDoseResponse)
  # )
  # drug_pairs$single_drug_model <- single_drug_model
  
  # 9. assemble output data
  data <- list(drug_pairs = drug_pairs, response = response)
  if (any(drug_pairs$replicate)){
    data$response_statistics <- replicate_response
  }
  return(data)
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
  noise <- matrix(
    stats::rnorm(
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
#' @return A list contains several data frames each of which contains 2 columns:
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
#' response <- data$response[data$response$block_id == 1,
#'                           c("conc1", "conc2", "response")]
#' single <- ExtractSingleDrug(response.df)
ExtractSingleDrug <- function(response) {
  concs <- grep("conc\\d", colnames(response), value = TRUE)
  single_drug <- vector("list", length(concs))
  names(single_drug) <- concs
  conc_sum <- rowSums(response[, concs])
  for (conc in concs) {
    index <- which(response[, conc] == conc_sum)
    single_drug[[conc]] <- data.frame(
      dose = unlist(response[index, conc]),
      response = response[index, "response"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  return(single_drug)
}
