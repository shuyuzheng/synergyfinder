# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, March 2021
#
# SynergyFinder
#
# Functions in this page:
#
# ReshapeData: pre-process the response data for further calculation and plot.
#
# Auxiliary functions:
# .AdjustColumnName: Adjust column names of input data table
# ExtractSingleDrug: Extract Single Drug Dose Response

#' Pre-process the Response Data for Further Calculation and Plot
#'
#' A function to transform the response data from data frame format to
#' dose-response matrices. Several processes could be chose to add noise, impute
#' missing values or correct base line to the dose-response matrix.
#'
#' @details The input data must contain the following columns:
#' (block_id/BlockId/PairIndex), (drug_row/DrugRow/Drug1),
#' (drug_col/DrugCol/Drug2), (response/Response/inhibition/Inhibition),
#' (conc_r/ConcRow/Conc1), (conc_c/ConcCol/Conc2), and
#' (ConcUnit/conc_r_unit, conc_c_unit/ConcUnit1, ConcUnit2, ConcUnit3)
#'
#' @param data drug combination response data in a data frame format
#' @param data_type a parameter to specify the response data type which can be
#'   either "viability" or "inhibition".
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
#' @param iteration An integer. It indicates the number of iterations for
#'   bootstrapping while calculating statistics for data with replicates.
#'
#' @return a list of the following components:
#'   \itemize{
#'     \item \strong{drug_pairs} A data frame contains the name of all the
#'     tested drugs, concentration unit, block IDs and a logical column 
#'     "replicate" to indicate whether there are replicates in the corresponding
#'     block.
#'     \item \strong{response} A data frame contains the columns: "concX" 
#'     concentrations for drugs from input data; "response_origin" response
#'     values from input data; "response" \% inhibition value for downstream
#'     analysis.
#'   \item \strong{response_statistics} A data frame. It will be output if
#'     there is block have replicated response values. It contains the
#'     block ID, the concentrations for all the tested drugs, and statistics for
#'     \% inhibition values across replicates (including mean, standard
#'     deviation, standard error of mean and 95\% confidence interval).
#'   }
#' 
#' @author
#' \itemize{
#'   \item Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'   \item Jing Tang \email{jing.tang@helsinki.fi}
#' }
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
ReshapeData <- function(data,
                        impute = TRUE,
                        impute_method = NULL,
                        noise = FALSE,
                        seed = NULL,
                        iteration = 10,
                        data_type = "viability") {
  data <- .AdjustColumnName(data)
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
      "(response/Response/inhibition/Inhibition),",
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
    dplyr::mutate(response_origin = response)# %>% 
    # unique()
  
  # 3. make sure the response values are % inhibition
  if (data_type == "viability") {
    response$response <- 100 - response$response
  } else if (data_type == "inhibition") {
    response$response <- response$response
  } else {
    stop(
      "Please specify the data type of response valuse: 'viability' or ",
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
    tidyr::nest(data = dplyr::all_of(concs)) %>% 
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
    if (!impute){
      stop(
        "The blocks: ", paste(non_complete_block, sep = ", "),
        " are not full combination matrices.",
        "Please complete them manually or run 'ReshapeData' with 'impute=TRUE'."
      )
    } else {
      warning(
        sum(is.na(combs$response)),
        " missing values are imputed."
      )
      imp <- suppressWarnings(
        mice::mice(combs, method = impute_method, printFlag = FALSE)
      )
      response <- suppressWarnings(mice::complete(imp)) %>% 
        dplyr::select(-response_origin) %>% 
        dplyr::left_join(
          response %>% 
            dplyr::select(-response),
          by = c("block_id", concs)
        )
    }
  }
  
  # 6. Dealing with replicates
  replicate_response <- response %>%
    dplyr::group_by(dplyr::across(c(-response, -response_origin)))%>%
    dplyr::summarise(
      response_sd = stats::sd(response),
      response_mean = mean(response),
      response_origin_sd = stats::sd(response_origin),
      response_origin_mean = mean(response_origin),
      n = dplyr::n(), .groups = "keep"
    ) %>%
    dplyr::filter(n > 1) %>% 
    dplyr::mutate(
      response_sem = response_sd / sqrt(n),
      response_CI95 = stats::qt(0.975, df = n - 1) * response_sem,
      response_ci_left = response_mean - response_CI95,
      response_ci_right = response_mean + response_CI95,
      response_origin_sem = response_origin_sd / sqrt(n),
      response_origin_CI95 = stats::qt(0.975, df = n - 1) * response_origin_sem,
      response_origin_ci_left = response_origin_mean - response_origin_CI95,
      response_origin_ci_right = response_origin_mean + response_origin_CI95,
    ) %>% 
    dplyr::select(-dplyr::ends_with("_sd"), -dplyr::ends_with("_CI95"), -dplyr::ends_with("_z"))
  dup_blocks <- replicate_response$block_id
  drug_pairs$replicate <- drug_pairs$block_id %in% dup_blocks
  # p value
  if (any(drug_pairs$replicate)){
    drug_pairs$response_p_value <- rep(NA, nrow(drug_pairs))
    drug_pairs$response_origin_p_value <- rep(NA, nrow(drug_pairs))
    blocks <- drug_pairs$block_id[drug_pairs$replicate]
    for (b in blocks) {
      index <- drug_pairs$block_id == b
      # Calculate synergy score. Different work flow for replicated data
      replicate <- drug_pairs$replicate[index]
        response_one_block <- response %>%
          dplyr::filter(block_id == b) %>%
          dplyr::select(-block_id) %>%
          dplyr::ungroup()
        concs <- grep("conc\\d", colnames(response_one_block), value = TRUE)
        iter_response <- pbapply::pblapply(seq(1, iteration), function(x){
          response_boot <- .Bootstrapping(response_one_block)
          s <- response_boot[, c("response", "response_origin")] %>% 
            colMeans()
          s <- as.data.frame(as.list(s))
          return(s)
        }) %>% 
          purrr::reduce(rbind.data.frame)
      p <- apply(iter_response, 2, function(x){
        z <- abs(mean(x)) / stats::sd(x)
        p <- exp(-0.717 * z - 0.416 * z ^2)
        p <- formatC(p, format = "e", digits = 2, zero.print = "< 2e-324")
        return(p)
      })
      names(p) <- paste0(names(p), "_p_value")
      drug_pairs$response_p_value[index] <- p["response_p_value"]
      drug_pairs$response_origin_p_value[index] <- p["response_origin_p_value"]
    }
  }
  
  # 7. assemble output data
  data <- list(
    drug_pairs = dplyr::ungroup(drug_pairs),
    response = dplyr::ungroup(response)
  )
  if (any(drug_pairs$replicate)){
    data$response_statistics <- dplyr::ungroup(replicate_response)
  }
  return(data)
}

# Auxiliary functions -----------------------------------------------------

#' Adjust Column Names of Input Data Table
#'
#' This function changes the column names in other format into the style:
#' block_id, drug1, drug2, conc1, conc2, response, conc_unit1, conc_unit2.
#'
#' @param data A data frame. It is the input data for function 
#'   \link{ReshapeData}
#'
#' @return The data frame with the changed column names.
#' 
#' @author
#' \itemize{
#'   \item Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'   \item Jing Tang \email{jing.tang@helsinki.fi}
#' }
#' 
#' @export
.AdjustColumnName <- function(data) {
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
  
  # if the column name for response is "inhibition
  colnames[which(colnames == "inhibition")] <- "response"
  
  colnames(data) <- colnames
  return(data)
}

#' Extract Single Drug Dose Response
#'
#' \code{ExtractSingleDrug} extracts the dose-response values of single drug
#'  from a drug combination dose-response matrix.
#'
#' @param response A data frame. It must contain the columns: "conc1", "conc2",
#' ..., for the concentration of the combined drugs and "response" for the
#' observed \%inhibition at certain combination.
#'
#' @return A list contains several data frames each of which contains 2 columns:
#'   \itemize{
#'     \item \strong{dose} The concertration of drug.
#'     \item \strong{response} The cell's response (inhibation rate) to
#'       corresponding drug concertration.
#' }
#'
#' @author
#' \itemize{
#'   \item Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'   \item Jing Tang \email{jing.tang@helsinki.fi}
#' }
#' 
#' @export
#'
#' @examples
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
#' response <- data$response[data$response$block_id == 1,
#'                           c("conc1", "conc2", "response")]
#' single <- ExtractSingleDrug(response)
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
