# TidyComb
# Functions for pipe all calculation togather
#
# Functions in this page:
#
# CalculateMat: Calculate scores from one dose-response matrix.
# CalculateTemplate: Calculate scores from template file which might contain
#                    several drug combination blocks.
# CalculateTemplateDebug: Debug version of CalculateTemplate.
# ParCalculateTemplate: Parallel calculate scores from template file. It using
#                    multi-core method to do parallel. (It might not work under
#                    Windows system)

#' Main function for synergy score calculation from one dose-response matrix
#'
#' \code{CalculateMat} is the main function for calculating synergy scores based
#' on model(ZIP, Bliss, Loewe, and HSA) fron one dose-response \strong{matrix}.
#'
#' The steps for calculation:
#' \enumerate{
#'   \item Pre-process Matrix
#'     \enumerate{
#'       \item Impute for missing values (with the average of values from the 
#'           nearest four cells) in original matrix by using function 
#'           \code{\link{ImputeNear}}.
#'       \item Add noise(A small random number ranging from 0 to 0.001) to 
#'           original matrix by using function \code{line{AddNoise}}.
#'       )
#'       \item Correct baseline to 0, if \code{correction} is \code{TRUE}.
#'     }
#'   \item Calculate synergy scores using the functions:
#'     \itemize{
#'       \item \code{\link{CalculateZIP}}
#'       \item \code{\link{CalculateBliss}}
#'       \item \code{\link{CalculateHSA}}
#'       \item \code{\link{CalculateLoewe}}
#'     }
#' }
#' 
#' @param response.mat A matrix contain the drug combination reaponse value.
#' Its column names are doses of drug added along columns. Its row name are
#' doses of drug added along rows.
#' 
#' @param method a parameter to specify which models to use to calculate the 
#'    synergy scores. Choices are "ZIP", "Bliss", "HSA" and "Loewe". Defaults to
#'    "ZIP".
#' 
#' @param noise a logical value. It indicates whether or not adding noise to
#' to the "response" values in the matrix. Default is \code{TRUE}.
#'
#' @param ... Other argumants required by nested functions. Some important
#' arguments are:
#' \itemize{
#'    \item \code{method} inherited from function \code{\link{CorrectBaseLine}};
#'    \item \code{Emin} and \code{Emax} inherited from function
#'          \code{FitDoseResponse}.
#' }
#'
#' @return a list. It contains 4 elements:
#'   \itemize{
#'     \item \strong{method} the method used to calculate the synergy scores.
#'     \item \strong{scores} It contains the modified response value and 4
#'     type of synergy scores of each drug dose response pair.
#'     \item \strong{dose.response.mats} The original input dose-response matrix.
#'     \item \strong{adjusted.response.mats} The dose response matrix adjusted by
#'       functions: \code{\link{AddNoise}}, \code{\link{ImputeNear}}, and 
#'       \code{\link{CorrectBaseLine}}.
#'  }
#'
#' @author Shuyu Zheng{shuyu.zheng@helsinki.fi}
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
CalculateMat <- function(response.mat, noise=TRUE, method="ZIP", 
                         data.type = "viability", correction = "non", ...) {

  options(scipen = 999)
  data <- list()
  
  data$dose.response.mats <- response.mat
  data$method <- method
  # 1. Pre-processing
  # 1.1 Impute data until there is no missing value in the response matrix
  while (sum(is.na(response.mat))) {
    response.mat <- ImputeNear(response.mat)
  }

  # 1.2. Add random noise to original matrix
  if (noise) {
    response.mat <- AddNoise(response.mat, method="random")
  }

  # 1.3. Correct baseline with corresponding "method". Available methods are
  #      "non", "part", "all".
    response.mat <- CorrectBaseLine(response.mat, method = correction, ...)
    
    
  ## 2 Calculate synergy scores
  if(!method %in% c("ZIP", "HSA", "Bliss", "Loewe")) {
    stop("The method parameter can only be one of the following: ZIP, HSA, Bliss
       and Loewe.")
  }
  scores <- switch(method,
                  ZIP = CalculateZIP(response.mat),
                  HSA = CalculateHSA(response.mat),
                  Bliss = CalculateBliss(response.mat),
                  Loewe = CalculateLoewe(response.mat))
  
  ## 3 Save data into the list
  data$adjusted.response.mats <- response.mat
  data$scores <- scores
  return(data)

  # clean up
  gc()
}

#' Calculate Drug Combination data in data frame format (template)
#'
#' @param template a dataframe in the format as template. Columns "block_id",
#' "drug_row", "drug_col", "response", "conc_r", "conc_c", "conc_r_unit",
#' "conc_c_unit","cell_line_name", "drug_row", "drug_col" are reqired.
#'
#' @param debug a logical value. If it is \code{TRUE}, block ID will be printed
#'    to console. The default setting is \code{FALSE}
#'    
#' @param data.type a parameter to specify the response data type which can be 
#' either "viability" or "inhibition".
#' 
#' @param method a parameter to specify which models to use to calculate the 
#'    synergy scores. Choices are "ZIP", "Bliss", "HSA" and "Loewe". Defaults to
#'    "ZIP".
#'    
#' @param ... Other arguments required by nested functions. Some important
#' arguments are:
#'  \itemize{
#'    \item \code{impute} and \code{noise} inherited from function
#'          \code{CalculateMat};
#'    \item \code{method} inherited from function \code{CorrectBaseLine};
#'    \item \code{Emin} and \code{Emax} inherited from function
#'          \code{FitDoseResponse}.
#' }
#' 
#' @return a list. It contains 5 elements:
#'   \itemize{
#'     \item \strong{method} the method used to calculate the synergy scores.
#'     \item \strong{scores} It contains the modified response value and 4
#'       type of synergy scores of each drug dose response pair.
#'     \item \strong{dose.response.mats} The original input dose-response matrix.
#'     \item \strong{adjusted.response.mats} The dose response matrix adjusted 
#'       by functions: \code{\link{AddNoise}}, \code{\link{ImputeNear}}, and 
#'       \code{\link{CorrectBaseLine}}.
#'     \item \strong{drug.pairs} the same as the input data component.
#'   }
#'
#' @author Shuyu Zheng{shuyu.zheng@helsinki.fi}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
CalculateDF <- function(template, debug=FALSE, method="ZIP", 
                              data.type = "viability", ...) {
  # 1. Check the input data
  # 1.1 check column names
  if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r", "conc_c",
             "conc_r_unit", "conc_c_unit") %in%
           colnames(template)))
    stop("The input data must contain the following columns: block_id, drug_row,
         drug_col, response, conc_r, conc_c, conc_r_unit, conc_c_unit")
  

  # 1.2 Check the data type
  if (data.type == "viability") {
    template$response <- 100 - template$response
  } else if (data.type == "inhibition") {
    template <- template
  }
  
  set.seed(1) # set random seed for AddNoise function.
  blocks <- unique(template$block_id)
  data <- list()
  data$method <- method
  data$drug.pairs <- dplyr::select(template, block_id, drug_row, drug_col, 
                                   conc_c_unit, conc_r_unit) %>% unique()
  # generate container
  scores <- vector(mode="list", length=length(blocks))
  names(scores) <- blocks
  dose.response.mats <- vector(mode="list", length=length(blocks))
  names(dose.response.mats) <- blocks
  adjusted.response.mats <- vector(mode="list", length=length(blocks))
  names(adjusted.response.mats) <- blocks

  for (block in blocks) {
    if (debug) {
      message(block)
      utils::flush.console()
    }
    # 1. Generate response matrix for each block
    response <- template %>%
      dplyr::filter(block_id == block)

    response.mat <- response %>%
      dplyr::select(conc_r, conc_c, response) %>%
      reshape2::acast(conc_r ~ conc_c, value.var = "response")

    # 2. Do calculation on matrix (with error control)
    tmp <- CalculateMat(response.mat = response.mat, ...)
    
    scores[[block]] <- tmp$scores
    dose.response.mats[[block]] <- tmp$dose.response.mats
    adjusted.response.mats[[block]] <- tmp$adjusted.response.mats

    # Clean temporary file
    tmp <- list()
    response <- data.frame()
    respons.mat <- matrix()
  }
  
  data$scores <- scores
  data$dose.response.mats <- dose.response.mats
  data$adjusted.response.mats <- adjusted.response.mats
  return(data)
  rm(.Random.seed)
}

# multiResultClass <- function(scores=NULL, dose.response.mat=NULL, 
#                              adjusted.response.mat = NULL) {
#   me <- list(
#     scores = scores,
#     dose.response.mat = dose.response.mat,
#     adjusted.response.mat = adjusted.response.mat
#   )
# 
#   ## Set the name for the class
#   class(me) <- append(class(me), "multiResultClass")
#   return(me)
# }

# ParCalculateTemplate <- function(template, cores = 1, ...) {
#   if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r", "conc_c",
#              "conc_r_unit", "conc_c_unit","cell_line_name", "drug_row",
#              "drug_col") %in%
#            colnames(template)))
#     stop("The input data must contain the following columns: ",
#          "block_id, drug_row, drug_col, response,\n",
#          "conc_r, conc_c, conc_r_unit, conc_c_unit, \n",
#          "cell_line_name.")
# 
#   blocks <- unique(template$block_id)
# 
#   cl <- parallel::makeForkCluster(cores)
#   doParallel::registerDoParallel(cl)
#   
#   res <- foreach::foreach (i = 1:length(blocks)) %dopar% {
#     set.seed(1)
#     result <- multiResultClass()
#     # 1. Generate response matrix for each block
#     response <- template %>%
#       dplyr::filter(block_id == block)
#     
#     response.mat <- response %>%
#       dplyr::select(conc_r, conc_c, response) %>%
#       reshape2::acast(conc_r ~ conc_c, value.var = "response")
#     
#     # 2. Do calculation on matrix (with error control)
#     tmp <- CalculateMat(response.mat = response.mat, ...)
#     
#     result$scores[[block]] <- tmp$scores
#     result$dose.response.mat[[block]] <- tmp$dose.response.mat
#     result$adjusted.response.mat[[block]] <- tmp$adjusted.response.mat
#     # Clean temporary file
#     tmp <- list()
#     response <- data.frame()
#     respons.mat <- matrix()
#     
#   return(data)
#     rm(.Random.seed)
#     result
#   }
# 
#   res2 <- list()
#   res2$synergy <- Reduce(function(x, y) {rbind.data.frame(x, y)},
#                          lapply(res, "[[" , "synergy"))
#   res2$surface <- Reduce(function(x, y) {rbind.data.frame(x, y)},
#                          lapply(res, "[[" , "surface"))
#   res2$summary <- Reduce(function(x, y) {rbind.data.frame(x, y)},
#                          lapply(res, "[[" , "summary"))
#   res2$curve <- Reduce(function(x, y) {rbind.data.frame(x, y)},
#                        lapply(res, "[[" , "curve"))
# 
#   return(res2)
# }
