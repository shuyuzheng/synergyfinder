# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, March 2021

#' A high-throughput 3 drug combination screening data
#'
#' A 3-drug combination screening data on Malaria. It is downloaded from 
#' [NCATS Matrix](https://matrix.ncats.nih.gov/) project 2321 "Malaria TACT".
#' It contains 2 blocks, one is synergistic and the other is antagonistic.
#' 
#' @format A data frame with the following columns: "block_id", "drug1",
#'   "drug2", "drug3", "conc1", "conc2", "conc3", "response", "conc_unit1",
#'   "conc_unit2", "conc_unit3"
#' @name NCATS_screening_data
NULL
