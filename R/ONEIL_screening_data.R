# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, March 2021

#' A high-throughput 2 drugs combination screening data with replication
#'
#' A 2-drug combination screening data on cancer cell lines. It is from
#' publication [(O'Neil, 2016)](https://pubmed.ncbi.nlm.nih.gov/26983881/).
#' It contains 2 blocks, one is synergistic and the other is antagonistic.
#' 
#' @format A data frame with the following columns: "block_id", "drug1",
#'   "drug2", "cell_line_name", "conc1", "conc2", "response", "conc_unit1",
#'   "conc_unit2"
#' @references O’Neil, J., Benita, Y., Feldman, I., Chenard, M., Roberts, B.,
#'   Liu, Y., Li, J., Kral, A., Lejnine, S., Loboda, A., Arthur, W.,
#'   Cristescu, R., Haines, B.B., Winter, C., Zhang, T., Bloecher, A.,
#'   Shumway, S.D., 2016. An Unbiased Oncology Compound Screen to Identify Novel
#'    Combination Strategies. Mol Cancer Ther 15, 1155–1162.
#' @name ONEIL_screening_data
NULL
