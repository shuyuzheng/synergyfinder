repResponse <- function(data){
  sum.rep <- data$response.df %>% 
    dplyr::group_by(.dots = colnames(data$response.df)[- grep("response", 
                                              colnames(data$response.df))]) %>% 
    dplyr::summarise(sd = sd(response),
                     response = mean(response),
                     n = dplyr::n(), .groups = "keep")
  rep_block <- unique(sum.rep$block_id[which(sum.rep$n > 1)])
  if (length(rep_block) > 1){
    data[['replicate.response']] <- sum.rep %>% 
      dplyr::mutate(sem = ifelse(n > 1, sd/sqrt(n), NA),
                    "95% confidence interval" = ifelse(n > 1, qt(0.975,df=n-1)*sem, NA))
  }
  data$drug.pairs <- data$drug.pairs %>% 
    dplyr::mutate(replicate = block_id %in% rep_block)
  return(data)
}

repSynergy <- function(response.df, method){
  pair <- unique(data$block_id)
  scores <- NULL
  for (p in pair){
    tmp <- dplyr::filter(data, PairIndex == p)
    repN <- tmp %>% 
      dplyr::group_by(PairIndex, Conc1, Conc2) %>% 
      dplyr::summarise(count = n(), .groups = "keep") %>%
      ungroup() %>% 
      dplyr::select(count) %>% 
      unique() %>% 
      unlist()
    if (length(unique(repN)) == 1){# the replicates in all wells are same
      rest <- tmp %>% 
        mutate(index = seq(1, n()))
      response <- NULL
      for (i in seq(1, (repN - 1))){
        t <- rest %>% 
          dplyr::group_by(PairIndex, Conc1, Conc2) %>% 
          dplyr::sample_n(1) %>% 
          dplyr::mutate(block_id = i)
        response <- rbind.data.frame(response, t)
        rest <- rest %>% 
          dplyr::filter(!index %in% t$index)
      } #else { # the replicates in each well are not same
        
      #}
      response <- rbind.data.frame(response, dplyr::mutate(rest, block_id = repN))
      # Calculate scoree
      blocks <- unique(response$block_id)
      score <- NULL
      for (i in blocks){
        response.mat <- reshape2::acast(Conc1 ~ Conc2, data = response[which(response$block_id == i),], 
                                        value.var = "Response")
        t <- switch(method,
                    ZIP = ZIP(response.mat),
                    HSA = HSA(response.mat),
                    Bliss = Bliss(response.mat),
                    Loewe = Loewe(response.mat))
        t <- reshape2::melt(t)
        colnames(t) <- c("Conc1", "Conc2", "synergy")
        score <- rbind.data.frame(score, t)
      }
      sum.score <- score %>% 
        dplyr::group_by(Conc1, Conc2) %>% 
        dplyr::summarise(mean = mean(synergy),
                         sd = sd(synergy),
                         n = n(), .groups = "keep") %>% 
        mutate(sem = sd/sqrt(n)) %>% 
        mutate(error = qt(0.975,df=n-1)*sem,
               lower_95CI = mean - error,
               upper_95CI = mean + error) %>% 
        mutate(PairIndex = p)
      scores <- rbind.data.frame(scores, sum.score)
    }
  }
  return(scores)
}