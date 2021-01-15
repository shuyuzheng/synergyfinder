repResponse <- function(data){
  sum.data <- data %>% 
    dplyr::group_by(PairIndex, Conc1, Conc2) %>% 
    dplyr::summarise(mean = mean(Response),
                     sd = sd(Response),
                     n = n(), .groups = "keep") %>% 
    mutate(sem = sd/sqrt(n)) %>% 
    mutate(error = qt(0.975,df=n-1)*sem,
           lower_95CI = mean - error,
           upper_95CI = mean + error) 
}

repSynergy <- function(data, method){
  pair <- unique(data$PairIndex)
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
    if (length(unique(repN)) == 1){
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
      }
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