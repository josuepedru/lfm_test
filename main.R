library(dplyr)
library(tidyr)
library(broom)
library(lmtest)
library(sandwich)

### time series test ####
ts_test <- function(testasset, factors, startdate = NULL, enddate = NULL, sd = c("standard", "neweywest")){
  
  sd <- match.arg(sd) # check if sd specified is in the list, if not provided "standard" is the default
  
  # Rename data frame columns 
  testasset <- testasset %>%
    rename_with(~ "Return", .cols = where(is.numeric)) %>% 
    rename_with(~ "Instrument", .cols = where(is.character)) %>% 
    rename_with(~ "Date", .cols = where(is.Date))
  if(!is.null(factors)) factors <- factors %>%rename_with(~ "Date", .cols = where(is.Date))
  
  # Filter the data based on startdate and enddate if provided
  if (!is.null(startdate) && !is.null(enddate)) {
    testasset <- testasset %>%
      filter(Date >= startdate, Date <= enddate)
    factors <- factors %>%
      filter(Date >= startdate, Date <= enddate)
  }
  
  # join testasset with the factors data frame
  data_regression <- testasset %>% 
    full_join(factors, by= "Date")
  
  # Generate the right-hand side of the formula string
  rhs_formula <- factors %>%
    ungroup() %>%
    dplyr::select(-Date) %>%
    summarise(rhs_formula = paste(sapply(names(.), function(x) if(grepl(" ", x)) paste0("`", x, "`") else x), collapse = " + ")) %>%
    pull(rhs_formula)
  
  # employ regression in order to extract models for each portfolio
  models <- data_regression %>%
    na.omit() %>% 
    group_by(Instrument) %>%
    do(model = lm(formula = paste("Return ~", rhs_formula), data = .))
  
  # employ regression in order to extract coefficients 
  coefficients <- models %>%
    rowwise() %>%
    mutate(tidy_model = list(
      if (sd == "neweywest") {
        coeftest_result <- coeftest(model, vcov. = NeweyWest(model))
        tibble(term = rownames(coeftest_result),
               estimate = coeftest_result[, "Estimate"],
               std.error = coeftest_result[, "Std. Error"],
               statistic = coeftest_result[, "t value"],
               p.value = coeftest_result[, "Pr(>|t|)"])
      } else if (sd == "standard") {
        broom::tidy(model)
      }
    )) %>%
    select(-model) %>%
    unnest(tidy_model) %>%
    select(-statistic)
  
  # loadings 
  loadings <- coefficients %>%
    filter(term != "(Intercept)")
  
  # intercept coefficients 
  intercepts <- coefficients %>%
    filter(term == "(Intercept)")
  
  # get R-squared and Adjusted R-squared from each model
  models <- models %>%
    rowwise(Instrument) %>%
    mutate(r_squared = summary(model)$r.squared,
           adj_r_squared = summary(model)$adj.r.squared) %>%
    select(Instrument, r_squared, adj_r_squared)
  
  # residuals
  residuals <- data_regression %>%
    na.omit() %>% 
    group_by(Instrument) %>%
    do(modelo= (lm(formula = paste("Return ~", rhs_formula), data = .))) %>%
    rowwise() %>%
    mutate(coeficientes = list(coef(modelo)),
           residuos = list(residuals(modelo))) %>%
    select(Instrument, residuos) %>%
    tidyr::unnest(residuos) %>%
    group_by(Instrument) %>%
    mutate(id = row_number())%>%
    pivot_wider(names_from = Instrument, values_from = residuos, id_cols = id) %>% 
    select(-id)
  
  ## GRS test
  # Sigma hat object
  Sigma_hat <- cov(residuals, use = "pairwise.complete.obs")
  
  # Calculate the number of assets, the number of factors, and the number of observations
  N <- ncol(testasset)
  K <- ncol(factors)
  T <- nrow(testasset)
  
  # Calculate the expected factor returns (also the risk premia)
  E_f <- factors %>%
    select(where(is.numeric)) %>%  
    colMeans(na.rm = TRUE)
  
  ## Calculate the Omega matrix
  # let the factors data frame in the matrix form
  factors_matrix <- factors %>% 
    select_if(is.numeric)
  factors_matrix <- as.matrix(factors_matrix)
  
  # compute the omega matrix
  Omega <- (1/T) * t(factors_matrix - E_f) %*% (factors_matrix - E_f)
  
  ## Calculate the test statistic and p-value
  # let the intercept data frame in the matrix form
  intercepts_matrix <- intercepts %>% 
    ungroup() %>%
    select(estimate)
  intercepts_matrix <- as.matrix(intercepts_matrix)
  
  # compute p-value
  test_statistic <- ((T - N - K)/N) * (1 + t(E_f) %*% solve(Omega) %*% E_f)^(-1) *t(intercepts_matrix) %*% Sigma_hat%*% (intercepts_matrix)
  p_value <- pf(test_statistic, N, T - N - K, lower.tail = FALSE)
  
  # output
  return(list(loadings = loadings, residuals = residuals, 
              intercepts = intercepts, p_value = p_value, 
              risk_premium = E_f, r_squared_adj = modeltests))
}


### Cross-sectional/Fama-Mcbeth test ####
cross_sec_test <- function(testasset, factors = NULL,
                           model = c("cross-section", "fama-macbeth"),
                           characteristic = NULL, char.only = FALSE,
                           startdate = NULL, enddate = NULL, sd = c("standard", "neweywest")) {
  
  # Rename data frame columns 
  testasset <- testasset %>%
    rename_with(~ "Return", .cols = where(is.numeric)) %>% 
    rename_with(~ "Instrument", .cols = where(is.character)) %>% 
    rename_with(~ "Date", .cols = where(is.Date))
  if(!is.null(characteristic)){ characteristic <- characteristic %>%
    rename_with(~ "Instrument", .cols = where(is.character)) %>% 
    rename_with(~ "Date", .cols = where(is.Date))}
  if(!is.null(factors)) factors <- factors %>%rename_with(~ "Date", .cols = where(is.Date))
  
  
  # Filter based on startdate and enddate
  if (!is.null(startdate) && !is.null(enddate)) {
    testasset <- testasset %>% filter(Date >= startdate, Date <= enddate)
    if (!is.null(factors)) factors <- factors %>% filter(Date >= startdate, Date <= enddate)
  }
  
  # First stage: time-series regression
  ts_test_results <- list()
  if (char.only == FALSE) {
    data_regression <- testasset %>% 
      left_join(factors, by = "Date")
    
    ts_test_results <- ts_test(testasset, factors, startdate, enddate, sd)
  }
  
  if (!is.null(ts_test_results$loadings)) {
    loadings_long <- ts_test_results$loadings %>%
      select(Instrument, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate)
    
    testasset <- testasset %>% 
      left_join(loadings_long, by = "Instrument")
  }
  
  if (char.only || !is.null(characteristic)) {
    famamcbeth_data <- testasset %>% 
      left_join(characteristic, by = c("Date", "Instrument"))
  }
  
  rhs_formula <- famamcbeth_data %>%
    ungroup() %>%
    dplyr::select(-Date, -Return, -Instrument) %>%
    {paste(names(.), collapse = " + ")}
  
  ### Second stage: cross-section or fama-macbeth
  ## Cross-section
  if (model == "cross-section"){
    famamcbeth_result <- famamcbeth_data %>%
      group_by(Date) %>%
      nest_by() %>%
      mutate(model = list(lm(formula = paste("Return ~", rhs_formula), data = data)),
             tidied = list(broom::tidy(model))) %>%
      select(-data, -model) %>%
      unnest(tidied)
    
    average_excess <- testasset %>% 
      pivot_longer(cols = -Date, 
                   names_to = "Instrument", 
                   values_to = "Return") %>% 
      group_by(Instrument) %>% 
      summarise(`Average Excess Return` = mean(Return))
    
    cross_data <- inner_join(average_excess, loadings_long, by = "Instrument") %>% 
      lm(formula = paste("`Average Excess Return` ~", rhs_formula), data = .)
    
    output <- list()
    # first stage results
    output$first_stage <- list(
      residuals = ts_test_results$residuals,
      loadings = ts_test_results$loadings,
      intercepts = ts_test_results$intercepts,
      p_value_GRS = ts_test_results$p_value,
      r_squared = ts_test_results$r_squared_adj)
    
    # second stage results  
    riskpremia <- summary(cross_data)
    output$second_stage <- list(
      riskpremia = riskpremia)  
    
  }
  
  ## Fama-Mcbeth
  if (model == "fama-macbeth"){
    
    famamcbeth_result <-  famamcbeth_data %>% 
      drop_na() %>% 
      nest(data=c(Instrument,names(famamcbeth_data%>%select(where(is.numeric))))) %>% 
      mutate(etimates = map(data, ~tidy(lm(formula = paste("Return ~", rhs_formula), 
                                           data = .x)))) %>% 
      select(-data) %>% 
      unnest()
    
    # Prepare output
    output <- list()
    
    # first stage results 
    if (!is.null(ts_test_results$loadings)) {
      output$first_stage <- list(
        residuals = ts_test_results$residuals,
        loadings = ts_test_results$loadings,
        intercepts = ts_test_results$intercepts,
        p_value_GRS = ts_test_results$p_value,
        r_squared = ts_test_results$r_squared_adj
      )
    }
    
    # second stage results
    riskpremia <- famamcbeth_result %>%
      filter(term != "(Intercept)") %>%
      select(Date, term, estimate)
    
    riskpremia_mean <- riskpremia %>%
      group_by(term) %>%
      summarise(mean_beta = mean(estimate))
    
    riskpremia_variance <- riskpremia %>%
      group_by(term) %>%
      summarise(variance_beta = sum((estimate-mean(estimate))^2) / (n()^2) )
    
    riskpremia_t_stat <- riskpremia_mean %>%
      left_join(riskpremia_variance, by = "term") %>%
      mutate(T = n(),
             t_stat = mean_beta / sqrt(variance_beta) ) %>%
      select(term, t_stat) 
    
    output$second_stage <- list(
      riskpremia = riskpremia,
      riskpremia_mean = riskpremia_mean,
      riskpremia_variance = riskpremia_variance,
      t_statistics = riskpremia_t_stat
    )}
  
  return(output)
}


