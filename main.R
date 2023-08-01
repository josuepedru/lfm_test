library(dplyr)
library(tidyr)
library(broom)

### time series test ####

timeseries_test <- function(testasset, factors){
  # identify variables with class Date
  date_cols <- names(testasset)[sapply(testasset, function(x) class(x) == "Date")]
  # rename those to "Date" in order to manipulate 
  testasset <- rename(testasset, Date = all_of(date_cols))
  
  # doing the same for the factors data frame
  date_cols <- names(factors)[sapply(factors, function(x) class(x) == "Date")]
  factors_renamed <- rename(factors, Date = all_of(date_cols))
  
  # pivot long the testasset data frame in order to join with the factors data frame
  data_regression <- testasset %>% 
    pivot_longer(cols = -Date, 
                 names_to = "Portfolio", 
                 values_to = "Return") %>% 
    full_join(factors_renamed, by= "Date")
  
  # employ regression in order to extract coefficients 
  coefficients <- data_regression %>%
    group_by(Portfolio) %>%
    do(broom::tidy(
      lm(Return ~ `Excess.Market.Return`+SMB, data = .),
    )) %>%
    select(-statistic)
  
  # loadings 
  loadings <- coefficients %>%
    filter(term != "(Intercept)")
  
  # intercept coefficients 
  intercepts <- coefficients %>%
    filter(term == "(Intercept)")
  
  # residuals
  residuals <- data_regression %>%
    group_by(Portfolio) %>%
    do(modelo= (lm(Return ~ `Excess.Market.Return`+SMB, data = .))) %>%
    rowwise() %>%
    mutate(coeficientes = list(coef(modelo)),
           residuos = list(residuals(modelo))) %>%
    select(Portfolio, residuos) %>%
    tidyr::unnest(residuos) %>%
    group_by(Portfolio) %>%
    mutate(id = row_number())%>%
    pivot_wider(names_from = Portfolio, values_from = residuos, id_cols = id) %>% 
    select(-id)
  
  ## GRS test
  # Sigma hat object
  Sigma_hat <- cov(residuals, use = "pairwise.complete.obs")
  
  # Calculate the number of assets, the number of factors, and the number of observations
  N <- ncol(testasset)
  K <- ncol(factors_renamed)
  T <- nrow(testasset)
  
  # Calculate the expected factor returns (also the risk premia)
  E_f <- factors_renamed %>%
    select(where(is.numeric)) %>%  
    colMeans(na.rm = TRUE)
  
  ## Calculate the Omega matrix
  # let the factors data frame in the matrix form
  factors_matrix <- factors_renamed %>% 
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
              intercepts= intercepts, p_value=p_value, risk_premium=E_f))
}


# import some factors 
ff_factors <- read.csv("F-F_Research_Data_Factors.csv") %>% 
  select(Date=X,`Excess Market Return`= Mkt.RF, SMB, HML) %>% 
  mutate(Date=as.Date(paste0(Date,"01"), format="%Y%m%d"))


# Applied the function 
testasset <- data.frame(date = ff_factors$Date,
                        Portfolio1=rnorm(nrow(ff_factors),mean = 0.05, sd = 0.2),
                        Portfolio2=rnorm(nrow(ff_factors),mean = 0.10, sd = 0.8),
                        Portfolio3=rnorm(nrow(ff_factors),mean = 0.15, sd = 0.10),
                        Portfolio4=rnorm(nrow(ff_factors),mean = 0.40, sd = 0.05))

factors <- data.frame(DATE = ff_factors$Date, 
                      `Excess Market Return`=ff_factors$`Excess Market Return`, 
                      SMB = ff_factors$SMB)

result <- timeseries_test(testasset, factors)









