library(dplyr)
library(tidyr)
library(broom)
library(lmtest)
library(sandwich)

### time series test ####

ts_test <- function(testasset, factors, startdate = NULL, enddate = NULL, sd = c("standard", "neweywest")){
  
  sd <- match.arg(sd) # check if sd specified is in the list, if not provided "standard" is the default
  
  # identify variables with class Date
  date_cols <- names(testasset)[sapply(testasset, function(x) class(x) == "Date")]
  
  # rename those to "Date" in order to manipulate 
  testasset <- rename(testasset, Date = all_of(date_cols))
  
  # doing the same for the factors data frame
  date_cols <- names(factors)[sapply(factors, function(x) class(x) == "Date")]
  factors <- rename(factors, Date = all_of(date_cols))
  
  # Filter the data based on startdate and enddate if provided
  if (!is.null(startdate) && !is.null(enddate)) {
    testasset <- testasset %>%
      filter(Date >= startdate, Date <= enddate)
    factors <- factors %>%
      filter(Date >= startdate, Date <= enddate)
  }
  
  # pivot long the testasset data frame in order to join with the factors data frame
  data_regression <- testasset %>% 
    pivot_longer(cols = -Date, 
                 names_to = "Portfolio", 
                 values_to = "Return") %>% 
    full_join(factors, by= "Date")
  
  # Generate the right-hand side of the formula string
  rhs_formula <- factors %>%
    dplyr::select(-Date) %>%
    summarise(rhs_formula = paste(names(.), collapse = " + ")) %>%
    pull(rhs_formula)
  
  # employ regression in order to extract models for each portfolio
  models <- data_regression %>%
    group_by(Portfolio) %>%
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
    rowwise(Portfolio) %>%
    mutate(r_squared = summary(model)$r.squared,
           adj_r_squared = summary(model)$adj.r.squared) %>%
    select(Portfolio, r_squared, adj_r_squared)
  
  # residuals
  residuals <- data_regression %>%
    group_by(Portfolio) %>%
    do(modelo= (lm(formula = paste("Return ~", rhs_formula), data = .))) %>%
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
              risk_premium = E_f, r_squared_adj = models))
}


# import some factors 
ff_factors <- read.csv("F-F_Research_Data_Factors.csv") %>% 
  select(Date=X,`Excess Market Return`= Mkt.RF, SMB, HML) %>% 
  mutate(Date=as.Date(paste0(Date,"01"), format="%Y%m%d"))

# generate random portfolios return 
testasset <- data.frame(date = ff_factors$Date,
                        Portfolio1=rnorm(nrow(ff_factors),mean = 0.05, sd = 0.2),
                        Portfolio2=rnorm(nrow(ff_factors),mean = 0.10, sd = 0.8),
                        Portfolio3=rnorm(nrow(ff_factors),mean = 0.15, sd = 0.10),
                        Portfolio4=rnorm(nrow(ff_factors),mean = 0.40, sd = 0.05))

# create factors data frame 
factors <- data.frame(DATE = ff_factors$Date, 
                      `Excess Market Return`=ff_factors$`Excess Market Return`, 
                      SMB = ff_factors$SMB,
                      HML = ff_factors$HML)

# apply the function for different windows 
result1 <- ts_test(testasset, factors, 
                  startdate = "2000-01-01", 
                  enddate = "2008-01-01")
result2 <- ts_test(testasset, factors, 
                   startdate = "2008-01-01", 
                   enddate = "2022-01-01")


### Cross-sectional/Fama-Mcbeth test ####
cross_sec_test <- function(testasset, factors = NULL,
                           model = c("cross-section", "fama-macbeth"),
                           characteristic = NULL, char.only = TRUE,
                           startdate = NULL, enddate = NULL, sd = c("standard", "neweywest")) {
  
  # Helper function to rename Date columns
  rename_date <- function(df) {
    date_cols <- names(df)[sapply(df, function(x) class(x) == "Date")]
    rename(df, Date = all_of(date_cols))
  }
  
  # Rename testasset date column
  testasset <- rename_date(testasset)
  
  # Rename characteristic and factors date column
  if (!is.null(factors)) factors <- rename_date(factors)
  if (!is.null(characteristic)) characteristic <- rename_date(characteristic)
  
  # Filter based on startdate and enddate
  if (!is.null(startdate) && !is.null(enddate)) {
    testasset <- testasset %>% filter(Date >= startdate, Date <= enddate)
    if (!is.null(factors)) factors <- factors %>% filter(Date >= startdate, Date <= enddate)
  }
  
  # First stage: time-series regression
  ts_test_results <- list()
  if (char.only == FALSE) {
    data_regression <- testasset %>% 
      pivot_longer(cols = -Date, names_to = "Portfolio", values_to = "Return") %>%
      left_join(factors, by = "Date")
    
    ts_test_results <- ts_test(testasset, factors, startdate, enddate, sd)
  }
  
  # Prepare data for second stage
  famamcbeth_data <- testasset %>%
    pivot_longer(cols = -Date, names_to = "Portfolio", values_to = "Return")
  
  if (!is.null(ts_test_results$loadings)) {
    loadings_long <- ts_test_results$loadings %>%
      select(Portfolio, term, estimate) %>%
      pivot_wider(names_from = term, values_from = estimate)
    
    famamcbeth_data <- famamcbeth_data %>% left_join(loadings_long, by = "Portfolio")
  }
  
  if (char.only || !is.null(characteristic)) {
    famamcbeth_data <- famamcbeth_data %>% left_join(characteristic, by = c("Date", "Portfolio"))
  }
  
  rhs_formula <- famamcbeth_data %>%
    select(-Date, -Return, -Portfolio) %>%
    summarise(rhs_formula = paste(names(.), collapse = " + ")) %>%
    pull(rhs_formula)
  
  ### Second stage: cross-section or fama-macbeth
  ## Cross-section
  if (model == "cross-section") {
    famamcbeth_result <- famamcbeth_data %>%
      group_by(Date) %>%
      nest_by() %>%
      mutate(model = list(lm(formula = paste("Return ~", rhs_formula), data = data)),
             tidied = list(broom::tidy(model))) %>%
      select(-data, -model) %>%
      unnest(tidied)
    
    average_excess <- testasset %>% 
      pivot_longer(cols = -Date, 
                   names_to = "Portfolio", 
                   values_to = "Return") %>% 
      group_by(Portfolio) %>% 
      summarise(`Average Excess Return` = mean(Return))
    
    cross_data <- inner_join(average_excess, loadings_long, by = "Portfolio") %>% 
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
    famamcbeth_result <- famamcbeth_data %>%
      group_by(Date) %>%
      nest_by() %>%
      mutate(model = list(lm(formula = paste("Return ~", rhs_formula), data = data)),
             tidied = list(broom::tidy(model))) %>%
      select(-data, -model) %>%
      unnest(tidied)
    
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
      summarise(mean_beta = mean(estimate, na.rm = TRUE))
    
    riskpremia_variance <- riskpremia %>%
      group_by(term) %>%
      summarise(variance_beta = var(estimate, na.rm = TRUE))
    
    output$second_stage <- list(
      riskpremia = riskpremia,
      riskpremia_mean = riskpremia_mean,
      riskpremia_variance = riskpremia_variance)}
  
  return(output)
}


# import some factors 
ff_factors <- read.csv("F-F_Research_Data_Factors.csv") %>% 
  select(Date=X,`Excess Market Return`= Mkt.RF, SMB, HML) %>% 
  mutate(Date=as.Date(paste0(Date,"01"), format="%Y%m%d"))

# generate random portfolios return 
# industry portfolio plus 40 portftolios anomalys
testasset <- read.csv("48_Industry_Portfolios.csv",
                      header = TRUE,skip=11)[c(1:1161),] %>% 
  mutate(Date = as.Date(paste0(X,"01"), format = "%Y%m%d")) %>% 
  mutate(across(-c(Date), as.numeric)) %>% select(-X) 

# create factors data frame 
factors <- data.frame(DATE = ff_factors$Date, 
                      `Excess Market Return`=ff_factors$`Excess Market Return`, 
                      SMB = ff_factors$SMB,
                      HML = ff_factors$HML)

# create characteristic data frame 
characteristic <- testasset %>%
  pivot_longer(cols = -Date, 
               names_to = "Portfolio", 
               values_to = "Return") %>% 
  select(-Return) %>% 
  mutate(Size = rnorm(n=55728, mean = 1.6, sd = 3.2), 
         B2M = rnorm(n=55728, mean = 2, sd = 4),
         Revenue = rnorm(n=55728, mean = 4, sd = 8))


# cross-sectional test
cross_section_opt <- cross_sec_test(testasset,
                                        factors,
                                        characteristic=NULL,
                                        model="cross-section",
                                        char.only = FALSE,
                                        sd = "standard")

## fama-mcbeth
# test both, factors and characteristics
both <- cross_sec_test(testasset,
                           factors,
                           characteristic,
                           model="fama-macbeth",
                           char.only = FALSE,
                           sd = "standard")
# test only characteristics
charac_only <- cross_sec_test(testasset,
                                  factors=NULL,
                                  characteristic,
                                  model="fama-macbeth",
                                  char.only = TRUE,
                                  sd = "standard")
# test only factors
factors_only <- cross_sec_test(testasset,
                         factors=factors,
                         characteristic=NULL,
                         model="fama-macbeth",
                         char.only = FALSE,
                         sd = "standard")
