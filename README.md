# Time Series and Cross-sectional/Fama-Macbeth Tests

This repository contains two primary functions for conducting time series tests and cross-sectional/Fama-Macbeth tests on financial data. These functions are designed to work with asset returns and factor data to evaluate the performance and risk of portfolios.

## Reference

The methodology and functions in this repository are primarily based on Chapter 12, "Regression-based tests of linear factor models", from John Cochrane's book "Asset Pricing". All equations and formulas implemented and referenced within this repository draw inspiration from this specific chapter. Users and readers are recommended to delve into this chapter for a more in-depth understanding of the underlying principles.

For further context and a comprehensive treatment of the topics discussed here, users are encouraged to consult Cochrane's "Asset Pricing" directly.


## Dependencies
- `dplyr`
- `tidyr`
- `broom`
- `lmtest`
- `sandwich`

## Functions
### 1. `ts_test()`

#### Description:
This function conducts a time series regression test on asset returns against specified factors.

#### Parameters:
- `testasset`: A dataframe containing asset returns with a Date column.
- `factors`: A dataframe containing factor returns with a Date column.
- `startdate`: Optional start date for the analysis.
- `enddate`: Optional end date for the analysis.
- `sd`: Standard deviation method. Options are "standard" and "neweywest". Default is "standard".

#### Returns:
A list containing:
- `loadings`: Factor loadings for each portfolio.
- `residuals`: Residuals from the time series regression.
- `intercepts`: Intercepts from the regression.
- `p_value_GRS`: P-value from the GRS test.
- `risk_premium`: Expected factor returns.
- `r_squared_adj`: Adjusted R-squared values.

#### Mathematical Background:

$$
T\\left[1+\\frac{E_T(f)}{\\sigma(f)}\\right]^{-1} \\alpha' \\Sigma^{-1} \\alpha \\sim \\chi_N^2
$$


### 2. `cross_sec_test()`

#### Description:
This function conducts a cross-sectional or Fama-Macbeth regression test on asset returns against specified factors and/or characteristics.

#### Parameters:
- `testasset`: A dataframe in long format that includes a Date class column, a column indicating the portfolio or asset names, and a column representing returns.
- `factors`: A dataframe containing factor returns with a Date column.
- `characteristic`: A dataframe in long format that includes a Date class column, a column indicating the portfolio or asset names, and columns representing signals/characteristics.
- `model`: Regression model to use. Options are "cross-section" and "fama-macbeth". Default is "cross-section".
- `char.only`: If TRUE, only uses characteristics in the regression. Default is FALSE.
- `startdate`: Optional start date for the analysis.
- `enddate`: Optional end date for the analysis.
- `sd`: Standard deviation method. Options are "standard" and "neweywest". Default is "standard".

#### Returns:
A list containing results from the first and second stages of the regression. The exact contents depend on the chosen model and input data.

P.s.: Please input data frames with variable names without spaces.

## Usage

After importing the required libraries and the functions, you can use them as follows:

```R
# Import data
ff_factors <- read.csv("F-F_Research_Data_Factors.csv")
testasset <- read.csv("48_Industry_Portfolios.csv")

# Conduct time series test
result1 <- ts_test(testasset, ff_factors)

# Conduct cross-sectional test
cross_section_opt <- cross_sec_test(testasset, factors)
```

## Contributing

If you find any bugs or have suggestions for improvements, please open an issue or submit a pull request.

---

You can adapt and expand this README as needed for your GitHub repository.
