# Value-at-Risk-Covid-19-Oil-Companies-MATLAB
# We used MATLAB to perform this analysis.
# This work aims to explain if semi-parametric model like Extreme Value Theory could be able to explain stress period like Covid-19 against traditional methods (Historical Simulation and Normal Distribution).
# The project is splitted in 3 sections:
# 1. --- Data Modelling ---
# In the first file we analyze the returns of four Oil Companies during the Covid-19 Crisis. (Form 1-Jan-2018 to 30-Aug-2022)
# The data has been downloaded from Datastream
# In order to work with standardized data we compute a GARCH(1,1) model.
# We conclude that any relevant AR process is useful and from ACF of Squared Residuals we confirm that we are working on I.I.D. Data
# 2. --- VaR estimation with traditional methods : Normal Distribution and Historical Simulation and the relative backtest procedures --- 
# We adopt a movable windows system starting from 1-Jan-2019 to 30-Aug-2022 with a WindowSize of 250 observations.
# We compute the analysis with 95% and 99% VaR confinance.
# For the backtesting a 95 % of the conditional distribution has been used.
# 3. --- VaR estimation with EVT: Extreme Value Theory ---
# The literature that has been used to build this model is the following:
#     - Estimation of tail-related risk measures for heteroscedastic financial time series: an extreme value approach (McNeil Frey, 2000)
#     - Using Extreme Value Theory and Copulas to Evaluate Market Risk (MATLAB Toolbox)
# We used a 10% rule in order to choose the threshold for the lower and the upper tail.
# We applied a movable window system in order to compare VaR with traditional methods.

## Results ##
# At 95% confidance level the traditional methods perform better in the VaR estimation looking at the backtesting.
# At 99% confidance level EVT method perform better than the 95% only for %TOT 
# At 99% confidance level all the methods seems to explain something different upon different companies.
# Extending the WindowSize to 500 for EVT doesn't improve it. 
