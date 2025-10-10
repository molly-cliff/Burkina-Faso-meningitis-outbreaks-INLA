

# Repository Contents
1. [Temporal Analysis](https://github.com/molly-cliff/Burkina-Faso-meningitis-outbreaks-INLA/blob/main/1.Temporal-analysis.R): Uses Bayesian spatiotemporal modeling with INLA (Integrated Nested Laplace Approximation) to investigate temporal trends in outbreak occurrence in Burkina Faso. The models explore how the probability of an outbreak changes over time (monthly from 2003-2022), using various temporal structures that capture smooth, correlated, or cyclical effects.

2. [Spatial Analysis](https://github.com/molly-cliff/Burkina-Faso-meningitis-outbreaks-INLA/blob/main/Spatial-analysis.R): Compares multiple spatial model structures(BYM, BYM2, Besag, and Besag Proper)to assess which best explains the occurrence of outbreaks across space and time. Time component used (AR1 cont-month) is from the best performing model in the previous Temporal Analysis. 

3. [Univariate analysis](https://github.com/molly-cliff/Burkina-Faso-meningitis-outbreaks-INLA/blob/main/3.%20Univariate%20analysis)): Univariate analysis of climatic factors, population density and vaccination status against outbreak occurance in Burkina Faso, with examination of variable collinearity. Space-time set up (AR1 cont-month and besag-proper) is from the best performing model in the previous Temporal/Spatial Analysis. 

4. [Multivariate Analysis](https://github.com/molly-cliff/Burkina-Faso-meningitis-outbreaks-INLA/blob/main/4.%20Multivariate%20analysis.R): Multivariate analysis of climatic factors against outbreak occurance, additionally testing the usage of non-linear efects of variables as well as logit v cloglog link as well as priors set up
  
5. [Cross validation] (https://github.com/molly-cliff/Burkina-Faso-meningitis-outbreaks-INLA/blob/main/5.%20Cross%20validation): Leave one out cross validation of best performing model against districts in Burkina Faso with >10 outbreak months over 20 year period.

