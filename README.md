Functions, simulation study and empirical data application for Regularized Least Square Latent Variable method (Regularized LS-LV).

## Authors
- Tra T. Le | Tilburg University
- Dr. Katrijn Van Deun (Supervisor)

## Description
The repository consists of accompanying codes for the Master's thesis of the Research Master in Social and Behavioral Science at Tilburg University. We developed the RLSLV method as an alternative for (exploratory) structural equation modeling in the high-dimensional data setting (small n, large p).
Here, you can find the **functions** developed in R to implement the method, the codes for the **Simulation Study** and **Empirical Data Analysis**.

## Functions
Folder [Functions] contains all the functions needed for the two methods: [CCLSLV.R](Functions/CCLSLV.R) for the cardinality constrained approach and [LSLVLASSO.R](Functions/LSLVLASSO.R) for the one with the LASSO penalty. Each file also includes functions to run **Model Selection** based on the Index of Spareness.

## Simulation Study
- Step 1: Run [corrmatrix.R](Simulation Study/corrmatrix.R) to generate the population correlation matrix
- Step 2: Run [datageneration_orthogonal.R](Simulation Study/datageneration_orthogonal.R) and [datageneration_correlated](Simulation Study/datageneration_correlated.R) to generate data for orthogonal and correlated latent variables scenarios, respectively. 
