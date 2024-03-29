Functions, simulation study and empirical data application for Regularized Least Square Latent Variable method (Regularized LS-LV).

## Authors
- Tra T. Le 
- Prof. dr. Jeroen Vermunt
- Prof. dr. Katrijn Van Deun

## Description
The repository consists of accompanying codes for the Master's thesis of the Research Master in Social and Behavioral Science at Tilburg University. We developed the RLSLV method as an alternative for (exploratory) structural equation modeling in the high-dimensional data setting (small n, large p).
Here, you can find the **functions** developed in R to implement the method, the codes for the **Simulation Study** and **Empirical Data Analysis**.

## Functions
Folder [Functions](Functions) contains all the functions needed for the two methods: 
- [CCLSLV.R](Functions/CCLSLV.R) for the cardinality constrained approach.
- [LSLVLASSO.R](Functions/LSLVLASSO.R) for the one with the LASSO penalty.
- [LSLVL_undoshrinkage.R](Functions/LSLVL_undoshrinkage.R) to undo the shrinkage for nonzero loadings.
The first two files also include functions to run **Model Selection** based on the Index of Spareness.

## Simulation Study
- Step 1: Run [corrmatrix.R](corrmatrix.R) to generate the population correlation matrix
- Step 2: Run [datageneration_orthogonal.R](datageneration_orthogonal.R) and [datageneration_correlated](datageneration_correlated.R) to generate data for orthogonal and correlated latent variables scenarios, respectively. 
- Step 3: Run one of the sim_*method*.R file to analyze the simulated data with one of the methods. For example, run [sim_CCLSLV.R](sim_CCLSLV.R) to run the analysis using the proposed method with cardinality constraint.

Additionally, the subfolder [RegSEM-failed](https://github.com/trale97/RegularizedLSLV/tree/main/Simulation%20Study/RegSEM-failed) contains codes for RegSEM to estimate the relationships among the factors/components, which produced problematic results. 

## Empirical Data Application
The folder [Application](Application) contains two files: [bigfive.R](bigfive.R) to analyze the Big Five Personality Data using all methods and [empirical_genetic.R](empirical_genetic.R) to analyze the Gene Expression data using the two proposed methods. 
- The Big Five data set can be obtained from the R package *qgraph* (Epskamp et al., 2012).
- The Gene Expression data set can be obtained from the NCBI GEO database using GEO accession number GSE7329 (or download directly from [here](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz)).
