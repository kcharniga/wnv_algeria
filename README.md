# West Nile virus in Algeria

This repository includes data and code used in the study “Patterns of Emergence and Circulation of West Nile Virus in Algeria” by Hachid et al. 

The file called “data.csv” includes: sample year, province, PRNT90 result (positive or negative), and age in 10-year groups. These data can be used to create the map of seroprevalence in Figure 1 and the seroprevalence with 95% confidence intervals by age group in Figure 2. More detailed data on patient ages have not be shared to maintain privacy.

The file called "severe_cases_data_table_s2.xlsx" contains the data on WNND cases from supplementary table 2. 

The file called “posterior_samples.RDS” is a list of lists containing the posterior samples for the best-fitting serocatalytic models for each province (the vectors for each parameter have length equal to the number of post-warmup iterations times the number of chains). 

We have included three R scripts and one Stan file. Note that these scripts will not run without the full dataset.
- descriptive_analysis.R: this script contains univariable and multivariable analyses with logistic regression used to identify risk factors for WNV infection
- serocatalytic_models.R: this script contains code for running the serocatalytic models with Rsero. Code for the “recent outbreak” model is located at the bottom of the script  
- compute_info_criteria2_function: this script makes a slight modification to the compute_info_criteria function in Rsero by allowing for moment matching (line 72) in the calculation of the PSIS-LOO
- simple_model.stan: this script contains the code for the “recent outbreak” model in the paper

