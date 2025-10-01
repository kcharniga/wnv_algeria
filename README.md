# West Nile virus in Algeria

This repository includes seroprevalence data on West Nile virus from Algeria used in the study “Patterns of Emergence and Circulation of West Nile Virus in Algeria” by Hachid et al. The data include: sample year, province, PRNT90 result (positive or negative), and age in 10-year groups. These data can be used to create the map of seroprevalence in Figure 1 and the seroprevalence with 95% confidence intervals by age group in Figure 2. More detailed data on patient ages have not be shared to maintain privacy.

We have included three R scripts and one Stan file. Note that these scripts will not run without the full dataset.
- descriptive_analysis.R: this script contains univariable and multivariable analyses with logistic regression used to identify risk factors for WNV infection
- serocatalytic_models.R: this script contains code for running the serocatalytic models with Rsero. Code for the “recent outbreak” model is located at the bottom of the script  
- compute_info_criteria2_function: this script makes a slight modification to the compute_info_criteria function in Rsero by allowing for moment matching (line 72) in the calculation of the PSIS-LOO
- simple_model.stan: this script contains the code for the “recent outbreak” model in the paper

