R code and detailed simulation results for the article "Logistic regression with covariate-dependent probability of misclassification"

List of the supplementary files:

estim09.r - fits the models M5, gamma1=0 (LZ), Sp=1 (M4), and beta1=0 - Y doesn't depend on X (M0)

estim09_M5_only.r - fits the model M5 only (this was used in the bias, SE and RMSE calculations)

analyse.r - calls estim09.r and summarizes the results (estimates, SEs, model comparisons)

simul_bias.r - simulation to assess bias, SE, RMSE of estimates

simul_power.r - simulation to assess power to detect that Y depends on X, Se depends on Z, and Sp < 1

simulation_results.zip - simulation results (10 rdata files, 5 for bias, 5 for power)

simulation_results_summary.pdf - summary of the simulation results

analysis_alcohol_drinking_driver.r - applies the method to the alcohol consumption data

d3_age_recoded.rdata - data file to the alcohol consumption analysis (Application example 1)

analysis_voting_Czechia.r - applies the method to the voting data

CzechRep.rdata - data file to the voting analysis (Application example 2)
