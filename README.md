# typhoidVIMC
VIMC typhoid transmission model

Data inputs (provided by VIMC) are not included.

The main file for running the model is "vax_sims_vimc_202111_withtables_bothscen.m". This file also includes the code for generating the burden tables and stochastic parameters table, which are saved as csv files. The differential equations for the model are defined in "fn_SIR_vaccine_vnv_strata.m" and "fn_SIR_vaccine_vnv_strata_h2o.m". The former is used for the burn-in period while the latter includes a decrease in the transmission rate over time associated with improvements in WASH.
