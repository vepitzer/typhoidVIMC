# typhoidVIMC
VIMC typhoid transmission model

Data inputs (provided by VIMC) are not included.

The main file for running the model is "vax_sims_vimc_202111_withtables_bothscen.m". This file also includes the code for generating the burden tables and stochastic parameters table, which are saved as csv files. The differential equations for the model are defined in "fn_SIR_vaccine_vnv_strata.m" and "fn_SIR_vaccine_vnv_strata_h2o.m". The former is used for the burn-in period while the latter includes a decrease in the transmission rate over time associated with improvements in WASH.

**December 2023 update**
Code for running the 2023 VIMC full model update is "vax_sims_vimc_2023.m". The main update to the code is that we have changed the way we estimate years of life lost to death from typhoid to use data on life expectancy at age of death (rather than life expectancy at birth); see lines 235-257. The output2023 data structure is converted to the VIMC tables using the code in "vimc_tables_2023.m".
