# Bayesian feedback in the framework of ecological sciences

In this repository, the code used to obtain the results in the article *Bayesian feedback in the framework of ecological sciences* is provided. The repository includes the main script to reproduce these results (*Main_code_feedback.R*) and a folder (*functions*) where the auxiliary scripts used in the main script are stored.

## Summary of the scripts

In particular, the different scripts are summarized below as a brief introduction to them.

  1. `Main_code_feedback.R`: This script contains the code to define the parameters that characterize the different scenarios, along with the function for executing the analyses and the feedback between models. In this process, functions defined in other scripts are called. These scripts are located in the functions folder.

  2. `dependent_model.R`: This script defines the function for analyzing dependent data using the preferential model.

  3. `fb_dependent_model.R`: This script defines the function for analyzing dependent data using a preferential model with feedback.

  4. `fb_independent_model.R`: This script defines the function for analyzing independent data using a geostatistical model with feedback.

  5. `independent_model.R`: This script defines the function for analyzing independent data using a geostatistical model.

  6. `sampling_functions.R`: This script defines the functions for the different sampling processes.

  7. `several_functions.R`: This script defines various general-purpose functions.

  8. `simulation_function.R`: This script defines the function for simulating the underlying biomass process, given the parameters that characterize each scenario.


