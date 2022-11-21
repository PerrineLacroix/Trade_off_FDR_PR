# Trade off between FDR and PR in model selection

This repository provides codes that generated the data and figures presented in the paper [Trade-off between prediction and FDR for high-dimensional Gaussian model selection](...).
In particular, it contains the implementation of the proposed algorithm.

## Summary

We propose a new calibration of the hyperparameter _K_ appearing in the penalty function during the penalized least-squares minimization in the ordered variable selection procedure. The goal is to avoid the selection of non active variables while maintaining predictive performances of the model selection procedure.

More precisely, the algorithm is tested in high-dimensional Gaussian linear regression. It is based on completely data-driven terms and realizes a trade-off between the predictive risk (PR) and the False Discovery Rate (FDR) controls.

## Composition of the Repository

The `script` folder contains tool functions to be launched.
The `Git_multiplicative_constant_K_creation_data_and_collection.R` function creates the data sets and the model collections. Then, its outputs are the inputs of the `Git_multiplicative_constant_K.R` function. This function:

- computes the PR and FDR values for each _K_ and each iteration,
- computes the lower and upper bounds of the FDR function,
- provides the implementation of the data-driven calibration algorithm,
- generates the figures included in the paper.

The `source` folder contains all the functions called by `Git_multiplicative_constant_K.R`. More precisely:

- `Git_cste_mult_variation_PR_FDP.R` : computes the least-squared, the PR and the FDR values for each model of the given collections, especially for the selected models.
- `Git_empirical_estimation.R` : applies the R function source `Git_cste_mult_variation_PR_FDP.R` on each available model collection.
- `Git_P_2r_estimation.R` : computes the _P\_2r_ terms estimation for each K of the considered grid.
- `Git_lower_bounds.R` : computes the theoretical lower bound on the FDR function with respect to the considered grid of K for all data sets.
- `Git_upper_bounds.R` : computes the theoretical upper bound on the FDR function with respect to the considered grid of K for all data sets.
- `Git_lower_bound_estimated_one_dataset.R` : evaluates the theoretical lower bound with the available dataset.
- `Git_upper_bound_estimated_one_dataset.R` : evaluates the theoretical upper bound with the available dataset.

The `data` folder is empty. Data generated by the `Git_multiplicative_constant_K_creation_data_and_collection.R` and `Git_multiplicative_constant_K.R` functions are saved in this folder.

The `results` folder is empty. Figures generated by the `Git_multiplicative_constant_K.R` function are saved in this folder.


## Versions and librairies

This code was developed and tested with Linux version 5.4.0 and `R version 3.6.3`.

It requires the following packages :
- `xtable`
- `ggplot2`
- `capushe`
