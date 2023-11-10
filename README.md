
# MScThesis

<!-- badges: start -->
<!-- badges: end -->

This project contains the R and python code for my M.Sc. Thesis with the topic "Estimation of higher-order species interactions from ecological time series".
It presents a collection of scripts and datasets that are used in the thesis.
The repository is structured into `explore`, `Python`, and `R` directories, each containing subdirectories with scripts and data related to the specific field.


## Repository Structure

### `explore` Directory

Contains all R code for the preparation of the datasets, including loading the data, cleaning, converting to `phyloseq` objects, aggregation, filtering, and saving pre-processed time series as CSV files.

- `01a-timeseries-BioTIME`: Includes all BioTIME studies.
- `01b-timeseries-CLVpaper`: Contains the C. Diff. Bucci dataset.
- `01c-timeseries-miaTIME`: Comprises the Silverman Artificial Gut dataset.
- `01d-timeseries-HumanGutData-Karwowska-paper`: Encompasses all microbial Human Gut datasets (for the subjects "donorA", "donorB", "male", "female").
- `01e-timeseries-miaSim_files`: Describes the simulation of the miaSimS4 time series, including adding noise.
- `01f-timeseries-NODEBNGMpaper`: Contains the 3DLV tri-trophic Lotka-Volterra time series and the Ushio fish population dataset.
- `01g-timeseries-simulation-VdP`: Describes the simulation of the Van der Pol equation plus adding noise.
- `02-filter_and_group_ts`: Structures the process of filtering and grouping the time series and saving them as CSV files.

### `Python` Directory

Contains scripts for alr transformation, the application of compositional Lotka Volterra and DeepMoD, and analysis plus creating plots.

- `ALR_transformation`: Code for ALR transformation of datasets, along with a folder for all ALR-transformed time series in CSV format.
- `CompLotkaVolterra`: Provides scripts for both generalized and compositional Lotka-Volterra applications, plus code from the CLV repository.
- `DeepMoD`: Contains code to run the DeePyMoD framework on time series data.
- `Analysis_and_Plots`: Includes code for post-modeling analyses and the construction of plots for the thesis.

### `R` Directory

Includes scripts for running the NODEBNGM model and some help functions.

- `NODEBNGM`: Contains the data and script from https://github.com/WillemBonnaffe/NODEBNGM/tree/main plus two additional scripts: `m_main_general.r` and `run_NODEBNGM.R` to run the NODEBNGM algorithm on our time series datasets.
- `Functions`: Houses helper functions "merge-methods-modified" and "fct_save_time_series".


## Getting Started

To use the code and datasets provided in this repository, clone the repo to your local machine using the following command:

```bash
git clone https://github.com/mariaproebstl/MScThesis.git
```


## Additional Material: Output Files

Additional output material from the application of generalized and compositional Lotka Volterra, as well as DeepMoD and NODEBGNM to the syntheic and real-world datasets can be found on

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10103346.svg)](https://doi.org/10.5281/zenodo.10103346)