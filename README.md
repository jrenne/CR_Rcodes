---
output:
  word_document: default
  html_document: default
---
# Replication package for "An Analytical Framework to Price Long-Dated Climate-Exposed Assets"

**Authors:** Pauline Chikhani and Jean-Paul Renne  
**Corresponding Author:** [Jean-Paul Renne](mailto:jean-paul.renne@unil.ch)  
**Version:** July 2025  


---

## Overview

This repository contains the code necessary to replicate the results of the paper titled "[An Analytical Framework to Price Long-Dated Climate-Exposed Assets](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3881262)."

---

## Runtime and Requirements


- **Software and Version Details for the Analysis**: The analysis published in the paper has been conducted using R version 4.3.1 (2023-06-16) -- "Beagle Scouts."

- **Libraries**: Ensure that all necessary libraries referenced in `main.R` are installed prior to running the code. Libraries needed: parallel, doParallel, colorspace, broman, optimx, MASS, expm.

- **Runtime**: On a setup with a 2.3 GHz 8-Core Intel Core i9 processor, all figures and tables (including those in the supplemental materials) can be generated in about 15 minutes. Depending on the setup, approximate time needed to reproduce the analyses on a standard (2025) desktop machine: 10-30 minutes

- **Storage**: Approximate storage space needed: < 15 MBytes

---

## Instructions to Replicators

### Replication using `main.R`

To replicate the results, the replicator can execute the script `main.R`.

### Step-by-step Execution Guide

This section outlines how to manually run the replication code. That is, it explains the different steps followed by `main.R`.

#### Set Working Directory

Begin by setting the working directory to the root folder containing the replication files.

```
setwd("~/Research/CR_Rcodes") # Adjust the path to your local setup
```

#### Libraries

To run the replication scripts, the following libraries are required:

* Parallel computing: `parallel`, `doParallel`
* Mathematical computation and optimization: `optimx`, `MASS`, `expm`
* Plotting and visualization: `colorspace`, `broman`

You can install missing packages with install.packages("package_name").

```
library(parallel)
library(doParallel)
library(colorspace)
library(broman)
library(optimx)
library(MASS)
library(expm)
```

#### Set Global Parameters

You should also specify how many CPU cores to use for parallel computations. Parallel processing is employed to generate certain figures that require multiple model runs, such as sensitivity analyses. By default, `number.of.cores` is set to 8. If the computer of the replicator has fewer than 8 cores, the maximum available number of cores will be used.

```
# Set number of cores used for parallel computing (sensitivity analysis) -------
number.of.cores <- 8 # Adjust based on your machine
```

#### Load Function Libraries

Load the custom R functions used throughout the code inside the folder `procedures`:

* The general library - `functions_general.R`: Core model functions used in the paper
* The figures library - `functions_figures.R`: Functions dedicated to figure generation
* The miscellaneous library - `functions_other_models.R`: Functions for alternative model comparisons

```
source("procedures/functions_general.R")
source("procedures/functions_figures.R")
source("procedures/functions_other_models.R")
```

#### Calibration the Baseline Model

This step initializes the model with the baseline parameter calibration used in the paper. The code can be found in the folder `estimations`.

```
# Calibrate baseline model -----------------------------------------------------
source("estimations/load_ini_model.R")
```

#### Generate Tables and Figures

Set binary indicators `indic_*_paper` to specify whether the code should generate plots and/or tables. Use `1` to activate generation, and `0` to skip.

The replicator can specify which figures to generate using the plots variable.

* To generate all figures, use: `plots <- 1:105`
* To generate a specific figure, replace the range with the corresponding number (e.g., `plots <- 4` for Figure 4).
* Figures in the appendix are numbered starting from 101.

Tables are generated automatically when `indic_tables_paper <- 1`.

```
# Determine whether tables and figures are generated ---------------------------
#Binary operators: 0 = NO, 1 = YES.
indic_plots_paper  <- 1 #produce paper's plots? see description below
indic_tables_paper <- 1 #produce paper's table? see description below
plots <- 1:105 # Set to desired figure numbers
if(indic_plots_paper==1){
  source("outputs/plots_paper.R")
}
if(indic_tables_paper==1){
  print("Preparing tables (couple of minutes)")
  source("outputs/tables_paper.R")
}

```

#### Figure Index Reference

The correspondence between figure numbers and their output can be found:
- At the bottom of the `main.R` file
- At the top of the `plots_paper` file in the folder `outputs`
- In the header of each individual figure script inside `outputs/make_figure_*.R`

This allows you to identify and selectively run individual figure scripts if desired.

We recommend that the replicator runs the figures through the `plots_paper.R` file, as some figure-specific parameters are defined at the beginning of that file (lines 35 to 57).


---

### Controlled Randomness

Random seed is set:

- at line 9 of program make_figure_gamma0_distri.R. This script produces the figure illustrating the gamma0 distribution (Figure 9).
- within function `compute.SCC.Lemoine`; this function utilizes Lemoine (2021) model to compute yields and Social Cost of Carbons (SCC). It is included in the script `procedures/functions_other_models.R`. It is used in `outputs/make_figures/make_figure_YC_RF.R` to produce Figure 6 (yield curves) and Figure S.4 (Term structure of temperature risk premiums).
- The other results of the paper are not based on pseudo-random number generators.


---

## Directory Structure

The repository is organized into several folders, detailed as follows:

1. **data/**: Contains data files used for calibration and plotting (see details in the section on data sources).
2. **estimations/**: Contains two R scripts used for calibrating the model.

    - `compute_alpha.R`: This script determines the value of the alpha parameter used in the  calibration exercise.
    - `load_ini_model.R`: This script constructs the baseline model. The calibration of some parameters is based on the literature; others are calculated so as to achieve specific targets.
3. **outputs/**: Includes scripts for producing figures and tables, along with the generated outputs (see table below). Specifically, it contains the following elements:

    - `results`: This folder is initially empty; it stores `.RData` files containing the results of relatively lengthy calculations (namely: `sensitivity_SCC_CRRA.Rdat`, `sensitivity_SCC_EZ.Rdat`, `SCC_vs_TRP_CR.Rdat`, and `m0_4_ShinyApps.Rdat`).
    - `make_figures` and `make_tables`: These folders contain, respectively, the scripts that generate the figures and tables shown in the paper and its appendices. These scripts are detailed in the next table of this README.
    - `Figures` and `Tables`: These folders store the generated figures (in `.pdf` format) and tables (in `.txt` format), respectively.
    - `plots_paper.R` and `tables_paper.R`: These are two scripts that initiate the construction of figures and tables, respectively. They call the scripts contained in `make_figures` and `make_tables`, respectively.
4. **procedures/**: Contains various procedures used to solve the model, price assets, and generate outputs. Specifically:

    - `functions_figures.R`: This script contains various procedures used to produce the figures.
    - `functions_general.R`: this script contains various procedures used to calibrate and solve the model, and to price various types of assets.
    - `functions_other_models.R`: This script contains procedures used to calibrate and solve alternative models. They are used to produce plots and tables involving outputs of these alternative models.

The main directory includes **main.R**, that is the primary script to source for replicating the paper's results.


### List of Tables and Figures Produced by the Codes

The provided code reproduces the tables and figures in the paper, including those in the supplemental appendix. The exception is Figure 1, which is a schema created using the TikZ package in LaTeX. Tables 3, 4, and 6 are also exceptions; they contain formulas written directly in LaTeX.

All scripts used to generate figures are located in the folder `outputs/make_figures/`. All scripts for producing tables are in `outputs/make_tables/`. The resulting figures are saved in PDF format within the folder `outputs/Figures/`, and all tables are written in LaTeX and stored as `.txt` files in `outputs/Tables/`.


| Figure/Table #    | Program                  | Output file                      |
|-------------|--------------------------|----------------------------------|
| Figure 2          | `make_figure_IRF1Gt.R` |  `Figure_IRF1GtC.pdf`  |
| Figure 3          | `make_figure_Damage_comparison.R` |  `Figure_Damage_comparison.pdf`  |
| Figure 4          | `make_figure_Tpdf.R` |  `Figure_Tat_P_and_Q_vector_CI.pdf`  |
| Figure 5          | `make_figure_Hpdf.R` |  `Figure_SL_P_and_Q_vector_CI.pdf`  |
| Figure 6          | `make_figure_YC_RF.R` |  `Figure_YC_RF.pdf`  |
| Figure 7          | `make_figure_breakeveninflation.R` |  `Figure_BreakEvenInflation.pdf`  |
| Figure 8          | `make_figure_options.R` |  `Figure_Option_Digital.pdf`  |
| Figure 9          | `make_figure_gamma0_distri.R` |  `Figure_gamma0.pdf`  |
| Figure 10          | `make_figure_gamma0_distri.R` |  `Figure_gamma0_Damages.pdf`  |
| Figure 11          | `make_figure_RCP_to_TAT.R` |  `Figure_RCP_to_TAT.pdf`  |
| Figure S.1          | `make_figure_calibration.R` |  `Figure_Calibration.pdf`  |
| Figure S.2          | `make_figure_mu.R` |  `Figure_Mitigation_comparison.pdf`  |
| Figure S.3          | `make_figure_RCP_to_TAT.R` |  `Figure_SCCvsTRP.pdf`  |
| Figure S.4          | `make_figure_YC_RF.R` |  `Figure_TRP_comparison.pdf`  |
| Figure S.5          | `make_figure_Merton.R` |  `Figure_Merton2.pdf`  |
| Table 1          | `make_table_Estimated_Param.R` |  `table_param_est.txt`  |
| Table 2          | `make_table_utility_solution.R` |  `table_utility_solution.txt`  |
| Table 5          | `make_table_SCC.R` |  `table_SCC_sensitiv.txt`  |
| Table 7          | `make_table_param.R` |  `table_param.txt`  |
| Table S.1          | `make_table_SCC.R` |  `table_TRP_sensitiv.txt`  |
| Table S.2          | `make_table_SCC.R` |  `table_LTR_sensitiv.txt`  |


---

## Data Availability and Provenance Statements

### Statement about Rights

We certify that the authors of the manuscript have legitimate access to and permission to use the data included in this manuscript. We also certify that the authors have permission to redistribute and publish the data contained within this replication package. Details are provided below. The remaining code (beyond the data mentioned above) is licensed by the authors under the [CC-BY-NC 4.0 Unported license](https://creativecommons.org/licenses/by-nc/4.0/). Further redistribution and use must adhere to the relevant license terms. Relevant permissions are documented in the [LICENSE.txt](LICENSE.txt) file.

Some data used to construct Figure 2 are obtained from the [USEPA/scghg](https://github.com/USEPA/scghg) repository, which provides code to replicate EPA (2023) results. These data are provided under the MIT License, allowing free use, modification, and distribution. We gratefully acknowledge Bryan Parthum for his support in helping us access and work with this data.

- Some data in the EPA repository (mentioned above) are outputs generated from [HECTOR](https://jgcri.github.io/hector/), an open-source model licensed under the [GPL v3.0](http://www.gnu.org/licenses/gpl-3.0.en.html). These outputs are used in Figure 2 and are subject to the terms of that license. Users should review GPL v3.0 to understand their rights and obligations.

- Some data from the EPA repository are outputs generated from the [FaIR (Finite-amplitude Impulse-Response)] model, licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/). These outputs, used in Figure 2, are governed by that license, which permits use, modification, and distribution according to its terms.

The present package also includes outputs generated with the [MAGICC6.0](https://www.ucl.ac.uk/~ucfbarc/MAGICC/) model (Meinshausen et al., 2011), originating from Christian Traeger's publicly available replication repository ([link](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view)), associated with Traeger (2023). The author of that package certifies he has legitimate access and permission to publish these data, which are governed by the [CC-BY-NC-SA 3.0 Unported license](https://creativecommons.org/licenses/by-nc-sa/3.0/). We used the MATLAB codes in Traeger's package (also licensed under CC-BY-NC 4.0) to load and process these data.

In our Figure 2, the line labeled "ACE-Joos" is based on data extracted from Christian Traeger's repository, which utilizes inputs originally obtained from Joos et al. (2013); the latter work is licensed under [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/).

Some data in this package originate from the replication package of Bauer and Rudebusch (2023), and are licensed under the [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) Public Domain Dedication. These data, used to produce Figure 6, were obtained from Michael Bauer's [personal website](https://www.michaeldbauer.com/files/sdr/sdrs.csv) and are used here in accordance with the license.

One csv file includes data derived from an Excel file created by William Nordhaus and Lint Barrage. (It is used in Figure S.2 of the Supplemental Appendix; the series is indicated by "DICE2023".) The data and associated program are @ DICE2023-Excel-b-4-3-10-v18.3, with copyrights held by William Nordhaus and Lint Barrage. The owners provide free and open access to the program. Commercial users must obtain permission for use in its original form or with changes from the copyright holder.

### Details on each Data Source

All data files (listed in the following table) can be found in the `~/data/` folder. The codes make use of the following `.csv` files:

| Data.Name  | Data.Files  | Citation |
| ------- | ------- | ------- |
| Social Discount Rates | `sdrs.csv` | Bauer and Rudebusch (2023) |
| Fair IRF from EPA report | `EPA_fair.csv` | Environmental Protection Agency (2023) |
| Hector IRF from EPA report | `EPA_hector.csv` | Environmental Protection Agency (2023) |
| Magicc IRF from EPA report | `EPA_magicc.csv` | Environmental Protection Agency (2023) |
| IRF from ACE | `IRF_Traeger_5y.csv` | Traeger (2023) |
| RCPs based on ACE model | `RCP_Mat_ACE.csv` | Traeger (2023) |
| RCPs based on MAGICC6.0 model | `RCP_Mat_MAGICC.csv` | Meinshausen et al. (2011), Traeger (2023) |
| DICE mitigation path | `mu_DICE.csv` | DICE2023; Barrage and Nordhaus (2024) |


- The file `sdrs.csv` (Bauer and Rudebusch, 2023), is available on Michael Bauer's website; it can be downloaded [here](https://www.michaeldbauer.com/files/sdr/sdrs.csv).
- The three files `EPA_fair.csv`, `EPA_hector.csv`, and `EPA_magicc.csv` contain the data of impulse response functions shown in Figure 2.2.3 of the EPA (2023) report. To generate these three files: (1) download the [EPA replication package](https://github.com/USEPA/scghg/tree/figures), (2) append the following three lines of code at the end of `plot_temperature_anomaly.R` (that file is available [here](https://github.com/USEPA/scghg/blob/figures/EPA/code/figures/plot_temperature_anomaly.R); i.e., in the folder `/EPA/code/figures/` of the EPA package), (3) set `/EPA/` as the working directory, (4) source `plot_temperature_anomaly.R`. This will create the three `.csv` files in the working directory.
```r
      write.csv(fair,   "EPA_fair.csv")
      write.csv(hector, "EPA_hector.csv")
      write.csv(magicc, "EPA_magicc.csv")
```
- The three `.csv` files: `IRF_Traeger_5y.csv`, `RCP_Mat_ACE.csv`, and `RCP_Mat_MAGICC.csv` are based on the [replication codes associated with Traeger (2023)](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view). (The file `IRF_Traeger_5y.csv` contains impulse response functions shown in Traeger, 2023; `RCP_Mat_ACE.csv` and `RCP_Mat_MAGICC.csv` contain trajectories of temperatures and carbon concentrations associated with Traeger's ACE 2023 model and the MAGICC6.0 model of Meinshausen et al., 2011, respectively.) These files can be generated by sourcing the script `make_csv_files.m` in Matlab; this script is in the folder `data/csv_from_Traeger/`. Before sourcing `make_csv_files.m`, one has to download and save (in `data/csv_from_Traeger/`) two data files from Traeger's package, namely: `impulse_timestep_5_logi_001.mat` directly available [here](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view?path=/openicpsr/154141/fcr:versions/V1/impulse_timestep_5_logi_001.mat&type=file), and `MagiccOcean_Traeger.m`, directly available  [here](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view?path=/openicpsr/154141/fcr:versions/V1/MagiccOcean.mat&type=file).
- The file `mu_DICE.csv` contains the trajectory of emission control rates extracted from the Excel file with outputs of the DICE2023 model. This Excel file is part of the [supplementary material](https://bit.ly/3TwJ5nO) associated with Barrage and Nordhaus (2024). It can be accessed directly via [this link](https://yale.app.box.com/s/whlqcr7gtzdm4nxnrfhvap2hlzebuvvm/file/1361579245945). The extracted series of emission control rates is located in row 71 of this Excel file, within the worksheet named `Opt`.

---

## References

[Barrage, L. and W. Nordhaus (2024)](https://doi.org/10.1073/pnas.2312030121),   Policies, projections, and the social cost of carbon: Results from the DICE-2023 model, Proc. Natl. Acad. Sci. U.S.A. 121 (13).

[Bauer, M. D. and G. D. Rudebusch (2023)](https://doi.org/10.1162/rest_a_01109). The Rising Cost of Climate Change: Evidence from the Bond Market. The Review of Economics and Statistics 105(5), 1255-1270.

[Environmental Protection Agency (2023)](https://www.epa.gov/environmental-economics/scghg). EPA Report on the Social Cost of Greenhouse Gases: Estimates Incorporating Recent Scientific Advances. EPA Report Docket ID No. EPA-HQ-OAR-2021-
0317, EPA.

[Joos, F., R. Roth, J. S. Fuglestvedt, G. P. Peters, I. G. Enting, W. von Bloh, V. Brovkin, E. J. Burke, M. Eby, N. R. Edwards, T. Friedrich, T. L. Fr^licher, P. R. Halloran, P. B. Holden, C. Jones, T. Kleinen, F. T. Mackenzie, K. Matsumoto, M. Meinshausen, G.-K. Plattner, A. Reisinger, J. Segschneider, G. Shaffer, M. Steinacher, K. Strassmann, K. Tanaka, A. Timmermann, and A. J. Weaver (2013)](https://acp.copernicus.org/articles/13/2793/2013/). Carbon Dioxide and Climate Impulse Response Functions for the Computation of Greenhouse Gas Metrics: a Multi-Model Analysis. Atmospheric Chemistry and Physics 13(5), 2793-2825.

[Lemoine, D. (2021)](https://www.journals.uchicago.edu/doi/abs/10.1086/710667?journalCode=jaere). The Climate Risk Premium: How Uncertainty Affects the Social Cost of Carbon.
Journal of the Association of Environmental and Resource Economists 8(1), 27-57.

[Meinshausen, M., S. C. B. Raper, and T. M. L. Wigley](https://acp.copernicus.org/articles/11/1417/2011/) (2011). Emulating Coupled Atmosphere-Ocean and Carbon Cycle Models with a Simpler Model, MAGICC6 - Part 1: Model Description and Calibration. Atmospheric Chemistry and Physics 11(4), 1417-1456.

[Traeger, C. P. (2023)](https://www.aeaweb.org/articles?id=10.1257/pol.20210297). ACE-Analytic Climate Economy. American Economic Journal: Economic Policy 15(3),
372-406.

