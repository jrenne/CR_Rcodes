<!-- --- -->
<!-- output: pdf_document -->
<!-- --- -->

# An Analytical Framework to Price Long-Dated Climate-Exposed Assets

**Authors:** Pauline Chikhani and Jean-Paul Renne  
**Corresponding Author:** [Jean-Paul Renne](mailto:jean-paul.renne@unil.ch)  
**Version:** July 2025  


---

## Overview

This repository contains the code necessary to replicate the results of the paper titled "[An Analytical Framework to Price Long-Dated Climate-Exposed Assets](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3881262)."

---

## Requirements

The analysis published in the paper has been conducted using R version 4.3.1 (2023-06-16) -- "Beagle Scouts."

### Instructions to Replicators

To replicate the results, execute the script `main.R`. 
   - The script begins by calibrating `mu_T` using the bisection method, which takes approximately 30 seconds.
   - Within `main.R`, you have the option to generate specific tables and figures, with detailed instructions included as commented code. By default, the code produces all the figures and tables in the paper.
   
### Libraries

Ensure that all necessary libraries referenced in `main.R` are installed prior to running the code. Libraries needed: parallel, doParallel, mgcv, colorspace, broman, optimx, MASS, expm.

### Runtime and Storage Requirements

- On a setup with a 2.3 GHz 8-Core Intel Core i9 processor, all figures and tables (including those in the supplemental materials) can be generated in about 15 minutes. Depending on the setup, approximate time needed to reproduce the analyses on a standard (2025) desktop machine: 10-30 minutes
- Approximate storage space needed: < 5 MBytes

### Controlled Randomness

Random seed is set:

- at line 9 of program make_figure_gamma0_distri.R. This script produces the figure illustrating the gamma0 distribution (Figure 9).
- within function `compute.SCC.Lemoine`; this function is included in the script `procedures/functions_other_models.R`. This function utilizes Lemoine (2021) model to compute yields and Social Cost of Carbons (SCC). It is used in `outputs/make_figures/make_figure_YC_RF.R` to produce Figure 6 (yield curves) and Figure VI.2 (Term structure of temperature risk premiums).
- The other results of the paper are not based on pseudo-random number generators (PRNGs).


---

## Directory Structure

The repository is organized into several folders, detailed as follows:

1. **data/**: Contains data files used for calibration and plotting.
2. **estimations/**: Contains R scripts for estimating and calibrating the model.
3. **outputs/**: Includes scripts for producing figures and tables, along with the generated outputs (see table below).
4. **procedures/**: Contains various procedures used to solve the model, price assets, and generate outputs.

The main directory includes **main.R**, that is the primary script to source for replicating the paper's results.


### List of Tables and Figures Produced by the Codes

The provided code reproduces tables and figures in the paper (as well as in the supplemental appendix). Figure 1 is an exception; it is a schema created using the tikz package of Latex. Tables 3, 4, 6 are exceptions; the first and third contain formulas and is directly written in Latex; the second contains SCC estimates taken from the literature.

All the scripts producing figures are in the folder `outputs/make_figures/`. All the scripts producing tables are in the folder `outputs/make_tables/`. All produced figures are in pdf format, and stored in the folder `output/Figures/`. All produced tables are written in Latex, and stored as txt files in `outputs/Tables/`.




| Figure/Table #    | Program                  | Output file                      |
|-------------------|--------------------------|----------------------------------|
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
| Figure III.1          | `make_figure_calibration.R` |  `Figure_Calibration.pdf`  |
| Figure V.1          | `make_figure_mu.R` |  `FFigure_Mitigation_comparison.pdf`  |
| Figure VI.1          | `make_figure_RCP_to_TAT.R` |  `Figure_SCCvsTRP.pdf`  |
| Figure VI.2          | `make_figure_YC_RF.R` |  `Figure_TRP_comparison.pdf`  |
| Figure VI.3          | `make_figure_Merton.R` |  `Figure_Merton2.pdf`  |
| Table 1          | `make_table_Estimated_Param.R` |  `table_param_est.txt`  |
| Table 2          | `make_table_utility_solution.R` |  `table_utility_solution.txt`  |
| Table 5          | `make_table_SCC.R` |  `table_SCC_sensitiv.txt`  |
| Table 7          | `make_table_param.R` |  `table_param.txt`  |
| Table VI.1          | `make_table_SCC.R` |  `table_TRP_sensitiv.txt`  |
| Table VI.2          | `make_table_SCC.R` |  `table_LTR_sensitiv.txt`  |



## Data Availability and Provenance Statements

### Statement about Rights

We certify that the authors of the manuscript have legitimate access to and permission to use the data used in this manuscript. We also certify that the authors of the manuscript have permission to redistribute/publish the data contained within this replication package. Details are given below. The remaining code (beyond the data aforementioned) is licensed by the authors under the [CC-BY-NC 4.0 Unported license](https://creativecommons.org/licenses/by-nc/4.0/). Further redistribution and use must adhere to the relevant license terms. Appropriate permission are documented in the [LICENSE.txt](LICENSE.txt) file.

Some data contained in this repository (to construct Figure 2 of the paper) are obtained from the [USEPA/scghg](https://github.com/USEPA/scghg) repository, which provides codes to replicate results displayed in EPA (2023). The data is provided under the MIT License, which permits free use, modification, and distribution. We gratefully acknowledge Bryan Parthum for his support in helping us access and work with the data.

- Some data contained in the EPA repository (aforementioned) are outputs generated from [HECTOR](https://jgcri.github.io/hector/), an open-source model licensed under the GNU General Public License v3.0. (These outputs are used in Figure 2 of the present paper.) As such, these outputs are also subject to the terms of the GPL v3.0. Users should review the GPL v3.0 license to understand their rights and obligations when using, modifying, or distributing these outputs. The original HECTOR license can be found [here](http://www.gnu.org/licenses/gpl-3.0.en.html).

- Some data contained in the EPA repository (aforementioned) are outputs generated from the [FaIR model (Finite-amplitude Impulse-Response)](https://docs.fairmodel.net/en/latest/); replication codes for this model are provided under the [Apache License, Version 2.0](http://www.apache.org/licenses/). (These outputs are used in Figure 2 of the present paper.) The resulting outputs are subject to the terms of the license, which allow use, modification, and distribution accordingly.

The present replication package includes outputs generated with the MAGICC6.0 model; these originate from [Christian Traeger´s publicly available replication repository](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view), associated with Traeger (2023). In the readme file associated with the latter package, the author of the package certifies that he has legitimate access and permission to publish the data generated with the MAGICC6.0 model, which is governed by the [CC-BY-NC-SA 3.0 Unported license](https://creativecommons.org/licenses/by-nc-sa/3.0/). To load these data from Traeger's package, we have made use of the (Matlab) codes present in Christian Traeger's package; these codes are licensed under the [CC-BY-NC 4.0 Unported license](https://creativecommons.org/licenses/by-nc/4.0/).

In our Figure 2, the line labeled "ACE-Joos" is based on data extracted from [Christian Traeger´s repository](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view), which utilizes inputs originally obtained from Joos et al. (2013); the latter work is licensed under [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/).

Some data used in this package, including the outputs in Figure 6, originate from the replication package of Bauer and Rudebusch (2023); these data are licensed under the [CC0 1.0 Universal Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/). This license effectively dedicates the data to the public domain, allowing free use, modification, and distribution without restrictions, as detailed in the full legal code. The data were extracted from Michael Bauer´s [personal website](https://www.michaeldbauer.com/files/sdr/sdrs.csv) and are used here in accordance with the license terms.


### Details on each Data Source


| Data.Name  | Data.Files | Location | Provided | Citation |
| -- | -- | -- | -- | -- | 
| Social Discount Rates | `Bauer_Rudebusch_sdrs.csv` | `data/` | TRUE | Bauer and Rudebusch (2023) |
| IRF from EPA report | `Figure223_EPA_data.Rdat` | `data/` | TRUE | Environmental Protection Agency (2023) |
| IRF from ACE | `IRF_Traeger_5y.csv` | `data/` | TRUE | Traeger (2023) |
| RCP temperature scenario and standard dev. | `mean_ssp.txt` | `data/` | TRUE | XXXX |
| RCPs based on ACE model | `RCP_Mat_ACE.csv` | `data/` | TRUE | Traeger (2023) |
| RCPs based on MAGICC6.0 model | `RCP_Mat_MAGICC.csv` | `data/` | TRUE | Meinshausen et al. (2011), Traeger (2023) |
| DICE mitigation path | `mu.csv` | `data/` | TRUE | DICE XXXX |


## References

[Bauer, M. D. and G. D. Rudebusch (2023)](https://doi.org/10.1162/rest_a_01109). The Rising Cost of Climate Change: Evidence from the Bond Market. The Review of Economics and Statistics 105(5), 1255-1270.

[Environmental Protection Agency (2023)](https://www.epa.gov/environmental-economics/scghg). EPA Report on the Social Cost of Greenhouse Gases: Estimates Incorporating Recent Scientific Advances. EPA Report Docket ID No. EPA-HQ-OAR-2021-
0317, EPA.

[Joos, F., R. Roth, J. S. Fuglestvedt, G. P. Peters, I. G. Enting, W. von Bloh, V. Brovkin, E. J. Burke, M. Eby, N. R. Edwards, T. Friedrich, T. L. Frölicher, P. R. Halloran, P. B. Holden, C. Jones, T. Kleinen, F. T. Mackenzie, K. Matsumoto, M. Meinshausen, G.-K. Plattner, A. Reisinger, J. Segschneider, G. Shaffer, M. Steinacher, K. Strassmann, K. Tanaka, A. Timmermann, and A. J. Weaver (2013)](https://acp.copernicus.org/articles/13/2793/2013/). Carbon Dioxide and Climate Impulse Response Functions for the Computation of Greenhouse Gas Metrics: a Multi-Model Analysis. Atmospheric Chemistry and Physics 13(5), 2793-2825.

[Lemoine, D. (2021)](https://www.journals.uchicago.edu/doi/abs/10.1086/710667?journalCode=jaere). The Climate Risk Premium: How Uncertainty Affects the Social Cost of Carbon.
Journal of the Association of Environmental and Resource Economists 8(1), 27-57.

[Meinshausen, M., S. C. B. Raper, and T. M. L. Wigley](https://acp.copernicus.org/articles/11/1417/2011/) (2011). Emulating Coupled Atmosphere-Ocean and Carbon Cycle Models with a Simpler Model, MAGICC6 - Part 1: Model Description and Calibration. Atmospheric Chemistry and Physics 11(4), 1417-1456.

[Traeger, C. P. (2023)](https://www.aeaweb.org/articles?id=10.1257/pol.20210297). ACE-Analytic Climate Economy. American Economic Journal: Economic Policy 15(3),
372-406.

---

If you have any questions or require further assistance, please reach out to the corresponding author.
