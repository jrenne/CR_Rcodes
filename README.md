
---

# An Analytical Framework to Price Long-Dated Climate-Exposed Assets

**Authors:** Pauline Chikhani and Jean-Paul Renne  
**Corresponding Author:** [Jean-Paul Renne](mailto:jean-paul.renne@unil.ch)  
**Version:** February 2025  

This repository contains the code necessary to replicate the results of the paper titled "[An Analytical Framework to Price Long-Dated Climate-Exposed Assets](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3881262)."

---

## A - Requirements

The analysis has been conducted using R version 4.3.1 (2023-06-16) -- "Beagle Scouts."

### Project Setup

1. **R Project**: The project is organized in the R project file named `CR_Rcodes.Rproj`.
2. **Execution**: To replicate the results, execute the script `main.R`. 
   - The script begins by calibrating `mu_T` using the bisection method, which takes approximately 30 seconds.
   - Within `main.R`, you have the option to generate specific tables and figures, with detailed instructions included as commented code.

### Libraries

Ensure that all necessary libraries referenced in `main.R` are installed prior to running the code.

### Computing Time

On a setup with a 2.3 GHz 8-Core Intel Core i9 processor, all figures and tables (including those in the supplemental materials) can be generated in about 15 minutes.

---

## B - Directory Structure

The repository is organized into several folders, detailed as follows:

1. **data/**: Contains data files used for calibration and plotting.
2. **estimations/**: Contains R scripts for estimating and calibrating the model.
3. **outputs/**: Includes scripts for producing figures and tables, along with the generated outputs.
4. **procedures/**: Contains various procedures used to solve the model, price assets, and generate outputs.

In addition to this README file, the main directory includes two essential files:
- **main.R**: The primary script to source for replicating the paper's results.

---

## C - Outputs

- **Figures**: Stored in `outputs/Figures/`.
- **Tables**: Located in `outputs/Tables/`.

---

If you have any questions or require further assistance, please reach out to the corresponding author.

--- 



