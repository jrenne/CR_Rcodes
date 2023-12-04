================================================================================
An Analytical Framework to Price Long-Dated Climate-Exposed Assets

Pauline Chikhani and Jean-Paul Renne

This version: December 2023
================================================================================

The present codes allow to replicate the results of the paper entitled "An Analytical Framework to Price Long-Dated Climate-Exposed Assets" that can be found at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3881262.


--------------------------------------------------------------------------------
A - Requirements

The paper's results have been generated under R version 4.3.1 (2023-06-16) -- "Beagle Scouts".

The name of the R project is "CR_Rcodes.Rproj." To use the codes, run "main.R." The code will start with calibrating mu_T using the bisection method; this takes about 30 seconds. In "main.R", you can specify if you want to create tables and figures. For the latter, you can choose which ones. Details are provided directly in "main.R", as commented code.

Make sure the libraries called in "main.R" are available before running the codes.

Computing time: Using a set up with 2.3 GHz 8-Core Intel Core i9, all figures and tables (including those appearing in the supplemental material) are generated in about 15 minutes.


--------------------------------------------------------------------------------
B - The different folders

The different folders are as follows:

1- "data" contains data files used for calibration and plots.

2- "estimations" contains R files used to estimate/calibrate the model.

3- "outputs" contains scripts producing figures and tables; it also contains these figures and tables.

4- "procedures" contains sets of procedures used to solve the model, price assets, and generate the various outputs.

In addition to the present README file, the main folder also contains two files "main.R" that has to be sourced in order to replicate the paper's results.


--------------------------------------------------------------------------------
C - Outputs

Figures are stored in "outputs/Figures/". Tables are stored in "outputs/Tables/".


