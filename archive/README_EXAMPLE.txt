# ====================================================
# FISCAL LIMITS AND THE PRICING OF EUROBONDS
# ====================================================
# Kevin PALLARA and Jean-Paul RENNE
# This version: July 2022.
# ====================================================


These codes replicate the results of the paper entitled
"Fiscal limits and the pricing of Eurobonds",
by Kevin Pallara and Jean-Paul Renne.

The paper us available at
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3891358


# ====================================================
A- How to run the codes?
# ====================================================

Run "main.R" to launch the estimation and produce the paper's charts & tables.

Before sourcing "main.R", modify the working directory (line 16 of main.R).

"main.R" launches the initial model simulation and, then, runs "run_estim.R". 

In "run_estim.R", by default, the baseline model is estimated by MLE and outputs are produced (charts & tables).

To launch the estimation of the full set of sensitivity exercises, the indicator "indic_estim_only_Baseline" should be set to 0.

In case results do not need to be saved and outputs do not need to be produced, the indicators "indic_save" and "indic_produce_outputs" should be set to 0.

(To launch the Bayesian estimation, "indic_run_Bayesian" should be set to 1.)

Note: using a laptop with a 2.3 GHz processor, estimating the baseline model and producing all tables and charts takes approximately 20 minutes.


# ====================================================
B- Organisation of codes
# ====================================================

The folders are as follows:

- "estimation" contains codes allowing for the estimation of the models.
- "procedures" contains various procedures, organised in different scripts according to their use.
- "data" contains the data used by the estimation (see next section).
- "results" contains R data files with the specification of previously saved models.
- "load.data.files" contains scripts that load and organise the data. 
- "outputs" contains output charts and tables, as well as the scripts producing these output files.
- "make.outputs" contains scripts generating the outputs.
- "initial.model" contains scripts that initialise the model (used as starting point of estimations).


# ====================================================
C- Data
# ====================================================

The data are gathered in the "data" folder. There is one folder per country. There is also a folder containing the prices of Eurobond "proxies" used to generate one of the paper's figures.


