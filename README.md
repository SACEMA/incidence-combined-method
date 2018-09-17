# Combining synthetic cohort and recency biomarker incidence estimates

Analytical code for estimating incidence using a combination of synthetic cohort
and recency biomarker-based approaches.

## Published article

The approach and methodology implemented in this code is described in the following article:

* Grebe E, Welte A, Johnson LF, Van Cutsem G, Puren A, Ellman T, Etard J-F,
CEPHIA, Huerga H. Population-level HIV incidence estimates using a combination
of synthetic cohort and recency biomarker approaches in KwaZulu-Natal, South
Africa. PLOS ONE. 2018;13(9):e0203638. PMID:30212513. doi:[10.1371/journal.pone.0203638](https://doi.org/10.1371/journal.pone.0203638). Preprint: bioRxiv. 2018. doi:[10.1101/306464](https://doi.org/10.1101/306464).

The article and appendices are also included in the directory `manuscript` for convenience.

## Data

A minimal, anonymised dataset is available with the article from the journal cite (doi:[10.1371/journal.pone.0203638.s003](https://doi.org/10.1371/journal.pone.0203638.s003)). For convenience, a copy of the anonymised dataset is also included in this repository. The dataset is contained in the file named `mbongolwane_eshowe_2013_anon_minimal.csv`. Note that the code was left as used in the primary analysis, and presumes the presence of several data files not included here.

The variables `id`, `cluster` and `ward` are randomised participant, cluster (primary sampling unit) and electoral ward (stratum) identifiers, respectively. To replicate the multistage sampling frame during boostrapping, wards were resampled with replacement, and within each sampled ward clusters were resampled with replacement. The `age_years` variable captures age in years, rounded to whole numbers as a safeguard against de-anonymisation. Final HIV status is captured in `hiv_status` and participants identified as HIV-infected using nucleic acid amplification testing have “True” in the `naat_yield` variable. Final LAg normalised optical density and Bio-Rad Avidity avidity index are captured in `LAgODn` and `BioRadAI`, respectively, and viral load in `viral_load`. For convenience, recency status according to the RITA utilised in this analysis is captured in recent.

In addition, the mortality table utilised in our analysis is provided (`mortality_table.feather`). These estimates are output from the [THEMBISA demographic model](https://www.thembisa.org/) (Johnson et al.) for KwaZulu-Natal and are no longer the most recent estimates. Updated estimates are available from the THEMBISA website.

The population table is based on the 2011 South African census is also included (`poptable.feather`), and provide sex and age-specific population estimates for Mbongolwane and Eshowe.

## Reproducing the results

If the code and data provided here are used to reproduce the results in the published article, small differences from the reported results are inevitable, given that the finely-grained age data has been rounded to the closest integer age.

## Using the code

This software is provided in the form used to produce the analyses described in the paper cited above. It has not been set up for easy re-use, for example, hard-coded paths and  for easy has not been code has not been set up for easy use with the minimal anonymised dataset or for application to alternative datasets.

The following scripts are provided:

* `incidence_functions.R`: Generic functions utilised in the main estimation scripts
* `estimate_incidence.R`: Produces the primary age-specific incidence estimates
* `average_incidence.R`: Estimates "average incidence" for age ranges, weighted by sampling density or susceptible population density
* `sensitivity_analysis.R`: Sensitivity analyses with respect to FRR, average incidence weighting scheme and possible change in the age-structured prevalence in secular time
* `graphics.R`: Produces the plots included in the manuscript.
