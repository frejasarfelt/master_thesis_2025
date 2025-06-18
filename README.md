# master_thesis_2025

All files are listed below with information about what their code contains and examples of figures, that were created.

**NB: Due to data protection regulations, the data is not published. This causes several of the files to fail running.**

## functions.R

All functions are gathered in this file. This keeps maintaining and updates easy.

## ABM_improved.R

Simulations of agent based models (ABMs).

![](Figures/ABM_cumulative_step_plot.png){width="430"}

![](Figures/ABM_facetted_scatterplot_scn3.png){width="431"}

![](Figures/ABM_infected_by_origin_scn2.png){width="300"}

## Binomial_Conditional_P.R

For age-specific data analysis. 

![](Figures/binom_p.png){width="431"}

## Censored_MLE_R.R

Computes maximum likelihood estimates of R for countries with at least 20 outbreaks of size minimum 5.

Also computes bootstrap CI and Fisher CI.

## descriptive_analysis_of_data.R

Produces all the numbers included in the descriptive analysis of TESSy data in the report.

## LLRT_Wald_test.R

For testing distributions against each other.

Also bootstrap test.

## outbreak_distributions_plots.R

Contains cumulative distribution plots of countries with at least 20 outbreaks of size minimum 5.

All plots assumes x = 1 is not in the data.

![](Figures/cum_prop_all_countries_q(x).png){width="430"}

## q(x).R

The function q(x) from Jansen & Stollenwerk is implemented. It's the probability of an outbreak size being exactly of size x.
Further, uncertainty associated with sample size its included.

![](Figures/q(x)_simulations.png){width="430"}

![](Figures/SE_vs_R_SH.png){width="300"} ![](Figures/Coverage_CI_vs_n.png){width="300"}

## q(x)\_cumulative.R

The function q(x) from Jansen & Stollenwerk is implemented cumulative. It's the probability of an outbreak size being of size x or greater.

![](Figures/q(x)_cumulative.png){width="430"}

## QQ.R

Produces QQ plots of outbreak size distributions from empirical data estimated R's.

![](Figures/QQ_plot.png){width="430"}

## Timeseries_outbreaks.R

Plots the outbreaks over time.

![](Figures/outbreaks_over_time_by_country.png){width="430"}
