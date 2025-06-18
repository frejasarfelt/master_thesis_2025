# master_thesis

Below is listed the different files. It lists what their code contains and examples of figures, that were created.

## functions.R

All functions are gathered in this file. This keeps maintaining and updates easy.

## ABM_improved.R

Simulations of agent based model, that divides people into n groups with possibility of different vaccination coverages. This investigates the meaning of heterogenity in subpopulations of a total population.

![Heterogenity: Two subpopulations - 9000 (97%) and 1000 (70%)](Figures/ABM_heterogenitet_2grupper.png){width="225"}

![Outbreak size distribution shown in a cumulative step plot for three different populations.](Figures/ABM_cumulative_step_plot.png){width="430"}

![Facetted scatterplot of scenario 3 showing the distribution of infections between the two subpopulations.](Figures/ABM_facetted_scatterplot_scn3.png){width="431"}

![Distribution of infected individuals in each group, when the outbreak starts in a specific population.](Figures/ABM_infected_by_origin_scn2.png){width="431"}

## q(x).R

The function q(x) from Jansen & Stollenwerk is implemented. It's the probability of an outbreak size being exactly of size x.

![q(x) plotted for different R values](Figures/q(x)_simulations.png){width="430"}

## q(x)\_cumulative.R

The function q(x) from Jansen & Stollenwerk is implemented cumulative. It's the probability of an outbreak size being of size x or greater.

![Cumulative q(x) plotted for different R values. Staircase: GW simulations. Dotted line: cumulative q(x)](Figures/q(x)_cumulative.png){width="430"}

## outbreak_distributions_plots.R

Contains cumulative distribution plots of countries with at least 20 outbreaks of size minimum 5.

All plots assumes x = 1 is not in the data. They are plotted both with and without q(x)

![](Figures/cum_prop_all_countries.png){width="496"}

## Censored_MLE_R.R

Computes maximum likelihood estimates of R for countries with at least 20 outbreaks of size minimum 5.

Also computes bootstrap CI and Fisher CI. And visualizes it.

![](images/Skærmbillede%202025-04-16%20kl.%2015.42.24.png)

## LLRT_Wald_test.R

For testing distributions against each other.

Also bootstrap test.

## Timeseries_outbreaks.R

Plots the outbreaks over time.

![](http://127.0.0.1:40829/graphics/deee78f7-cc2b-48ef-8195-bc64174f1c52.png)

## descriptive_analysis_of_data.R

This script produces all the numbers included in the descriptive analysis of TESSy data in the report.

## Binomial_Conditional_P.R

For age-specific data analysis.
