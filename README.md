# JW_2024_Clinch-River-Mussles

The intent of this project is to investigate the population dynamics of adjacent mussel beds in the Clinch River whose populations interact non-symmetrically via the movement of glochidia during spread of a pathogen and, separately, while absorbing pollutants from the river.  The population is assumed to reproduce for a single season each year and only undergo aging and mortality otherwise, so population sizes for each bed in the population follow a 4-periodic, nonautonomous, recurrence relation.

For the purposes of the initial manuscript, we investigated populations of 5 beds isolated beds and 5 beds embedded in a larger population of 9 beds.  For simplicity, these two scenarios were given separate m-files and are denoted throughout with either a "5" of "9" suffix.  Each file is self-contained and all parameters are set within the file rather than being passed as parameters.  

**DemographOnly5** and **DemographyOnly9**: The files generate baseline simulations for the population in the absence increased mortality caused by pathogens or pollutants.  The number of time-steps for the recurrence relation is given at the top of the file.  The next several lines establish paramters for plotting.  The initial state of the population is set on line 17 with several values of potentiall interesting initial states commented out.  The "SubX" populations require that the State_SubX_Prop matrices be loaded and be appropriately named in the variabled space.


