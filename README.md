# JW_2024_Clinch-River-Mussles

The intent of this project is to investigate the population dynamics of adjacent mussel beds in the Clinch River whose populations interact non-symmetrically via the movement of glochidia during spread of a pathogen and, separately, while absorbing pollutants from the river.  The population is assumed to reproduce for a single season each year and only undergo aging and mortality otherwise, so population sizes for each bed in the population follow a 4-periodic, nonautonomous, recurrence relation.  The system allows for emigration of glochidia but not immigration.

For the purposes of the initial manuscript, we investigated populations of 5 beds isolated beds and 5 beds embedded in a larger population of 9 beds.  For simplicity, these two scenarios were given separate m-files and are denoted throughout with either a "5" of "9" suffix.  Each file is self-contained and all parameters are set within the file rather than being passed as parameters.  

**DemographOnly5.m** and **DemographyOnly9.m**: The files generate baseline simulations for the population in the absence increased mortality caused by pathogens or pollutants.  The number of time-steps for the recurrence relation is given at the top of the file.  The next several lines establish paramters for plotting.  The variables nBeds and nAge should not be changed.  The initial state of the population is set on line 17 with several values of potentiall interesting initial states commented out.  The "SubX" populations require that the State_SubX_Prop matrices be loaded and be appropriately named in the variabled space.  Sruvival probabilities, p_surv, and fecundity, fec, are set in lines 31 and 32 and used to build the appropriate projection matrices.  Parameters b_m1, b_0, b_p1, and b_p2 determine what proportion of glochidia produceed by bed i move to bed i-1 (b_m1), stay at i (b_0), move to i+1 (b_p1), and move to i+2 (b_p2).

**mature_demX.mat**: These files store matrices of the mature population size (excludes glochidia) once the demography only model has reached a stable geometric cycle (depends only on the periodicity of the system).  These are used to set initial coniditions for variations on the model.  
**PathogenX.m**:  These files are used to model the population of 5 or 9 beds in the presence of a pathogen for which the probability of infection for an individual is fixed (assuming the infection is present in the apprpriate bed) and movement of infection requires relocation of infected glochidia.

**DensityPathogenX.m**: These files model the population assuming the presence of a density-dependent infection process.  

**ContaminationX.m**: These files model the population in the presence of a pollutant or contimant.  Appropriate parameter descriptions can be found in the initial paper draft.


