The code in this folder reproduces the analysis in STAT8810 - Final Report: 

"COMPARING PSEUDO-INPUT SELECTION METHODS IN SPARSE GAUSSIAN PROCESS REGRESSION 
 WITH APPLICATIONS TO HEART RATE DATA"

by Hengrui Luo and Giovanni Nattino. 



The code is organized as follows:

- "simulateGP.R" contains functions to draw random realizations of GPs.
- "fitPseudoInputsGP.R" contains functions to fit GP models with pseudo-inputs. The 
  five approaches described in the report are implemented with this code.
- "summaryFitPseudoInputsGP.R" contains functions to summarize the result of the GP fit
  (producing plots and mean prediction error)

- "simulations_effect_m_1.R" implements the simulation study investigating the effect of 
  m (Simulation 1 in the report) for the "smooth" function. 
- "simulations_effect_m_2.R" implements the simulation study investigating the effect of 
  m (Simulation 1 in the report) for the "wiggly" function.
- "simulation_effect_draws_1.R" implements the simulation study investigating the average 
  performance of the methods (Simulation 2 in the report) for the "smooth" function.
- "simulation_effect_draws_2.R" implements the simulation study investigating the average 
  performance of the methods (Simulation 2 in the report) for the "wiggly" function.

- "apply gp to hr data.R" applies the implemented code to the real HR.