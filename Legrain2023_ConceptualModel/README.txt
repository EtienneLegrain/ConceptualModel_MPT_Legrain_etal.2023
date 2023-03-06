# Conceptual model used from the study : Legrain et al., (2023). A gradual mechanism is more likely to have caused the Mid-Pleistocene Transition than an abrupt event.

Please contact Etienne Legrain at etienne.legrain@univ-grenoble-alpes.fr for questions about running or development.

The conceptual model is coded in Python 3.9

Each simulation mentionned in the manuscript is make up of two .py files : 
XXX_simulation.py compute the best parameters vector using the emcee python hammer.
XXX_simulation_plot.py plot three graphs using the best parameters vector determined with the XXX_simulation.py file. This vector should be past line 17 of the code.

The Orbital_Parameters_2Ma_1ka, Orbital_Parameters_2Ma_1ka_synthetic, Sea-level_Berends_2Ma_1ka files are the input data used for each simulations.

The BestParam-simulations.txt file contains the best parameter vector for each simulations.
The BestRuns_output.csv file contains the modelled ice volume from the best parameter vector of each simulations.
