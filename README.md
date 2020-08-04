# CellCulture-Simulator (CC-Sim)
CC-Sim (CellCulture-Simulator) aims to digitally model the cell culture process, 
combining a cell model, control strategies, bioreactors, actuation devices, 
assays, and media definitions in order to simulate the defined platform and 
increase process understanding. 

For more information, visit: https://portfolio.dustinsands.com/cc-sim/

Currently implemented:
  Various assays, all with their own forms of error (systematic, drift, random)
    and calibration schemes
  Modular control strategies and actuation devices
  pH calc from bicarb buffer equations (no weak acids / bases yet) 
  Customizable bioreactor setup
  Macroscopic CHO cell model
    Randomly generated cell lines with 30+ varying parameters
    Metabolism with limiting reactants
    Includes PQ data (correlated (incorrectly) to environmental conditions)
    Delayed death
    Size changes depending on phase of growth / environmental conditions
    Glucose, O2, and AA (lumped) consumption reflect real data
  Units to safeguard whether implementation is correct
  Compartamentalized simulation and observation data to keep controls from using
    real instead of observed data
  Randomized sample times, mixture errors, assay errors, cell parameters, etc.
  Visualization / Plotting Functions

  
To run:
  main.py is configured to generate a random cell line, and run a 2L scale 
  fed-batch reactor for 14 days and plot the results.  Use it to get an idea
  of how experiments can be configured, run, and finally have the results 
  visualized.
  
  For an example run: 
  https://portfolio.dustinsands.com/cc-sim/cc-sim-example-run/
  
Best practices:
  Use units.  E.G. quantities.Quantity(1, 'kg') for one kilogram
    Don't include units for temperature (relational units are not supported)
  It is the class' responsibility to convert to SI units / float depending on mode

Sources: https://portfolio.dustinsands.com/cc-sim/appendix/
