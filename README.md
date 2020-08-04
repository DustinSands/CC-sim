# CellCulture-Simulator (CC-Sim)
An upstream bioprocess simulator that comes with a CHO-cell model.
Allows for the generation of cell-culture data.

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


Best practices:
  Use units.  E.G. quantities.Quantity(1, 'kg') for one kilogram
    Don't include units for temperature (relational units are not supported)
  It is the class' responsibility to convert to SI units / float depending on mode
  
To run:
  main.py is configured to generate a random cell line, and run a 2L scale 
  fed-batch reactor for 14 days and plot the results.  Use it to get an idea
  of how experiments can be configured, run, and finally have the results 
  visualized.


CHO Cell stages, size increase, composition
  Pan X, Dalm C, Wijffels RH, Martens DE. Metabolic characterization of a CHO 
  cell size increase phase in fed-batch cultures. Appl Microbiol Biotechnol. 
  2017;101(22):8101-8113. doi:10.1007/s00253-017-8531-y

Shear equations
  R. Bowen, Unraveling the mysteries of shear-sensitive mixing systems,
  Chem. Eng. 9 (June) (1986) 55–63.
      
kLa equations    
  Liu, K.; Phillips, J.R.; Sun, X.; Mohammad, S.; Huhnke, R.L.; Atiyeh, H.K. 
  Investigation and Modeling of Gas-Liquid Mass Transfer in a Sparged and 
  Non-Sparged Continuous Stirred Tank Reactor with Potential Application in 
  Syngas Fermentation. Fermentation 2019, 5, 75.
  
CO2 Equilibria:
https://folk.ntnu.no/skoge/prost/proceedings/dycops2007-and-cab2007/CAB/Monday/Estimation/M32_89_Final_ManuscriptStampedno.pdf

OUR and CO2UR
https://aiche.onlinelibrary.wiley.com/doi/abs/10.1002/btpr.646

Chemically Defined Media:
  Pan X, Streefland M, Dalm C, Wijffels RH, Martens DE. Selection of chemically 
  defined media for CHO cell fed-batch culture processes. Cytotechnology. 
  2017;69(1):39-56. doi:10.1007/s10616-016-0036-5