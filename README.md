# Cell-Culture
Allows for the generation of cell-culture data.

Whenever possible, data is generated from core mechanics so that an arbitrary
change to experimental design will create data that reflects this.  E.G. a 
change to feeding strategy or agitation rate will accurately be reflected in
things like VCD.

This tool could, in theory, be used to:
  Look at how important accuracy and precision in assays are
  Generate DoE data and assess accuracy of results
  Train and test ML algorithms*
  Teach students what real data would look like and how to derive parameters 
    from it.  
  
*obviously, don't use this as training data to predict on real world data.

Currently implemented:
  Various assays, all with their own forms of error (systematic, drift, random)
    and calibration schemes
  Modular control strategies and actuation devices
  pH calc from bicarb buffer equations (no weak acids / bases yet) 
  Customizable bioreactor setup**
  Macroscopic cell model
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
    
    
A note on the cell model: This code uses a macroscopic cell simulation to approximate
how a cell could react in a given situation.  


Best practices:
  Use units.  E.G. quantities.Quantity(1, 'kg') for one kilogram
    Don't include units for temperature (relational units are not supported)
  It is the class' responsibility to convert to SI units / float depending on mode
  
To run:
  main.py is configured to generate a random cell line, and run 2x2L scale 
  fed-batch reactors for 14 days and plot the results.  Use it to get an idea
  of how experiments can be configured, run, and finally have the results 
  visualized.
  
  
**kLa, shear equations will be most accurate at 2L scale


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