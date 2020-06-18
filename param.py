# -*- coding: utf-8 -*-
"""
Holds global parameters for the simulation
"""
import numpy as np
import quantities as q


#The simulation resolution
resolution = np.timedelta64(1, 'm')


# Same as simulation resolution, but in quantities package
q_intermediate = resolution / np.timedelta64(1, 's')
q_res = q.Quantity(q_intermediate, 's').simplified

cells = q.unitquantity.IrreducibleUnit('cells', symbol = 'ce')
cm = q.UnitQuantity('e5c', 1e5*cells, 'e5c')
hcm = q.UnitQuantity('e6c', 1e6*cells, 'e6c')
q.CD = q.CompoundUnit('e5c/ml')
q.HCD = q.CompoundUnit('e6c/ml')
# cell_density = q.UnitQuantity('cell density', cells/q.mL, symbol = 'e5c/ml')
# cells = q.unitquantity.IrreducibleUnit('cell', symbol = 'cell')
# cells = q.UnitQuantity('e5 cells', cells*1e5, symbol='e5c')

tracked_components = ['CO2', 'O2', 'Glucose', 'H2CO3']
gas_components = ['CO2', 'O2', 'air']

# ERROR CONSTANTS
instrumentation={
  'BGA':{
    # Slope
    'O2_random_error_CV' : 0.003,
    'O2_systematic_error_CV' : 0.005,    # Low confidence
    'CO2_random_error_CV' : 0.005,       # Low confidence
    'CO2_systematic_error_CV' : 0.01,    # Low confidence
    # Net
    'pH_random_error_sigma' : 0.005,
    'pH_systematic_error_sigma' : 0.002, # Low confidence
    },
  'Probe':{
    'O2':{  
      'drift_CV' : 0.1,  #s
      't98' : np.timedelta64(60, 's'), #s
      'random_CV' : 0.001,     # relative
      },
    'pH':{
      'drift_sigma' : 1.73, #s
      't98' : np.timedelta64(30, 's'), #s
      'random_sigma' : 0.001,     # 'net'
      },
    'temperature':{  
      't98' : np.timedelta64(30, 's'),
      'random_sigma' : 0.01,      # net
      'systematic_sigma': 0.05    # net
      },
    },
  'Cell_Counter':{
    'density_random_error_CV':0.05,      # s
    'density_systematic_error_CV':0.03,  # s
    'viability_random_error_sigma':1.5, #net, s
    'viability_systematic_error_sigma':2, #net, s
    'size_random_error_sigma':0.1,
    'size_systematic_error_sigma':0.2,
    }
  }

# Constants
kla_ratio = {
  'O2':1,
  'CO2':0.95,
  'CO':1.03,
  }