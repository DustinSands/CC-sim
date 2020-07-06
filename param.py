# -*- coding: utf-8 -*-
"""
Holds global parameters for the simulation
"""
import numpy as np
import quantities as q
from quantities import Quantity as Q


#The simulation resolution
resolution = np.timedelta64(1, 'm')


# Same as simulation resolution, but in quantities package
q_intermediate = resolution / np.timedelta64(1, 's')
q_res = Q(q_intermediate, 's').simplified

cells = q.unitquantity.IrreducibleUnit('cells', symbol = 'ce')
pg = q.UnitQuantity('pg', 1e-15*q.kg, 'picogram')
ce = q.UnitQuantity('ce', cells, 'ce')
cm = q.UnitQuantity('e5c', 1e5*cells, 'e5c')
hcm = q.UnitQuantity('e6c', 1e6*cells, 'e6c')
q.CD = q.CompoundUnit('e5c/ml')
q.HCD = q.CompoundUnit('e6c/ml')
# cell_density = q.UnitQuantity('cell density', cells/q.mL, symbol = 'e5c/ml')
# cells = q.unitquantity.IrreducibleUnit('cell', symbol = 'cell')
# cells = q.UnitQuantity('e5 cells', cells*1e5, symbol='e5c')

liquid_components = ['dCO2', 'dO2', 'glucose', 'H2CO3', 
                      'iron', 'LDH', 'lactate', 'amino_acids', 'component_A',
                      'component_B', 'component_C', 'component_D', 'IGG_a',
                      'IGG_b', 'IGG_n']
gas_components = ['CO2', 'O2', 'air']
cell_components = ['dO2', 'glucose', 'iron', 'amino_acids', 'dCO2', 'IGG_a',
                      'IGG_b', 'IGG_n']

gravity = Q(9.81, 'm/s**2')
actual_cc_density = Q(980, 'g/L') #Assumed constant.  Can be deleted when calculated density encoded
expected_cc_density = Q(1000, 'g/L')
viscosity = Q(0.001, 'Pa*s')

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
  'bioHT':{
    'random_error_CV':0.04,     #LC
    'systematic_error_CV':0.03, #LC
    },
  'Cell_Counter':{
    'density_random_error_CV':0.05,      # s
    'density_systematic_error_CV':0.03,  # s
    'viability_random_error_sigma':1.5, #net, s
    'viability_systematic_error_sigma':2, #net, s
    'size_random_error_sigma':0.1,
    'size_systematic_error_sigma':0.2,
    },
  'Scale':{
    'random_error_sigma':0.02, #grams
    }
  }

cells = {
  'component_A_production_rate':(Q(1, 'pg/ce/day')*q_res).simplified,
  'extinction_coeff': 200,
  'mass_transfer_rate': Q(1e-15, 'L/m**2/min'), # 8e-7 m**3 min, corresponds to 2/min
  'O2_consumption':Q(5,'g/min/L'),
  'glucose_consumption':Q(3e-11,'g/ce/hour'),
  'aa_consumption':Q(5,'g/day/L'),
  }

actuation = {
  'MFC':{
    'systematic_error_CV': 0.005,
    'break_chance': Q(1, '1/yr'),
    },
  'peristaltic':{
    'systematic_error_CV': 0.03,
    'break_chance': Q(1, '1/yr'),
    },
  'agitator':{
    'systematic_error_CV': 0.01,
    'break_chance': Q(0, '1/yr'),
    },
  'mixtures':{
    'component_CV': 0.0005}
  }

# Constants
kla_ratio = {
  'O2':1,
  'CO2':0.95,
  'CO':1.03,
  }