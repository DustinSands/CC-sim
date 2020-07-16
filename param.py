# -*- coding: utf-8 -*-
"""
Holds global parameters for the simulation
"""
import numpy as np
import quantities as q
from quantities import Quantity as Q


#The simulation resolution
resolution = np.timedelta64(1, 'm')

# Safe mode; off checks that all units are compatible when operating on them
# on is much faster, but should not be used when implementing new features
skip_units = 0


# Same as simulation resolution, but in quantities package
q_intermediate = resolution / np.timedelta64(1, 'm')
q_res = Q(q_intermediate, 'min').simplified

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

cell_components = ['dO2', 'glucose', 'iron', 'amino_acids', 'dCO2', 'IGG_a',
                      'IGG_b', 'IGG_n', 'component_A', 'Na', 'K']
liquid_components = [*cell_components, 'LDH', 'lactate', 
                      'component_B', 'component_C', 'component_D', 'Cl',
                      ]
gas_components = ['CO2', 'O2', 'air']
positively_charged = ['Na', 'K']
negatively_charged = ['Cl']

molecular_weight = {
  'glucose':Q(180.156, 'g/mol'),
  'dO2':Q(32., 'g/mol'),
  'dCO2':Q(44.01, 'g/mol'),
  'acetic_acid':Q(60.052, 'g/mol'),
  'amino_acids':Q(137., 'g/mol'),
  'butyric_acid':Q(88., 'g/mol'),
  'citric_acid':Q(192.12, 'g/mol'),
  'formic_acid':Q(46.03, 'g/mol'),
  'isovaleric_acid':Q(102.13, 'g/mol'),
  'lactate':Q(90.08, 'g/mol'),
  'adenine':Q(135.13, 'g/mol'),
  'iron':Q(200., 'g/mol'),   #forms of iron are generally in the 150-250 g/mol range
  'LDH':Q(140., 'kg/mol'),
  'IGG_a':Q(150., 'kg/mol'),
  'IGG_b':Q(150., 'kg/mol'),
  'IGG_n':Q(150., 'kg/mol'),
  'IGG':Q(150, 'kg/mol'),
  'component_A':Q(100, 'kg/mol'),
  'component_B':Q(50, 'kg/mol'),
  'component_C':Q(500, 'g/mol'),
  'component_D':Q(200, 'g/mol'),
  'Na':Q(22.99,'g/mol'),
  'K':Q(39.1, 'g/mol'),
  'Cl':Q(35.453, 'g/mol'),
  }

gravity = Q(9.81, 'm/s**2')
actual_cc_density = Q(980, 'g/L') #Assumed constant.  Can be deleted when calculated density encoded
expected_cc_density = Q(1000, 'g/L')
viscosity = Q(0.001, 'Pa*s')
environment_temperature = 26
volumetric_heat_capacity = Q(4180, 'J/L') # per K 

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
  'component_A_production_rate':Q(1.e-17, 'mol/ce/day'),
  'extinction_coeff': 200.,
  'mass_transfer_rate': Q(1.1e-11, 'm/s'), # 1.1e-11 is 0.2/min for 20um cell
  'O2_consumption':Q(0.15,'M/min'),  # 5 g/min/L
  'glucose_consumption':Q(1.6667e-13,'mol/hour'),
  # 'glucose_consumption':Q(1.6667e-13,'mol/hour'), #3e-11 g/h per cell
  'aa_consumption':Q(1.,'M/day'), # 150 g/day/L
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
  'dO2':1,
  'dCO2':0.95,
  'dCO':1.03,
  }