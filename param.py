# -*- coding: utf-8 -*-
"""
Holds global parameters for the simulation
"""
import pdb

import numpy as np
import quantities as q
from quantities import Quantity as Q



#The simulation resolution
resolution = np.timedelta64(1, 'm')

#Restrict observed quantities to real assays (no cheating)
realistic_mode = 0

# Safe mode; off checks that all units are compatible when operating on them
# on is much faster, but should not be used when implementing new features
skip_units = 1


# Allows things like peristaltic pumps to break (generally 1 / yr)
allow_breaking = 0


# Same as simulation resolution, but in quantities package
q_intermediate = resolution / np.timedelta64(1, 'm')
q_res = Q(q_intermediate, 'min').simplified
if skip_units:
  step_size = float(q_res)
else: step_size = q_res

cells = q.unitquantity.IrreducibleUnit('cells', symbol = 'ce')
pg = q.UnitQuantity('pg', 1e-15*q.kg, 'picogram')
ce = q.UnitQuantity('ce', cells, 'ce')
cm = q.UnitQuantity('e5c', 1e5*cells, 'e5c')
hcm = q.UnitQuantity('e6c', 1e6*cells, 'e6c')
q.CD = q.CompoundUnit('e5c/ml')
q.HCD = q.CompoundUnit('e6c/ml')

#Constants
gravity = Q(9.81, 'm/s**2').simplified
actual_cc_density = Q(980, 'g/L').simplified #Assumed constant.  Can be deleted when calculated density encoded
expected_cc_density = Q(1000, 'g/L').simplified
viscosity = Q(0.001, 'Pa*s').simplified
environment_temperature = 26
volumetric_heat_capacity = Q(4180, 'J/L').simplified # per K 
# if skip_units:
#   gravity = float(gravity)
#   actual_cc_density =float(actual_cc_density)
#   expected_cc_density = float(expected_cc_density)
#   viscosity = float(viscosity) 
#   environment_temperature = float(environment_temperature)
#   volumetric_heat_capacity = float(volumetric_heat_capacity)



cell_components = ['dO2', 'glucose', 'iron', 'amino_acids', 'dCO2', 'IGG_a',
                      'IGG_b', 'IGG_n', 'component_A', 'Na', 'K']
liquid_components = [*cell_components, 'LDH', 'lactate', 
                      'component_B', 'component_C', 'component_D', 'Cl',
                      ]
gas_components = ['CO2', 'O2', 'air']
positively_charged = ['Na', 'K']
negatively_charged = ['Cl', 'lactate']

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
    'viability_random_error_sigma':0.015, #net, s
    'viability_systematic_error_sigma':0.02, #net, s
    'size_random_error_sigma':Q(0.1, 'um'),
    'size_systematic_error_sigma':Q(0.2, 'um'),
    },
  'Scale':{
    'random_error_sigma':Q(0.02, 'g'), #grams
    },
  'Osmo':{
    'random_sigma':2,
    'systematic_sigma':1,
    },
  }

cells = {
  'component_A_production_rate':Q(4e-16, 'mol/ce/day'),
  'extinction_coeff': 200.,
  # 'mass_transfer_rate': Q(1.8e-5, 'm/s'), # "kLa" of 12/min for cells
  'dO2_consumption':Q(0.050,'M/min'),  # 5 g/min/L
  'glucose_consumption':Q(1.6667e-13,'mol/hour/ce'),
  'aa_consumption':Q(2.e-12,'mol/day/ce'), # 150 g/day/L
  'aa_limiting_concentration':Q(0.5, 'mM'),
  'O2_limiting_concentration':Q(0.001, 'mM'),
  'glucose_limiting_concentration':Q(0.1, 'mM'),
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
"""CO2 kLa is actually 0.95 of O2 kLa.  However, there occurs significant 
depletation / saturation of CO2 in gaseous phase that causes C* to change.  
Temporary fix is dividing kLa by 10."""
kLa_ratio = {
  'dO2':1,
  'dCO2':0.95/15,
  'dCO':1.03,
  }


for variable in [molecular_weight, cells, actuation, instrumentation]:
  import helper_functions
  helper_functions.simplify(variable)
  if skip_units:
    helper_functions.remove_units(variable)
