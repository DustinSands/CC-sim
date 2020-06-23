# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:14:22 2020

@author: Racehorse
"""

# cell_size [um]
# VCD [e5c/ml]
# viability [%]
# uptake: {all tracked components}

import random

from quantities import Quantity as Q

import tests, param

"""These are NOT intended to exactly emulate real cells.  Parameters are 
intended to roughly represent how cells might behave, as well as provide
a variance between randomly generated cell lines."""
global_cell_line_mean_parameters = \
   {'growth_osmo': 320, 
    'growth_osmo_tolerance':20,
    'production_osmo':380,
    'production_osmo_tolerance':40,
    'viable_osmo':350,
    'viable_osmo_tolerance':60,
    'mean_production_rate':Q(0.1, 'g/e6c/d'),
    'titer_variance':,
    'shear_tolerance':,
    'shear_requirement':,
    'temperature':36,
    'temperature_tolerance':3,
    'min_glucose':,
    'dO2': 60,
    'dO2_tolerance':,
    'pH': 7,
    'pH_tolerance':,
    'metabolism_efficiency':1,
    'LDH_density':
    
    }
# The standard deviation of each randomly generated variable
global_cell_line_STD_parameters = \
   {'growth_osmo': 15, 
    'growth_osmo_tolerance':10,
    'production_osmo':20,
    'production_osmo_tolerance':15,
    'viable_osmo':15,
    'viable_osmo_tolerance':20,
    'mean_production_rate':,
    'titer_variance':,
    'shear_tolerance':,
    'shear_requirement':,
    'temperature':36,
    'temperature_tolerance':3,
    'min_glucose':,
    'dO2': 60,
    'dO2_tolerance':,
    'pH': 7,
    'pH_tolerance':,
    'metabolism_efficiency':1,
    'LDH_density':,
      }

def gen_cell_line(global_cell_params):
  cell_param = {}
  for variable in global_cell_line_mean_parameters:
    cell_param[variable] = random.gauss(global_cell_line_mean_parameters[variable],
                                        global_cell_line_STD_parameters[variable])
  
  return cell_param
  
class cell_wrapper:
  def __init__(self, cell_line, seed_density):
    self.VCD = seed_density
    self.cell_param = cell_line
    
  def step(self, env):
    cells = {}
    cells['uptake'] = {}
    for component in param.tracked_components:
      cells['uptake'][component] = Q(0, 'g/L/min')
    
    #Growth calculations
    growth_preference = math.exp(-(env['osmo']-self.cp['growth_osmo'])**2/
                                 (2*self.cp['growth_osmo_tolerance']**2))
    production_preference = math.exp(-(env['osmo']-self.cp['production_osmo'])**2/
                                 (2*self.cp['production_osmo_tolerance']**2))
    stagnate_weight = 1
    total_weight = growth_preference + production_preference + stagnate_weight
    growth_energy_ratio = growth_preference / total_weight
    production_energy_ratio = production_preference / total_weight
    
    
    return cells
  
if __name__ == '__main__':
  tests.cell_tests()