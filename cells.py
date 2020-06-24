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
a variance between randomly generated cell lines.

Tolerance represents the point at which 10% of the cells die a day."""
global_cell_line_mean_parameters = \
   {'growth_osmo': 320, 
    'growth_osmo_range':20,
    'production_osmo':380,
    'production_osmo_range':40,
    'osmo':350,
    'osmo_tolerance':60,
    'max_production_rate':Q(50, 'pg/ce/d'),
    'max_shear' : Q(3500, '1/s'),
    'shear_tolerance':Q(1000, '1/s'),
    'shear_requirement':Q(300, '1/s'),
    'temperature':36,
    'temperature_tolerance':2,
    'dO2': 60,
    'dO2_tolerance':30,
    'pH': 7,
    'pH_tolerance':0.2,
    'metabolism_efficiency':1,
    'metabolism_rate':1,
    'LDH_density':1,  #scale unknown
    'death_delay': 44, #hours
    'delayed_death_transition': 4, #hours
    'growth_diameter': Q(12, 'um'),
    'production_diameter_increase': Q(2, 'um'),
    
    }
# The standard deviation of each randomly generated variable
global_cell_line_STD_parameters = \
   {'growth_osmo': 15, 
    'growth_osmo_range':10,
    'production_osmo':20,
    'production_osmo_range':15,
    'osmo':15,
    'osmo_tolerance':20,
    'max_production_rate':Q(15, 'pg/ce/d'),
    'max_shear':Q(500, '1/s'),
    'shear_tolerance':Q(500, '1/s'),
    'shear_requirement':Q(100, '1/s'),
    'temperature':1,
    'temperature_tolerance':1,
    'dO2': 5,
    'dO2_tolerance':10,
    'pH': 0.05,
    'pH_tolerance':0.05,
    'metabolism_efficiency':0.2,
    'metabolism_rate':0.2,
    'LDH_density':0.2,
    'death_delay': 6,
    'delayed_death_transition': 1,
    'growth_diameter': Q(0.3, 'um'),
    'production_diameter': Q(0.4, 'um'),
      }

def random_cell_line():
  cell_param = {}
  for variable in global_cell_line_mean_parameters:
    cell_param[variable] = random.gauss(global_cell_line_mean_parameters[variable],
                                        global_cell_line_STD_parameters[variable])
  
  return cell_param

def gen_cell_line():
  while true:
    cell_param = random_cell_line()
    # Critic; check that parameters are within reasonable levels
    param_wrong = 0
    param_wrong += cell_param['growth_osmo_range'] < 10
    param_wrong += cell_param['production_osmo_range'] < 15
    param_wrong += cell_param['osmo_tolerance'] < 30
    param_wrong += cell_param['shear_tolerance'] < 300
    param_wrong += cell_param['temperature_tolerance'] < 0.5
    param_wrong += cell_param['dO2_tolerance'] < 15
    param_wrong += cell_param['pH_tolerance'] < 0.1
    param_wrong += cell_param['delayed_death_transition'] < 2
    
    if param_wrong==0:
      return cell_param
  
class cell_wrapper:
  def __init__(self, cell_line, seed_density):
    self.VCD = seed_density
    self.cp = cell_line
    self.delayed_death_buckets = [0]*round(self.cp['death_delay'])
    self.hour_countdown = np.timedelta64(1, 'h')/param.resolution
    self.dying_cells = Q(0, 'e6c/ml')
    self.dead_cells = Q(0, 'e6c/ml')
    
    
  def step(self, env):
    cells = {}
    cells['uptake'] = {}
    for component in param.tracked_components:
      cells['uptake'][component] = Q(0, 'g/L/min')

      
    #Growth calculations
    growth_preference = math.exp(-(env['osmo']-self.cp['growth_osmo'])**2/
                                 (2*self.cp['growth_osmo_range']**2))
    production_preference = math.exp(-(env['osmo']-self.cp['production_osmo'])**2/
                                 (2*self.cp['production_osmo_range']**2))
    stagnate_weight = 1
    total_weight = growth_preference + production_preference + stagnate_weight
    growth_energy_ratio = growth_preference / total_weight
    production_energy_ratio = production_preference / total_weight
    
    
    #Death Calculations
    
    #Once every hour
    self.hour_countdown -= 1
    if self.hour_countdown <= 0:
      self.hour_countdown += np.timedelta64(1, 'h') / param.resolution
      self.dying_cells += self.delayed_death_buckets[-1]
      self.delayed_death_buckets = [0]+self.delayed_death_buckets[:-1]
    extinction_rate = 0
    # If environmental variables are out of bounds, the cells start to die
    for condition in ['osmo', 'temperature', 'd02', 'pH']:
      tolerance = condition + '_tolerance'
      # diff = max(0, abs(env[condition]-self.cp[condition])-self.cp[tolerance])
      diff = env[condition]-self.cp[condition]
      extinction_rate += 1/(200* math.exp(-(diff/self.cp[tolerance])**2))-0.005
    #Convert extinction rate from per day to per step
    extinction_rate /= np.timedelta64(1, 'D') / param.resolution
    #If death rate is low enough, the death is delayed.  Using a logistic
    #curve to calculate ratio
    death_ratio = 1/(1+math.exp(-20*(extinction_rate - 0.3)))
    self.delayed_death_buckets[0] += (1-death_ratio)*extinction_rate*self.VCD
    self.dead_cells += extinction_rate*self.VCD
    
    
    
    VCD = sum(self.delayed_death_buckets)+self.dying_cells + self.VCD
    cells['viability'] = VCD / (VCD+self.dead_cells)*100
    cells['VCD'] = VCD
    

    #shock - not implemented
    

    
    return cells
  
if __name__ == '__main__':
  tests.cell_tests()