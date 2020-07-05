# -*- coding: utf-8 -*-
"""
cell model.  Current purpose is to reflect roughly what a CHO cell producing
IGG would look like.  Greatly simplified model that does not cover many cases.
"""

# cell_size [um]
# VCD [e5c/ml]
# viability [%]
# mass_transfer: {all tracked components}

import random
import math

from quantities import Quantity as Q
import numpy as np

import tests, param, main

"""These are NOT intended to exactly emulate real cells.  Parameters are 
intended to roughly represent how cells might behave, as well as provide
a variance between randomly generated cell lines.

Tolerance represents the point at which 10% of the cells die a day."""
global_cell_line_mean_parameters = \
   {'growth_osmo': 320, 
    'growth_osmo_range':20,
    'osmo':350,
    'osmo_tolerance':60,
    'max_production_rate':Q(50, 'pg/ce/d'),
    'shear_tolerance' : Q(3500, '1/s'),
    'shear_sensitivity':Q(0.001, 's'),
    'shear_requirement':Q(300, '1/s'),
    'temperature':36,
    'temperature_tolerance':3,
    'dO2': 60*Q(6.72e-5, 'g/L'),
    'dO2_tolerance':30*Q(6.72e-5, 'g/L'),
    'pH': 7,
    'pH_tolerance':0.2,
    'metabolism_rate':1,
    'LDH_density':1,  #scale unknown
    'death_delay': 44, #hours
    'delayed_death_transition': 4, #hours
    'growth_diameter': Q(14, 'um'),
    'production_diameter_increase': Q(5, 'um'),
    'component_A_sensitivity':Q(1/10, 'L/mg'),
    'component_A_tolerance':Q(10, 'mg/L'),
    'component_B_sensitivity':Q(1/5, 'L/mg'),
    'component_B_tolerance':Q(20, 'mg/L'),
    'component_C_sensitivity':Q(1/15, 'L/mg'),
    'component_C_tolerance':Q(30, 'mg/L'),
    'component_D_sensitivity':Q(0, 'L/mg'),
    'component_D_tolerance':Q(40, 'mg/L'),
    'acidic_variants':0.24,
    'basic_variants':0.11,
    'respiratory_quotient':0.97,
    'growth_rate': Q(26, 'h'),
    }
# The standard deviation of each randomly generated variable
global_cell_line_STD_parameters = \
   {'growth_osmo': 15, 
    'growth_osmo_range':10,
    'osmo':15,
    'osmo_tolerance':20,
    'max_production_rate':Q(15, 'pg/ce/d'),
    'shear_tolerance':Q(500, '1/s'),
    'shear_sensitivity':Q(0.0005, 's'),
    'shear_requirement':Q(100, '1/s'),
    'temperature':1,
    'temperature_tolerance':1,
    'dO2': 5*Q(6.72e-5, 'g/L'),
    'dO2_tolerance':10*Q(6.72e-5, 'g/L'),
    'pH': 0.05,
    'pH_tolerance':0.05,
    'metabolism_rate':0.2,
    'LDH_density':0.2,
    'death_delay': 6,
    'delayed_death_transition': 1,
    'growth_diameter': Q(1, 'um'),
    'production_diameter_increase': Q(2, 'um'),
    'component_A_sensitivity':Q(1/5, 'L/mg'),
    'component_A_tolerance':Q(2, 'mg/L'),
    'component_B_sensitivity':Q(1/1, 'L/mg'),
    'component_B_tolerance':Q(2, 'mg/L'),
    'component_C_sensitivity':Q(1/5, 'L/mg'),
    'component_C_tolerance':Q(10, 'mg/L'),
    'component_D_sensitivity':Q(1/10, 'L/mg'),
    'component_D_tolerance':Q(5, 'mg/L'),
    'acidic_variants':0.08,
    'basic_variants':0.04,
    'respiratory_quotient':0.01,
    'growth_rate': Q(1, 'h'),
      }

def random_cell_line():
  cell_param = {}
  for variable in global_cell_line_mean_parameters:
    cell_param[variable] = random.gauss(global_cell_line_mean_parameters[variable],
                                        global_cell_line_STD_parameters[variable])
  
  return cell_param

def gen_cell_line():
  while True:
    cell_param = random_cell_line()
    # Critic; check that parameters are within reasonable levels
    param_wrong = 0
    param_wrong += cell_param['growth_osmo_range'] < 10
    param_wrong += cell_param['osmo_tolerance'] < 30
    param_wrong += cell_param['shear_tolerance'] < 300
    param_wrong += cell_param['temperature_tolerance'] < 1
    param_wrong += cell_param['dO2_tolerance'] < 15*Q(6.72e-5, 'g/L')
    param_wrong += cell_param['pH_tolerance'] < 0.1
    param_wrong += cell_param['delayed_death_transition'] < 2
    param_wrong += cell_param['acidic_variants'] < 0.06
    param_wrong += cell_param['basic_variants'] < 0.04
    param_wrong += cell_param['production_diameter_increase'] < 1
    
    if param_wrong==0:
      return cell_param
  
class cell_wrapper:
  p = param.cells
  def __init__(self, cell_line, seed_cells):
    self.viable_cells = seed_cells
    self.cp = cell_line
    self.delayed_death_buckets = Q([0.]*round(self.cp['death_delay']),'ce')
    self.hour_countdown = np.timedelta64(1, 'h')/param.resolution
    self.dying_cells = Q(0., 'ce')
    self.dead_cells = Q(0., 'ce')
    self.extinction_conversion = np.timedelta64(1, 'D') / param.resolution
    self.diameter = self.cp['growth_diameter']
    self.volume = math.pi/6*self.diameter**3
    self.growth_per_tick = (param.q_res/self.cp['growth_rate']).simplified*(
      self.volume)
    self.out_of_range = {}
    self.concentration = {component:Q(0., 'g/L') for component in param.cell_components}
    
  def calc_mass_transfer(self, env):
    """
    Resource intake
    Modelled as passive first-order single-layer with high diffusion value
    Future updates:
      active mass_transfer
      boundary layers
      """
    mass_transfer = {}
    for component in param.liquid_components:
      mass_transfer[component] = Q(0, 'g/min')
    for component in param.cell_components:
      mass_transfer[component] = math.pi*self.diameter**2*self.p['mass_transfer_rate']*\
        (env[component]-self.concentration[component]) #[g/min]
    return mass_transfer
  
  def calc_growth(self, env, metabolism_rate):
    """Calculate where the cell is expending its energy and what stage it is 
    in.  Currently based purely on present conditions, can include internal
    state in future updates.
    
    (Although it is unclear to what extent cells do this)
    
    """
    self.volume += self.growth_per_tick
      
    sensitivity = self.cp['component_A_sensitivity']
    tolerance = self.cp['component_A_tolerance']
    growth_inhibition_coeff= 1/(1+math.exp(-sensitivity*(env['component_A']-
                                                         tolerance-3/sensitivity)))
    growth_inhibition_coeff *= math.exp(-(env['osmo']-self.cp['growth_osmo'])**2/
                                 (2*self.cp['growth_osmo_range']**2))
    #Convert some biomass to VCD if target diameter is smaller than current
    #Otherwise, just grow diameter
    target_diameter = self.diameter + (1-growth_inhibition_coeff)*self.cp[
      'production_diameter_increase']
    target_volume = math.pi/6*target_diameter**3
    volume_diff = self.volume - target_volume
    if volume_diff > 0:
      if volume_diff > 2*self.growth_per_tick:
        volume_diff = 2*self.growth_per_tick
      self.volume -= volume_diff
      self.VCD *= volume_diff / self.volume + 1
        
    

  
  def calc_death(self, env, limiting_ratio, metabolism_rate):
    """
    Death Calculations
    Some cell death is delayed.  This is done by bucketizing the deaths by the 
    hour.  Once an hour we start killing the cells in the last bucket and
    create a new one.
    """
    self.hour_countdown -= 1
    if self.hour_countdown <= 0:
      self.hour_countdown += np.timedelta64(1, 'h') / param.resolution
      self.dying_cells += self.delayed_death_buckets[-1]
      self.delayed_death_buckets = [Q(0., 'ce')]+self.delayed_death_buckets[:-1]
    extinction_rate = 0
    # If environmental variables are out of bounds, the cells start to die
    for condition in ['osmo', 'temperature', 'dO2', 'pH']:
      tolerance = condition + '_tolerance'
      # diff = max(0, abs(env[condition]-self.cp[condition])-self.cp[tolerance])
      diff = env[condition]-self.cp[condition]
      extinction_rate += (math.exp((diff/self.cp[tolerance])**2)-1)/\
        self.p['extinction_coeff']
    #Convert extinction rate from per day to per step
    extinction_rate /= self.extinction_conversion
    #If death rate is low enough, the death is delayed.  Using a logistic
    #curve to calculate ratio
    death_ratio = 1/(1+math.exp(-20*(extinction_rate - 0.3)))
    self.delayed_death_buckets[0] += (1-death_ratio)*extinction_rate*self.viable_cells
    self.dead_cells += extinction_rate*self.viable_cells

  def metabolism(self, env, mass_transfer):
    """Calculate usage of resources by first finding limiting reactant, then
    removing resources used and creating biomass / product.
    
    Future update: include iron for biomass only."""
    metabolism_rate = self.cp['metabolism_rate']*0.19*math.exp(0.01458*env['temperature'])

    aa_vol_consumption_rate = self.p['aa_consumption']*metabolism_rate
    O2_vol_consumption_rate = self.p['O2_consumption']*metabolism_rate
    glucose_consumption_rate = self.p['glucose_consumption']*metabolism_rate*Q(1, 'ce')
    glucose_mass = self.concentration['glucose']*self.volume
    limiting_ratio = min(glucose_mass/(mass_transfer['glucose']-glucose_consumption_rate),
          self.concentration['amino_acids']/(mass_transfer['amino_acids']/self.volume-aa_vol_consumption_rate),
          self.concentration['dO2']/(mass_transfer['dO2']/self.volume-O2_vol_consumption_rate),
          param.q_res)
    metabolic_scalar = metabolism_rate * limiting_ratio 
    self.concentration['glucose'] -= metabolic_scalar*glucose_consumption_rate/self.volume
    self.concentration['dO2'] -= metabolic_scalar*O2_vol_consumption_rate
    self.concentration['CO2'] += metabolic_scalar*O2_vol_consumption_rate
    self.concentration['amino_acids'] -= metabolic_scalar*aa_vol_consumption_rate
    
    self.calc_death(env, limiting_ratio, metabolism_rate)
    self.calc_growth(env, metabolism_rate)
    
  def update_in_range(self, env):
    """
    Calculate roughly how close vars are in being in range
    0 = in range
    1 = out of range
    """
    
    for variable in ['component_A', 'component_B', 'component_C', 'component_D', 
                     'shear']:
      sensitivity = self.cp[variable+'_sensitivity']
      tolerance = self.cp[variable+'_tolerance']
      self.out_of_range[variable] = 1/(1+math.exp((
        sensitivity*(env[variable]-1/sensitivity*3-tolerance)).simplified))
    for variable in ['pH', 'dO2', 'osmo', 'temperature']:
      tolerance = variable+'_tolerance'
      self.out_of_range[variable] = \
        math.exp((env[variable]-self.cp[variable])**2/self.cp[tolerance]**2)
    
  def step(self, env):
    cells = {}
    self.update_in_range(env)

    """
    Calculate mass_transfer first, but do not change concentration.  Then calculate 
    metabolism and perform both intake and consumption back-to-back.
    """
    mass_transfer = self.calc_mass_transfer(env)
    self.metabolism(env, mass_transfer)
    
    living_cells = self.delayed_death_buckets.sum()+self.dying_cells + self.viable_cells
    cells['viability'] = living_cells / (living_cells+self.dead_cells)*100
    cells['living_cells'] = living_cells
    cells['total_cells'] = living_cells+self.dead_cells
    cells['diameter'] = self.diameter
    cells['volume'] = self.volume
    
    """Calculate PQ.  Mostly intended to create deviations in PQ if any of the 
    variables known to affect PQ go out of range for the cell line.  Might be 
    worth making the affect condependent on multiple variables if aiming to 
    optimize PQ in the future (and make it more difficult so things like DOE 
    are useful?)"""
    IGG = {}
    acidic_fraction = self.out_of_range['component_C']*-0.1+self.out_of_range['shear']*0.1+\
      self.out_of_range['pH']*0.4
    basic_fraction = self.out_of_range['component_B']*0.25+self.out_of_range['dO2']*0.5+\
      self.out_of_range['temperature']*-0.1
    
    #shock - not implemented
    cells['mass_transfer'] = {component:value*self.viable_cells for 
                              component, value in mass_transfer.items()}
    cells['mass_transfer']['component_A'] = -self.viable_cells * self.p['component_A_production_rate']

    
    return cells
  
if __name__ == '__main__':
  main.run_sim()
  tests.cell_tests()