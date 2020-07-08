# -*- coding: utf-8 -*-
"""
Models properties of the intercellular fluid in bioreactor.

Takes control device outputs (like oxygen flow rates, feed flow rates, etc.)
and outputs the conditions within the bioreactor.

Currently, modeled as a completely homogenous solution.  At some point I'd like
to leave room for heterogeneity.

"""
import math
import pdb

from quantities import Quantity as Q

import param, tests, main

"""Defines solubility equation for components that vary with temperature."""
solubility = {
  'dO2': lambda t: Q(float(0.0005 - t*0.00000475), 'g/L/kPa'),
  'dCO2': lambda t: Q(float(0.021 - t*0.000275), 'g/L/kPa')
  }


class hollow_fiber:
  """Perfusion hollow fiber."""
  
class agitator:
  """Holds impeller parameters.
  """
  def __init__(self, impeller_type = 'rushton', 
               number = 1, 
               diameter = Q(6, 'cm'),
               width = Q(2, 'cm')
               ):
    self.type = impeller_type
    self.number = number
    self.diameter = diameter
    self.width = width
    if impeller_type == 'rushton':
      self.power_number = 5.5
    elif impeller_type == 'marine':
      self.power_number = 2.2
    else: raise ValueError('impeller not recognized!')

    
    self.ungassed_power_coeff = (self.power_number*Q(1000, 'kg/m**3')*self.diameter**5).simplified
    
  def ungassed_power(self, RPS):
    return RPS**3*self.ungassed_power_coeff
  

class bioreactor:
  """Set up bioreactor configuation, then perform step to update environment
  in the bioreactor.
  
  Assumed to be a perfectly cylindrical vessel for volume calculations."""
  def __init__(self, start_time,
                 initial_components,
                 seed_volume,
                 initial_temperature,
                 volume = Q(3, 'L'),
                 agitator = agitator(),
                 diameter = Q(13, 'cm'),     
                 sparger_height = Q(2, 'cm'), # Height from bottom of vessel
                 sparger_pore_size = Q(20, 'um'),
                 # num_pores = 10000,
                 cell_separation_device = None,     # Perfusion only
                 head_pressure = Q(760, 'mmHg'),
                 heat_transfer_coeff = Q(50, 'W/m**2'), #per degree (not supported)
                 ):
    self.volume = volume
    self.working_volume = seed_volume
    self.agitator = agitator
    self.CSA = diameter**2/4*math.pi
    self.diameter = diameter
    self.current_time = start_time
    self.sparger_pore_size = sparger_pore_size
    # self.sparger_num_pores = num_pores
    self.kla_func = self.create_kla_function()
    self.old = {}
    self.pressure = (sparger_height/2*param.actual_cc_density*param.gravity+head_pressure).simplified
    self.mass = {component:Q(0., 'g') for component in param.liquid_components}
    self.mass.update(initial_components)
    self.sparger_height = sparger_height
    self.temperature = initial_temperature
    self.overall_heat_transfer_coeff = heat_transfer_coeff*\
      (2*self.CSA+self.volume/self.CSA*math.pi*self.diameter)   #external area


    
  def update_shear(self, RPS):
    """Updates the shear rates with new RPS.
    
    Equations from:
      R. Bowen, Unraveling the mysteries of shear-sensitive mixing systems,
      Chem. Eng. 9 (June) (1986) 55–63.
      """
    ratio = (self.agitator.diameter / self.diameter)**0.3
    self.mean_shear = Q((4.2*RPS*ratio*self.agitator.diameter /\
                         self.agitator.width).simplified, '1/s')
    self.max_shear = self.mean_shear / 4.2 * 9.7
  
  def check_and_update(self, actuation):
    self.current_time += param.resolution
    if not (actuation == self.old and actuation['liquid_volume'] == Q(0, 'L/min')):
      self.old = actuation
      # Is this in per hour or per minute?
      self.kla = self.kla_func(actuation['RPS'], actuation['gas_volume'], self.working_volume)
      self.gas_percentages = self.calc_gas_percentages(actuation)
      self.update_shear(actuation['RPS'])
      
  def calc_gas_percentages(self,actuation):
    total = Q(0., 'L/min')
    percent = {}
    for component in param.gas_components:
      total += actuation[component]
    for component in param.gas_components:
      percent['d'+component] = actuation[component]/total
    # Assume N2 is irrelevant
    percent['dO2'] += percent['dair']*0.21
    return percent
    
  def step(self, actuation, cells):
    """Calculate environmental changes."""
    self.check_and_update(actuation)
    cell_fraction = (cells['total_cells']*cells['volume']/Q(1, 'ce')/self.working_volume).simplified

    for component in param.liquid_components:
      if component[1:] in param.gas_components:
        C_star = self.gas_percentages[component]*self.pressure*solubility[component](self.temperature)
        
        transfer_rate = (self.kla*param.kla_ratio[component])*\
          (C_star - self.mass[component]/(self.working_volume*(1-cell_fraction)))
        transfer_rate = (transfer_rate*self.working_volume).simplified
      else: 
        transfer_rate = actuation[component]
      # print('BR', component, transfer_rate, self.mass[component])
      self.mass[component] += (transfer_rate - cells['mass_transfer'][component]) * param.q_res
    self.working_volume += actuation['liquid_volume']*param.q_res
    
    
    #Temperature
    heat_transfer = actuation['heat'] + \
      (param.environment_temperature-self.temperature)*\
        self.overall_heat_transfer_coeff
    self.temperature += (heat_transfer*param.q_res /\
      (self.working_volume*param.volumetric_heat_capacity)).simplified
    
    osmo = (self.total_moles()/Q(1, 'mmol/L')/\
            (self.working_volume*(1-cell_fraction))).simplified
    
    #Build output
    environment={'shear':self.mean_shear, 
                 'max_shear':self.max_shear,
                 'volume':self.working_volume,
                 'mOsm':osmo,
                 'pH': self.pH(),
                 'temperature':self.temperature,
                 'time': self.current_time
                        }
    for component in param.liquid_components:
      environment.update({component:self.mass[component]/self.working_volume})
      
    return environment
    
  def total_moles(self):
    total = Q(0., 'mmol')
    for component in self.mass:
      total += self.mass[component]/param.molecular_weight[component]
      if component in ['NaCL', 'KCl']:  #Counts twice
        total += self.mass[component]/param.molecular_weight[component]
    return total
  
  def pH(self):
    """Calculates pH of system.  Needs future update."""
    return 7.3 - self.mass['dCO2']/Q(1, 'g/L')/self.working_volume
  
  def create_kla_function(self):
    """Calculate and return a function for oxygen transfer rate.
    
    Equations used based on: 
    Liu, K.; Phillips, J.R.; Sun, X.; Mohammad, S.; Huhnke, R.L.; Atiyeh, H.K. 
    Investigation and Modeling of Gas-Liquid Mass Transfer in a Sparged and 
    Non-Sparged Continuous Stirred Tank Reactor with Potential Application in 
    Syngas Fermentation. Fermentation 2019, 5, 75."""
    diam_ratio = (self.agitator.diameter / self.diameter).simplified
    A = 5.3 * math.exp(-5.4*diam_ratio)
    B = 0.47 * diam_ratio**1.3
    C = 0.64 - 1.1 * diam_ratio
    froude_coeff = (self.agitator.diameter / Q(9.81, 'm/s**2')).simplified
    CSA = (math.pi*self.diameter**2/4).simplified
      
    def kla_func(RPS, gas_flow, working_volume):
      froude = (froude_coeff*(RPS)**2)
      aeration = (gas_flow / (RPS*self.agitator.diameter**3)).simplified
      ungassed_power = self.agitator.ungassed_power(RPS)
      lower_gassed_power = ungassed_power *\
        (1-(B-A*param.viscosity/Q(1, 'Pa*s'))*froude**0.25*\
        math.tanh(C*aeration))
      upper_gassed_power = (self.agitator.number - 1)* ungassed_power* \
        (1-(A+B*froude)*aeration**(C+0.04*froude))
      superficial_velocity = float((gas_flow / CSA).simplified)
      kla = 1080*float(((lower_gassed_power+upper_gassed_power) / working_volume).simplified)** \
        0.39 * superficial_velocity**0.79
      kla = Q(kla, '1/min')
      return kla
    return kla_func


if __name__ == '__main__':
  main.run_sim()
  tests.bioreactor_test()