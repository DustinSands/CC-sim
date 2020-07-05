# -*- coding: utf-8 -*-
"""
Models properties of the intercellular fluid in bioreactor.

Takes control device outputs (like oxygen flow rates, feed flow rates, etc.)
and outputs the conditions within the bioreactor.

Currently, modeled as a completely homogenous solution.  At some point I'd like
to leave room for heterogeneity.

"""
import math

from quantities import Quantity as Q

import param, tests

"""Defines solubility equation for components that vary with temperature."""
solubility = {
  'dO2': lambda t: Q(0.05 - t*0.000475, 'g/L/kPa'),
  'CO2': lambda t: Q(2.1 - t*0.0275, 'g/L/kPa')
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

    
    self.ungassed_power_coeff = self.power_number*Q(1000, 'kg/m**3')*self.diameter**5
    
  def ungassed_power(self, RPS):
    return RPS**3*self.ungassed_power_coeff
  

class bioreactor:
  """Set up bioreactor configuation, then perform step to update environment
  in the bioreactor.
  
  Assumed to be a perfectly cylindrical vessel for volume calculations."""
  def __init__(self, start_time,
                 initial_components,
                 volume = Q(3, 'L'),
                 agitator = agitator(),
                 diameter = Q(13, 'cm'),     
                 sparger_height = Q(2, 'cm'), # Height from bottom of vessel
                 sparger_pore_size = Q(20, 'um'),
                 cell_separation_device = None,     # Perfusion only
                 head_pressure = Q(760, 'mmHg'),
                 ):
    self.volume = volume
    self.agitator = agitator
    self.CSA = diameter**2/4*math.pi
    self.diameter = diameter
    self.check_list = self.build_checklist()
    self.current_time = start_time
    self.sparger_pore_size = sparger_pore_size
    self.kla_func = self.create_kla_function()
    self.old = {}
    self.pressure = sparger_height/2*param.actual_cc_density*param.gravity+head_pressure
    self.mass = initial_components
    self.sparger_height = sparger_height

    
  def build_checklist(self):
    """The checklist for items that affect:
      gas percentages
      kla (RPS, gas flow, volume)
      working volume
      """
    checklist = ['liquid_volume', '']
      
    return checklist
    
  def update_shear(self, RPS):
    """Updates the shear rates with new RPS.
    
    Equations from:
      R. Bowen, Unraveling the mysteries of shear-sensitive mixing systems,
      Chem. Eng. 9 (June) (1986) 55–63.
      """
    ratio = (self.impeller.diameter / self.diameter)**0.3
    self.mean_shear = 4.2*RPS*ratio*self.impeller.diameter / self.impeller.width
    self.max_shear = mean_shear / 4.2 * 9.7
  
  def check_and_update(self, actuation):
    if not actuation == self.old:
      self.old = actuation
      self.working_volume += actuation['liquid_volume']
      # Is this in per hour or per minute?
      self.kla = self.kla_func(actuation['RPS'], actuation['gas_volume'], self.working_volume)
      self.gas_percentages = self.calc_gas_percentages(actuation)
      self.update_shear(actuation['RPS'])
      
  def calc_gas_percentages(actuation):
    total = 0
    percent = {}
    for component in param.gas_components:
      total += actuation[component]
    for component in param.liquid_components:
      percent[component] = actuation[component] / total
    # Assume N2 is irrelevant
    percent['O2'] += percent['air']*0.21
    return percent
    
  def step(self, actuation, cells):
    """Calculate environmental changes."""
    self.check_and_update(actuation)
    self.current_time += param.resolution
    cell_fraction = (cells['total_cells']*cells['cell_diameter']**3*math.pi/6).simplified
    for component in param.liquid_components:
      if component in param.gas_components:
        C_star = self.gas_percentages[component]*self.pressure*solubility[component](self.temperature)
        transfer_rate = (self.kla*param.kla_ratio[component])*\
          (C_star - self.mass[component]/(self.volume*cell_fraction))
        transfer_rate *= self.volume
      else:
        transfer_rate = actuation[component]
      self.mass[component] += (transfer_rate - cells['mass_transfer'][component]*\
                               self.volume) * param.resolution
      self.mass[component] -= cells['mass_transfer'][component]
    
    self.volume += actuation['liquid_volume']
    
    environment={'shear':self.mean_shear, 
                        'max_shear':self.max_shear,
                        'volume':self.volume,
                        
                        }
    for component in param.liquid_components:
      environment.update({component:self.concentration[component]})
      
    environment['ph'] = 7.3 - self.concentration['CO2']
    
  def create_kla_function(self):
    """Calculate and return a function for oxygen transfer rate.
    
    Equations used based on: 
    Liu, K.; Phillips, J.R.; Sun, X.; Mohammad, S.; Huhnke, R.L.; Atiyeh, H.K. 
    Investigation and Modeling of Gas-Liquid Mass Transfer in a Sparged and 
    Non-Sparged Continuous Stirred Tank Reactor with Potential Application in 
    Syngas Fermentation. Fermentation 2019, 5, 75."""
    diam_ratio = self.agitator.diameter / self.diameter
    A = 5.3 * math.exp(-5.4*diam_ratio)
    B = 0.47 * diam_ratio**1.3
    C = 0.64 - 1.1 * diam_ratio
    froude_coeff = self.agitator.diameter / Q(9.81, 'm/s**2')
    velocity_coeff = 1/(math.pi * self.sparger_pore_size**2/4)
      
    def kla_func(RPS, gas_flow, working_volume):
      froude = froud_coeff*RPS
      aeration = gas_flow / (RPS+self.imp_diam**3)
      ungassed_power = self.agitator.ungassed_power(RPS)
      lower_gassed_power = ungassed_power * (1-(B-A*viscosity)*froude**0.25*\
                                             math.tanh(C*aeration))
      upper_gassed_power = (self.agitator.number - 1)* ungassed_power* \
        (1-(A+B*froude)*aeration**(C+0.04*froude))
      superficial_velocity = gas_flow * velocity_coeff
      kla = 1080*((lower_gassed_power+upper_gassed_power) / working_volume)** \
        0.39 * superficial_velocity**0.79
      kla = Q(kla, '1/min')
      return kla
    return kla_func









if __name__ == '__main__':
  tests.bioreactor_test()