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
get_solubility = {
  'dO2': lambda t: 0.0000158 - t*0.00000015,
  'dCO2': lambda t: 0.0004839 - t*0.0000061,
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
    
    self.diameter = diameter.simplified
    self.width = width
    self.type = impeller_type
    self.number = number
    
    if impeller_type == 'rushton':
      self.power_number = 5.5
    elif impeller_type == 'marine':
      self.power_number = 2.2
    else: raise ValueError('impeller not recognized!')
    
    self.ungassed_power_coeff = \
      (self.power_number*param.actual_cc_density*self.diameter**5).simplified
      
    if param.skip_units == 1:
      self.diameter = float(self.diameter)
      self.width = float(self.width)
      self.ungassed_power_coeff = float(self.ungassed_power_coeff)

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
                 cell_separation_device = None,     # Perfusion only
                 head_pressure = Q(760, 'mmHg'),
                 heat_transfer_coeff = Q(50, 'W/m**2'), #per degree (not supported)
                 ):
    
    self.volume = volume.simplified
    self.working_volume = seed_volume.simplified
    self.diameter = diameter.simplified
    self.sparger_pore_size = sparger_pore_size.simplified
    self.pressure = (sparger_height/2*param.actual_cc_density*param.gravity+head_pressure).simplified
    self.mole = {component:Q(0., 'millimole').simplified for component in param.liquid_components}
    self.mole.update(initial_components)
    for component in self.mole:
      self.mole[component] = self.mole[component].simplified
      if param.skip_units:
        self.mole[component] = float(self.mole[component])
    self.sparger_height = sparger_height.simplified
    self.CSA = self.diameter**2/4*math.pi
    self.overall_heat_transfer_coeff = heat_transfer_coeff.simplified*\
      (2*self.CSA+self.volume/self.CSA*math.pi*self.diameter)   #external area
    self.one_cell = Q(1, 'ce')
    self.solubility_units = Q(1, 's**2*mol/(kg*m**2)') # Solubility / pressure
    self.volumetric_heat_capacity = param.volumetric_heat_capacity
    self.molar = Q(1, 'M').simplified
    self.kLa = Q(1e-10, '1/s') #Needs initial nonzero value
    self.old = {component:Q(0, 'mol/s') for component in param.liquid_components}
    
    if param.skip_units == 1:
      self.volume = float(self.volume)
      self.working_volume = float(self.working_volume)
      self.diameter = float(self.diameter)
      self.sparger_pore_size = float(self.sparger_pore_size)
      self.pressure = float(self.pressure)
      self.sparger_height = float(self.sparger_height)
      self.CSA = float(self.CSA)
      self.overall_heat_transfer_coeff = float(self.overall_heat_transfer_coeff)
      self.one_cell = 1
      self.solubility_units = float(self.solubility_units)
      self.volumetric_heat_capacity = float(self.volumetric_heat_capacity)
      self.molar = float(self.molar)
      self.kLa = float(self.kLa)
      self.old = {component: 0 for component in param.liquid_components}

    
    self.agitator = agitator
    self.current_time = start_time
    # self.sparger_num_pores = num_pores
    self.kla_func = self.create_kla_function()
    self.temperature = initial_temperature
    self.gas_percentages = {'dO2':0.21,
                            'dCO2':0.00045}
    
    

    
  def update_shear(self, RPS):
    """Updates the shear rates with new RPS.
    
    Equations from:
      R. Bowen, Unraveling the mysteries of shear-sensitive mixing systems,
      Chem. Eng. 9 (June) (1986) 55–63.
      """
    ratio = (self.agitator.diameter / self.diameter)**0.3
    self.mean_shear = 4.2*RPS*ratio*self.agitator.diameter /\
                         self.agitator.width
    self.max_shear = self.mean_shear / 4.2 * 9.7
  
  def check_and_update(self, actuation):
    self.current_time += param.resolution
    if not (actuation == self.old and actuation['liquid_volumetric_rate'] == 0):
      self.old = actuation.copy()
      # Is this in per hour or per minute?
      self.kLa = self.kla_func(actuation['RPS'], actuation['gas_volumetric_rate'], self.working_volume)
      self.gas_percentages = self.calc_gas_percentages(actuation)
      self.update_shear(actuation['RPS'])
      
  def calc_gas_percentages(self,actuation):
    """Takes all of the gas flowrates being input and calculates the 
    percentages of each component in it.  Must also convert air to O2 and CO2.
    """
    if param.skip_units: total = 0
    else: total = Q(0., 'L/min').simplified
    percent = {}
    for component in param.gas_components:
      total += actuation[component]
    for component in param.gas_components:
      percent['d'+component] = actuation[component]/total
    # Assume N2 is irrelevant
    percent['dO2'] += percent['dair']*0.21
    percent['dCO2'] += percent['dair']*0.00045
    return percent
  
  def calc_max_consumption(self, actuation):
    """Calculates the maximum average consumption that can be used by cells 
    before the concentration will reach negative values next step."""
    max_consumption = {}
    for component in param.liquid_components:
      if component[1:] in param.gas_components:
        solubility = get_solubility[component](self.temperature)*\
          self.pressure*self.solubility_units
        C_star = self.gas_percentages[component]*solubility
        kLa = self.kLa*param.kLa_ratio[component]
        k_t = math.exp(-kLa*param.step_size)
        max_consumption[component] = self.working_volume * kLa * C_star + \
          self.mole[component]*kLa*k_t/(1-k_t)
      else:
        max_consumption[component] = actuation[component] + \
          self.mole[component]/param.step_size
    return max_consumption
    
  def step(self, actuation, cells):
    """Calculate environmental changes."""

    #Calculate updates from previous step (before updating params)
    for component in param.liquid_components:
      if component[1:] in param.gas_components:
        solubility = get_solubility[component](self.temperature)*\
          self.pressure*self.solubility_units
        C_star = self.gas_percentages[component]*solubility
        kLa = self.kLa*param.kLa_ratio[component]
        k_t = math.exp(-kLa*param.step_size)
        # if component == 'dCO2':
        mol = k_t*self.mole[component] + (1-k_t)*\
          (C_star*self.working_volume - cells['mass_transfer'][component]/kLa)
        if mol < -1e-10:
          pdb.set_trace()
        self.mole[component] = mol
        
      else: 
        transfer_rate = self.old[component]
        # if component =='glucose':
        #   print('BR', component, transfer_rate, self.mole[component])
        self.mole[component] += (transfer_rate - cells['mass_transfer'][component]) * param.step_size

    self.check_and_update(actuation)
    self.working_volume += actuation['liquid_volumetric_rate']*param.step_size
    # print(cells['mass_transfer']['dO2'], cells['mass_transfer']['dCO2'])
    #Temperature
    k_t = math.exp(-self.overall_heat_transfer_coeff*param.step_size/\
                   self.working_volume/self.volumetric_heat_capacity)
    self.temperature = (1-k_t)*(param.environment_temperature+\
                        actuation['heat']/self.overall_heat_transfer_coeff)+\
      self.temperature*k_t

    cell_fraction = cells['total_cells']*cells['volume']/self.one_cell/self.working_volume
    if cell_fraction > 1:
      print('Wayyyyy too many cells')
    osmo = self.total_moles()/(self.working_volume*(1-cell_fraction))
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
      environment.update({component:self.mole[component]/self.working_volume})
    environment['max_consumption'] = self.calc_max_consumption(actuation)
    return environment
    
  def total_moles(self):
    if param.skip_units:
      total = 0
    else: total = Q(0., 'mol')
    for component in self.mole:
      total += self.mole[component]
    return total
  
  def pH(self):
    """Calculates pH of system.  Uses Henderson–Hasselbalch equation for bicarb
    (pKa = 6.36)
    
    net_charge refers to net charge ignoring [H+] and [OH-]
    All species except bicarb are considered strong bases / acids (no equilibria)
    """
    net_charge = 0
    
    for species in self.mole:
      if species in param.positively_charged:
        net_charge += self.mole[species]/self.working_volume/self.molar
      if species in param.negatively_charged:
        net_charge -= self.mole[species]/self.working_volume/self.molar
    pH = -math.log10((-net_charge/2+math.sqrt((net_charge/2)**2+\
      10**-14+10**-6.36*1.0017*self.mole['dCO2']/self.working_volume/self.molar)))
    return pH
  
  def create_kla_function(self):
    """Calculate and return a function for oxygen transfer rate.
    """
    if param.skip_units:
      viscosity = float(param.viscosity)
      gravity = float(param.gravity)
    else: 
      viscosity = param.viscosity / Q(1, 'Pa*s')
      gravity = param.gravity
    diam_ratio = self.agitator.diameter / self.diameter
    A = 5.3 * math.exp(-5.4*diam_ratio)
    B = 0.47 * diam_ratio**1.3
    C = 0.64 - 1.1 * diam_ratio
    froude_coeff = self.agitator.diameter / gravity
    CSA = math.pi*self.diameter**2/4
    
    
      
    def kla_func(RPS, gas_flow, working_volume):
      froude = (froude_coeff*(RPS)**2)
      aeration = gas_flow / (RPS*self.agitator.diameter**3)
      ungassed_power = self.agitator.ungassed_power(RPS)
      
      lower_gassed_power = ungassed_power *\
        (1-(B-A*viscosity)*froude**0.25*\
        math.tanh(C*aeration))
      upper_gassed_power = (self.agitator.number - 1)* ungassed_power* \
        (1-(A+B*froude)*aeration**(C+0.04*froude))
      superficial_velocity = float(gas_flow / CSA)
      kla = 1080*float((lower_gassed_power+upper_gassed_power) / working_volume)** \
        0.39 * superficial_velocity**0.79
      if param.skip_units:
        return kla/60
      return Q(kla, '1/min').simplified
    return kla_func


if __name__ == '__main__':
  main.run_sim()
  tests.bioreactor_test()