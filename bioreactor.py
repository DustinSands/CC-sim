# -*- coding: utf-8 -*-
"""
Models properties of the intercellular fluid in bioreactor.

Takes control device outputs (like oxygen flow rates, feed flow rates, etc.)
and outputs the conditions within the bioreactor.

Modeled as a completely homogenous solution.

"""
import math
import pdb

from scipy.optimize import root_scalar
from quantities import Quantity as Q

import param, helper_functions

"""Defines solubility equation for components that vary with temperature."""
get_solubility = {
  'dO2': lambda t: 0.0000158 - t*0.00000015,
  'dCO2': lambda t: 0.0004839 - t*0.0000061,
  }

class hollow_fiber:
  """Perfusion hollow fiber master class.
  Every step, takes flowrate out, recirculation rate, and molarity of components
  Returns a dict of sieving rate for each component."""
  def __init__(self, lumen_ID = Q(1, 'mm'), length = Q(0.2, 'm'), n = 75,
               pore_size = Q(0.2, 'um')):
    self.number = n
    self.length = length.simplified
    self.pore_size = pore_size.simplified
    self.ID = lumen_ID.simplified
    self.TMP_resistance = Q(0.01, '1/m**3')  #Initial Value, will be updated
    
    if param.skip_units:
      self.length = float(self.length)
      self.pore_size = float(self.pore_size)
      self.ID = float(self.ID)
      self.TMP_resistance = float(self.TMP_resistance)
      
    self.CSA = self.number * math.pi * self.ID**2 / 4
    self.SA = self.number*math.pi*self.ID * self.length
    
  def shear_rate(self, recirc_rate):
    return 32* recirc_rate / (math.pi * self.ID**3 * self.number)
  
  def flux(self, perfusion_flowrate):
    return perfusion_flowrate / self.SA
  
  def calc_resistances(self, viscosity, cell_fraction):
    # print(f'Viscosity: {viscosity}')
    R = 128*viscosity*self.length/(math.pi * self.ID**4 * self.number) 
    #arbitrary modifier for high cell fraction. Could use research and update
    cell_modifier = (1-cell_fraction**2)
    recirc_pressure = lambda q: R*q
    flux_pressure = lambda q: self.TMP_resistance*q*viscosity/cell_modifier
    return recirc_pressure, flux_pressure
  
  def step(self, perfusion_flowrate, recirc_rate, molarity):
    """calcs shear, flux, and passes to subclassed sieving function."""
    shear = self.shear_rate(recirc_rate)
    flux = self.flux(perfusion_flowrate)
    sieving = self.calc_sieving(shear, flux, molarity)
    return shear, sieving
    
class perfect_filter(hollow_fiber):
  """Doesn't suffer from sieving."""
  membrane_type = 'perfect'
  def calc_sieving(self, shear, flux, molarity):
    """No sieving issues. Perfectly filters cells."""
    sieving = {component:1 for component in param.liquid_components}
    return sieving
  
class PES(hollow_fiber):
  """Sieving decreases due to sheared cells."""
  membrane_type = 'PES'
  p = param.filters['PES']
  def __init__(self, sieving_efficiency = 1, **kwargs):
    """sieving_efficiency is a modifier to make siever better or worse."""
    super().__init__(**kwargs)
    # comp B ~0.001 mM is 50% sieiving
    self.sieving_efficiency = sieving_efficiency * self.p['sieving_efficiency']
    if param.skip_units:
      self.sieving_efficiency= float(self.sieving_efficiency)
      
  def calc_sieving(self, shear, flux, molarity):
    if flux == 0:
      #edge case
      if param.skip_units:
        efficiency = 1000
      else:
        efficiency = Q(1000, 'mol/m**3').simplified
    else:
      efficiency = self.sieving_efficiency*shear/flux
    ratio = 1-molarity['component_B']/\
      (molarity['component_B']+efficiency)
    self.TMP_resistance = self.p['flux_resistance']/ratio
    sieving = {component:1 for component in param.liquid_components}
    sieving['IGG_a'] = ratio
    sieving['IGG_b'] = ratio
    sieving['IGG_n'] = ratio
    sieving['component_B'] = 0.05
    # sieving['component_A'] = 0.9
    return sieving
  

  
class agitator:
  """Impeller class.  Holds impeller parameters and calculates power."""
  p = param.environment
  def __init__(self, impeller_type = 'rushton', 
               number = 1, 
               diameter = Q(6, 'cm'),
               width = Q(2, 'cm')
               ):
    
    self.diameter = diameter.simplified
    self.width = width.simplified
    self.type = impeller_type
    self.number = number
    
    if impeller_type == 'rushton':
      self.power_number = 5.5
    elif impeller_type == 'marine':
      self.power_number = 2.2
    else: raise ValueError('impeller not recognized!')
    
    self.ungassed_power_coeff = \
      (self.power_number*self.p['actual_cc_density']*self.diameter**5).simplified
      
    if param.skip_units == 1:
      self.diameter = float(self.diameter)
      self.width = float(self.width)
      self.ungassed_power_coeff = float(self.ungassed_power_coeff)

  def ungassed_power(self, RPS):
    return RPS**3*self.ungassed_power_coeff
  

class bioreactor:
  """Standard cylindrical bioreactor environment.  Refer to technical doc
  for more details.
  Supports:
    Temperature
    Mass transfer (aeration, cells, liquid input)
    Perfusion
    pH (using carbonate)
    Shear
    Cell-free volume
  """
  p = param.environment
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
    head_pressure = head_pressure.simplified
    sparger_height = sparger_height.simplified
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
    self.molar = Q(1, 'M').simplified
    self.kLa = Q(1e-10, '1/s') #Needs initial nonzero value
    self.old = {component:Q(0, 'mol/s') for component in param.liquid_components}
    self.viscosity = self.p['viscosity']
    self.fiber_shear = Q(0, '1/s')
    self.perfusion_flowrate = Q(0, 'm**3/s')
    self.sieving = {component:1 for component in param.liquid_components}
    self.recirc_flowrate= Q(0, 'm**3/s') #initial guess
    self.flowrate_units = Q(1, 'm**3/s')
    if cell_separation_device == None:
      self.perfusion = False
    else:
      self.hf = cell_separation_device
      self.perfusion = True
      
    if param.skip_units:
      #Convert to float to save processing power
      self.volume = float(self.volume)
      self.working_volume = float(self.working_volume)
      self.diameter = float(self.diameter)
      self.sparger_pore_size = float(self.sparger_pore_size)
      head_pressure = float(head_pressure)
      sparger_height = float(sparger_height)
      self.sparger_height = float(self.sparger_height)
      self.CSA = float(self.CSA)
      self.overall_heat_transfer_coeff = float(self.overall_heat_transfer_coeff)
      self.one_cell = 1
      self.solubility_units = float(self.solubility_units)
      self.molar = float(self.molar)
      self.kLa = float(self.kLa)
      self.old = {component: 0 for component in param.liquid_components}
      self.viscosity = float(self.viscosity)
      self.fiber_shear = float(self.fiber_shear)
      self.perfusion_flowrate = float(self.perfusion_flowrate)
      self.recirc_flowrate = float(self.recirc_flowrate)
      self.flowrate_units = float(self.flowrate_units)

    self.cellfree_volume = self.working_volume #only for first step
    self.pressure = sparger_height/2*self.p['actual_cc_density']*self.p['gravity']+head_pressure
    self.base_viscosity = self.viscosity
    self.agitator = agitator
    self.current_time = start_time
    self.kla_func = self.create_kla_function()
    self.temperature = initial_temperature
    self.gas_percentages = {'dO2':0.21,
                            'dCO2':0.00045}
    
  def update(self, actuation, cells):
    """Updates the state of the bioreactor (sans mass transfer).  Includes things 
    like temperature, shear, osmo, etc."""
    self.current_time += param.resolution
    self.old = actuation.copy()
    #Aeration
    self.kLa = self.kla_func(actuation['RPS'], actuation['gas_volumetric_rate'], self.working_volume)
    self.gas_percentages = self.calc_gas_percentages(actuation)
    volume_added = actuation['liquid_volumetric_rate']*param.step_size
    #Temperature from conduction
    k_t = math.exp(-self.overall_heat_transfer_coeff*param.step_size/\
                   self.working_volume/self.p['volumetric_heat_capacity'])
    self.temperature = (1-k_t)*(self.p['environment_temperature']+\
                        actuation['heat']/self.overall_heat_transfer_coeff)+\
                        self.temperature*k_t
    #Temperature from convection
    self.working_volume += volume_added
    alpha = volume_added / self.working_volume
    self.temperature += alpha*(self.p['environment_temperature']-self.temperature)
    #calc osmo, viscosity
    self.cell_fraction = cells['total_cells']*cells['volume']/self.one_cell/self.working_volume
    if self.cell_fraction > 1:
      print('Wayyyyy too many cells.  Something is wrong.')
    self.cellfree_volume = (1-self.cell_fraction)*self.working_volume
    #viscosity and volumetric cell density is roughly linearly correlated
    self.viscosity = cells['dry_weight'] / self.working_volume * self.p['viscosity_effect'] + self.base_viscosity
    self.molarity = {component:self.mole[component]/self.cellfree_volume for 
                     component in param.liquid_components}
    if self.perfusion:
      self.update_perfusion(actuation)
    self.update_shear(actuation)
    # don't add 'mOsm' until after updating perfusion
    self.update_osmo(self.molarity)
  
  def update_shear(self, actuation):
    """Updates the shear rate in BR from agitation rate."""
    ratio = (self.agitator.diameter / self.diameter)**0.3
    self.mean_shear =  4.2*actuation['RPS']*ratio*self.agitator.diameter /\
                         self.agitator.width
    if self.perfusion:
      self.max_shear = max(self.mean_shear / 4.2 * 9.7, #shear from impeller
                           self.fiber_shear, #shear from hollow fiber
                           actuation['recirc_shear'], #shear from recirc pump
                           )
    else:
      self.max_shear = self.mean_shear / 4.2 * 9.7
  
  def calc_flowrate(self, function_A, function_B, initial_guess):
    """Given two functions, find where their intersection lies.  
    
    Input Functions take flowrate and should return resistance.  
    
    Returns the volumetric flowrate and pressure difference."""
    helper_functions.timer['root'].start()
    cost_function = lambda flowrate:function_A(flowrate)-function_B(flowrate)
    solution = root_scalar(cost_function, x0=initial_guess, method='brentq', 
                           bracket = [0, 0.003], xtol = 1e-10)
    flowrate = solution.root*self.flowrate_units
    assert solution.converged==True
    pressure = function_A(flowrate)
    helper_functions.timer['root'].stop()
    return flowrate, pressure
    
  def update_perfusion(self, actuation):
    """Calculates resistances to recirculation and shear from the hollow-fiber.
    Then, uses the functions passed by actuation to calculate flowrates."""
    recirc_resistance, flux_resistance = self.hf.calc_resistances(self.viscosity, self.cell_fraction)
    self.recirc_flowrate, self.recirc_pressure = self.calc_flowrate(
      recirc_resistance, actuation['recirc_func'], self.recirc_flowrate)
    self.perfusion_flowrate, self.perfusion_pressure = self.calc_flowrate(
      flux_resistance, actuation['flux_func'], self.perfusion_flowrate)
    #calculate hollow fiber parameters
    self.fiber_shear, self.sieving = \
      self.hf.step(self.perfusion_flowrate,
                   self.recirc_flowrate, 
                   self.molarity)
    self.permeate_molarity = {component:value*self.sieving[component] for 
                              component, value in self.molarity.items()}
    self.update_osmo(self.permeate_molarity)
    self.working_volume -= self.perfusion_flowrate*param.step_size
      
  def calc_gas_percentages(self, actuation):
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
    if self.perfusion:
      #Use perfusion mass balances
      for component in param.liquid_components:
        if component[1:] in param.gas_components:
          #Mass balance includes gas mass transfer
          solubility = get_solubility[component](self.temperature)*\
            self.pressure*self.solubility_units
          C_star = self.gas_percentages[component]*solubility
          kLa = self.kLa*param.kLa_ratio[component]
          k = (kLa*self.working_volume+self.perfusion_flowrate*self.sieving[component])/self.cellfree_volume
          if k*param.step_size > 50:
            #prevent exponential overflow
            max_consumption[component] = kLa*self.working_volume*C_star
          else: 
            max_consumption[component] = k*self.mole[component]/(
              math.exp(k*param.step_size)-1)+kLa*self.working_volume*C_star
        else:
          k = self.perfusion_flowrate*self.sieving[component]/\
              self.cellfree_volume
          if k == 0:
            max_consumption[component] = actuation[component] + \
              self.mole[component]/param.step_size
            #Use batch mode equation
          else:
            max_consumption[component] = actuation[component]+self.mole[component]*k/\
              (math.exp(k*param.step_size)-1)
          if max_consumption[component] < -1e-10:
            pdb.set_trace()
    else:
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
    
  def update_osmo(self, molarity):
    """Calculates and updates the osmo for the given sample."""
    if param.skip_units:
      total = 0
    else: total = Q(0., 'mol/m**3')
    for component in molarity:
      total += molarity[component]
    molarity['mOsm'] = total
  
  def pH(self):
    """Calculates pH of system.  Uses Henderson–Hasselbalch equation for bicarb
    (pKa = 6.36)
    
    net_charge refers to net charge ignoring [H+] and [OH-]
    All species except bicarb are considered strong bases / acids (no equilibria)
    """
    net_charge = 0
    
    for species in self.mole:
      if species in param.positively_charged:
        net_charge += self.mole[species]/self.cellfree_volume/self.molar
      if species in param.negatively_charged:
        net_charge -= self.mole[species]/self.cellfree_volume/self.molar
    pH = -math.log10((-net_charge/2+math.sqrt((net_charge/2)**2+\
      10**-14+10**-6.36*1.0017*self.mole['dCO2']/self.cellfree_volume/self.molar)))
    return pH

  def create_kla_function(self):
    """Calculate and return a function for oxygen transfer rate.
    """
    if param.skip_units:
      vu = 1 #viscosity units
      gravity = float(self.p['gravity'])
    else: 
      vu = Q(1, 'Pa*s')
      gravity = self.p['gravity']
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
        (1-(B-A*self.viscosity/vu)*froude**0.25*\
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
  
  def update_moles(self, cells):
    """Calculate updates from previous step (before updating params)"""
    if self.perfusion:
      #Expects there to be perfusion-related variables defined
      for component in param.liquid_components:
        if component[1:] in param.gas_components:
          #Mass balance includes gas transfer
          solubility = get_solubility[component](self.temperature)*\
            self.pressure*self.solubility_units
          C_star = self.gas_percentages[component]*solubility
          kLa = self.kLa*param.kLa_ratio[component]
          k = kLa*self.working_volume + self.perfusion_flowrate*self.sieving[component]
          self.mole[component] = (kLa*self.working_volume*C_star - cells['mass_transfer'][component])*\
            (1-math.exp(-k*param.step_size / self.cellfree_volume))*self.cellfree_volume/k+\
              self.mole[component]*math.exp(-k*param.step_size/self.cellfree_volume)
          
        else:
          #Mass balance comes from actuation (constant)
          k = self.perfusion_flowrate*self.sieving[component]/\
              self.cellfree_volume
          if k == 0:
            # Use batch mode equation
            self.mole[component] += (self.old[component] - cells['mass_transfer'][component]) * param.step_size
          else:
            self.mole[component] = (self.old[component]-cells['mass_transfer'][component])*\
              (1-math.exp(-k*param.step_size))/k+self.mole[component]*math.exp(-k*param.step_size)
          # print(f'{component}:{self.mole[component]}, {self.old[component]}, {cells["mass_transfer"][component]}')
        if self.mole[component] < -1e-15:
          pdb.set_trace()
    else:
      #Batch mode
      for component in param.liquid_components:
        if component[1:] in param.gas_components:
          solubility = get_solubility[component](self.temperature)*\
            self.pressure*self.solubility_units
          C_star = self.gas_percentages[component]*solubility
          kLa = self.kLa*param.kLa_ratio[component]
          k_t = math.exp(-kLa*param.step_size)
          self.mole[component] = k_t*self.mole[component] + (1-k_t)*\
            (C_star*self.working_volume - cells['mass_transfer'][component]/kLa)
        else:
          transfer_rate = self.old[component]
          self.mole[component] += (transfer_rate - cells['mass_transfer'][component]) * param.step_size
          
  def step(self, actuation, cells):
    """Calculate environmental changes."""
    #Update molarity from previous step
    self.update_moles(cells)

    #Update conditions
    self.update(actuation, cells)

    #Build output
    environment={'shear':self.mean_shear, 
                 'max_shear':self.max_shear,
                 'volume':self.working_volume,
                 'pH': self.pH(),
                 'temperature':self.temperature,
                 'time': self.current_time,
                 'viscosity': self.viscosity,
                 'mass':self.working_volume*self.p['actual_cc_density'],
                 'br_molarity':self.molarity,
                 'cell_fraction':self.cell_fraction
                 }
    if self.perfusion:
      environment.update({'recirc_rate': self.recirc_flowrate,
                          'fiber_shear': self.fiber_shear,
                          'recirc_RPM': actuation['recirc_RPM'],
                          'recirc_shear': actuation['recirc_shear'],
                          'permeate_molarity': self.permeate_molarity,
                          'perfusion_flowrate': self.perfusion_flowrate,
                          'TMP': self.perfusion_pressure,
                          })
    environment['max_consumption'] = self.calc_max_consumption(actuation)
    return environment
    
if __name__ == '__main__':
  tests.bioreactor_test()