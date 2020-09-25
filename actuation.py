# -*- coding: utf-8 -*-
"""
Holds all physical actuation equipment.  Takes a setpoint and computes actual
quantities.

Wrapper gets all quantities added and returns an actuation dictionary of rates.
"""
import random, pdb, math

from quantities import Quantity as Q
import param, helper_functions

class mixture:
  """Holds components in a mixture.  Introduces error in the exact amount of 
  components in the mixture as well as handles things like changeouts.
  
  Changeouts happen whenever empty.  Can be updated later.
  
  Max added is the maximum mixture available.  None is no limit.
  
  Also responsible for calculating and updating osmo for each mixture."""
  p2 = param.actuation['mixtures']
  def __init__(self, mixture_components, reservoir_size, max_added = None):
    """mixture_components:
      dict of component along with target concentration (g/L)
    """
    if max_added == None:
      max_added = Q(1000000, 'm**3')
    self.max = max_added.simplified
    self.mixture_definition = mixture_components.copy()
    helper_functions.simplify(self.mixture_definition)
    self.empty = False
    self.total = Q(0., 'm**3')
    if reservoir_size ==None:
      reservoir_size = Q(1., 'L').simplified
    self.size = reservoir_size.simplified
    for component, concentration in self.mixture_definition.items():
      if concentration.dimensionality == \
        Q(1, 'kg/m**3').dimensionality:
        self.mixture_definition[component] = (concentration/\
                                  param.molecular_weight[component]).simplified
    
    if param.skip_units:
      self.size = float(self.size)
      helper_functions.remove_units(self.mixture_definition)
      self.total = float(self.total)
      self.max = float(self.max)

    self.replace_source()
      
  def dispense(self, liquid_rate):
    """Calculates number of moles of each component and total liquid volume
    dispensed.  Removes from total."""
    output = {}
    if not self.empty:
      self.remaining -= liquid_rate*param.step_size
      if self.remaining < 0:    #Update at some point to happen during offline assays?
        self.replace_source()
      for component, molarity in self.molarity.items():
        output[component] = molarity*liquid_rate
      output['liquid_volumetric_rate'] = liquid_rate
    return output
    
  def replace_source(self):
    """Generates new mixture with fresh errors."""
    
    if self.total + self.size < self.max:
      # Full topup
      self.molarity = helper_functions.create_media(self.mixture_definition)
      # Osmo is measured for every batch and is part of the definition
      self.mixture_definition['mOsm'] = self.calc_osmo()
      self.remaining = self.size
      self.total += self.size
    elif self.total != self.max:
      # Last topup
      print('Last topup used.')
      self.molarity = helper_functions.create_media(self.mixture_definition)
      # Osmo is measured for every batch and is part of the definition
      self.mixture_definition['mOsm'] = self.calc_osmo()
      self.remaining = self.max - self.total
      self.total = self.max
    else:
      # Empty
      print('Empty mixture!')
      self.empty = True

    
  def calc_osmo(self):
    """Technically this is using actual osmo instead of osmo from media definiton.
    However, since osmo machines are fairly accurate, this won't matter much."""
    if param.skip_units:
      total = 0
    else:
      total = Q(0., 'mol/m**3')
    for component, molarity in self.molarity.items():
      total += molarity
    return total

class MFC:
  """Takes setpoint and returns actual amount dispensed.  Includes constant 
  systematic error and random error.  Both are minimal.
  
  In: Volumetric flowrate (as setpoint)
  Out: Volumetric flowrate"""
  p = param.actuation['MFC']
  def __init__(self, component, self_correcting = False, 
               error_CV = p['systematic_error_CV'], break_chance = 
               p['break_chance']):
    """test doc"""
    self.systematic_error = helper_functions.gauss(1, error_CV)
    # 1 = working fine, 0 = broken
    self.broken_mult = 1
    self.set_point = Q(0, 'L/min').simplified
    self.component = component
    self.break_chance = float((break_chance*param.q_res).simplified)*param.allow_breaking
    
    if param.skip_units:
      self.set_point = 0

  def step(self):
    if self.set_point != 0:
      # Only called if pump is on.  Constant chance of breaking whenever
      # the pump is on.
      if random.random()<self.break_chance:
        self.broken_mult = 0
    addition_rate = self.systematic_error * self.set_point * self.broken_mult
    return {self.component: addition_rate}

class peristaltic():
  """Takes the setpoint and returns actual amount dispensed.  Includes error
  in pump calibration (Default: sigma 3%) as well as any self-correcting measures. 
  
  Includes chance to break. (Default: 1/yr)
  
  Can also be used for perfusion by setting mode to 'recirc' or 'perfusion'
  """
  p = param.actuation['peristaltic']
  def __init__(self, mixture_components = None,
               max_added = None,
               source_size = None, 
               self_correcting = False, 
               error_CV = p['systematic_error_CV'], 
               break_chance = p['break_chance'],
               mode = 'add',
               max_flowrate = Q(300, 'ml/min'),):
    """
    Add mode:
      mixture_components is passed to mixture.  Should be dict of target
      molarities of the mixture.
      
      source_size affects how often it is changed out with a new batch (new errors)
    permeate / recirc mode:
      Used for perfusion.  
    """
    self.systematic_error = helper_functions.gauss(1, error_CV)
    # 1 = working fine, 0 = broken
    self.broken_mult = 1
    self.set_point = Q(0, 'ml/min').simplified
    self.break_chance = float((break_chance*param.q_res).simplified)*param.allow_breaking
    self.shear_units = Q(1, '1/s')
    self.max_flowrate = max_flowrate.simplified
    self.max_pressure = Q(1, 'atm').simplified
    self.flowrate_units = Q(1, 'm**3/s')
    
    if mode == 'add':
      self.output_func = self.add
      if mixture_components == None:
        raise ValueError('Peristaltic Pump needs mixture components in add mode!')
      if type(mixture_components) == mixture:
        self.source = mixture_components
      else:
        self.source = mixture(mixture_components, source_size, max_added = max_added)
    elif mode == 'recirc':
      self.output_func = self.recirc
      if mixture_components != None:
        raise ValueError('Mixture components supplied in recirc mode!')
    elif mode == 'permeate':
      self.output_func = self.perfusion
      if mixture_components != None:
        raise ValueError('Mixture components supplied in permeate mode!')
    else:
      raise ValueError('Unrecognized mode!')
    
    if param.skip_units:
      self.set_point = float(self.set_point)
      self.shear_units = float(self.shear_units)
      self.max_flowrate = float(self.max_flowrate)
      self.max_pressure = float(self.max_pressure)
      self.flowrate_units = float(self.flowrate_units)
    
  def add(self, liquid_rate):
    # Normal mode
    return self.source.dispense(liquid_rate)
    
  def recirc(self, liquid_rate):
    """ Mode for perfusion recirculation. If you knew what the shear generated 
    by a peristaltic pump was, it would go here.  This is an example that
    loosely correlates to 2L shear."""
    shear_rate = (1000+liquid_rate*3000) * self.shear_units*100
    self.ideal_flow_rate = liquid_rate
    return {'recirc_func':self.calc_flow,
            'recirc_shear': shear_rate}
  
  def perfusion(self, liquid_rate):
    """Same as recirc, except for permeate mode (different keyword).  
    Doesn't need shear since there are no cells."""
    self.ideal_flow_rate = liquid_rate
    return {'flux_func':self.calc_flow}
  
  def calc_flow(self, flowrate):
    """Funciton to be passed to environment.  Calculates and returns pressure
    given supplied flowrate."""
    q = float(flowrate) * self.flowrate_units
    if q >= self.ideal_flow_rate:
      #Invalid flowrate.  Create negative pressure so there's a slope in correct direction
      if self.ideal_flow_rate > 0:
        pressure = -(q-self.ideal_flow_rate)/self.ideal_flow_rate*self.max_pressure
      else:
        pressure = -(q)/self.flowrate_units*self.max_pressure
    else:
      pressure = self.max_pressure*(1-1/50*math.log((2*self.ideal_flow_rate)/(self.ideal_flow_rate-q)-1))
    if not type(flowrate)==Q:
      pressure /= self.flowrate_units
    return pressure
    
  def step(self):
    """Called every timestep by actuation.  Calculates whether it breaks, and
    returns dict based on what mode it is in."""
    if self.set_point != 0:
        # Only called if pump is on.  Constant chance of breaking whenever
        # the pump is on.
      if random.random()<self.break_chance:
        self.broken_mult = 0
    #if setpoint > max_flowrate, reduce it
    if self.max_flowrate < self.set_point:
      self.set_point = self.max_flowrate
    liquid_rate = self.systematic_error * self.set_point * self.broken_mult
    return self.output_func(liquid_rate)


class heating_jacket:
  """Adds energy to bioreactor.  Doesn't break."""
  def __init__(self, max_wattage = Q(100, 'W')):
    self.max_wattage = max_wattage.simplified
    self.set_point = 0
    
    if param.skip_units:
      self.max_wattage = float(self.max_wattage)
    
  def step(self):
    return {'heat':self.set_point/100*self.max_wattage}
    
class agitator:
  """ Controls the motor for the agitator."""
  p = param.actuation['agitator']
  def __init__(self, RPS = Q(250, '1/min'),
               error_CV = p['systematic_error_CV'], 
               break_chance = p['break_chance']):

    self.systematic_error = helper_functions.gauss(1, error_CV)
    # 1 = working fine, 0 = broken
    self.broken_mult = 1
    self.set_point = RPS.simplified
    self.break_chance = float((break_chance*param.q_res).simplified)*param.allow_breaking
    
    if param.skip_units:
      self.set_point = float(self.set_point)
    
  def step(self):
    if self.set_point != 0:
      # Only called if pump is on.  Constant chance of breaking whenever
      # the pump is on.
      if random.random()<self.break_chance:
        self.broken_mult = 0
    rate = self.systematic_error * self.set_point * self.broken_mult
    return {'RPS':rate}
  
class levitronix:
  """600-SU model.  Affects recirculation flowrate by adding a pressure."""
  p = param.actuation['levitronix']
  def __init__(self, error_CV = p['systematic_error_CV'], 
               break_chance = p['break_chance'],
               max_RPM = Q(150, '1/s')):
    self.max_RPM = max_RPM.simplified
    self.set_point = Q(0, '1/s')
    self.systematic_error = helper_functions.gauss(1, error_CV)
    # 1 = working fine, 0 = broken
    self.broken_mult = 1
    self.set_point = Q(0, '1/s')
    self.break_chance = float((break_chance*param.q_res).simplified)*param.allow_breaking
    self.pressure_conversion = Q(1, 'kg/(m**4*s)').simplified
    
    if param.skip_units:
      self.max_RPM = float(self.max_RPM)
      self.set_point = float(self.set_point)
      self.pressure_conversion = float(self.pressure_conversion)
      self.set_point = float(self.set_point)
    
  def step(self):
    if self.set_point != 0:
      # Only called if pump is on.  Constant chance of breaking whenever
      # the pump is on.
      if random.random()<self.break_chance:
        self.broken_mult = 0
    if self.set_point > self.max_RPM:
      self.set_point = self.max_RPM
    self.rate = self.systematic_error * self.set_point * self.broken_mult
    shear = self.rate
    return {'recirc_RPM':self.rate, 
            'recirc_func': self.calc_pressure,
            'recirc_shear':shear}
  
  def calc_pressure(self, flowrate):
    """Funciton to be passed to environment.  Calculates and returns flowrate
    based on supplied resistance."""
    q = float(flowrate)
    rps = float(self.rate)
    #constants from fitted curve
    pressure = self.p['a']*rps**self.p['b']+\
      self.p['c']*q**self.p['d']*rps**self.p['e']+\
      self.p['f']*q**self.p['g']
    if type(flowrate)==Q:
      pressure *= flowrate.units
    pressure *= self.pressure_conversion
    return pressure

class wrapper:
  """Actuation wrapper.  Holds the entire actuation list, and cycles through
  to call each step method in order to obtain all actuation additions."""
  def __init__(self, actuation_list):
    """Takes the actuation list that must be cycled through every step."""
    self.actuation_list = actuation_list
    
  def step(self):
    #Create a value of 0 for everything
    if param.skip_units:
      actuation = {component:0 for component in [*param.liquid_components, 
                *param.gas_components, 'liquid_volumetric_rate',
                'gas_volumetric_rate', 'RPS', 'heat', 'recirc_pressure',
                'perfusion_rate',]}
    else:
      actuation =  {}
      actuation['liquid_volumetric_rate']=Q(0., 'm**3/s')
      actuation['gas_volumetric_rate']=Q(0., 'm**3/s')
      actuation['RPS'] = Q(0., '1/s')
      actuation['heat'] = Q(0., 'kg*m**2/s**3')
      for item in param.liquid_components:
        actuation[item]=Q(0., 'mol/s')
      for item in param.gas_components:
        actuation[item] = Q(0., 'm**3/s')
    
    # Cycle through all actuation components and sum results
    for item in self.actuation_list:
      components_added = item.step()
      for key, value in components_added.items():
        if key in ['liquid_volumetric_rate', *param.liquid_components, 
                   *param.gas_components]:
          actuation[key] += value
        else:
          actuation[key] = value
    
    #Computer total gas volumetric rate
    for component in param.gas_components:
      actuation['gas_volumetric_rate'] += actuation[component]
      
    return actuation
  
if __name__ == '__main__':
  tests.actuation_test()