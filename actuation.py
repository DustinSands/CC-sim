# -*- coding: utf-8 -*-
"""
Holds all physical actuation equipment.  Takes a setpoint and computes actual
quantities.

Wrapper gets all quantities added and returns an actuation dictionary of rates.
"""
import random, pdb

from quantities import Quantity as Q
import param, helper_functions

class mixture:
  """Holds components in a mixture.  Introduces error in the exact amount of 
  components in the mixture as well as handles things like changeouts.
  
  Changeouts happen whenever empty.  Can be updated later."""
  p2 = param.actuation['mixtures']
  def __init__(self, mixture_components, reservoir_size):
    """mixture_components:
      dict of component along with target concentration (g/L)
    """
    if reservoir_size ==None:
      reservoir_size = Q(1., 'L').simplified
    
    self.target_molarity = {}
    self.mixture_definition = mixture_components.copy()
    self.size = reservoir_size.simplified
    
    for component, concentration in mixture_components.items():
      if concentration.simplified.dimensionality == \
        Q(1, 'kg/m**3').dimensionality:
        self.target_molarity[component] = (concentration/\
                                           param.molecular_weight[component]).simplified
      else: 
        self.target_molarity[component] = concentration.simplified
    if param.skip_units:
      self.size = float(self.size)
      helper_functions.remove_units(self.target_molarity)
    
    self.replace_source()
      
  def dispense(self, liquid_rate):
    output = {}
    self.remaining -= liquid_rate*param.step_size
    if self.remaining < 0:    #Update at some point to happen during offline assays?
      self.replace_source()
    for component, molarity in self.molarity.items():
      output[component] = molarity*liquid_rate
    output['liquid_volumetric_rate'] = liquid_rate
    return output
    
  def replace_source(self):
    """Generates new mixture with fresh errors."""
    self.molarity = helper_functions.create_media(self.target_molarity)
    self.remaining = self.size
    
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
  
  In: Volumetric flowrate
  Out: Volumetric flowrate"""
  p = param.actuation['MFC']
  def __init__(self, component, self_correcting = False, 
               error_CV = p['systematic_error_CV'], break_chance = 
               p['break_chance']):
    """test doc"""
    self.systematic_error = random.gauss(1, error_CV)
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
  
  Includes chance to break. (Default: 1/yr)"""
  p = param.actuation['peristaltic']
  def __init__(self, mixture_components,
               source_size = None, 
               self_correcting = False, 
               error_CV = p['systematic_error_CV'], 
               break_chance = p['break_chance']):
    """mixture_components is passed to mixture.  Should be dict of target
    molarities of the mixture.
    
    source_size affects how often it is changed out with a new batch (new errors)
    """
    self.systematic_error = random.gauss(1, error_CV)
    # 1 = working fine, 0 = broken
    self.broken_mult = 1
    self.set_point = Q(0, 'ml/min').simplified
    self.source = mixture(mixture_components, source_size)
    self.break_chance = float((break_chance*param.q_res).simplified)*param.allow_breaking
    
    if param.skip_units:
      self.set_point = float(self.set_point)
    
  def step(self):
    if self.set_point != 0:
      # Only called if pump is on.  Constant chance of breaking whenever
      # the pump is on.
      if random.random()<self.break_chance:
        self.broken_mult = 0
    liquid_rate = self.systematic_error * self.set_point * self.broken_mult
    return self.source.dispense(liquid_rate)

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
  def __init__(self, RPS,
               error_CV = p['systematic_error_CV'], 
               break_chance = p['break_chance']):

    self.systematic_error = random.gauss(1, error_CV)
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
# class scale:
#   """A scale.  Assumes you're using it within range.  These things are pretty
#   reliable and accurate, so no sources of error introduced here."""

class wrapper:
  def __init__(self, actuation_list):
    self.actuation_list = actuation_list
    # self.initial_actuation =  {}
    # self.initial_actuation['liquid_volume']=Q(0., 'L/min').simplified
    # self.initial_actuation['gas_volume']=Q(0., 'L/min').simplified
    # self.initial_actuation['RPS'] = Q(0., '1/s').simplified
    # self.initial_actuation['heat'] = Q(0., 'W').simplified
    # for item in param.liquid_components:
    #   self.initial_actuation[item]=Q(0., 'mol/min').simplified
    # for item in param.gas_components:
    #   self.initial_actuation[item] = Q(0., 'L/min').simplified
    # if param.skip_units:
    #   helper_functions.remove_units(self.initial_actuation)
    
  def step(self):

    # actuation = self.initial_actuation.copy()
    actuation =  {}
    actuation['liquid_volume']=Q(0., 'L/min').simplified
    actuation['gas_volume']=Q(0., 'L/min').simplified
    actuation['RPS'] = Q(0., '1/s').simplified
    actuation['heat'] = Q(0., 'W').simplified
    for item in param.liquid_components:
      actuation[item]=Q(0., 'mol/min').simplified
    for item in param.gas_components:
      actuation[item] = Q(0., 'L/min').simplified
    if param.skip_units:
      actuation = {component:0 for component in [*param.liquid_components, 
                *param.gas_components, 'liquid_volumetric_rate',
                'gas_volumetric_rate', 'RPS', 'heat']}
    
    
    for item in self.actuation_list:
      components_added = item.step()
      for key, value in components_added.items():
        # print(item, key, value)
        actuation[key] += value
        
    gas_volume = Q(0., 'L/min')
    for component in param.gas_components:
      actuation['gas_volumetric_rate'] += actuation[component]
    return actuation
    # for component in param.liquid_components:
    #   actuation.update({component:    })
    #   total +=
    # {'liquid_volume_rate':total}
    pass    