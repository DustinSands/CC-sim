# -*- coding: utf-8 -*-
"""
Holds all physical actuation equipment.  Takes a setpoint and computes actual
quantities.

Wrapper gets all quantities added and returns an actuation dictionary of rates.
"""
import random

from quantities import Quantity as Q
import param

class mixture:
  """Holds components in a mixture.  Introduces error in the exact amount of 
  components in the mixture as well as handles things like changeouts.
  
  Changeouts happen whenever empty.  Can be updated later."""
  p2 = param.actuation['mixtures']
  def __init__(self, mixture_components):
    """mixture_components:
      dict of component along with target concentration
    """
    if reservoir_size ==None:
      reservoir_size = Q(1, 'L')
    self.concentration = {}
    self.target_concentrations = mixture_components
    self.size = reservoir_size
    self.replace_source()

  def dispense(self, rate):
    output = {}
    total_liquid = rate*param.q_res
    self.remaining -= total_liquid
    if self.remaining < 0:    #Update at some point to happen during offline assays?
      self.replace_source()
    for component, concentration in self.concentration:
      output[component] = concentration*total_liquid
    output['liquid_volume'] = total_liquid
    return output
    
  def replace_source(self):
    """Generates new mixture with fresh errors."""
    for component, concentration in self.target_concentrations.items():
      error = random.gauss(0, self.p2['component_CV'])
      self.concentration[component] = concentration*(1+error)
    self.remaining = self.size
      
    

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
    self.set_point = Q(0, 'L/min')
    self.component = component

  
  def step(self):
    if self.set_point != 0:
      # Only called if pump is on.  Constant chance of breaking whenever
      # the pump is on.
      if random.random()<float(break_chance*param.resolution):
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
    concentrations of the mixture.
    
    source_size affects how often it is changed out with a new batch (new errors)
    """
    self.systematic_error = random.gauss(1, error_CV)
    # 1 = working fine, 0 = broken
    self.broken_mult = 1
    self.set_point = Q(0, 'ml/min')
    self.source = mixture(mixture_components, source_size)
    
  def step(self):
    # Only called if pump is on.  Constant chance of breaking whenever
    # the pump is on.
    if random.random()<float(break_chance*param.resolution):
      self.broken_mult = 0
    total_liquid = self.error * self.set_point * self.broken_mult
    return self.dispense(total_liquid)
    
  
# class scale:
#   """A scale.  Assumes you're using it within range.  These things are pretty
#   reliable and accurate, so no sources of error introduced here."""

class wrapper:
  def __init__(self, actuation_list):
    self.actuation_list = actuation_list
    
  def step(self):
    actuation = {}
    for item in param.liquid_components:
      actuation[item]=0
    actuation['liquid_volume']=0
    
    for item in self.actuation_list:
      components_added = item.step()
      for key, value in components_added.items():
        actuation[key] += value
        
    gas_volume = Q(0, 'L/min')
    for component in param.gas_components:
      total_gas += actuation[component]
    return actuation
    # for component in param.liquid_components:
    #   actuation.update({component:    })
    #   total +=
    # {'liquid_volume_rate':total}
    pass    