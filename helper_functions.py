# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 18:44:54 2020

@author: Racehorse
"""
import random, pdb, time

from quantities import Quantity as Q
import quantities as q

import param

def gauss(mean, std):
  """Mock random function to ensure consistent results in debug."""
  if param.debug:
    return mean
  else:
    return random.gauss(mean, std)

def create_media(target_molarity):
  """Takes the media definition and returns a created media.
  Includes: 
    Errors from component weighing
    Only (and all) tracked liquid components
    Dissociates molecules.
  Outputs with same units from input.
    """
  zero_molar = next(iter(target_molarity.values()))*0
  end_molarity = {component:zero_molar for component in param.liquid_components}
  
  for component, molarity in target_molarity.items():
    #If component dissociates, add species separately
    if component =='NaHCO3':
      #HCO3 is in equilibrium with CO2 and pH.  Charge is accounted for in Na
      component = ['Na', 'dCO2']
    elif component == 'NaCl':
      component = ['Na', 'Cl']
    elif component == 'KCl':
      component = ['K', 'Cl']
    elif component == 'mOsm':
      #mOsm isn't a component; skip it
      continue
    else: component = [component]
    
    error_ratio = gauss(1, param.actuation['mixtures']['component_CV'])
    
    for species in component:
      if species in param.liquid_components:
        end_molarity[species] = end_molarity[species]+molarity*error_ratio
  
  return end_molarity

def create_media_as_mole(media_definition, volume):
  """Creates media from definition (in molarity) and volume.
  Also simplifies.
  Slow."""
  actual_molarity = create_media(media_definition)
  moles = {}
  volume = volume.simplified
  for species in actual_molarity:
    moles[species] = (volume*actual_molarity[species]).simplified
  return moles

def remove_units(entity):
  if type(entity) == dict:
    for key in entity:
      entity[key] = remove_units(entity[key])
  elif type(entity) == list:
    for index in range(len(entity)):
      entity[index] = remove_units(entity[index])
  elif type(entity) == Q:
    entity = float(entity.simplified)
  return entity
        
def simplify(entity):
  if type(entity) == dict:
    for key in entity:
      entity[key] = simplify(entity[key])
  elif type(entity) == list:
    for index in range(len(entity)):
      entity[index] = simplify(entity[index])
  elif type(entity) == Q:
    entity = entity.simplified
  return entity

#Units of these variables, in SI units
base_unit_lookup = {
  'VCD':'ce/m**3',
  'cell_diameter':'m',
  'mass': 'kg', 
  'BGA_dO2': 'mol/m**3', 
  'BGA_dCO2': 'mol/m**3', 
  'glucose':  'kg/m**3', 
  'IGG': 'kg/m**3', 
  'viability':'dimensionless',
  'glucose addition rate':'kg/s',
  'glucose feed sp':'m**3/s',
  'mOsm addition rate':'mol/s',
  'mOsm feed sp':'m**3/s',
  'air':'m**3/s',
  'CO2':'m**3/s',
  'O2':'m**3/s',
  'target_diameter':'m',
  'rVCD':'ce/m**3',
  'total_cell_volume':'m**3',
  'recirc_rate':'m**3/s',
  'recirc_RPM':'1/s',
  'VVD':'1/s',
  'permeate rate':'m**3/s',
  'topoff':'m**3/s',
  'fiber_shear':'1/s',
  'rglucose':'mol/m**3',
  'amino_acids':'mol/m**3',
  'sieving':'',
  'rperfusion':'m**3/s',
  }

#Units that we actually view them at (for rescaling)
rescale_unit_lookup = {
  'VCD':q.CD,
  'cell_size':'um',
  'mass': 'kg', 
  'BGA_dO2': 'mM', 
  'BGA_dCO2': 'mM', 
  'glucose':  'g/L', 
  'IGG': 'g/L', 
  'cell_diameter': 'um',
  'viability':'percent',
  'glucose addition rate':'g/h',
  'glucose feed sp':'ml/h',
  'mOsm addition rate':'mmol/h',
  'mOsm feed sp':'ml/h',
  'air':'ml/min',
  'CO2':'ml/min',
  'O2':'ml/min',
  'target_diameter':'um',
  'rVCD':'e5c/ml',
  'total_cell_volume':'L',
  'recirc_rate':'L/min',
  'recirc_RPM':'1/min',
  'VVD':'1/d',
  'permeate rate':'ml/min',
  'topoff':'ml/min',
  'fiber_shear':'1/s',
  'rglucose':'mM',
  'amino_acids':'mM',
  'sieving':'%',
  'rperfusion':'ml/min',
  }

timer = {}

class time_tracker():
  """Tracks total time spent between start and stop.  
  """
  def __init__(self):
    self.on = False
    self.total = 0
    if False:
      self.start = self.debug_start
      self.stop = self.debug_stop
  
  def start(self):
    self.start_time = time.perf_counter()
  
  def stop(self):
    self.stop_time = time.perf_counter()
    elapsed = self.stop_time - self.start_time
    self.total += elapsed
    
  def debug_start(self):
    assert self.on ==False
    self.on = True
    self.start_time = time.perf_counter()
  
  def debug_stop(self):
    assert self.on == True
    self.on = False
    self.stop_time = time.perf_counter()
    elapsed = self.stop_time - self.start_time
    self.total += elapsed
    
  def print(self):
    print(f'{self.name} total elapsed:{self.total} s')
  
  def reset(self):
    self.on = False
    self.total = 0
    
def reset_times():
  for name in timer_list:
    timer[name].reset()

def print_times():
  """Helper function to print how long the agent has spent doing various tasks.
  Resets all timers at end. """
  percent = {}
  for name in timer:
    percent[name] = round(timer[name].total / timer['total'].total *100,1)
  print(f'Total time:{timer["total"].total}')
  print(f'Env:{percent["env"]}%')
  print(f'-Rootfinding:{percent["root"]}%')
  print(f'Cells:{percent["cells"]}%')
  print(f'Control, Actuation:{percent["con"]}%')
  print(f'Assays:{percent["assays"]}%')
  print(f'Metrics:{percent["metrics"]}%')
  for name in timer:
    timer[name].reset()

def get_plotfunc(ax, param):
  """Decides whether the appropriate plotting function is a step or line plot.
  
  Also adds markers for offline variables.
  """
  steps = ['glucose addition rate', 'glucose feed sp', 'mOsm feed sp', 
           'mOsm addition rate', 'Concentrated Feed', 'VVD', 'permeate rate',
           ]
  offline = ['IGG', 'VCD', 'glucose', 'mOsm', 'viability', 'cell_diameter']
  if param in steps:
    func = lambda x, y, **kwargs: ax.step(x, y, where='post', **kwargs)
  elif param in offline:
    func = lambda *args, **kwargs: ax.plot(*args, **kwargs, marker = 'x')
  else:
    func = lambda *args, **kwargs: ax.plot(*args, **kwargs)
  return func

def scale_units(assays):
  """Adds units and rescales for human interaction / visualization."""
  if type(assays) == list:
    for entry in assays:
      entry = scale_units(entry)
  else:
    for key in assays:
      if type(assays[key]) == dict:
        scale_units(assays[key])
      else:
        if param.skip_units:
          if key in base_unit_lookup:
            units = base_unit_lookup[key]
            assays[key] = Q(assays[key], units)
        if key in rescale_unit_lookup:
          assays[key] = assays[key].rescale(rescale_unit_lookup[key])
  

  