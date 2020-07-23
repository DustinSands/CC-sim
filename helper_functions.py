# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 18:44:54 2020

@author: Racehorse
"""
import random, pdb

from quantities import Quantity as Q
import quantities as q

import param

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
    else: component = [component]
    error_ratio = random.gauss(1, param.actuation['mixtures']['component_CV'])
    
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
  'glucose feed rate':'m**3/s',
  'mOsm addition rate':'mol/s',
  'mOsm feed rate':'m**3/s',
  'air':'m**3/s',
  'CO2':'m**3/s',
  'O2':'m**3/s',
  'target_diameter':'m',
  }

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
  'glucose feed rate':'ml/h',
  'mOsm addition rate':'mmol/h',
  'mOsm feed rate':'ml/h',
  'air':'ml/min',
  'CO2':'ml/min',
  'O2':'ml/min',
  'target_diameter':'um',
  }

def get_plotfunc(ax, param):
  """Decides whether the appropriate plotting function is a step or line plot.
  """
  steps = ['glucose addition rate', 'glucose feed rate', 'mOsm feed rate', 
           'mOsm addition rate']
  if param in steps:
    func = lambda x, y, **kwargs: ax.step(x, y, where='post', **kwargs)
  else:
    func = ax.plot
  return func

def scale_units(assays):
  """Adds units and rescales for human interaction / visualization."""
  assays = assays.copy()
  for key in assays:
    if param.skip_units:
      if key in base_unit_lookup:
        units = base_unit_lookup[key]
        assays[key] = Q(assays[key], units)
    if key in rescale_unit_lookup:
      assays[key] = assays[key].rescale(rescale_unit_lookup[key])
  return assays
  

  