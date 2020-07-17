# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 18:44:54 2020

@author: Racehorse
"""
import random

from quantities import Quantity as Q
import quantities as q

import param

def create_media_as_mass(media_definition, volume):
  """Creates media from definition (in molarity) and volume."""

  media = {component:Q(0., 'g').simplified for component in param.liquid_components}
  
  for component, molarity in media_definition.items():
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
        media[species] += \
          (param.molecular_weight[species]*molarity*volume*error_ratio).simplified
        
  return media

def create_media_as_mole(media_definition, volume):
  """Creates media from definition (in molarity) and volume."""
  

  media = {component:Q(0., 'mol').simplified for component in param.liquid_components}
  
  for component, molarity in media_definition.items():
    #If component dissociates, add species separately
    if component =='NaHCO3':
      #HCO3 is in equilibrium with CO2 and pH.  Charge is accounted for in Na
      component = ['Na', 'dCO2']
    elif component == 'NaCl':
      component = ['Na', 'Cl']
    elif component == 'KCl':
      component = ['K', 'Cl']
    else: component = [component]
    error_ratio = random.gauss(1, 0.005)
    
    for species in component:
      if species in param.liquid_components:
        media[species] += \
          (molarity*volume*error_ratio).simplified
        
  return media

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
  'glucose_solution':'m**3/s',
  'mass': 'kg', 
  'BGA_dO2': 'mol/m**3', 
  'BGA_dCO2': 'mol/m**3', 
  'glucose':  'kg/m**3', 
  'IGG': 'kg/m**3', 
  'viability':'percent',
  'glucose_feed':'kg/s',
  }

rescale_unit_lookup = {
  'VCD':q.CD,
  'cell_size':'um',
  'glucose_solution':'ml/h',
  'mass': 'kg', 
  'BGA_dO2': 'mM', 
  'BGA_dCO2': 'mM', 
  'glucose':  'g/L', 
  'IGG': 'g/L', 
  'cell_diameter': 'um',
  'viability':'percent',
  'glucose_feed':'g/hour',
  }

def scale_assays(assays):
  """Adds units and rescales assays for human interaction / visualization."""
  assays = assays.copy()
  for key in assays:
    if param.skip_units:
      if key in base_unit_lookup:
        units = base_unit_lookup[key]
        assays[key] = Q(assays[key], units)
    if key in rescale_unit_lookup:
      assays[key] = assays[key].rescale(rescale_unit_lookup[key])
  return assays
  

  