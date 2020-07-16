# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 18:44:54 2020

@author: Racehorse
"""
import random

from quantities import Quantity as Q

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
      if type(entity[key]) == Q:
        entity[key] = float(entity[key].simplified)
  if type(entity) == list:
    for index in range(len(entity)):
      if type(entity[index]) == Q:
        entity[index] = float(entity[index].simplified)