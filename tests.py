# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 06:18:03 2020

@author: sckuz
"""
import random

import numpy as np
import quantities as q
from quantities import Quantity as Q

import controls, assays, bioreactor, cell_sim
import param

def set_seeds(seed):
  np.random.seed(seed)
  random.seed(seed)

def controls_tests():
  controls.fed_batch_feed({'glucose':Q(500, 'g/L')},
                             initial_volume = 2, sample_interval = 1440, 
                             cpp = 'glucose', set_point = 2,
                             initial_time = np.datetime64('2020-01-01'), target_seeding_density = 10)
  aeration = controls.DH_aeration(60)
  # print(aeration.step(obs1))
  # print(aeration.actuation[0].set_point)
  # print(aeration.step(obs2))
  # print(aeration.actuation[0].set_point)
  assert aeration.step(obs1)['aeration_PID']==5
  assert aeration.step(obs2)['aeration_PID']==11.083333333333332
  
  print('passed feeding strategy')
  
def quantities_test():
  assert float((Q(2, 'h')/Q(1, 'min')).simplified) == 120
  assert Q(2, 'kg/min')*Q(1, 'h') == Q(120, 'kg')
  assert np.timedelta64(1, 'D')/np.timedelta64(1, 'h') == 24
  
  VCD = Q(10, 'e5c/ml')
  sVCD = VCD.simplified
  sVCD = sVCD.rescale(q.CD)
  assert float(VCD) == float(sVCD)
  print('passed quantities')

def assays_test():
  # set_seeds(0)
  BGA = assays.BGA()
  oxygen = assays.O2_probe(np.datetime64('2018-01-01T00:00'), cal_reference = BGA)
  pH = assays.pH_probe(np.datetime64('2018-01-01T00:00'), cal_reference = BGA)

  oxygen.one_point(np.datetime64('2018-01-01T00:00'))
  
  print(oxygen.read_value(env1, cells))
  # assert 61.30524904989422 == oxygen.read_value(state, cells)['O2']
  print(oxygen.read_value(env2, cells))
  # assert 61.331671360837824==oxygen.read_value(state2, cells)['O2']
  print(pH.read_value(env1, cells))
  # assert 6.9982197908112544 == pH.read_value(state, cells)['pH']
  print(pH.read_value(env2, cells))
  # assert 7.107234728864488==pH.read_value(state2, cells)['pH']
  
  vicell = assays.cell_counter()
  print(vicell.read_value(env1, cell_state))
  # assert Q(10.72273548515264, 'e5c/ml') == vicell.read_value(state, cells)
  
  print('passed assays')
  
def cell_tests():
  cell_line = cell_sim.gen_cell_line()
  cell_culture = cell_sim.cell_wrapper(cell_line, Q(1, 'e6c/ml'))
  pass

def bioreactor_test():
  pass

def all_tests():
  pass
env1 = {'time':np.datetime64('2018-01-05T00:00'),
         'dO2': 60,
         'pH':7,
         }
env2 = {'time':np.datetime64('2018-05-16T00:00'),
         'dO2': 60,
         'pH':7
         }
obs1 = {'dO2': 60}
obs2 = {'dO2':55}
cell_state = {'VCD': Q(10, 'e5c/ml'),
         'cell_diameter': Q(14, 'um'),
         'viability': 92}

# def graph_cells(parameter1, parameter2):
#   cell_line = cells.gen_cell_line()
#   instance = cell_sim.cell_wrapper(cell_line, Q(1, 'e6c/ml'))
#   for step in range(14*24*60):
#     cells_out = instance.step(env1)
  

def run_tests():
  controls_tests()
  quantities_test()
  assays_test()
  bioreactor_test()
  cell_tests()
  all_tests()
  


  
if __name__ =='__main__':
  run_tests()