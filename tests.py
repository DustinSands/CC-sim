# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 06:18:03 2020

@author: sckuz
"""
import random

import numpy as np
import quantities as q
from quantities import Quantity as Q

import control_strategy, assays
import param

def set_seeds(seed):
  np.random.seed(seed)
  random.seed(seed)

def feeding_strategy_tests():
  control_strategy.fed_batch(control_variable = 'glucose',
                             initial_volume = 2, sample_interval = 1440, 
                             cpp = 'glucose', set_point = 2,
                             time = np.datetime64('2020-01-01'), target_seeding_density = 10)
  
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
  set_seeds(0)
  BGA = assays.BGA()
  oxygen = assays.O2_probe(np.datetime64('2018-01-01T00:00'), cal_reference = BGA)
  pH = assays.pH_probe(np.datetime64('2018-01-01T00:00'), cal_reference = BGA)
  state = {'time':np.datetime64('2018-01-05T00:00'),
           'O2': 60,
           'pH':7
           }
  state2 = {'time':np.datetime64('2018-05-16T00:00'),
           'O2': 60,
           'pH':7
           }
  cells = {'VCD': Q(10, 'e5c/ml')}
  oxygen.one_point(np.datetime64('2018-01-01T00:00'))
  
  print(oxygen.read_value(state, cells))
  # assert 61.192538759595784 == oxygen.read_value(state, cells)
  print(oxygen.read_value(state2, cells))
  # assert 59.07032195548049==oxygen.read_value(state2, cells)
  print(pH.read_value(state, cells))
  # assert 7.000070542379912 == pH.read_value(state, cells)
  print(pH.read_value(state2, cells))
  # assert 7.228386164133978==pH.read_value(state2, cells)
  
  vicell = assays.cell_counter()
  # print(vicell.read_value(state))
  assert Q(10.72273548515264, 'e5c/ml') == vicell.read_value(cells)
  
  print('passed assays')

  

def run_tests():
  feeding_strategy_tests()
  quantities_test()
  assays_test()

  
if __name__ =='__main__':
  run_tests()