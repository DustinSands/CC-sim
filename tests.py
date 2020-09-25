# -*- coding: utf-8 -*-
"""
Regression Testing
"""
import bioreactor, cell_sim, controls, actuation, assays, helper_functions, simulator, param
example_env =

example_cells = 

example_act = 

example_obs = 


def env_tests():
  
def cell_tests():
  
def controls_tests():
  
def actuation_tests():
  
def assays_tests():
  
def sim_tests():
  
def helper_tests():
  
def run_all():
  env_test()
  cell_tests()
  controls_tests()
  actuation_tests()
  assays_tests()
  sim_tests()
  helper_tests()
  
if __name__ == '__main__':
  param.debug = 1
  param.skip_units = 1
  run_all()
  print('All skip units tests passed!')
  param.skip_units = 0
  run_all()
  print('Use units tests passed!')
  
