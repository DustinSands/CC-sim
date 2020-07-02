# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:10:00 2020

@author: Racehorse
"""


#Initialize Global Cell Line Params


#Define experimental setup

#Initialize current states by experimental setup
# Need: environment, cells, actuation

import random

import numpy as np
from quantities import Quantity as Q

import actuation, assays, cells, controls, bioreactor, param

def create_config(num_experiments):
  start_time = np.datetime64('2020-01-01')
  initial_volume = :Q(0.5, 'L')
  seed_density = Q(1, 'e6c/ml')
  starting_cells = (initial_volume*seed_density).simplified
  cell_line = cells.gen_cell_line()
  assay_setup = [assays.BGA(), start_time]  #Same BGA instance for all
  fed_batch_setup = {'feed_mixture':{'glucose': Q(500, 'g/L')}, 
                   'initial_volume':initial_volume,
                   'sample_interval':np.timedelta64(24, 'h'), 
                   'cpp':'glucose', 
                   'set_point':Q(2, 'g/L'), 
                   'initial_time': start_time, 
                   'target_seeding_density':Q(1, 'e5c/ml')}
  aeration_setup = {'setpoint':60, 'max_air':Q(0.2, 'L/min'), 
               'max_O2':Q(0.1, 'L/min')}
  control_setup = [(controls.fed_batch_feed,[],fed_batch_setup),
                   (controls.DH_aeration,[],aeration_setup)]
  
  br_setup = [start_time]
  cell_setup = [cell_line, starting_cells]
  config = (assay_setup, control_setup, br_setup, cell_setup)
  return [config]*num_experiments

def run_experiments(config, days):
  assay_wrappers = []
  control_wrappers = []
  actuation_wrappers = []
  env_wrappers = []
  cell_wrappers = []
  for assay_setup, control_setup, br_setup, cell_setup in config:
    instance = assays.wrapper(*assay_setup)
    assay_wrappers.append(instance)
    next_control_wrapper = controls.wrapper(control_setup)
    control_wrappers.append(next_control_wrapper)
    actuation_wrappers.append(actuation.wrapper(next_control_wrapper.actuation_list))
    env_wrappers.append(bioreactor.bioreactor(*br_setup))
    cell_wrappers.append(cells.cell_wrapper(*cell_setup))
  
  steps_per_day = np.timedelta64(1, 'D') / param.resolution
  total_steps = days*steps_per_day+random.gauss(0, steps_per_day*0.05)
  
  
  day = 0
  next_offline_step = 0

  for step in range(total_steps):
    for br in range(len(env_wrapper)):
      if step == next_offline_step:
        offline = True
        day += 1
        print(f'Day {day})!')
        if day == days:
          next_offline_step = total_steps - 1
        else:
          next_offline_step = steps_per_day*day +\
            random.gauss(0,0.05*steps_per_day)
      else:
        offline = False
      
      # Run simulation for a step
      obs = assay_wrapper[br].step(environment, cells, offline)
      #This will update setpoints for actuation, so actuation doesn't need an 
      #input.
      control_metrics = control_wrapper[br](obs, offline)
      actuation_out = actuation_wrapper[br]()
      # Environment will increment timestep
      environment = env_wrapper[br](actuation_out, cells_output)
      cells_output = cell_wrapper[br](environment)
      
      #Record important data
      
if __name__ == '__main__':
  config = create_config(2)
  run_experiments(config, np.timedelta64(14, 'D'))