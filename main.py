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

import actuation, assays, cells, controls, bioreactor, param

def create_config(num_experiments):
  start_time = np.datetime64('2020-01-01')
  assay_setup = [assays.BGA(), start_time]
  fed_batch_setup = [{'glucose', Q(500, 'g/L')}, np.timedelta64(24, 'h'), 
                   'glucose', Q(2, 'g/L'), start_time, Q(1, 'e5c/ml')]
  aeration_setup = {'setpoint':60, 'max_air':Q(0.2, 'L/min'), 
               'max_O2':Q(0.1, 'L/min')}]
  control_setup = [controls.fed_batch_feed(*fed_batch_setup),
                   controls.DH_aeration(**aeration_setup)]
  
  br_setup = [start_time]
  cell_setup = []

def run_experiments(config, duration):
  assay_wrappers = []
  control_wrappers = []
  actuation_wrappers = []
  env_wrappers = []
  cell_wrappers = []
  for assay_setup, control_setup, br_setup, cell_setup in config:
    assay_wrappers.append(assays.wrapper(**assay_setup))
    next_control_wrapper = controls.wrapper(**control_setup)
    control_wrappers.append(next_control_wrapper)
    actuation_wrappers.append(actuation.wrapper(next_control_wrapper.actuation_list))
    env_wrappers.append(bioreactor.bioreactor(**br_setup))
    cell_wrappers.append(cells.cell_instance(**cell_setup))
  
  total_steps = duration / param.resolution
  
  
  day = 0
  next_offline_step = 0

  for step in range(total_steps):
    for br in range(len(env_wrapper)):
      if step == next_offline_step:
        offline = True
        day += 1
        next_offline_step = day / param.resolution +\
          random.gauss()
      else:
        offline = False
      
      # Run simulation for a step
      assays = assay_wrapper[br].step(environment, cells, offline)
      #This will update setpoints for actuation, so actuation doesn't need an 
      #input.
      control_metrics = control_wrapper[br](assays, offline)
      actuation = actuation_wrapper[br]()
      # Environment will increment timestep
      environment = env_wrapper[br](actuation, cells)
      cells = cell_wrapper[br](environment)
      
      #Record important data