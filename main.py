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

import actuation, assays, cells, control_strategy, bioreactor, param


def run_experiments(config, duration):
  assay_wrapper = []
  control_wrapper = []
  actuation_wrapper = []
  env_wrapper = []
  cell_wrapper = []
  for assay_setup, control_setup, actuation_setup, br_setup, cell_setup in config:
    assay_wrapper.append(assays.wrapper(**assay_setup))
    control_wrapper.append(control_strategy.wrapper(**control_setup))
    actuation_wrapper.append(actuation.wrapper(**actuation_setup))
    env_wrapper.append(bioreactor.bioreactor(**br_setup))
    cell_wrapper.append(cells.cell_instance(**cell_setup))
  
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
      controls = control_wrapper[br](assays, actuation)
      actuation = actuation_wrapper[br](controls)
      # Environment will increment timestep
      environment = env_wrapper[br](actuation, cells)
      cells = cell_wrapper[br](environment)
      
      #Record important data