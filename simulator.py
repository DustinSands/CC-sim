# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 19:31:43 2020

@author: Racehorse
"""
import math

import numpy as np
from quantities import Quantity as Q

import param
import assays, cell_sim, controls, helper_functions, actuation, bioreactor


"""The simulation needs initial starting conditions for actual and cells to pass to 
modules for the first step."""

initial_actuation = {parameter:Q(0, 'mol/min') for parameter in param.liquid_components}
initial_actuation.update({
  'heat':Q(50, 'W').simplified,
  'RPS':Q(5, '1/s'),
  'air':Q(0.01, 'L/min').simplified,
  'O2':Q(0, 'm**3/s'),
  'CO2':Q(0, 'm**3/s'),
  'gas_volumetric_rate':Q(0.01,'L/min').simplified,
  'liquid_volumetric_rate':Q(0,'m**3/s'),
  })

no_cells = {'mass_transfer':{parameter:Q(0, 'mol/min').simplified for parameter in param.liquid_components},
                'total_cells':Q(0, 'ce').simplified,
                'diameter':Q(15, 'um').simplified,
                'volume':Q(15,'um').simplified**3*math.pi/6,
  }

if param.skip_units:
  helper_functions.remove_units(initial_actuation)
  helper_functions.remove_units(no_cells)

def run_experiments(config, duration = np.timedelta64(14, 'D')):
  assay_wrappers = []
  control_wrappers = []
  actuation_wrappers = []
  env_wrappers = []
  cell_wrappers = []
  for assay_setup, control_setup, br_setup, cell_setup in config:
    instance = assays.wrapper(**assay_setup)
    assay_wrappers.append(instance)
    next_control_wrapper = controls.wrapper(control_setup)
    control_wrappers.append(next_control_wrapper)
    actuation_wrappers.append(actuation.wrapper(next_control_wrapper.actuation_list))
    env_wrappers.append(bioreactor.bioreactor(**br_setup))
    cell_wrappers.append(cell_sim.cell_wrapper(**cell_setup))
  
  steps_per_day = np.timedelta64(1, 'D') / param.resolution
  total_steps = round(duration/param.resolution+helper_functions.gauss(0, steps_per_day*0.05))
  
  
  #Define initial values to pass

  
  actuation_out = [initial_actuation for x in range(len(env_wrappers))]
  cells_output = [no_cells]*len(env_wrappers)
  environment = [[]]*len(env_wrappers)
  obs = [[]]*len(env_wrappers)
  
  online_metric_interval = np.timedelta64(3, 'm')/param.resolution
  day = [-1]*len(env_wrappers)
  days = int(round(duration/np.timedelta64(1, 'D')))
  next_offline_step = [0]*len(env_wrappers)
  metrics = [[]for y in range(len(env_wrappers))]
  
  
  #1.5 hours of equilibration - cell wrapper is skipped
  for step in range(int(round(steps_per_day / 16))):
    for br in range(len(env_wrappers)):
      environment[br] = env_wrappers[br].step(actuation_out[br], cells_output[br])
      # Run simulation for a step
      obs[br] = assay_wrappers[br].step(environment[br], cells_output[br], False)
      #This will update setpoints for actuation, so actuation doesn't need an 
      #input.
      control_metrics = control_wrappers[br].step(obs[br], False)
      actuation_out[br] = actuation_wrappers[br].step()
    
  
  for step in range(int(round(total_steps))):
    for br in range(len(env_wrappers)):
      # print(f'step {step}!')
      
      if step == next_offline_step[br]:
        offline = True
        day[br] += 1
        print(f'Day {day[br]}!')
        if day[br]+1 == days:
          next_offline_step[br] = total_steps - 1
        else:
          next_offline_step[br] = round(steps_per_day*(day[br]+1) +\
            helper_functions.gauss(0,0.05*steps_per_day))
      else:
        offline = False
      environment[br] = env_wrappers[br].step(actuation_out[br], cells_output[br])
      cells_output[br], cheater_metrics = cell_wrappers[br].step(environment[br])
      # Run simulation for a step
      obs[br] = assay_wrappers[br].step(environment[br], cells_output[br], offline)
      #This will update setpoints for actuation, so actuation doesn't need an 
      #input.
      control_metrics = control_wrappers[br].step(obs[br], offline)
      actuation_out[br] = actuation_wrappers[br].step()
      # Environment will increment timestep
      
      # cells_output = cell_wrapper[br](environment)
      
      if offline == True or step%online_metric_interval < 1:
        step_metrics = {**control_metrics,
                        **obs[br],
                        }
        metrics[br].append(helper_functions.scale_units(step_metrics))
        #CHEATER METRICS
        if not param.realistic_mode:
          cheater_metrics.update({
                             'amino_acids':environment[br]['amino_acids'],
                             'dCO2':environment[br]['dCO2'],
                             'rVCD':cells_output[br]['living_cells']/environment[br]['volume'],
                             'rglucose':environment[br]['glucose'],
                             'compA':environment[br]['component_A'],
                             'total_cell_volume':cells_output[br]['living_cells']*cells_output[br]['volume'],
                             'total_living_cells':cells_output[br]['living_cells'],
                             })
          metrics[br][-1].update(helper_functions.scale_units(cheater_metrics))
      if offline == True:
        # dual_plot('pH', 'pH_PID', metrics[br])
        # dual_plot('dO2', 'aeration_PID', metrics[br])
        pass

  
  return metrics