# -*- coding: utf-8 -*-
"""
Module that runs the experiments.  

-Creates wrapper for each module for each experiment
-Initializes experiment before adding cells
-Runs experiment and periodically stores metrics
--Optionally stores cheater metrics as well
-Decides when to take offline samples
"""
import math, pdb

import numpy as np
from quantities import Quantity as Q

import param
import assays, cell_sim, controls, helper_functions, actuation, bioreactor


timers = {tracker:helper_functions.time_tracker() for \
          tracker in ['total', 'env', 'cells', 'assays', 'con', 'metrics', 'root']}
helper_functions.timer = timers


"""The simulation needs initial starting conditions for actuation and cells to pass to 
modules for the first step."""

if param.skip_units:
  flow_func = lambda R: 100.*R
else: flow_func = lambda q: Q(100., 'kg/(m**4*s)')*q
initial_actuation = {parameter:Q(0, 'mol/min') for parameter in param.liquid_components}
initial_actuation.update({
  'heat':Q(50, 'W').simplified,
  'RPS':Q(5, '1/s'),
  'air':Q(0.01, 'L/min').simplified,
  'O2':Q(0, 'm**3/s'),
  'CO2':Q(0, 'm**3/s'),
  'gas_volumetric_rate':Q(0.01,'L/min').simplified,
  'liquid_volumetric_rate':Q(0,'m**3/s'),
  'recirc_func': flow_func,
  'flux_func': flow_func,
  'recirc_shear': Q(0, '1/s'),
  'recirc_RPM':Q(0, '1/s'),
  })
no_cells = {'mass_transfer':{parameter:Q(0, 'mol/min').simplified for parameter in param.liquid_components},
            'total_cells':Q(0, 'ce').simplified,
            'diameter':Q(15, 'um').simplified,
            'volume':Q(15,'um').simplified**3*math.pi/6,
            'dry_weight':Q(0, 'g').simplified,
  }

if param.skip_units:
  helper_functions.remove_units(initial_actuation)
  helper_functions.remove_units(no_cells)

def run_experiments(config, duration = np.timedelta64(14, 'D')):
  """Takes the experimental setup and (optionally) run duration, runs the 
  experiments, and returns metrics.  Metrics is a list (by time), each point 
  being a dict of metrics."""
  assay_wrappers = []
  control_wrappers = []
  actuation_wrappers = []
  env_wrappers = []
  cell_wrappers = []
  
  #Initialize wrappers for each module
  for assay_setup, control_setup, br_setup, cell_setup in config:
    instance = assays.wrapper(**assay_setup)
    assay_wrappers.append(instance)
    next_control_wrapper = controls.wrapper(control_setup)
    control_wrappers.append(next_control_wrapper)
    actuation_wrappers.append(actuation.wrapper(next_control_wrapper.actuation_list))
    env_wrappers.append(bioreactor.bioreactor(**br_setup))
    cell_wrappers.append(cell_sim.cell_wrapper(**cell_setup))
  
  #Calculate time information
  steps_per_day = np.timedelta64(1, 'D') / param.resolution
  total_steps = round(duration/param.resolution+helper_functions.gauss(0, steps_per_day*0.05))
  
  #Define initial values to pass
  actuation_out = [initial_actuation for x in range(len(env_wrappers))]
  cells_output = [no_cells]*len(env_wrappers)
  environment = [[]]*len(env_wrappers)
  obs = [[]]*len(env_wrappers)
  
  online_metric_interval = np.timedelta64(30, 'm')/param.resolution
  day = [-1]*len(env_wrappers)
  days = int(round(duration/np.timedelta64(1, 'D')))
  next_offline_step = [0]*len(env_wrappers)
  metrics = [[]for y in range(len(env_wrappers))]
  if param.skip_units:
    cell = 1
  else:
    cell = param.ce
  
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
    
  timers['total'].start()
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
      # Run simulation for a step
      timers['env'].start()
      environment[br] = env_wrappers[br].step(actuation_out[br], cells_output[br])
      timers['env'].stop()
      timers['cells'].start()
      cells_output[br], cheater_metrics = cell_wrappers[br].step(environment[br])
      timers['cells'].stop()
      timers['assays'].start()
      obs[br] = assay_wrappers[br].step(environment[br], cells_output[br], offline)
      timers['assays'].stop()
      timers['con'].start()
      #This will update setpoints for actuation, so actuation doesn't need an 
      #input.
      control_metrics = control_wrappers[br].step(obs[br], offline)
      actuation_out[br] = actuation_wrappers[br].step()
      timers['con'].stop()
      timers['metrics'].start()
      if offline == True or step%online_metric_interval < 1:
        step_metrics = {**control_metrics,
                        **obs[br],
                        }
        metrics[br].append(step_metrics)
        #Cheater metrics; not possible in real life
        if not param.realistic_mode:
          cheater_metrics.update({
                              'amino_acids':environment[br]['br_molarity']['amino_acids'],
                              'dCO2':environment[br]['br_molarity']['dCO2'],
                              'rVCD':cheater_metrics['viable_cells']/environment[br]['volume'],
                              'rglucose':environment[br]['br_molarity']['glucose'],
                              'rOsmo':environment[br]['br_molarity']['mOsm'],
                              'compA':environment[br]['br_molarity']['component_A'],
                              'compB':environment[br]['br_molarity']['component_B'],
                              'total_cell_volume':cells_output[br]['total_cells']*cells_output[br]['volume']/cell,
                              'cell_fraction':environment[br]['cell_fraction'],
                              })
          if 'recirc_rate' in environment[br]:
            cheater_metrics.update({'fiber_shear': environment[br]['fiber_shear'],
                                    'TMP': environment[br]['TMP'],
                                    'rperfusion':environment[br]['perfusion_flowrate'],
                                    })
          metrics[br][-1].update(cheater_metrics)
      timers['metrics'].stop()
  timers['total'].stop()
  helper_functions.print_times()
  global latest
  latest = metrics
  return metrics
