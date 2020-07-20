# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:10:00 2020

@author: Racehorse
"""




import random
import math
import pdb

import numpy as np
from quantities import Quantity as Q
import quantities as q
from matplotlib import pyplot as plt

import param
import actuation, assays, cell_sim, controls, bioreactor, helper_functions

example_media = {
            'NaHCO3':Q(22., 'mM'), 
            'dO2':Q(0.21, 'mM'), 
           # 'H2CO3':Q(, 'g/L'), 
           # 'iron':Q(, 'g/L'), 
           # 'LDH':Q(, 'g/L'), 
           # 'lactate':Q(, 'g/L'), 
           'glucose':Q(8.39, 'mM'), 
           'amino_acids':Q(150, 'mM'),
           'acetic_acid':Q(3.15, 'mM'),
           'butyric_acid':Q(0.76, 'mM'),
           'citric_acid':Q(0.64, 'mM'),
           'formic_acid':Q(7.6, 'mM'),
           'isovaleric_acid':Q(1.53, 'mM'),
           'lactate':Q(4.99, 'mM'),
           'adenine':Q(0.68, 'mM'),
           'NaCl':Q(40., 'mM'),
           'KCl':Q(50., 'mM')}

"""The simulation needs initial starting conditions for actual and cells to pass to 
modules for the first step."""

initial_actuation = {parameter:Q(0, 'mol/min') for parameter in param.liquid_components}
initial_actuation.update({
  'heat':Q(50, 'W').simplified,
  'RPS':Q(5, '1/s').simplified,
  'air':Q(0.01, 'L/min').simplified,
  'O2':Q(0, 'L/min').simplified,
  'CO2':Q(0, 'L/min').simplified,
  'gas_volume':Q(0.01,'L/min').simplified,
  'liquid_volume':Q(0,'L/min').simplified,
  })

no_cells = {'mass_transfer':{parameter:Q(0, 'mol/min').simplified for parameter in param.liquid_components},
                'total_cells':Q(0, 'ce').simplified,
                'diameter':Q(15, 'um').simplified,
                'volume':Q(15,'um').simplified**3*math.pi/6,
  }

if param.skip_units:
  helper_functions.remove_units(initial_actuation)
  helper_functions.remove_units(no_cells)


        

def create_config(num_experiments):
  start_time = np.datetime64('2020-01-01')
  initial_volume = Q(0.5, 'L')
  seed_density = Q(1, 'e6c/ml')
  starting_cells = (initial_volume*seed_density).simplified
  random.seed(0)
  cell_line = cell_sim.gen_cell_line()
  assay_setup = [assays.osmo(),
                 assays.BGA(), 
                 assays.cell_counter(),
                 start_time,
                 assays.bioHT(),
                 ['glucose', 'IGG']]  #Same BGA, cell_counter, bioHT instance for all
  fed_batch_setup = {'feed_mixture':{'glucose': Q(500, 'g/L')}, 
                   'initial_volume':initial_volume,
                   'sample_interval':Q(24, 'h'), 
                   'cpp':'glucose', 
                   'set_point':Q(1.5, 'g/L'), 
                   'initial_time': start_time, 
                   'target_seeding_density':Q(10, 'e5c/ml')}
  concentrated_media = {component:2*concentration for component, concentration
                        in example_media.items()}
  secondary_feed_setup = {'feed_mixture':concentrated_media,
                          'initial_volume':initial_volume,
                          'sample_interval':Q(24, 'h'), 
                          'cpp':'mOsm', 
                          'set_point':Q(300, 'mM'), 
                          'initial_time': start_time, 
                          'target_seeding_density':Q(10, 'e5c/ml')}
  aeration_setup = {'setpoint':60, 'max_air':Q(0.2, 'L/min'), 
               'max_O2':Q(0.1, 'L/min')}
  control_setup = [(controls.basic_fed_batch_feed,[],fed_batch_setup),
                   (controls.aeration,[],aeration_setup),
                   (controls.temperature, [36], {}),
                   (controls.pH, [7.15], {}),
                   (controls.basic_fed_batch_feed, [], secondary_feed_setup),
                   ]

  media = helper_functions.create_media_as_mole(example_media, initial_volume)
  br_setup = [start_time, media, initial_volume, 36]
  cell_setup = [cell_line, starting_cells]
  config = (assay_setup, control_setup, br_setup, cell_setup)
  return [config]*num_experiments

def run_experiments(config, duration):
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
    cell_wrappers.append(cell_sim.cell_wrapper(*cell_setup))
  
  steps_per_day = np.timedelta64(1, 'D') / param.resolution
  total_steps = round(duration/param.resolution+random.gauss(0, steps_per_day*0.05))
  
  
  #Define initial values to pass

  
  actuation_out = [initial_actuation for x in range(len(env_wrappers))]
  cells_output = [no_cells]*len(env_wrappers)
  environment = [[]]*len(env_wrappers)
  obs = [[]]*len(env_wrappers)
  
  online_metric_interval = np.timedelta64(10, 'm')/param.resolution
  day = [-1]*len(env_wrappers)
  days = int(round(duration/np.timedelta64(1, 'D')))
  next_offline_step = [5]*len(env_wrappers)
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
            random.gauss(0,0.05*steps_per_day))
      else:
        offline = False
      environment[br] = env_wrappers[br].step(actuation_out[br], cells_output[br])
      cells_output[br] = cell_wrappers[br].step(environment[br])
      # Run simulation for a step
      obs[br] = assay_wrappers[br].step(environment[br], cells_output[br], offline)
      #This will update setpoints for actuation, so actuation doesn't need an 
      #input.
      control_metrics = control_wrappers[br].step(obs[br], offline)
      actuation_out[br] = actuation_wrappers[br].step()
      # Environment will increment timestep
      
      # cells_output = cell_wrapper[br](environment)
      
      if offline == True or step%online_metric_interval < 1:
        metrics[br].append(helper_functions.scale_units(control_metrics))
        metrics[br][-1].update(helper_functions.scale_units(obs[br]))
        #CHEATER METRICS
        metrics[br][-1].update({'dCO2':environment[br]['dCO2']})
        metrics[br][-1].update({'amino_acids':environment[br]['amino_acids']})
      if offline == True:
        # dual_plot('pH', 'pH_PID', metrics[br])
        # dual_plot('dO2', 'aeration_PID', metrics[br])
        pass

  
  return metrics
       
        
      

def dual_plot( y1_param, y2_param, metrics=None, total_days=14):
  """Plot two different metrics on one graph."""
  if metrics == None:
    metrics = default_metrics[0]
  fig, ax1 = plt.subplots()
  tuples = [((metrics[point]['time']-metrics[0]['time'])/np.timedelta64(1,'D'),
         metrics[point][y1_param])
        for point in range(len(metrics)) if y1_param in metrics[point]]
  if len(tuples) > 0:
    x1, y1 = list(zip(*tuples))
  else:
    x1, y1 = [], []
  ax1.set_xlabel('Day')
  if len(y1)>0 and type(y1[0]) == Q:
    ax1.set_ylabel(y1_param+f' {y1[0].dimensionality}', color = 'red')
  else:
    ax1.set_ylabel(y1_param, color = 'red')
  ax1.plot(x1, y1, color = 'red')
  ax1.set_xlim([0, total_days])
  if max(y1)/2> min(y1):
    ax1.set_ylim(-max(y1)/100)

  tuples = [((metrics[point]['time']-metrics[0]['time'])/np.timedelta64(1,'D'),
          metrics[point][y2_param])
        for point in range(len(metrics)) if y2_param in metrics[point]]
  if len(tuples) > 0:
    x2, y2 = list(zip(*tuples))
  else:
    x2, y2 = [], []
  ax2 = ax1.twinx()
  if len(y2)>0 and type(y2[0]) ==Q:
    ax2.set_ylabel(y2_param+f' {y2[0].dimensionality}', color = 'blue')
  else:
    ax2.set_ylabel(y2_param, color = 'blue')
  ax2.plot(x2, y2, color='blue')
  if max(y2)/2>min(y2):
    ax2.set_ylim(-max(y2)/100)
  fig.tight_layout()
  fig.dpi = 200
  plt.show()
  plt.close()
  
def multireactor_plot(param, metrics = None, total_days = 14):
  """Plots all experiments on the same plot for a single parameter"""
  if metrics == None:
    metrics = default_metrics
  fig, ax1 = plt.subplots()
  for experiment in range(len(metrics)):
    tuples = [((metrics[experiment][point]['time']-metrics[0]['time'])/np.timedelta64(1,'D'),
           metrics[experiment][point][y1_param])
          for point in range(len(metrics)) if y1_param in metrics[experiment][point]]
    if len(tuples) > 0:
      x, y = list(zip(*tuples))
    else:
      x, y = [], []
    ax1.plot(x1, y1, color = 'red')
  ax1.set_xlabel('Day')
  if len(y1)>0 and type(y1[0]) == Q:
    ax1.set_ylabel(y1_param+f' {y1[0].dimensionality}', color = 'red')
  else:
    ax1.set_ylabel(y1_param, color = 'red')
 
  ax1.set_xlim([0, total_days])

def run_sim():
  config = create_config(1)
  metrics = run_experiments(config, np.timedelta64(14, 'D'))
  return metrics
  
if __name__ == '__main__':
  default_metrics = run_sim()
  dual_plot('viability','VCD' )
  dual_plot('mOsm', 'cell_diameter')
  dual_plot('pH', 'pH_PID')
  dual_plot('dO2', 'aeration_PID')
  dual_plot('glucose', 'glucose addition rate')
  dual_plot('mass', 'mOsm feed rate')
  dual_plot('mOsm', 'mOsm feed rate')
  dual_plot('dCO2', 'amino_acids')
