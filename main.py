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
            'NaHCO3':Q(24., 'mM'), 
            'dO2':Q(0.21, 'mM'), 
           # 'H2CO3':Q(, 'g/L'), 
           # 'iron':Q(, 'g/L'), 
           # 'LDH':Q(, 'g/L'), 
           # 'lactate':Q(, 'g/L'), 
           'glucose':Q(8.39, 'mM'), 
           'amino_acids':Q(61.55, 'mM'),
           'amino_acids':Q(150, 'mM'),
           'acetic_acid':Q(3.15, 'mM'),
           'butyric_acid':Q(0.76, 'mM'),
           'citric_acid':Q(0.64, 'mM'),
           'formic_acid':Q(7.6, 'mM'),
           'isovaleric_acid':Q(1.53, 'mM'),
           'lactate':Q(4.99, 'mM'),
           'adenine':Q(0.68, 'mM'),
           'NaCl':Q(50., 'mM'),
           'KCl':Q(70., 'mM')}

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
  random.seed(1)
  cell_line = cell_sim.gen_cell_line()
  assay_setup = [assays.BGA(), 
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
                   'target_seeding_density':Q(1, 'e5c/ml')}
  aeration_setup = {'setpoint':60, 'max_air':Q(0.2, 'L/min'), 
               'max_O2':Q(0.1, 'L/min')}
  control_setup = [(controls.basic_fed_batch_feed,[],fed_batch_setup),
                   (controls.aeration,[],aeration_setup),
                   (controls.temperature, [36], {}),
                   (controls.pH, [7.15], {})]

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
  
  
  day = [-1]*len(env_wrappers)
  days = int(round(duration/np.timedelta64(1, 'D')))
  next_offline_step = [5]*len(env_wrappers)
  metrics = [[{} for x in range(1+days)]for y in range(len(env_wrappers))]
  
  #1 hour of equilibration - cell wrapper is skipped
  for step in range(int(round(steps_per_day / 24))):
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
      
      if offline == True:
        metrics[br][day[br]].update(helper_functions.scale_assays(control_metrics))
        metrics[br][day[br]].update(helper_functions.scale_assays(obs[br]))
        dual_plot(metrics[br], 'viability', 'VCD')
        dual_plot(metrics[br], 'glucose', 'glucose_feed')
  
  return metrics
       
        
      
      #Record important data
def dual_plot(metrics, y1_param, y2_param, total_days=14):
  print(metrics)
  x = [(metrics[day]['time']-metrics[0]['time'])/np.timedelta64(1,'D')
       for day in range(len(metrics)) if 'time' in metrics[day]]
  y1 = [metrics[day][y1_param] for day in range(len(metrics)) if y1_param in metrics[day]]
  y2 = [metrics[day][y2_param] for day in range(len(metrics)) if y2_param in metrics[day]]
  print(x, y1, y2)
  #Plot two sets of data
  fig, ax1 = plt.subplots()
  ax1.set_xlabel('Day')
  if type(y1[0]) ==Q:
    ax1.set_ylabel(y1_param+f' {y1[0].dimensionality}', color = 'red')
  else:
    ax1.set_ylabel(y1_param, color = 'red')

  ax1.plot(x, y1, color = 'red')
  ax1.set_xlim([0, total_days])
  # ax1.set_ylim(top=80)

  ax2 = ax1.twinx()
  # ax2.set_ylim(bottom=-1)
  if type(y2[0]) ==Q:
    ax2.set_ylabel(y2_param+f' {y2[0].dimensionality}', color = 'blue')
  else:
    ax2.set_ylabel(y2_param, color = 'blue')
  ax2.plot(x, y2, color='blue')
  fig.tight_layout()
  fig.dpi = 200
  plt.show()
  plt.close()


def run_sim():
  config = create_config(1)
  metrics = run_experiments(config, np.timedelta64(14, 'D'))
  return metrics
  
if __name__ == '__main__':
  metrics = run_sim()
