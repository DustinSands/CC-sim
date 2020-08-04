# -*- coding: utf-8 -*-
"""
Main script to set up desired config, run simulation, and plot results from

@author: Dustin Sands
"""

import numpy as np
from quantities import Quantity as Q
from matplotlib import pyplot as plt


import assays, cell_sim, controls, helper_functions
from simulator import run_experiments

#Define the media to be used
example_media = {
            'NaHCO3':Q(22, 'mM'), 
            'dO2':Q(0.21, 'mM'), 
           # 'H2CO3':Q(, 'g/L'), 
           # 'iron':Q(, 'g/L'), 
           # 'LDH':Q(, 'g/L'), 
           # 'lactate':Q(, 'g/L'), 
           'glucose':Q(8.39, 'mM'), 
           'amino_acids':Q(70, 'mM'),
           'acetic_acid':Q(3.15, 'mM'),
           'butyric_acid':Q(0.76, 'mM'),
           'citric_acid':Q(0.64, 'mM'),
           'formic_acid':Q(7.6, 'mM'),
           'isovaleric_acid':Q(1.53, 'mM'),
           'lactate':Q(4.99, 'mM'),
           'adenine':Q(0.68, 'mM'),
           'NaCl':Q(70., 'mM'),
           'KCl':Q(40., 'mM')}

# Example concentrated media, e.g. for CFB or Perfusion
example_concentrated_media = {'NaHCO3':Q(22, 'mM'), 
                              'dO2':Q(0.21, 'mM'), 
                              # 'glucose':Q(100, 'mM'), 
                              'amino_acids':Q(700, 'mM'),
                              'lactate':Q(4.99, 'mM'),
                              'NaCl':Q(100., 'mM'),
                              'KCl':Q(80., 'mM')
                              }




        

def create_config(num_experiments, 
                  cell_line = None,
                  glucose_sp = Q(1.5, 'g/L'),
                  seed_density = Q(10, 'e5c/ml'),
                  equipment_instances = None):
  """Create a configuration for an experiment.
  """
  
  start_time = np.datetime64('2020-01-01')
  initial_volume = Q(0.5, 'L')
  seed_density = seed_density
  starting_cells = initial_volume*seed_density
  if cell_line == None:
    cell_line = cell_sim.gen_cell_line()
  if equipment_instances == None:
    equipment_instances = {'osmo':assays.osmo(),
                           'BGA_instance':assays.BGA(), 
                           'cell_counter':assays.cell_counter(),
                           'bioHT':assays.bioHT(),}

  assay_setup = {**equipment_instances,
                 'start_time':start_time,
                 'bioHT_list':['glucose', 'IGG'],  #Same BGA, cell_counter, bioHT instance for all
                 }
  
  # Set up controls.  Each control has its own setup
  glucose_feed = (controls.basic_fed_batch_feed,[],
                    {'feed_mixture':{'glucose': Q(500, 'g/L')}, 
                    'initial_volume':initial_volume,
                    'sample_interval':Q(24, 'h'), 
                    'cpp':'glucose', 
                    'set_point':glucose_sp, 
                    'initial_time': start_time, 
                    'target_seeding_density':seed_density}
                  )
  concentrated_feed = (controls.basic_fed_batch_feed, [], 
                        {'feed_mixture':example_concentrated_media,
                        'initial_volume':initial_volume,
                        'sample_interval':Q(24, 'h'), 
                        'cpp':'mOsm', 
                        'set_point':Q(330, 'mM'), 
                        'initial_time': start_time, 
                        'target_seeding_density':seed_density,
                        'max_added':Q(1., 'L'),}
                      )
  aeration_setup = (controls.aeration,[],
                    {'setpoint':60, 'max_air':Q(0.2, 'L/min'), 
                    'max_O2':Q(0.1, 'L/min'),}
                    )
  control_setup = [glucose_feed,
                   aeration_setup,
                   (controls.temperature, [36], {}),
                   (controls.pH, [7.15], {}),
                   concentrated_feed,
                   ]
  media = helper_functions.create_media_as_mole(example_media, initial_volume)
  br_setup = {'start_time':start_time,
              'initial_components':media,
              'seed_volume':initial_volume,
              'initial_temperature':36,
              }
  cell_setup = {'cell_line':cell_line, 
                'seed_cells':starting_cells}
  config = (assay_setup, control_setup, br_setup, cell_setup)
  return [config]*num_experiments


       
        
def dual_subplot(y11_param, y12_param, y21_param, y22_param, metrics = None,
                 big_title = None, left_title = None, right_title = None):
  """Plots two variables on two each of plots, then puts them next to each
  other."""
  fig = plt.figure(figsize = (10, 4))
  # spec = fig.add_gridspec(2, 1)
  ax1 = fig.add_subplot(121)
  plot_var(ax1, y11_param, metrics, color = 'red')
  
  ax2 = ax1.twinx()
  plot_var(ax2, y12_param, metrics, color = 'blue')
  
  if left_title != None:
    ax1.set_title(left_title)
  
  ax3 = fig.add_subplot(122)
  plot_var(ax3, y21_param, metrics, color = 'red')
  
  ax4 = ax3.twinx()
  plot_var(ax4, y22_param, metrics, color = 'blue')
  
  if right_title != None:
    ax3.set_title(right_title)
  
  if big_title != None:
    fig.suptitle(big_title, y = 1.04)
  fig.tight_layout()
  fig.dpi = 200
  plt.show()
  plt.close()

def dual_plot( y1_param, y2_param, metrics=None):
  """Plot two different metrics on one graph."""
  fig, ax1 = plt.subplots()
  plot_var(ax1, y1_param, metrics = metrics, color = 'red')
  
  ax2 = ax1.twinx()
  plot_var(ax2, y2_param, metrics = metrics, color = 'blue')
  
  fig.tight_layout()
  fig.dpi = 200
  plt.show()
  plt.close()
  
def plot_var(axis, variable, metrics, color = 'black'):
  """Takes an axis to plot on, the metrics, and a variable of interest.  Then
  plots that variable on that axis."""
  if metrics == None:
    metrics = default_metrics
  plot_func = helper_functions.get_plotfunc(axis, variable)
  for metric in metrics:
    tuples = [((metric[point]['time']-metric[0]['time'])/np.timedelta64(1,'D'),
           metric[point][variable])
          for point in range(len(metric)) if variable in metric[point]]
    if len(tuples) > 0:
      x, y = list(zip(*tuples))
      plot_func(x, y, color = color)

  if len(tuples) > 0:
    if max(y)/3> min(y):
      axis.set_ylim(bottom = -max(y)/100, top = None)
    if type(y[0]) == Q:
      axis.set_ylabel(variable+f' {y[0].dimensionality}', color = color)
    else:
      axis.set_ylabel(variable, color = color)

  axis.set_xlabel('Day')
  axis.set_xlim(0, (metric[-1]['time']-metric[0]['time'])/np.timedelta64(1,'D'))
  
    
  
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
  metrics = run_experiments(config)
  return metrics
  
if __name__ == '__main__':
  shared_equipment = {'osmo':assays.osmo(),
                      'BGA_instance':assays.BGA(), 
                      'cell_counter':assays.cell_counter(),
                      'bioHT':assays.bioHT(),}
  
  cell_line = cell_sim.gen_cell_line()
  low_config = create_config(2, 
                             cell_line = cell_line, 
                             equipment_instances = shared_equipment,
                             seed_density = Q(2, 'e5c/ml'),
                             )
  medium_config = create_config(2, 
                             cell_line = cell_line, 
                             equipment_instances = shared_equipment,
                             seed_density = Q(5, 'e5c/ml'),
                             )
  high_config = create_config(2, 
                             cell_line = cell_line, 
                             equipment_instances = shared_equipment,
                             seed_density = Q(10, 'e5c/ml'),
                             )
  low_metrics = run_experiments(low_config)
  
  dual_subplot('viability', 'VCD', 'mOsm', 'IGG', metrics = low_metrics,
               big_title = r'Seed Density: 2 e5c/ml', left_title = 'Viability vs VCD',
               right_title = 'Osmo vs IGG')
  # dual_plot('mOsm', 'mOsm feed rate', metrics = low_metrics)

  # dual_plot('pH', 'pH_PID', metrics = low_metrics)
  # dual_plot('dO2', 'aeration_PID', metrics = low_metrics)
  
  medium_metrics = run_experiments(medium_config)
  dual_subplot('viability', 'VCD', 'mOsm', 'IGG', metrics = medium_metrics,
               big_title = r'Seed Density: 5 e5c/ml', left_title = 'Viability vs VCD',
               right_title = 'Osmo vs IGG')
  high_metrics = run_experiments(high_config)
  dual_subplot('viability', 'VCD', 'mOsm', 'IGG', metrics = high_metrics,
               big_title = r'Seed Density: 10 e5c/ml', left_title = 'Viability vs VCD',
               right_title = 'Osmo vs IGG')
  
  dual_plot('mass', 'mOsm feed rate', metrics = low_metrics)

  dual_plot('mass', 'mOsm feed rate', metrics = high_metrics)
  
  


  
