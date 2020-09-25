# -*- coding: utf-8 -*-
"""
Main script to set up desired config, run simulation, and plot results from

@author: Dustin Sands
"""
import pdb

import numpy as np
from quantities import Quantity as Q
from matplotlib import pyplot as plt

import assays, cell_sim, controls, helper_functions, bioreactor
from simulator import run_experiments
from analysis import dual_plot, dual_subplot, analyze

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
example_perfusion_media = {'NaHCO3':Q(22, 'mM'), 
                           'dO2':Q(0.21, 'mM'), 
                           'glucose':Q(110, 'mM'), 
                           'amino_acids':Q(180, 'mM'),
                           'lactate':Q(4.99, 'mM'),
                           'NaCl':Q(70., 'mM'),
                           'KCl':Q(40., 'mM'),
                           }
# example_media_concentrate = {'NaHCO3':Q(88, 'mM'), 
#                               'dO2':Q(0.21, 'mM'), 
#                               'glucose':Q(50, 'mM'), 
#                               'amino_acids':Q(700, 'mM'),
#                               'lactate':Q(20, 'mM'),
#                               'NaCl':Q(200., 'mM'),
#                               'KCl':Q(150., 'mM')
#                               }

def create_perfusion_config(num_experiments, 
                            cell_line = None,
                            glucose_sp = Q(3, 'g/L'),
                            seed_density = Q(10, 'e6c/ml'),
                            equipment_instances = None):
  """Config creator for perfusion experiments."""
  start_time = np.datetime64('2020-01-01')
  initial_volume = Q(2, 'L')
  starting_cells = initial_volume*seed_density
  if cell_line == None:
    cell_line = cell_sim.gen_cell_line()
  if equipment_instances == None:
    #Same BGA, cell_counter, bioHT instance for all
    equipment_instances = {'osmo':assays.osmo(),
                           'BGA_instance':assays.BGA(), 
                           'cell_counter':assays.cell_counter(),
                           'bioHT':assays.bioHT(),}

  #Set up assays to be run
  assay_setup = {**equipment_instances,
                 'start_time':start_time,
                 'bioHT_list':['glucose', 'IGG'],  
                 'recirc_flowmeter':True, #flowmeter
                 'levitronix_RPM': True, #levitronix RPM
                 }
  
  # Set up controls.  Each control has its own setup
  perfusion_feed = (controls.basic_perfusion_feed, [], 
                    {'feed_mixture':example_perfusion_media,
                    'initial_volume':initial_volume,
                    'sample_interval':Q(24, 'h'), 
                    'cpp':'mOsm', 
                    'set_point':Q(350, 'mM'), 
                    'initial_time': start_time, 
                    'target_seeding_density':seed_density,
                    'max_added':Q(50, 'L'),
                    'initial_VVD':Q(0.4, '1/d'),
                    'glucose_setpoint':glucose_sp,}
                    )
  aeration_setup = (controls.aeration,[],
                    {'setpoint':60, 'max_air':Q(0.2, 'L/min'), 
                    'max_O2':Q(0.3, 'L/min'),}
                    )
  control_setup = [perfusion_feed,
                   aeration_setup,
                   (controls.temperature, [36], {}),
                   (controls.pH, [7.15], {}),
                   (controls.recirculation,[Q(1, 'L/min')],{}),
                   ]
  media = helper_functions.create_media_as_mole(example_media, initial_volume)
  br_setup = {'start_time':start_time,
              'initial_components':media,
              'seed_volume':initial_volume,
              'initial_temperature':36,
              'cell_separation_device':bioreactor.PES()
              }
  cell_setup = {'cell_line':cell_line, 
                'seed_cells':starting_cells}
  config = (assay_setup, control_setup, br_setup, cell_setup)
  return [config]*num_experiments

def create_batch_config(num_experiments, 
                  cell_line = None,
                  glucose_sp = Q(1.5, 'g/L'),
                  seed_density = Q(10, 'e5c/ml'),
                  equipment_instances = None):
  """Config creator for fed-batch experiments."""
  start_time = np.datetime64('2020-01-01')
  initial_volume = Q(0.5, 'L')
  seed_density = seed_density
  starting_cells = initial_volume*seed_density
  if cell_line == None:
    cell_line = cell_sim.gen_cell_line()
  if equipment_instances == None:
    #Same BGA, cell_counter, bioHT instance for all
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

def run_seed_density_comparison():
  """As an example run, different seed densities in fed batch were compared.
  This function creats and runs that example."""
  #Each config uses the same shared equipment
  shared_equipment = {'osmo':assays.osmo(),
                      'BGA_instance':assays.BGA(), 
                      'cell_counter':assays.cell_counter(),
                      'bioHT':assays.bioHT(),}
  #Each experiment needs to use the same cell line
  cell_line = cell_sim.gen_cell_line()
  #Create configs for each seed density (two runs each)
  low_config = create_batch_config(2, 
                             cell_line = cell_line, 
                             equipment_instances = shared_equipment,
                             seed_density = Q(2, 'e5c/ml'),
                             )
  medium_config = create_batch_config(2, 
                             cell_line = cell_line, 
                             equipment_instances = shared_equipment,
                             seed_density = Q(5, 'e5c/ml'),
                             )
  high_config = create_batch_config(2, 
                             cell_line = cell_line, 
                             equipment_instances = shared_equipment,
                             seed_density = Q(10, 'e5c/ml'),
                             )
  #Run and plot each configuration
  low_metrics = run_experiments(low_config)
  analyze(low_metrics)
  dual_subplot('viability', 'VCD', 'mOsm', 'IGG', metrics = low_metrics,
               big_title = r'Seed Density: 2 e5c/ml', left_title = 'Viability vs VCD',
               right_title = 'Osmo vs IGG')
  
  medium_metrics = run_experiments(medium_config)
  analyze(medium_metrics)
  dual_subplot('viability', 'VCD', 'mOsm', 'IGG', metrics = medium_metrics,
                big_title = r'Seed Density: 5 e5c/ml', left_title = 'Viability vs VCD',
                right_title = 'Osmo vs IGG')
  
  high_metrics = run_experiments(high_config)
  analyze(high_metrics)
  dual_subplot('viability', 'VCD', 'mOsm', 'IGG', metrics = high_metrics,
                big_title = r'Seed Density: 10 e5c/ml', left_title = 'Viability vs VCD',
                right_title = 'Osmo vs IGG')
  
  #Compare final mass of highest and lowest config
  dual_plot('mass', 'mOsm feed sp', metrics = low_metrics)
  dual_plot('mass', 'mOsm feed sp', metrics = high_metrics)


def run_perfusion():
  """Creates a perfusion configuration and runs it, then plots a couple results."""
  config = create_perfusion_config(1, glucose_sp = Q(3, 'g/L'))
  metrics = run_experiments(config)
  analyze(metrics)

  dual_subplot('VCD', 'viability', 'glucose', 'IGG')
  dual_plot('sieving', 'permeate rate')
  return metrics

if __name__ == '__main__':
  metrics = run_perfusion()


  
  


  
