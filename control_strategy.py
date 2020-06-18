# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:30:40 2020

@author: Racehorse
"""
import random

import tests



class feed_strategy:
  """Base class for feed strategy.
  
  cpp is the critical process parameter that the feed is attemping to control.
  
  volume is bioreactor volume
  
  sample_interval is the expected time between offline samples
  
  set_point is the setpoint to maintain the cpp at
  
  addition_rate is the desired setpoint of addition of feed media
  """
  def __init__(self, control_variable, initial_volume, sample_interval, cpp, set_point, time, target_seeding_density):
    self.initial_volume = initial_volume
    self.seeding_density = target_seeding_density
    self.seed_time = time
    self.VCD = self.seeding_density
    self.interval = sample_interval
    self.sp = set_point
    self.cpp = cpp
    self.control = control_variable
  
  def step(self, assays, actuation):
    #Check if there is a new reading for the control parameter
    if self.cpp in assays:
      self.update_control(assays)
    return {self.control: self.addition_rate}
      
  
class fed_batch(feed_strategy):
  """Fed-batch reactor.  Adds a constant amount of specified feed in order to 
  adjust cpp to setpoint at next sample interval.
  
  Assumes constant:
    cell-specific consumption rate (since last obs)
    Volume
    VCD (VCD of last time period is averaged)
  
  cpp must be a concentration.
  
  """
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.addition_rate = 0
    self.last_time = self.seed_time
    
    
  def update_control(self, obs):
    seconds_since_last_obs = (obs['time']-self.last_time)/np.timedelta64(1, 's')
    average_VCD = (obs['VCD']+self.last_VCD)/2
    mass = obs[self.cpp]*obs['Volume']
    consumption_rate = \
      (self.addition_rate - (mass-self.last_mass)/seconds_since_last_obs) /\
      average_VCD
    self.addition_rate = (self.sp - obs[self.cpp])*obs['Volume']/self.interval +\
      obs['VCD']*consumption_rate
    self.last_mass = mass
    self.last_VCD = obs['VCD']
    self.last_time = obs['time']
    

  

class dynamic_perfusion(feed_strategy):
  """Perfusion reactor.  Adds a constant amount of nutrient in order to maintain
  nutrient at set point based on consumption rate since last observation.
  
  Assumes cell-specific consumption rate is constant.
  Assumes cells grow according to cell growth model."""
  def __init__(self):
    raise ValueError('Not implemented!')

  
  
class continuous_perfusion(feed_strategy):
  """Perfusion reactor.  Adds a constant amount of nutrient in order to maintain
  nutrient at set point based on consumption rate since last observation.
  
  Assumes constant consumption rate."""
  def __init__(self):
    raise ValueError('Not implemented!')




if __name__ == '__main__':
  tests.feeding_strategy_tests()