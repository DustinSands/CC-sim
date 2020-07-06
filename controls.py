# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:30:40 2020

@author: Racehorse
"""
import random

from quantities import Quantity as Q
import numpy as np

import tests, actuation, param




class feed_strategy:
  """Base class for feed strategy.
  
  cpp is the critical process parameter that the feed is attemping to control.
  
  volume is bioreactor volume
  
  sample_interval is the expected time between offline samples
  
  set_point is the setpoint to maintain the cpp at
  
  addition_rate is the desired setpoint of addition of feed media
  
  Future work: should be generalized to anything that updates offline
  """
  def __init__(self, initial_volume, sample_interval, cpp, 
               set_point, initial_time, target_seeding_density):
    self.initial_volume = initial_volume
    self.seeding_density = target_seeding_density
    self.seed_time = initial_time
    self.VCD = self.seeding_density
    self.last_VCD = target_seeding_density
    self.last_mass = self.initial_volume*param.expected_cc_density
    self.interval = sample_interval
    self.sp = set_point
    self.cpp = cpp
  
  def step(self, obs, offline):
    #Check if there is a new reading for the control parameter
    if offline:
      metric = self.update_control(obs)
    return metric
  
class fed_batch_feed(feed_strategy):
  """Fed-batch reactor.  Adds a constant amount of specified feed in order to 
  adjust cpp to setpoint at next sample interval.
  
  Assumes constant:
    cell-specific consumption rate (since last obs)
    
    Volume
    
    VCD (VCD of last time period is averaged)
  
  cpp must be a concentration.
  
  """
  def __init__(self, feed_mixture = None, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.addition_rate = Q(0, 'g/min')
    self.last_time = self.seed_time-np.timedelta64(1, 'D')
    
    if feed_mixture == None:
      feed_mixture = {'glucose':Q(500, 'g/L')}
    self.feed_mixture = feed_mixture
    self.actuation = [actuation.peristaltic(feed_mixture)]

    
  # def update_control(self, obs):
  #   time_since_last_obs = Q((obs['time']-self.last_time)/np.timedelta64(1, 'm'), 'min')
  #   average_VCD = (obs['VCD']+self.last_VCD)/2
  #   consumption_rate = \
  #     (self.addition_rate - (obs['mass']-self.last_mass)/time_since_last_obs) /\
  #     average_VCD
  #   self.predicted_VCD = obs['VCD'] #Maybe update this later....
  #   self.addition_rate = (self.sp - obs[self.cpp])*obs['Volume']/self.interval +\
  #     self.predicted_VCD*consumption_rate
  #   self.last_mass = obs['mass']
  #   self.last_VCD = obs['VCD']
  #   self.last_time = obs['time']
  #   self.actuation.set_point = self.addition_rate
    
  def update_control(self, obs):
    if not self.cpp in obs:
      raise ValueError('CPP not included in assays!')
    time_since_last_obs = Q((obs['time']-self.last_time)/np.timedelta64(1, 'm'), 'min')
    average_VCD = (obs['VCD']+self.last_VCD)/2
    new_mass = obs[self.cpp]*obs['mass']/param.expected_cc_density
    consumption_rate = (self.addition_rate - \
      (new_mass-self.last_mass)/time_since_last_obs)/average_VCD
    self.predicted_VCD = 2*obs['VCD']-average_VCD #Maybe update this later....
    self.addition_rate = (self.sp - obs[self.cpp])*obs['mass']/param.expected_cc_density/\
      self.interval +\
      self.predicted_VCD*consumption_rate
    self.last_mass = new_mass
    self.last_VCD = obs['VCD']
    self.last_time = obs['time']
    self.actuation[0].set_point = self.addition_rate / self.feed_mixture[self.cpp]
    return {'glucose_solution':self.actuation[0].set_point}

class dynamic_perfusion_feed(feed_strategy):
  """Perfusion reactor.  Adds a constant amount of nutrient in order to maintain
  nutrient at set point based on consumption rate since last observation.
  
  Assumes cell-specific consumption rate is constant.
  Assumes cells grow according to cell growth model."""
  def __init__(self):
    raise ValueError('Not implemented!')
  
class continuous_perfusion_feed(feed_strategy):
  """Perfusion reactor.  Adds a constant amount of nutrient in order to maintain
  nutrient at set point based on consumption rate since last observation.
  
  Assumes constant consumption rate."""
  def __init__(self):
    raise ValueError('Not implemented!')

class PID:
  """Simple PID controller.  Can be used to control things like aeration.
  
  Assumes constant sample interval.
  
  response(% of max) = Kp * D + int(Ki * D) dt - Kd * d value/dt
  where D is (value - setpoint)
  """
  alpha = 0.35 #smoothing for derivative
  def __init__(self, Kp, Ki, Kd, set_point, initial_value):
    self.kp = Kp
    self.ki = float((Ki * param.q_res).simplified)
    self.kd = float((Kd / param.q_res).simplified)
    self.error_int = initial_value
    self.sp = set_point
    self.deriv = 0
    self.last_value = 60
    
  def step(self, value):
    D = self.sp - value
    self.error_int += self.ki*D
    self.deriv = (1-self.alpha) * self.deriv + self.alpha * (value - self.last_value)
    self.last_value = value
    response = self.kp * D + self.error_int - self.deriv*self.kd
    if not 0 <= response <= 100:
      self.error_int -= self.ki*D     #Undo int error if not within range to prevent windup
      if response < 0: response = 0
      else: response = 100    
    return response
  
class aeration:
  """aeration strategy.  Uses a PID to first add air to max air, then
  adds oxygen to max oxygen.
  """
  def __init__(self, setpoint, max_air = Q(0.2, 'L/min'), 
               max_O2 = Q(0.1, 'L/min')):
    """max_air and max_O2 should be in volumetric flow rates."""
    self.PID = PID(0.1, Q(1/15, '1/min'), Q(3, 'min'), 
                   setpoint, 5)
    self.actuation = [actuation.MFC('air'), actuation.MFC('O2'), actuation.agitator(300/60)]
    self.max_air = max_air
    self.max_O2 = max_O2
    
  def step(self, obs, offline):
    PID_out = self.PID.step(obs['dO2'])
    if PID_out < 80:
      self.actuation[0].set_point = PID_out/80 * self.max_air
      self.actuation[1].set_point =  Q(0,'L/min')
    else:
      self.actuation[0].set_point = self.max_air
      self.actuation[1].set_point = (PID_out-80)/20 * self.max_O2
    return {'aeration_PID':PID_out}

  
class wrapper:
  def __init__(self, control_setup):
    """control_setup is list of controls that are already initialized.
    Setup creates a list of actuation.
      """
    self.actuation_list = []
    self.control_list = []
    for control_type, args, kwargs in control_setup:
      control = control_type(*args, **kwargs)
      self.control_list.append(control)
      self.actuation_list.extend(control.actuation)
      
  def step(self, assays, offline):
    metrics = {}
    for control in self.control_list:
      metric = control.step(assays, offline)
      metrics.update(metric)
    
      
    return metrics



if __name__ == '__main__':
  tests.controls_tests()