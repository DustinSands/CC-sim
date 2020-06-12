# -*- coding: utf-8 -*-
"""
Holds all instrumentation devices.

TO DO:
  Allow probes to catastrophically fail
  Program fail behavior
  
  
drift_sigma refers to drift per year variance

"""
import random

import numpy as np

import tests, param


class machine:
  """All machines have a random calibration error.  They then drift each time 
  they are used until recalibrated."""
  # def __init__(self):
  #   self.systematic_error = random.gauss(1, self.sys_err_sigma)


class probe(machine):
  """Moves towards correct value each step.
  """
  sys_err_sigma = 0.03
  def __init__(self, cal_reference):
    # Timestep propagation delay ratio
    self.ratio = 1-(1-0.98)**(param.resolution/self.p['t98'])
    self.reference = cal_reference

    
  
class O2(probe):
  """All errors (except one-point) are introduced as a slope, not offset."""
  # Slope drift.  <5% every 2 months
  p = param.instrumentation['Probe']['O2']
  
  def __init__(self, **kwargs):
    super(O2, self).__init__(**kwargs)
    self.drift_slope = random.gauss(0, self.p['drift_sigma'])
    #Large offset as probe isn't yet calibrated
    self.sys_error = 5
    self.value = 100
    
    
  def one_point(self, time, value = 100):
    self.cal_time = time
    self.sys_error = self.reference.read_O2_value(value) - \
      value + random.gauss(0, self.p['random_sigma'])    
  
  def read_value(self, state):
    # time_delta is time since last calibration
    
    time_delta = state['time'] - self.cal_time
    drift_error = time_delta / np.timedelta64(365, 'D')*self.drift_slope
    random_error = random.gauss(0, self.p['random_sigma'])
    value = state['O2']*(1+drift_error+random_error)+self.sys_error
    self.value = self.ratio*value+(1-self.ratio)*self.value
    return self.value
    
  
class pH(probe):
  """pH is a log scale; any error in the slope expresses itself as a net error.
  """
  # Net drift: < 0.1 pH/week
  p = param.instrumentation['Probe']['pH']
  
  def __init__(self, **kwargs):
    super(pH, self).__init__(**kwargs)
    self.drift_slope = random.gauss(0, self.p['drift_sigma'])
    #Large offset as probe isn't yet calibrated
    self.sys_error = 5
    self.value = 30
    
    
  def one_point(self, time, value = 7):
    self.cal_time = time
    self.sys_error = self.reference.read_pH_value(value) - \
      value + random.gauss(0, self.p['random_sigma'])    
  
  def read_value(self, state):
    # time_delta is time since last calibration
    
    time_delta = state['time'] - self.cal_time
    drift_error = time_delta / np.timedelta64(365, 'D')*self.drift_slope
    random_error = random.gauss(0, self.p['random_sigma'])
    value = state['pH']+drift_error+random_error+self.sys_error
    self.value = self.ratio*value+(1-self.ratio)*self.value
    return self.value
    
class BGA(machine):
  """BGA have very high drift, and as a result are designed to constantly 
  self-calibrate.  Therefor, we assume there is no drift error and instead have 
  only random error."""
  p = param.instrumentation['BGA']
  
  def __init__(self):
    self.O2_systematic_error = random.gauss(0, self.p['O2_systematic_error_sigma'])
    self.pH_systematic_error = random.gauss(0, self.p['pH_systematic_error_sigma'])
    
  def read_O2_value(self, value):
    random_error = random.gauss(0, self.p['O2_random_error_sigma'])
    return value * (1+ random_error + self.O2_systematic_error)
    
    
  def read_pH_value(self, value):
    random_error = random.gauss(0, self.p['pH_random_error_sigma'])
    return value + random_error + self.pH_systematic_error

    
def bioHT(machine):
  """Performs various assays.  Takes tests to perform as argument."""
  
  def __init__(self):

class cell_counter(machine):
  """No drift implemented, as was not able to find any figures for this."""
  p = param.instrumentation['Cell_Counter']
  
  def __init__(self, **kwargs):
    super(cell_counter, self).__init__(**kwargs)
    self.update_calibration()

  def update_calibration(self):
    self.sys_error = random.gauss(0, self.p['systematic_error_sigma'])

  def read_value(self, state):
    # time_delta is time since last calibration
    random_error = random.gauss(0, self.p['random_error_sigma'])
    value = state['VCD']*(1+random_error+self.sys_error)
    return value
  
if __name__ == '__main__':
  tests.assays_test()
  
