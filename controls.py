# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:30:40 2020

@author: Racehorse
"""
import random
import pdb

from quantities import Quantity as Q
import numpy as np

import assays, actuation, param, helper_functions




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
               set_point, initial_time, target_seeding_density,
               max_added = None, feed_mixture = None,
               expected_cc_density = param.expected_cc_density):
    self.initial_volume = initial_volume.simplified
    self.seeding_density = target_seeding_density.simplified
    self.sp = set_point.simplified
    # if self.sp.dimensionality == Q(1,'kg/m**3').dimensionality:
    #   self.sp = self.sp / param.molecular_weight[cpp]
    self.interval = sample_interval.simplified
    self.addition_rate = 0*self.sp.units*Q(1,'m**3/s')
    self.total_added = Q(0, 'm**3')
    self.max_added = max_added 
    self.expected_density = expected_cc_density.simplified
    self.one_second = Q(1, 's')
    if feed_mixture == None:
      self.feed_mixture = {'glucose':Q(500, 'g/L').simplified}
    else:
      self.feed_mixture = feed_mixture.copy()
      helper_functions.simplify(self.feed_mixture)
    
    if param.skip_units:
      self.initial_volume = float(self.initial_volume)
      self.seeding_density = float(self.seeding_density)
      self.sp = float(self.sp)
      self.interval = float(self.interval)
      self.addition_rate = float(self.addition_rate)
      self.total_added = float(self.total_added)
      self.expected_density = float(self.expected_density)
      self.one_second = float(self.one_second)
    self.seed_time = initial_time
    self.VCD = self.seeding_density
    self.last_VCD = self.seeding_density
    self.last_volume = self.initial_volume
    self.ignored_first = 0
    self.cpp = cpp
    
  def step(self, obs, offline):
    #Check if there is a new reading for the control parameter
    metric = {}
    if offline:
      metric = self.update_control(obs)
    return metric
  
  
class constant_feed(feed_strategy):
  """As the name implies, adds a constant feed.  The name argument specifies
  what to call the returned metric."""
  def __init__(self, feed_rate, name, start_day = 0,
               feed_mixture = None, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.name = name
    self.start_day = start_day
    
    self.actuation = [actuation.peristaltic(self.feed_mixture)]
    self.feed_rate = feed_rate.simplified
    
    if param.skip_units: 
      helper_functions.remove_units(self.feed_mixture)
      self.feed_rate = float(self.feed_rate)
    self.feed_mixture['mOsm'] = self.actuation[0].source.calc_osmo()    
    if param.skip_units:
      for component in self.feed_mixture:
        self.feed_mixture[component] = float(self.feed_mixture[component])
      
  def update_control(self, obs):
    if (obs['time']-self.seed_time)/np.timedelta64(1, 'D')>self.start_day:
      self.actuation[0].set_point = self.feed_rate 
    return {f'{self.name}':self.actuation[0].set_point}


class basic_fed_batch_feed(feed_strategy):
  """Fed-batch reactor.  Adds a constant amount of specified feed in order to 
  adjust cpp to setpoint at next sample interval.
  
  Assumes mostly constant:
    cell-specific consumption rate (since last obs)
    
    Volume
    
    VCD (VCD of last time period is averaged)
  
  cpp must be a concentration.
  
  """
  def __init__(self, 
               *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.last_time = self.seed_time-np.timedelta64(1, 'D')
    self.actuation = [actuation.peristaltic(self.feed_mixture, 
                                            max_added = self.max_added)]
    if self.cpp == 'mOsm':
      self.conversion = 1
    else:
      #Assuming it is controlling it off a concentration (probably glucose)
      self.conversion = param.molecular_weight[self.cpp]
  

  def update_control(self, obs):
    """Calculates required flowrate and updates actuation."""
    volume = obs['mass']/self.expected_density
    time_since_last_obs = (obs['time']-self.last_time)/np.timedelta64(1, 's')
    self.total_added = self.actuation[0].set_point * time_since_last_obs   

    if self.ignored_first:
      if not self.cpp in obs['bioreactor']:
        raise ValueError('CPP not included in assays!')
      if not param.skip_units:
        time_since_last_obs *= Q(1, 's')
      average_VCD = (obs['bioreactor']['VCD']+self.last_VCD)/2
      average_volume = (self.last_volume+volume)/2
      CSCR = -(obs['bioreactor'][self.cpp]*volume - self.last_concentration*self.last_volume - \
              self.addition_rate * time_since_last_obs) /\
        (average_VCD * average_volume * time_since_last_obs)
      predicted_VCD = 2*obs['bioreactor']['VCD']-average_VCD #Maybe update this later....
      predicted_volume = 2*volume - average_volume
      self.addition_rate = ((self.sp*predicted_volume - obs['bioreactor'][self.cpp]*volume)/\
        self.interval + predicted_VCD*CSCR*predicted_volume)
      if self.addition_rate < 0:
        self.addition_rate *= 0
      
    self.actuation[0].set_point = self.addition_rate /\
      self.actuation[0].source.mixture_definition[self.cpp]/self.conversion
    self.last_VCD = obs['bioreactor']['VCD']
    self.last_time = obs['time']
    self.last_volume = volume
    self.last_concentration = obs['bioreactor'][self.cpp]
    self.ignored_first = 1
    
    return {f'{self.cpp} feed sp':self.actuation[0].set_point,
            f'{self.cpp} addition rate':self.addition_rate,}
  


class dynamic_perfusion_feed(feed_strategy):
  """Perfusion reactor.  Adds a constant amount of nutrient in order to maintain
  nutrient at set point based on consumption rate since last observation.
  
  Assumes cell-specific consumption rate is constant.
  Assumes cells grow according to cell growth model."""
  def __init__(self):
    raise ValueError('Not implemented!')
  
class basic_perfusion_feed(feed_strategy):
  """Perfusion reactor.  Adds a constant amount of feed in order to maintain
  nutrient at set point based on consumption rate since last observation.
  
  Assumes constant consumption rate."""
  def __init__(self, max_topoff_pump_speed = Q(300, 'ml/min'), max_VVD = Q(2., '1/d'),
               min_VVD = Q(0.5, '1/d'), volume_SP = Q(2, 'L'), 
               initial_VVD = Q(1, '1/d'), glucose_setpoint = None,
               *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.volume_SP = volume_SP.simplified
    self.mass_SP = (self.volume_SP * param.expected_cc_density).simplified
    self.max_flowrate = max_topoff_pump_speed.simplified
    self.actuation = [actuation.peristaltic(self.feed_mixture, self.max_added,
                                            max_flowrate = self.max_flowrate),
                      actuation.peristaltic(mode='permeate'), #permeate pump
                      ]
    self.VVD =initial_VVD.simplified
    self.max_VVD = max_VVD.simplified
    self.min_VVD = min_VVD.simplified
    self.glucose_control = False
    if glucose_setpoint != None:
      self.glucose_control = True
      self.actuation.append(actuation.peristaltic({'glucose': Q(500, 'g/L')}))
      self.glucose_setpoint = glucose_setpoint
    if self.cpp == 'mOsm':
      self.conversion = 1
    else:
      #Assuming it is controlling it off a concentration (probably glucose)
      self.conversion = param.molecular_weight[self.cpp]
    
    if param.skip_units:
      self.volume_SP = float(self.volume_SP)
      self.mass_SP = float(self.mass_SP)
      self.max_flowrate = float(self.max_flowrate)
      self.VVD = float(self.VVD)
      self.max_VVD = float(self.max_VVD)
      self.min_VVD = float(self.min_VVD)
      self.glucose_setpoint = float(self.glucose_setpoint) #0 if false; irrelevant
    
    #topup PID
    self.PID = PID(Q(100, '1/kg'), Q(10, '1/min/kg'), Q(10, 'min/kg'), 
                   self.mass_SP, 0)

  def update_control(self, obs):
    metrics = {}
    volume = obs['mass']/self.expected_density
    
    # Update every time in case of changeout
    feed_concentration = self.actuation[0].source.mixture_definition[self.cpp]
    if self.ignored_first:
      time_since_last_obs = (obs['time']-self.last_time)/np.timedelta64(1, 's')
      if not self.cpp in obs['bioreactor']:
        pdb.set_trace()
        raise ValueError('CPP not included in assays!')
      if not param.skip_units:
        time_since_last_obs *= Q(1, 's')
      self.total_added = self.actuation[0].set_point * time_since_last_obs
      average_VCD = (obs['bioreactor']['VCD']+self.last_VCD)/2
      average_volume = (self.last_volume+volume)/2
      predicted_VCD = 2*obs['bioreactor']['VCD']-average_VCD #Maybe update this later....
      predicted_volume = 2*volume - average_volume #shouldn't be chaning in perfusion mode
      if self.glucose_control:
        #if using glucose control section...
        glucose_feed_concentration = \
          self.actuation[2].source.mixture_definition['glucose']*param.molecular_weight['glucose']
        glucose_rate = self.actuation[2].set_point*glucose_feed_concentration+\
          self.actuation[0].source.mixture_definition['glucose']*self.VVD*\
          average_volume*param.molecular_weight['glucose']
        glucose_diff = (obs['bioreactor']['glucose']-self.last_glucose_concentration)/\
          (time_since_last_obs)
        average_glucose = (obs['bioreactor']['glucose']+self.last_glucose_concentration)/2
        GCR = (glucose_rate/average_volume - glucose_diff - self.VVD*average_glucose)/average_VCD
        metrics['GCR'] = GCR
      #cpp section  
      rate_change = (obs['bioreactor'][self.cpp]-self.last_concentration)/time_since_last_obs*average_volume
      CSCR = (self.VVD*average_volume*(self.last_feed_concentration-
              obs['bioreactor'][self.cpp])-rate_change)/(average_VCD*average_volume)
      target_VVD = (CSCR*predicted_VCD+(self.sp-obs['bioreactor'][self.cpp])/\
                  self.interval)/(feed_concentration-self.sp)
      self.VVD = 0.*self.VVD+1*target_VVD
      if self.VVD < self.min_VVD:
        self.VVD = self.min_VVD
      elif self.VVD > self.max_VVD:
        self.VVD = self.max_VVD
      metrics['CSCR'] = CSCR
      
      #Finish glucose control after calculating new VVD
      if self.glucose_control:
        target_glucose_rate = (GCR*predicted_VCD+self.VVD*self.glucose_setpoint+\
          (self.glucose_setpoint-obs['bioreactor']['glucose'])/self.interval)*predicted_volume          
        media_glucose_rate = self.actuation[0].source.mixture_definition['glucose']*\
          self.VVD*predicted_volume*param.molecular_weight['glucose']
        glucose_addition_rate = target_glucose_rate - media_glucose_rate
        if glucose_addition_rate < 0:
          glucose_addition_rate *= 0
        self.actuation[2].set_point = glucose_addition_rate/glucose_feed_concentration
        metrics['glucose addition rate'] = glucose_addition_rate
    
    metrics['VVD'] = self.VVD
    self.last_feed_concentration = feed_concentration
    self.actuation[1].set_point = self.VVD*self.volume_SP
    self.last_VCD = obs['bioreactor']['VCD']
    self.last_time = obs['time']
    self.last_volume = volume
    self.last_concentration = obs['bioreactor'][self.cpp]
    self.last_glucose_concentration = obs['bioreactor']['glucose']
    self.ignored_first = 1
    metrics['permeate rate']=self.actuation[1].set_point
    if self.glucose_control:
      metrics['glucose feed sp'] = self.actuation[2].set_point
    return metrics 
  
  def step(self, obs, offline):
    metrics = {}
    PID_out = self.PID.step(obs['mass'])
    self.actuation[0].set_point = PID_out/100*self.max_flowrate
    metrics['topoff'] = PID_out/100*self.max_flowrate

    if offline:
      metrics.update(self.update_control(obs))
    
    return metrics   
  
class PID:
  """Simple PID controller.  Can be used to control things like aeration.
  
  Assumes constant sample interval.
  
  Response is on scale of 0 to 100
  
  response(% of max) = Kp * D + int(Ki * D) dt - Kd * d value/dt
  where D is (value - setpoint)
  
  Needs a controlling strategy.
  """
  alpha = 0.35 #smoothing for derivative
  def __init__(self, Kp, Ki, Kd, set_point, initial_value):
    self.kp = helper_functions.remove_units(Kp)
    self.ki = float((Ki * param.q_res).simplified)
    self.kd = float((Kd / param.q_res).simplified)
    self.error_int = initial_value
    self.sp = helper_functions.remove_units(set_point)
    self.deriv = 0
    self.last_value = self.sp
   
  def step(self, value):
    value = float(value)
    D = self.sp - value
    self.error_int += self.ki*D
    self.deriv = (1-self.alpha) * self.deriv + self.alpha * (value - self.last_value)
    self.last_value = value
    response = self.kp * D + self.error_int - self.deriv*self.kd
    if (response < 0 and self.ki*D < 0) or \
      (response > 100 and self.ki*D > 0):
      self.error_int -= self.ki*D     #Undo int error if not within range to prevent windup
      if response < 0: response = 0
      else: response = 100    
    return response
  
class recirculation:
  """Controls the recirculation flowrate by changing pump RPM."""
  def __init__(self, setpoint,
               K_p = Q(0.5, 'min/L'), 
               K_i = Q(0.05, '1/L'),
               K_d = Q(0.2, 'min**2/L'),):
    # self.recirc_PID = PID(20000, Q(5000, 's/min/m**3'), Q(1000, '1/m**3'), 
    #                       self.recirc_setpoint, 10)
    self.sp = setpoint.simplified
    self.PID = PID(K_p, K_i, K_d, setpoint, 5)
    self.actuation = [actuation.levitronix()]
    if param.skip_units:
      self.sp = float(self.sp)
    
  def step(self, obs, offline):
    """Doesn't need to return metrics; obs has it."""
    PID_out = self.PID.step(obs['recirc_rate'])
    self.actuation[0].set_point = PID_out/100*self.actuation[0].max_RPM
    return {'recirc_PID':PID_out,}
    
  
class aeration:
  """aeration strategy.  Uses a PID to first add air to max air, then
  adds oxygen to max oxygen.
  """
  def __init__(self, setpoint, 
               min_air = Q(0.01, 'L/min'), 
               max_air = Q(0.2, 'L/min'), 
               max_O2 = Q(0.1, 'L/min'),
               K_p = 0.05, 
               K_i = Q(0.02, '1/min'),
               K_d = Q(0.05, 'min'),):
    """max_air and max_O2 should be in volumetric flow rates."""
    
    self.sp = setpoint
    
    self.PID = PID(K_p, K_i, K_d, setpoint, 5)
    self.actuation = [actuation.MFC('air'), actuation.MFC('O2'), actuation.agitator(Q(300/60., '1/s'))]
    self.max_air = max_air.simplified
    self.max_O2 = max_O2.simplified
    self.min_air = min_air.simplified
    self.zero = Q(0, 'L/min').simplified
    
    if param.skip_units:
      self.max_air = float(self.max_air.simplified)
      self.max_O2 = float(self.max_O2)
      self.min_air = float(self.min_air)
      self.zero = float(self.zero)
    
  def step(self, obs, offline):
    
    PID_out = self.PID.step(obs['dO2'])
    if PID_out < 50:
      self.actuation[0].set_point = PID_out/50 * (self.max_air-self.min_air)+self.min_air
      self.actuation[1].set_point =  self.zero
    else:
      self.actuation[0].set_point = self.max_air
      self.actuation[1].set_point = (PID_out-50)/50 * self.max_O2
    return {'aeration_PID':PID_out,
            'O2':self.actuation[1].set_point,
            'air':self.actuation[0].set_point,}

class pH:
  """CO2 Aeration control strategy to control pH within set range.
  
  Currently implemented: max pH controlled with CO2
  Future: minimum pH copntrolled by adding base."""
  def __init__(self, max_pH, 
               min_CO2 = Q(0, 'L/min'), 
               max_CO2 = Q(0.1, 'L/min'),):
    """max_air and max_O2 should be in volumetric flow rates."""
    # self.PID = PID(-3/5, -Q(5, '1/min')/20, -Q(2, 'min')/10, 
    self.PID = PID(-3, -Q(5, '1/min'), -Q(2, 'min'), 
                   max_pH, 5)
    self.actuation = [actuation.MFC('CO2')]
    self.max_CO2 = max_CO2.simplified
    self.min_CO2 = min_CO2.simplified
    self.sp = max_pH
    
    if param.skip_units:
      self.max_CO2 = float(self.max_CO2)
      self.min_CO2 = float(self.min_CO2)

    
  def step(self, obs, offline):
    PID_out = self.PID.step(obs['pH'])
    self.actuation[0].set_point = PID_out/100 * (self.max_CO2-self.min_CO2)+self.min_CO2
    return {'pH_PID':PID_out,
            'CO2':self.actuation[0].set_point,}

class temperature:
  """Adds heat to control temperature."""
  def __init__(self, setpoint, max_heating = Q(100, 'W')):
    self.PID = PID(10, Q(1, '1/min'), Q(3, 'min'), 
                   setpoint, 50)
    self.actuation = [actuation.heating_jacket(max_heating)]
    
  def step(self, obs, offline):
    PID_out = self.PID.step(obs['temperature'])
    self.actuation[0].set_point = PID_out
    return {'heating_PID':PID_out}
  
class wrapper:
  def __init__(self, control_setup):
    """control_setup is list of tuples for each control.
    Tuple is in format (control, *args, **kwargs)
    
    Setup initializes each control and creates an actuation list.
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