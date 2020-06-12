# -*- coding: utf-8 -*-
"""
**Depricated.  Replaced with library 'quantities'

Holds various datatypes and defines how those data types can be manipulated

E.G. addition rate * time should equal a quantity added
"""
import numpy as np

import tests


class addition_rate:
  def __init__(self, amount, time_period):
    """Takes a quantity and unit time.  Can pass str ('s' for 1/s) or a numpy
    timedelta for unit time. """
    self.amount = amount
    if type(time_period) == str:
      time_period = np.timedelta64(1, time_period)
    self.time_period = time_period
    
  def __mul__(self, other):
    #Number of time periods
    multiplier = other / self.time_period
    #The total quantity added over time
    quantity = self.amount * multiplier
    return quantity
  
  def __add__(self, other):
    """Converts other addition rate to right units, and returns new addition
    rate. """
    assert type(other) == addition_rate, "can't add to non-addition rate"
    time_multiplier = self.time_period / other.time_period
    return addition_rate(self.amount+other.amount*multiplier, self.time_period)
  
  def __sub__(self, other):
    """Converts other addition rate to right units, and returns new addition
    rate. """
    assert type(other) == addition_rate, "can't add to non-addition rate"
    multiplier = self.time_period / other.time_period
    return addition_rate(self.amount-other.amount*multiplier, self.time_period)
  
  def __str__(self):
    return(f'{self.amount}/{str(self.time_period)}')
  
class physical_quantity:  
  def __init__(self, amount):
    self.amount = amount

class volume(physical_quantity):
  def __init__(self, amount, units):
    super().__init__(amount)
    
  
class mass(physical_quantity):
  
if __name__ == '__main__':
  tests.dtypes_tests()