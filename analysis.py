# -*- coding: utf-8 -*-
"""
Holds all of the analysis that people do (not necessary for experiments).

Things like plotting variables, calculating sieving, etc.
"""
import pdb

import numpy as np
import matplotlib.pyplot as plt
from quantities import Quantity as Q

import simulator, helper_functions


def dual_subplot(y11_param, y12_param, y21_param=None, y22_param=None, metrics = None,
                 big_title = None, left_title = None, right_title = None,
                 y11_source = 'bioreactor', y12_source = 'bioreactor',
                 y21_source = 'bioreactor', y22_source = 'bioreactor'):
  """Plots two variables on each of two plots, then puts them next to each
  other."""
  fig = plt.figure(figsize = (10, 4))
  ax1 = fig.add_subplot(121)
  plot_var(ax1, y11_param, metrics, color = 'red', source = y11_source)
  
  if y21_param != None:
    ax2 = ax1.twinx()
    plot_var(ax2, y12_param, metrics, color = 'blue', source = y12_source)
  
  if left_title != None:
    ax1.set_title(left_title)
  
  
  ax3 = fig.add_subplot(122)
  if y21_param == None:
    plot_var(ax3, y12_param, metrics, color = 'red', source = y12_source)
  else:
    plot_var(ax3, y21_param, metrics, color = 'red', source = y21_source)
  
  if y22_param != None:
    ax4 = ax3.twinx()
    plot_var(ax4, y22_param, metrics, color = 'blue', source = y22_source)
  
  if right_title != None:
    ax3.set_title(right_title)
  
  if big_title != None:
    fig.suptitle(big_title, y = 1.04)
  fig.tight_layout()
  fig.dpi = 200
  plt.show()
  plt.close()

def dual_plot( y1_param, y2_param=None, metrics=None, title = None, 
              y1_source = 'bioreactor', y2_source = 'bioreactor'):
  """Plots (up to) two different metrics on one graph."""
  fig, ax1 = plt.subplots()
  plot_var(ax1, y1_param, metrics = metrics, color = 'red', source = y1_source)
  if y2_param != None:
    ax2 = ax1.twinx()
    plot_var(ax2, y2_param, metrics = metrics, color = 'blue', source = y2_source)
  if title != None:
    ax1.set_title(title)
  fig.tight_layout()
  fig.dpi = 200
  plt.show()
  plt.close()
  
def plot_var(axis, variable, metrics, color = 'black', source = 'bioreactor'):
  """Takes an axis to plot on, the metrics, and a variable of interest.  Then
  plots that variable on that axis."""
  if metrics == None:
    metrics = simulator.latest
  plot_func = helper_functions.get_plotfunc(axis, variable)
  for metric in metrics:
    tuples = [((metric[point]['time']-metric[0]['time'])/np.timedelta64(1,'D'),
           metric[point][variable])
          for point in range(len(metric)) if variable in metric[point]]
    if len(tuples) == 0:
      tuples = [((metric[point]['time']-metric[0]['time'])/np.timedelta64(1,'D'),
             metric[point][source][variable])
            for point in range(len(metric)) if source in metric[point]]
    if len(tuples) > 0:
      x, y = list(zip(*tuples))
      plot_func(x, y, color = color)
  if len(tuples) > 0:
    if max(y)/3> min(y):
      axis.set_ylim(bottom = -float(max(y))/100, top = None)
    if type(y[0]) == Q:
      axis.set_ylabel(variable+f' {y[0].dimensionality}', color = color)
    else:
      axis.set_ylabel(variable, color = color)

  axis.set_xlabel('Day')
  axis.set_xlim(0, (metric[-1]['time']-metric[0]['time'])/np.timedelta64(1,'D'))
  
def analyze(metrics):
  """Takes observed metrics and performs calculations on them to get "human" 
  analysis results such as sieving.  Also updates metrics with quantities."""
  for time_series in metrics: #unpack bioreactor
    for timepoint in time_series:
      try:
        sieving = timepoint['permeate']['IGG']/timepoint['bioreactor']['IGG']
        timepoint['sieving'] = sieving
      except KeyError:  #no permeate / not perfusion
        pass
      except ZeroDivisionError:
        #No IGG in bioreactor
        pass
  helper_functions.scale_units(metrics)
      
if __name__ == '__main__':
  import main
  main.run_perfusion()

        
    
  