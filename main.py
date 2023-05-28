import numpy as np

def Veff(z):
  '''
  Effective velosity in SI unit
  '''
  v = 3E4 * min(1, (1+z)/1000)
  return v
