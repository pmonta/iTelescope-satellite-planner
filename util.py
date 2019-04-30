import numpy as np

# radians to degrees

def rad2deg(r):
  return r*(180/np.pi)

# degrees to radians

def deg2rad(x):
  return x*(np.pi/180)

# right ascension (mixed radix) to radians

def ra2rad(h,m,s):
  ra = h + m/60.0 + s/3600.0
  ra = (15*ra)*(np.pi/180)
  return ra

# declination (mixed radix) to radians

def dec2rad(d,m,s):
  dec = d + m/60.0 + s/3600.0
  dec = dec*(np.pi/180)
  return dec
