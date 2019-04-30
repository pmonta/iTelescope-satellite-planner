import numpy as np
import spiceypy as spice
import random
from util import *

spice.furnsh("kernels.txt")

au = 149597870.700       # astronomical unit
mu = 1.327124400e11      # GM_sun
e = (np.pi/180)*(23+(26.0/60)+(21.406/3600))   # obliquity of the ecliptic
m = np.array([[1,0,0],[0,np.cos(e),np.sin(e)],[0,-np.sin(e),np.cos(e)]])  # rotation matrix, equatorial to ecliptic

def comet(et,p,v):
  state = [p[0],p[1],p[2],v[0],v[1],v[2]]   # km, km/s
  elts = spice.oscelt(state,et,mu)          # compute osculating Keplerian parameters from a state vector
  RP,ECC,INC,LNODE,ARGP,M0,T0,MU = elts
  a = np.abs(RP/(1-ECC))
  n = np.sqrt(MU/(a*a*a))    # rad/s
  TP = T0 - M0/n  #fixme: need a case distinction for the sign of this correction?
  RP = RP/au
  INC = rad2deg(INC)
  LNODE = rad2deg(LNODE)
  ARGP = rad2deg(ARGP)
  TP = float(spice.timout(TP,"JULIAND.######## ::TDB"))  #fixme
  return RP,ECC,ARGP,LNODE,INC,TP

# quantize cometary orbital parameters to precision available in the MPC format

def quantize_comet(c):
  RP,ECC,ARGP,LNODE,INC,TP = c
  RP = np.around(RP,6)
  ECC = np.around(ECC,6)
  ARGP = np.around(ARGP,4)
  LNODE = np.around(LNODE,4)
  INC = np.around(INC,4)
  TP = np.around(TP,4)
  return RP,ECC,ARGP,LNODE,INC,TP

# compute error in position between initial state vector and quanitzed cometary orbit

def comet_position_error(et,p,elts):
  RP,ECC,ARGP,LNODE,INC,TP = elts
  qRP = RP*au
  qECC = ECC
  qINC = deg2rad(INC)
  qLNODE = deg2rad(LNODE)
  qARGP = deg2rad(ARGP)
  T0 = et
  qTP = spice.str2et('%.14f JD'%TP)
  qTP = qTP - 69.184  #fixme
  qa = np.abs(qRP/(1-qECC))
  qn = np.sqrt(mu/(qa*qa*qa))
  qM0 = qn*(T0-qTP)
  qelts = qRP,qECC,qINC,qLNODE,qARGP,qM0,T0,mu
  state = spice.conics(qelts,et)
  p_c = state[0:3]
  e = p - p_c
  return np.sqrt(np.sum(e*e))

# print the cometary orbint in Stellarium format

def print_comet_stellarium(c,name):
  RP,ECC,ARGP,LNODE,INC,tp = c
  print('# %s' % name)
  print('[%s]' % name)
  print('absolute_magnitude=19')
  print('albedo=0.1')
  print('coord_func=comet_orbit')
  print('name=%s' % name)
  print('orbit_ArgOfPericenter=%.4f' % ARGP)
  print('orbit_AscendingNode=%.4f' % LNODE)
  print('orbit_Eccentricity=%.6f' % ECC)
  print('orbit_Inclination=%.4f' % INC)
  print('orbit_PericenterDistance=%.6f' % RP)
  print('orbit_TimeAtPericenter=%.4f' % tp)
  print('type=comet')
  print('')

# print the cometary orbint in Horizons format

def print_comet_horizons(c,name,et):
  RP,ECC,ARGP,LNODE,INC,TP = c
  T0 = float(spice.timout(et,"JULIAND.#### ::TDB"))
  print("EPOCH=%.4f EC=%.6f QR=%.6f TP=%.4f OM=%.4f W=%.4f IN=%.4f" % (T0,ECC,RP,TP,LNODE,ARGP,INC))

# print the cometary orbint in MPC format

def print_comet_mpc(c,name):
  RP,ECC,ARGP,LNODE,INC,tp = c
  t_perihelion = spice.unitim(tp,"JDTDB","ET")
  s = spice.timout(t_perihelion,"YYYY MM DD.#### ::RND ::TDB",16)
  year = int(s[0:4])
  month = int(s[5:7])
  day = float(s[8:])
  mag = 8.0
  slope = 4.0
  return '    P%s    %04d %02d %7.4f %9.6f  %8.6f  %8.4f  %8.4f  %8.4f  20171201  %4.1f %4.1f  P/%s' % (name,year,month,day,RP,ECC,ARGP,LNODE,INC,mag,slope,name)

# Earth's state vector (heliocentric)

def earth_state(et):
  state, ltime = spice.spkezr('399',et,'ECLIPJ2000','NONE','10')
  p = state[0:3]
  v = state[3:6]
  return p,v

# heliocentric unit vector in the direction of the given equatorial celestial coordinates

def dir_ra_dec(ra,dec):
  v = np.array([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])
  return m.dot(v)

# compute a cometary orbit using a (geocentric) direction, speed, and range

def mcomet(et,ra,dec,dra,ddec,c_range):
  et_earth = et
  et_orbit = et
  p_earth0,v_earth0 = earth_state(et_earth)
  p_earth1,v_earth1 = earth_state(et_earth+5)  # fixme: something better than a finite-difference derivative estimate across a 5-second span
  ra0 = ra2rad(ra[0],ra[1],ra[2])
  dec0 = dec2rad(dec[0],dec[1],dec[2])
  dir0 = dir_ra_dec(ra0,dec0)
  p0 = p_earth0 + c_range*dir0                 # extend outward by c_range
  dir1 = dir_ra_dec(ra0+5*dra,dec0+5*ddec)
  p1 = p_earth1 + c_range*dir1
  v = (p1-p0)/5
  c = comet(et_orbit,p0,v)
  c = quantize_comet(c)
  err0 = comet_position_error(et_orbit,p0,c)     # error in km at nominal time
  err1 = comet_position_error(et_orbit+5,p1,c)   # error in km at nominal + 5 seconds
  return err0,err1,c

# Search for a cometary orbit that's a good match to a celestial object at (ra,dec) moving at (dra,ddec)

def comet_orbit_search(et,ra,dec,dra,ddec):
  err0_l = []
  err1_l = []
  elts_l = []
  for i in range(500):
    c_range = 250000.0 + random.uniform(-30000.0,30000.0)  # 500 random trials near 250000 km
    err0,err1,elts = mcomet(et,ra,dec,dra,ddec,c_range)    # compute the heliocentric cometary orbit
    err0_l.append(err0)
    err1_l.append(err1)
    elts_l.append(elts)
  j = np.argmin(np.array(err0_l))                          # pick the best
  return err0_l[j],err1_l[j],elts_l[j]
