#!/usr/bin/env python

import numpy as np
import spiceypy as spice
import orbit_comet
import util
from skyfield.api import Topos, load

from argparse import ArgumentParser

def read_objects(filename):
  f = open(filename)
  t = []
  for x in f.readlines():
    t.append(x.rstrip())
  return t

parser = ArgumentParser()
parser.add_argument("--start", dest="time_start", help="start time in UTC, e.g. 2019-04-28T18:20:00", required=True)
parser.add_argument("--telescope", dest="telescope", help="telescope name, e.g. T24", required=True)
parser.add_argument("--objects", dest="objects_filename", help="filename of object list", required=True)
parser.add_argument("--tle", dest="tle_filename", help="filename of TLEs", required=True)
args = parser.parse_args()

filter_broadband_table = {
  'T8': 'Luminance',
  'T9': 'Luminance',
  'T12': 'Luminance',
  'T17': 'Luminance',
  'T30': 'Luminance',
  'T31': 'Luminance',
  'T32': 'Luminance',
  'T33': 'Luminance',
  'T63': 'Luminance',
  'T68': None,
  'T24': 'Luminance',
  'T40': 'Luminance',
  'T5': 'Clear',
  'T11': 'Luminance',
  'T14': 'Luminance',
  'T20': 'Luminance',
  'T21': 'Luminance',
  'T7': 'Luminance',
  'T16': 'Luminance',
  'T18': 'Luminance',
}  

def filter_broadband(t):
  return filter_broadband_table[t]

if args.telescope in ['T8','T9','T12','T17','T30','T31','T32','T33']:      # Australia (Siding Spring Observatory)
  site = Topos('31.2733055 S', '149.0697472 E', elevation_m=1165)
elif args.telescope in ['T63','T64','T68']:                                # Australia (Bathurst Observatory)
  site = Topos('33.3966667 S', '149.4916667 E', elevation_m=650)
elif args.telescope in ['T24']:                                            # California (Sierra Remote Observatory)
  site = Topos('37.070278 N','119.412778 W',elevation_m=1405)
elif args.telescope in ['T40']:                                            # Chile
  site = Topos('30.4533333 S','70.7533333 W',elevation_m=1560)
elif args.telescope in ['T5','T11','T14','T20','T21']:                     # New Mexico
  site = Topos('32.9036111 N','105.5283333 W',elevation_m=2225)
elif args.telescope in ['T7','T16','T18']:                                 # Spain
  site = Topos('38.1655555 N','2.3263333 W',elevation_m=1650)
else:
  raise 'unknown telescope'

filter = filter_broadband(args.telescope)

# telescope timing model

t_init = 8.0            # time to initialize
t_large_slew = 70.0     # average time to execute a large slew
t_small_slew = 7.0      # small slew, less than a few degrees
t_save = 25.0           # time to save an image
t_small_exposure = 2    # exposure time, initial image
t_large_exposure = 10   # exposure time, main image

# ACP (telescope control system) apparent errors

acp_time_bias = -146.0  # bias in time when propagating comet orbits
acp_lag_bias = -16.0

# print script header

print("""; telescope %s
#NOADAPTIVEFOCUS
#EXPRESS
#BINNING 1""" % (args.telescope))
if filter:
  print("#FILTER %s" % filter)

ts = load.timescale()

et = spice.str2et(args.time_start)
et = et + t_init

satellites = load.tle(args.tle_filename)
objects = read_objects(args.objects_filename)

for sat_name in objects:
  satellite = satellites[sat_name]
  difference = satellite - site

  t = spice.unitim(et,"ET","JDTDB")
  t = ts.tdb(jd=t)

  topocentric = difference.at(t)
  ra0,dec0,distance = topocentric.radec()
  alt, az, distance = topocentric.altaz()
  if alt.degrees<30.0:
    continue

  t = spice.unitim(et+5,"ET","JDTDB")
  t = ts.tdb(jd=t)

  topocentric = difference.at(t)
  ra1,dec1,distance = topocentric.radec()

  ra = ra0.hms()
  dec = dec0.dms()
  dra = (ra1._degrees-ra0._degrees)/5 # fixme: take modulo 24 hours
  ddec = (dec1.degrees-dec0.degrees)/5
  dra = dra/(180/np.pi)
  ddec = ddec/(180/np.pi)

  number = satellite.model.satnum

  err0,err1,elts = orbit_comet.comet_orbit_search(et+acp_time_bias+acp_lag_bias,ra,dec,dra,ddec)

  print('; satellite: <%s> alt: %.0f az: %.0f' % (satellite.name,alt.degrees,az.degrees))
  print('; time: %s' % spice.timout(et,"YYYY-MM-DDTHR:MN:SC.### ::RND ::UTC",24))
  print('; ra %d %d %.2f   dec %d %d %.2f   dra %.3f ddec %.3f' % (ra[0],ra[1],ra[2],dec[0],dec[1],dec[2],dra*(180/np.pi)*3600,ddec*(180/np.pi)*3600))
  x = orbit_comet.print_comet_mpc(elts,number)
  print('#INTERVAL %d' % t_small_exposure)     # take a short exposure after a (potentially) large slew
  print('#TRACKOFF')                           # no orbital tracking: round stars (possibly useful for an initial astrometric estimate), trailed satellite (if visible)
  print(x)
  print('#INTERVAL %d' % t_large_exposure)     # take the main exposure after a smaller slew
  print('#TRACKON')                            # orbital tracking: round satellite, trailed stars
  print(x)

  et = et + t_large_slew + t_small_exposure + t_save + t_small_slew + t_large_exposure + t_save
