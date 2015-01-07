# -*- coding: utf-8 -*-
"""
Created on Tue Dec 09 18:15:45 2014

@author: Drew Leonard
"""

from matplotlib import use
use('agg')
import numpy as np
import matplotlib.pyplot as plt
from sunpy import wcs
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import datetime as dt
from os.path import join
from sys import path
path.append('/imaps/holly/home/ajl7/CoronaTemps/')
from temperature import TemperatureMap as tmap
from astropy import units


def flareclass_to_flux(flareclass):
    """
    Converts a GOES flare class into the corresponding X-ray flux.
    
    Parameters
    ----------
    flareclass : string
        The flare class to convert into X-ray flux, as a string. 
        E.g.: 'X3.2', 'm1.5', 'A9.6'.
    
    Returns
    -------
    flux : astropy.units.Quantity
        X-ray flux between 1 and 8 Angstroms as measured near Earth in W/m^2
    
    Examples
    --------
    >>> flareclass_to_flux('A1.0')
    1e-08
    >>> flareclass_to_flux('c4.7')
    4.7e-06
    >>> flareclass_to_flux('X2.4')
    0.00024

    """
    assert isinstance(flareclass, str)
    flareclass = flareclass.upper()
    conversions = {'A': 1.0e-8, 'B': 1.0e-7, 'C': 1.0e-6, 'M': 1.0e-5,
                   'X': 1.0e-4}
    fluxval = float(flareclass[1:]) * conversions[flareclass[0]]
    flux = units.Quantity(fluxval, "W/m^2")
    
    return flux


start = parse('2011-02-01')
end = parse('2011-04-01')

client = hek.HEKClient()
flares = client.query(hek.attrs.Time(start, end),
                      hek.attrs.EventType('FL'))

flares = [fl for fl in flares if (fl['ar_noaanum'] > 11137 and 
                                  fl['ar_noaanum'] < 11184)]

ar_rad = 75
ar_temps = []
fl_classes = []

for flare in flares:
  try:
    flaretime = parse(flare['event_starttime'])
    starttime = flaretime-dt.timedelta(hours=0.5)
    timerange = tr(starttime, flaretime)

    region = client.query(hek.attrs.EventType('AR'),
                          hek.attrs.Time(flaretime-dt.timedelta(minutes=5), 
                                         flaretime))
    if flare['ar_noaanum'] == 11149:
        print 'HEK is fucked, using 11147'
        flare['ar_noaanum'] = 11147
    region = [r for r in region if r['ar_noaanum'] == flare['ar_noaanum']]
    if isinstance(region, list):
        try:
            region = region[0]
        except:# IndexError:
            print "An error occured for active region AR{}".format(flare['ar_noaanum'])
            continue
    
    # Define times for maps
    delta = dt.timedelta(minutes=1)
    ntimes = int(timerange.seconds()/delta.total_seconds())
    times = [time.start() for time in timerange.split(ntimes)]
    
    """home = '/imaps/sspfs/archive/sdo/aia'
    data_dir = join(home, 'activeregions/AR11153/data/')
    maps_dir = join(home, 'activeregions/AR11153/images/')"""
    data_dir = '/imaps/sspfs/archive/sdo/aia/activeregions/AR{}/data/'.format(flare['ar_noaanum'])
    maps_dir = '/imaps/holly/home/ajl7/tempmaps'
    
    means = []
    for time in times:
        # Load/calculate temperature map data
        try:
            thismap = tmap(time, data_dir=data_dir, maps_dir=maps_dir)
            thismap.save()
        
            # Crop temperature map to active region
            x, y = wcs.convert_hg_hpc(region['hgc_x'], region['hgc_y'],
                                      b0_deg=thismap.heliographic_latitude,
                                      l0_deg=thismap.carrington_longitude)
            thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])
            
            # Append appropriate temperature values to list
            means.append(np.nanmean(thismap.data))
        except:
            print "Failed", time
            means.append(np.nan)
            raise
    
    # Convert time values to time before flare
    times = [(t - flaretime).total_seconds()/60 for t in times]
    # Plot temperature values of AR with time for that flare
    fig = plt.figure()
    plt.plot(times, means)
    fname = 'AR{}/{}__{}-class'.format(flare['ar_noaanum'], flare['event_starttime'], flare['fl_goescls'])
    fname = fname.replace('.', '_')
    plt.savefig(join('/imaps/holly/home/ajl7/tempplots/', fname))
    plt.close()
    # Append  temperature values for final temperature map to list
    ar_temps.append(means[-1])
    # Append class of flare to list
    fl_classes.append(flareclass_to_flux(str(flare['fl_goescls'])).value)
  except:
    print 'Failed for {} flare at {}'.format(flare['fl_goescls'], flare['event_starttime'])
    raise

# Plot instantaneous temperatures of active regions for all flares against flare class
fig = plt.figure()
plt.plot(ar_temps, fl_classes)
plt.savefig("/imaps/holly/home/ajl7/tempplots/allflares")
plt.close()
