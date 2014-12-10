# -*- coding: utf-8 -*-
"""
Created on Tue Dec 09 18:15:45 2014

@author: Drew Leonard
"""

from matplotlib import use
use('agg')
from sunpy import wcs
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import datetime as dt
from os.path import join
from sys import path
path.append('/imaps/holly/home/ajl7/CoronaTemps/')
from temperature import TemperatureMap as tmap

start = parse('2011-02-10')
end = parse('2011-04-01')

client = hek.HEKClient()
flares = client.query(hek.attrs.Time(start, end),
                      hek.attrs.EventType('FL'))

flares = [fl for fl in flares if (fl['ar_noaanum'] > 11137 and 
                                  fl['ar_noaanum'] < 11184)]

fluxes = []
temps = []
coords = []
datefile = open('dates.txt', 'w')
datefile.write('')
datefile.close()

ar_rad = 75

for flare in flares:
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
        except IndexError:
            continue
    
    # Define times for maps
    delta = dt.timedelta(minutes=1)
    ntimes = int(timerange.seconds()/delta.total_seconds())
    times = [time.start() for time in timerange.split(ntimes)]
    
    home = '/imaps/sspfs/archive/sdo/aia'
    data_dir = join(home, 'activeregions/AR11153/data/')
    maps_dir = join(home, 'activeregions/AR11153/images/')
    
    for time in times:
        # Load/calculate temperature map data
        thismap = tmap(time, data_dir=data_dir, maps_dir=maps_dir)
        thismap.save()
        
        # Crop temperature map to active region
        x, y = wcs.convert_hg_hpc(region['hgc_x'], region['hgc_y'],
                                  b0_deg=thismap.heliographic_latitude,
                                  l0_deg=thismap.carrington_longitude)
        thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])
        
        # Append appropriate temperature values to list
    
    # Append  temperature values for final temperature map to list
    # Plot temperature values of AR with time for that flare
    
# Plot instantaneous temperatures of active regions for all flares