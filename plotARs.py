# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 15:41:43 2015

@author: Drew Leonard
"""

from matplotlib import use
use('agg')
from sunpy import wcs
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import datetime as dt
from sys import path
path.append('/home/drew/CoronaTemps/')
from temperature import TemperatureMap as tmap


start = parse('2011-02-15 01:45')
end = parse('2011-02-16 02:15')

client = hek.HEKClient()
flares = client.query(hek.attrs.Time(start, end),
                      hek.attrs.EventType('FL'))

flare = [fl for fl in flares if (fl['ar_noaanum'] == 11158)][0]

flaretime = parse(flare['event_starttime'])
starttime = flaretime-dt.timedelta(hours=6)
timerange = tr(starttime, flaretime)

region = client.query(hek.attrs.EventType('AR'),
                      hek.attrs.Time(flaretime-dt.timedelta(minutes=5), 
                                     flaretime),
                      hek.attrs.AR.NOAANum == 11158)

if isinstance(region, list): region = region[1]

print region['ar_noaanum']

# Define times for maps
delta = dt.timedelta(hours=2)
ntimes = int(timerange.seconds()/delta.total_seconds())
times = [time.start() for time in timerange.split(ntimes)]

for time in times:
    # Load/calculate temperature map data
    thismap = tmap(time)
    thismap.save()

    # Crop temperature map to active region
    x, y = wcs.convert_hg_hpc(region['hgc_x'], region['hgc_y'],
                              b0_deg=thismap.heliographic_latitude,
                              l0_deg=thismap.carrington_longitude)
    #thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])
    thismap = thismap.region_map(region, mapradius=75)
    
    # Plot temperature map with it's 171A counterpart and a context image
    thismap.compare('171', context_wlen='171', 
                    savedir='/home/drew/statistical_challenges_in_SIP/')
