# -*- coding: utf-8 -*-
"""
Created on Tue Dec 09 18:15:45 2014

@author: Drew Leonard
"""

from matplotlib import use
use('agg')
from sunpy import wcs
from sunpy.map import Map
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import datetime as dt
from repeat_tempmaps import repeat
from os.path import join

start = parse('2011-02-10')
end = parse('2011-02-16')
#end = parse('2011-01-04')

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

for flare in flares:
    flaretime = parse(flare['event_starttime'])
    starttime = flaretime-dt.timedelta(hours=0.5)
    timerange = tr(starttime, flaretime)

    fits_dir = '/imaps/sspfs/archive/sdo/aia/activeregions/'.format(starttime)
    fname = 'aia*171*{0:%Y?%m?%d}?{0:%H?%M?%S}*lev1?fits'.format(flaretime)
    filepath = join(fits_dir, fname)
    temp_im = Map(filepath)
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
    coords = wcs.convert_hg_hpc(region['hgc_x'], region['hgc_y'],
                                b0_deg=temp_im.heliographic_latitude,
                                l0_deg=temp_im.carrington_longitude)
    
    datefile = open('dates.txt', 'a')
    datefile.write('\nAR{}; {}\n'.format(flare['ar_noaanum'],
                                         flare['fl_goescls']))
    datefile.close()
    means, fluxes, maxes, times = repeat(starttime, flaretime, timeres=1.0/60.0,
                                         coords=coords, ar=flare['ar_noaanum'], 
                                         plotminmax=True)
