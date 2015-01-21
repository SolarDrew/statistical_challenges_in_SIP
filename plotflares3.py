# -*- coding: utf-8 -*-
"""
Created on Tue Dec 09 18:15:45 2014

@author: Drew Leonard
"""

from matplotlib import use
use('agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colours
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
end = parse('2011-03-01')

client = hek.HEKClient()
flares = client.query(hek.attrs.Time(start, end),
                      hek.attrs.EventType('FL'))

flares = [fl for fl in flares if (fl['ar_noaanum'] > 11137 and 
                                  fl['ar_noaanum'] < 11184)]

ar_rad = 75
ar_temps_fltime = []
ar_temps_1 = []
ar_temps_10 = []
ar_temps_30 = []
fl_classes = []

flarelist = open("/imaps/holly/home/ajl7/tempplots/flarelist.txt", "w")

fig = plt.figure(figsize=(24, 8))
ax1 = fig.add_subplot(1, 3, 1)
ax1.set_title('Temperature value over time')
ax2 = fig.add_subplot(1, 3, 2)
ax2.set_title('Running difference')
ax3 = fig.add_subplot(1, 3, 3)
ax3.set_title('Difference from temperature at flare time')

# Set up some colourmap stuff for line-plotting later
cmap = cm = plt.get_cmap('afmhot')
cNorm  = colours.Normalize(vmin=0, vmax=1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

for flare in flares:
  try:
    flaretime = parse(flare['event_starttime'])
    starttime = flaretime-dt.timedelta(hours=0.5)
    timerange = tr(starttime, flaretime)

    flarelist.write("{} {} & {} & {} \\\\ \n".format(flaretime.date(), flaretime.time(), flare['fl_goescls'], flare['ar_noaanum']))
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
    maps_root = '/imaps/sspfs/archive/sdo/aia/fulldisk/images/' #'/imaps/holly/home/ajl7/tempmaps'
    
    # Create empty lists to store temperature values and running differences
    # Called means because I did mean first, but used for max, percentiles,
    # etc sometimes
    means = []
    runningdiffT = []
    for time in times:
        # Load/calculate temperature map data
        print time
        try:
            maps_dir = join(maps_root, "{:%Y/%m/%d}/temperature/".format(time))
            thismap = tmap(time, data_dir=data_dir, maps_dir=maps_dir)
            #thismap.save()
        
            # Crop temperature map to active region
            x, y = wcs.convert_hg_hpc(region['hgc_x'], region['hgc_y'],
                                      b0_deg=thismap.heliographic_latitude,
                                      l0_deg=thismap.carrington_longitude)
            thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])
            
            # Append appropriate temperature values to list
            means.append(np.nanmean(thismap.data))
            #means.append(np.nanmax(thismap.data))
            #means.append(np.percentile(thismap.data, 95))
            if len(means) > 1:
                runningdiffT.append(means[-1]-means[-2])
            else:
                runningdiffT.append(0)
        except:
            print "Failed", time
            means.append(np.nan)
            #raise
    
    print means
    # Convert time values to time before flare
    times = [(t - flaretime).total_seconds()/60 for t in times]
    # Decide line colour based on flare class
    flarecolours = {'A': 0.3, 'B': 0.4, 'C': 0.5, 'M': 0.6, 'X': 0.7}
    colourVal = scalarMap.to_rgba(flarecolours[str(flare['fl_goescls'])[0].upper()])
    # Plot temperature values of AR with time for that flare
    ##fig = plt.figure()
    ax1.plot(times, means, color=colourVal)
    ax2.plot(times, runningdiffT, color=colourVal)
    ax3.plot(times, [mean-means[-1] for mean in means], color=colourVal) # Absolute difference
    #ax3.plot(times, [((mean-means[-1])/means[-1])*100 for mean in means], color=colorVal) # Percentage difference
    ##fname = 'AR{}/{}__{}-class'.format(flare['ar_noaanum'], flare['event_starttime'], flare['fl_goescls'])
    ##fname = fname.replace('.', '_')
    ##plt.savefig(join('/imaps/holly/home/ajl7/tempplots/', fname))
    ##plt.close()
    # Append  temperature values for final temperature map to list
    ar_temps_fltime.append(means[-1])
    ar_temps_1.append(means[-2])
    ar_temps_10.append(means[-11])
    ar_temps_30.append(means[0])
    # Append class of flare to list
    #fl_classes.append(flareclass_to_flux(str(flare['fl_goescls'])).value)
    fl_classes.append(np.log10(flareclass_to_flux(str(flare['fl_goescls'])).value))
  except:
    print 'Failed for {} flare at {}'.format(flare['fl_goescls'], flare['event_starttime'])
    #raise

flarelist.close()

#ax1.ylim(5.9, 6.3)
plt.savefig("/imaps/holly/home/ajl7/tempplots/allars")
plt.close()

"""# Plot ratio of temperatures to given temperature over time
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(2, 2, 1)
plt.title("T(t)/T(t=x)")
ratio = [mean/means[-1] for mean in means]
print min(ratio), max(ratio)
plt.plot(times, [mean/means[-1] for mean in means], label='x = 0')
plt.ylabel("Temperature ratio")
ax2 = fig.add_subplot(2, 2, 2)
plt.title("T(t)/T(t=-1)")
ratio =	[mean/means[-2] for mean in means]
print min(ratio), max(ratio)
plt.plot(times, [mean/means[-2] for mean in means], label='x = -1')
ax3 = fig.add_subplot(2, 2, 3)
plt.title("T(t)/T(t=-10)")
ratio =	[mean/means[-11] for mean in means]
print min(ratio), max(ratio)
plt.plot(times, [mean/means[-11] for mean in means], label='x = -10')
plt.ylabel("Temperature ratio")
plt.xlabel("Time (minutes)")
ax4 = fig.add_subplot(2, 2, 4)
plt.title("T(t)/T(t=-30)")
ratio =	[mean/means[0] for mean in means]
print min(ratio), max(ratio)
plt.plot(times, [mean/means[0] for mean in means], label='x = -30')
plt.xlabel("Time (minutes)")
#plt.ylim(0, 1)
#plt.legend()
plt.savefig("/imaps/holly/home/ajl7/tempplots/tempratios")
plt.close()"""

# Plot instantaneous temperatures of active regions for all flares against flare class
fig = plt.figure(figsize=(16, 16))
ax1 = fig.add_subplot(2, 2, 1)
plt.title("At time of flare")
plt.scatter(ar_temps_fltime, fl_classes)
plt.ylabel("GOES flux of flare")
#ax1.set_yscale('log')
ax2 = fig.add_subplot(2, 2, 2)
plt.title("1 minute before flare")
plt.scatter(ar_temps_1, fl_classes)
#ax2.set_yscale('log')
ax3 = fig.add_subplot(2, 2, 3)
plt.title("10 minutes before flare")
plt.scatter(ar_temps_10, fl_classes)
plt.ylabel("GOES flux of flare")
plt.xlabel("Mean temperature of active region")
#ax3.set_yscale('log')
ax4 = fig.add_subplot(2, 2, 4)
plt.title("30 minutes before flare")
plt.scatter(ar_temps_30, fl_classes)
plt.xlabel("Mean temperature of active region")
#ax4.set_yscale('log')
#plt.ylim(0.0, 0.0005)
plt.savefig("/imaps/holly/home/ajl7/tempplots/allflares")
plt.close()
