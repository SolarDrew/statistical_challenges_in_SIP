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
from matplotlib import patches
from sunpy import wcs
from sunpy.map import Map
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import datetime as dt
import os
from os.path import join, exists, dirname, basename, expanduser
import sys
from sys import path, argv
path.append('/imaps/holly/home/ajl7/CoronaTemps/')
from temperature import TemperatureMap as tmap
from astropy import units
import glob
from scipy.io.idl import readsav


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


density = None#'211'

parameter = 'all'#'min'
allpars = ['min', '5th %-ile', '10th %-ile', 'mean', '90th %-ile', '95th %-ile', 'max', 'stdev',
           'n over 6.1', 'n over 6.3', 'n over 6.5']
if parameter == 'all':
    pars = allpars[:7]
else:
    pars = [parameter]

home = expanduser('~')
savedir = join(home, "tempplots/{}/".format(parameter.replace(' ', '_').replace('.', '_')))
if density:
    savedir = savedir.replace('tempplots', 'densplots_{}'.format(density))
if not exists(savedir):
    os.makedirs(savedir)
data_root = '/imaps/sspfs/archive/sdo/aia/drew_ar11166/'

ar_rad = 75

starttime = parse('2011-03-05')
endtime = parse('2011-03-13')
timerange = tr(starttime, endtime)

# Get Carrington coords of active region
client = hek.HEKClient()
ars = client.query(hek.attrs.Time(starttime, endtime),
                hek.attrs.EventType('AR'),
                hek.attrs.AR.NOAANum == 11166)
coords = ars[0]['hgc_x'], ars[0]['hgc_y']
print coords

# Define times for maps
delta = dt.timedelta(minutes=2)
ntimes = int(timerange.seconds()/delta.total_seconds())
alltimes = [time.start() for time in timerange.split(ntimes)]
windows = [[]]
print len(alltimes)
for t in alltimes:
    if glob.glob(join(data_root, '{:%Y/%m/%d/*/*_%H%M_*.fits}'.format(t))) != []:
        windows[-1].append(t)
        diff = (t - windows[-1][0]).total_seconds()/60
        if diff >= 30:
            windows.append([])
print [len(times) for times in windows]

absfig, axa1 = plt.subplots(5, 2, figsize=(10, 20))
axa1 = axa1.flatten()
absfig.suptitle(parameter)

for ax, times in zip(axa1, windows):
    for t in times:
        print str(t)
    tmapfig, axt = plt.subplots(figsize=(16, 16))
    ntimes = len(times)
    thisendtime = times[-1]
    try:
        paramvals_fname = join(savedir, 'params_ar11166_{:%Y%m%d-%H%M}'.format(times[0]))
        if density:
            paramvals_fname += '_density{}'.format(density)
        print 'Loading ', paramvals_fname
        
        if not exists(paramvals_fname):
          paramvals = np.zeros((11, ntimes))
          for t, time in enumerate(times):
            # Load/calculate temperature map data
            print time
            try:
                data_dir = data_root#join(data_root, '{:%Y/%m/%d}'.format(time))
                maps_dir = join(home, "AR-tmaps/{:%Y/%m/%d}".format(time))

                thismap = tmap(time, data_dir=data_dir, maps_dir=maps_dir)#, verbose=True)#,
                #               submap=([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad]))
                thismap.save()
                # Crop temperature map to active region
                x, y = wcs.convert_hg_hpc(coords[0], coords[1],
                                          b0_deg=thismap.heliographic_latitude,
                                          l0_deg=thismap.carrington_longitude)
                print x, y
                thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])

                if density:
                    thismap = thismap.calculate_density(wlen=density)
                # Crop temperature map to active region
                #x, y = flare['xpos'], flare['ypos']
                #largemap = thismap.submap([x-200, x+200], [y-200, y+200])
                #thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])
                plt.sca(axt)
                plt.cla()
                """if density:
                    largemap.plot(vmin=np.nanmean(largemap.data) - (2*(np.nanstd(largemap.data))),
                                  vmax=np.nanmean(largemap.data) + (2*(np.nanstd(largemap.data))))
                else:
                    largemap.plot()
                rect = patches.Rectangle([x-ar_rad, y-ar_rad], ar_rad*2, ar_rad*2, color='white',
                                         fill=False)
                axt.add_artist(rect)"""
                thismap.plot()
                plt.title('Max: {:.2f}, Mean: {:.2f}, Min: {:.2f}'.format(thismap.max(), thismap.mean(), thismap.min()))
                output_dir = join('/imaps/holly/home/ajl7/AR-tmaps/AR11166_b/{:%Y/%m/%d}'.format(time))
                if density:
                    output_dir = output_dir.replace('tmaps', 'nmaps_{}'.format(density))
                if not exists(output_dir):
                    os.makedirs(output_dir)
                plt.savefig(join(output_dir, '{:%Y-%m-%dT%H%M}'.format(time)))
            
                data = thismap.data
                paramvals[:, t] = [np.nanmin(data),
                                   np.percentile(data, 5),
                                   np.percentile(data, 10),
                                   np.nanmean(data),
                                   np.percentile(data, 90),
                                   np.percentile(data, 95),
                                   np.nanmax(data),
                                   np.std(data),
                                   np.count_nonzero(data > 6.1),
                                   np.count_nonzero(data > 6.3),
                                   np.count_nonzero(data > 6.5)]
            except:
                print "Failed", time
                paramvals[:, t] = np.nan
                #np.savetxt(paramvals_fname, paramvals)
                #raise
          np.savetxt(paramvals_fname, paramvals)
        else:
          paramvals = np.loadtxt(paramvals_fname)
    
        for par in pars:
            # Rename a thing so I don't have to change a load of code.
            # I'm a bad man
            parind = allpars.index(par)
            means = paramvals[parind, :]
            print means
            if 'n over' not in parameter and 0 in means:
                print 'Anomalous data - skipping'
            # Convert time values to time before flare
            times2 = [(ti - thisendtime).total_seconds()/60 for ti in times]
            # Set current axis and plot
            plt.sca(ax)
            plt.plot(times2, means)
        print times[0], times[-1]
        print times2
    except:
        print 'Failed'
        raise

    print 'plot 1'
    plt.ylim(5.9, 6.7)

absfig.savefig(join(savedir, "ar11166"))
plt.close('all')
