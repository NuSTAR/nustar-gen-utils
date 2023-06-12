
from nustar_gen import info
ns = info.NuSTAR()
import os

def goes_lightcurve(obs, show_sun=False, show_sky=False, show_impact=True):
    """
    Module for producing a GOES lightcurve and the NuSTAR in/out of Sun times
    
    Parameters
    -----------
    obs : Class
        Output from nustar_gen.info.Observation() class
    
        
    Other Parameters
    ----------------
    show_sun : bool
        Show just the Sun period (default is False)
    show_sky : bool
        Show the unocculted periods (default is false)
    show_impact : bool
        Show the periods where you are in Sun and the source is unocculted.
        Default is True
        Note that setting either show_sun or show_sky prevents show_impact

    Returns
    -------
    None

    """
    
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.patches import Rectangle
    from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
    import numpy as np
    from astropy.visualization import time_support
    from astropy.io.fits import getheader, getdata

    try:
        from sunpy import timeseries as ts
        from sunpy.net import Fido
        from sunpy.net import attrs as a
    except:
        print('This module requires SunPy. Please install this first.')
        return
        
    # Get the header info
    evf = obs.evtfiles['A'][0]
    hdr = getheader(evf)
    tstart = (ns.met_to_time(hdr['TSTART']))
    tend = (ns.met_to_time(hdr['TSTOP']))
    
    result = Fido.search(a.Time(tstart.fits, tend.fits), a.Instrument("XRS"))
    files = Fido.fetch(result, progress=False)
    goes_all = ts.TimeSeries(files, concatenate=True)
    goes = goes_all.truncate(tstart.iso, tend.iso)



    fig, ax = plt.subplots(figsize=(12, 8))
    goes.plot(columns=["xrsb"])
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d %H:%M')) 
    fig.autofmt_xdate() 



    # Okay, now we want to add on some NuSTAR info.
    
    attorb_file = os.path.join(obs.evdir, f'nu{obs.seqid}A.attorb')
    assert os.path.isfile(attorb_file), f'make_goes_json_image.py: No attorb file {attorb_file}'

    attorb, hdr = getdata(attorb_file, 1, header=True)

    obs_start = ns.met_to_time(hdr['TSTART'])
    obs_end = ns.met_to_time(hdr['TSTOP'])

    ax.set_ylim(1e-7, 1e-3)
    # Do the sunshine part
    if show_sun:
        show_impact=False
        in_sun = (np.where(attorb['SUNSHINE'] ==1))[0]
        diff = (np.diff(in_sun))
        transition = (np.where(diff > 1))[0]
        # Check to see if you started in sun:

        lims = ax.get_ylim()
        height = lims[1] - lims[0]
        for ind, edge in enumerate(transition):
            if (ind == 0):
                if (in_sun[0] == 0):        
                    left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][0]).datetime)
                    right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
                else:
                    left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[0]]).datetime)
                    right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
            else:
        
                right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
                left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][left_edge_ind]).datetime)
    
            prev_right_edge = in_sun[edge]
            left_edge_ind = prev_right_edge + diff[edge]
        
            width = right_edge - left_edge
            rect = Rectangle((left_edge, lims[0]), width, 1, color='yellow', alpha = 0.5)
            ax.add_patch(rect)
    
        # Now get the last one:
        left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][left_edge_ind]).datetime)
        right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun].max()).datetime)
        rect = Rectangle((left_edge, lims[0]), width, 1, color='yellow', alpha = 0.5)
        ax.add_patch(rect)


    if show_sky:
        show_impact=False

        # Now do the occulted / unocculted part
    
        ## Add the unocculted periods ##
        in_sun = (np.where(attorb['ELV'] > 3))[0]
        diff = (np.diff(in_sun))
        transition = (np.where(diff > 1))[0]
        # Check to see if you started in sun:

        lims = ax.get_ylim()
        height = lims[1] - lims[0]
        for ind, edge in enumerate(transition):
            if (ind == 0):
                if (in_sun[0] == 0):        
                    left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][0]).datetime)
                    right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
                else:
                    left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[0]]).datetime)
                    right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
            else:
        
                right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
                left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][left_edge_ind]).datetime)
    
            prev_right_edge = in_sun[edge]
            left_edge_ind = prev_right_edge + diff[edge]
        
            width = right_edge - left_edge
            rect = Rectangle((left_edge, lims[0]), width, 1, color='blue', alpha = 0.2)
            ax.add_patch(rect)
    
        # Now get the last one:
        left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][left_edge_ind]).datetime)
        right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun].max()).datetime)
        rect = Rectangle((left_edge, lims[0]), width, 1, color='blue', alpha = 0.2)
        ax.add_patch(rect)

    if show_impact:
        ## Show both in-sun and unocculted
        in_sun = (np.where((attorb['ELV'] > 3)&(attorb['SUNSHINE'] ==1)))[0]
        diff = (np.diff(in_sun))
        transition = (np.where(diff > 1))[0]
        # Check to see if you started in sun:

        lims = ax.get_ylim()
        height = lims[1] - lims[0]
        for ind, edge in enumerate(transition):
            if (ind == 0):
                if (in_sun[0] == 0):        
                    left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][0]).datetime)
                    right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
                else:
                    left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[0]]).datetime)
                    right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
            else:
        
                right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun[edge]]).datetime)
                left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][left_edge_ind]).datetime)
    
            prev_right_edge = in_sun[edge]
            left_edge_ind = prev_right_edge + diff[edge]
        
            width = right_edge - left_edge
            rect = Rectangle((left_edge, lims[0]), width, 1, color='orange', alpha = 0.2)
            ax.add_patch(rect)
    
        # Now get the last one:
        left_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][left_edge_ind]).datetime)
        right_edge = mdates.date2num(ns.met_to_time(attorb['TIME'][in_sun].max()).datetime)
        rect = Rectangle((left_edge, lims[0]), width, 1, color='orange', alpha = 0.2)
        ax.add_patch(rect)

    plt.show()
            
        
    