
from nustar_gen import info
ns = info.NuSTAR()
import os
import astropy.units as u

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
    from astropy.time import Time

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
    
    if tstart < Time('2019-01-01T01:00:00'):
        gind = 15
    else:
        gind = 17
    
    result = Fido.search(a.Time(tstart.fits, tend.fits), a.Instrument("XRS"),
                        a.goes.SatelliteNumber(gind))
    
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
    return

def compare_sun_spec(obs,mod='A',src_rad = 2*u.arcmin,
    en_range = (3, 20), en_bins = 34, ratio = False):
    """ 
    Uses the attorb file to generate a spectrum when NuSTAR is in/out
    of sunlight to check for solar activity.
    
    Parameters
    -----------
    obs : Class
        Output from nustar_gen.info.Observation() class
    
        
    Other Parameters
    ----------------
    mod : str
        Which module, should 'A' or 'B'. Default is 'A'
    src_rad : astropy unit
        Radius to throw away source counts. Default is 2-arcmin
    en_range : list
        Energy range to use. Default is (3, 20).
    en_bins : int
        Number of bins to use for histogram (default is 34)
    ratio : bool 
        Whether to plot the spectra or their ratio. Default is False.
   
    Returns
    -------
    None
    """
    import os
    import warnings
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
    warnings.simplefilter('ignore', category=RuntimeWarning)
    from astropy.table import Table
    from numpy import interp, histogram, sqrt
    from astropy.io.fits import getheader
    from nustar_gen.utils import chan_to_energy
    
    
    import matplotlib.pyplot as plt
    
    
    # Load the 01 event file:
    evf = obs.science_files[mod][0]
    
    
    # In the future, maybe come back and do the thing where we get rid of the
    # source and just look at the background. But we typically get rid of that info
    # 
    
    limit_pix = ((src_rad / ns.pixel).cgs).value

    # Check and see if nuskytodet has been run
    sky2det = os.path.join(obs.evdir, f'nu{obs.seqid}_sky2det{mod}.fits')
    if not os.path.isfile(sky2det):
        # Spawn this in the background:
        attfile = os.path.join(obs.evdir, f'nu{obs.seqid}_att.fits')
        mastfile = os.path.join(obs.evdir, f'nu{obs.seqid}_mast.fits')        
        hdr = getheader(obs.evtfiles[mod][0])
        cmdstring = f"nuskytodet pntra={hdr['RA_OBJ']} pntdec={hdr['DEC_OBJ']} "
        cmdstring += f"attfile={attfile} instrument=FPM{mod} skydetfile={sky2det}"
        cmdstring += f" mastaspectfile={mastfile}"
        print('No sky2det file found, running nuskytodet')
        
        sky2det_log = os.path.join(obs.evdir, f'nu{obs.seqid}_sky2det{mod}.log')
        cmdstring += f'> {sky2det_log}'
        print(cmdstring)   
        os.system(cmdstring)        
    detsrc = Table.read(sky2det, hdu =1)
    
    attorbf = os.path.join(obs.evdir, f'nu{obs.seqid}{mod}.attorb')
    assert os.path.isfile(attorbf), f'Attorb file not found: {attorbf}. Run nupipeline first.'
    
    attorb = Table.read(attorbf, hdu=1, memmap=True)
    evt = Table.read(evf, hdu='EVENTS', memmap=True)
     
    evt['SUN'] = interp(evt['TIME'], attorb['TIME'], attorb['SUNSHINE'])
    xr = interp(evt['TIME'], detsrc['TIME'], detsrc['DET1X'])
    yr = interp(evt['TIME'], detsrc['TIME'], detsrc['DET1Y'])
    evt['RAD'] = sqrt((evt['DET1X'] - xr)**2 + (evt['DET1Y'] - yr)**2)
    
    
    no_sun = chan_to_energy(evt[(evt['SUN'] == 0)&(evt['RAD']>limit_pix)]['PI'])
    sun = chan_to_energy(evt[(evt['SUN'] ==1)&(evt['RAD']>limit_pix)]['PI'])
    
    scale_no_sun = (len(no_sun) / len(evt))
    scale_sun = (len(sun) / len(evt))
    
    hist_sun, edges = histogram(sun, range = en_range, bins = en_bins)
    hist_nosun, edges = histogram(no_sun, range = en_range, bins = en_bins)
    center = (edges[:-1] + edges[1:])/2
    
    if ratio is False:
        ax = plt.figure(figsize=(12, 8)).subplots()
        ax.step(center, hist_nosun / scale_no_sun, label = 'No Sun')    
        ax.step(center, hist_sun / scale_sun, label = 'Sun')    
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        ax.set_xlabel('Energy (keV)', fontsize=12)
        ax.set_ylabel('Scaled counts', fontsize=12)
        plt.show()
    else:
        ax = plt.figure(figsize=(12, 8)).subplots()
        ax.step(center, (hist_sun / scale_sun) / (hist_nosun / scale_no_sun),
                label = 'Ratio (Sun / No Sun)')    
        ax.set_xscale('log')
        ax.set_xlabel('Energy (keV)', fontsize=12)
        ax.set_ylabel('Ratio', fontsize=12)
        ax.legend()
        ax.axhline(1.0, color = 'green')
        plt.show()
    
    return
    

    