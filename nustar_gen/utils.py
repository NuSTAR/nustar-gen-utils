import os
import warnings
import numpy as np



def energy_to_chan(keV):
    '''
    Convert keV to NuSTAR PI channels: PI = (keV - 1.6) / 0.04
    
    Parameters
    ----------
    
    keV: float
        Unitless float giving the energy you want to convert to channels.
    
    Example
    ---------
    >>> from nustar_gen.utils import energy_to_chan
    >>> chan = energy_to_chan(10.)
    >>> np.isclose(chan, 210)
    True

    '''
    return int((keV - 1.6) / 0.04)

def chan_to_energy(chan):
    '''
    Convert NuSTAR PI channels to keV: keV = PI * 0.04 + 1.6
    
    
    Parameters
    ----------
    
    channel: int
        Integer PI channel to be converted into keV    
    
    Example
    -------
    >>> from nustar_gen.utils import chan_to_energy
    >>> chan = chan_to_energy(210)
    >>> np.isclose(chan, 10.)
    True

    '''
    try:
        en = [float(x) *  0.04 + 1.6 for x in chan]
    except:
        en = chan * 0.04 + 1.6
    return np.array(en)

def make_usr_gti(input_gtis, outfile='usrgti.fits', **kwargs):
    '''
    Utility script to make a usrgti file. Uses the template distributed with the
    repo.
    
    Parameters
    ----------
    gtis :  dict
        Dictionary assumed to have a 'START' and 'STOP' keys.
        Can be an array of dicts.
    
    Other parameters
    ----------
    outfile: str
        Full path to desired outfile location. Default is 'usrgti.fits' in current
        working directory


    '''
    import astropy.io.fits as fits
    from astropy.time import Time
    
    # Check if the keys exists:
    assert('START' in input_gtis[0]), "utils.make_usr_gti: gtis missing tstart key"
    assert('STOP' in input_gtis[0]), "utils.make_usr_gti: gtis missing tend key"
    
    curdir = os.path.dirname(__file__)
    tempdir = os.path.join(curdir, 'templates')
    tempfile = os.path.join(tempdir, 'usrgti.fits')
    gti, hdr = fits.getdata(tempfile, 1, header=True)
    
    
    # Time header hdr:
    new_hdr = hdr[0:17]
    isot = Time.now().fits.split('.')[0]
    new_hdr.set('DATE', isot, 'USRGTI creation date')
    
    for copy_key in ['MJDREFI', 'MJDREFF', 'HDUCLASS', 'HDUCLAS1', 'HDUCLAS2']:
        new_hdr.set(copy_key, hdr[copy_key])
    
    # Now actually update the GTIs
    gti_stub = gti[0:1]
    for ind, pair in enumerate(input_gtis):
        gti_stub['START'] = pair['START']
        gti_stub['STOP'] = pair['STOP']
        if ind == 0:
            new_gti = gti_stub.copy()
        else:
            new_gti = append_fits_entry(new_gti, gti_stub)
        
    
    fits.writeto(outfile, new_gti, new_hdr, **kwargs)

    return


def append_fits_entry(base_rec, new_entry):
    '''
    Helper script to append an entry to a FITS_rec object.
    Based on astropy documentation here:
    https://docs.astropy.org/en/stable/io/fits/usage/table.html

    '''
    from astropy.io.fits import BinTableHDU
    
    old_rows = base_rec.shape[0]
    new_rows = old_rows + 1
    new_rec = hdu = BinTableHDU.from_columns(base_rec.columns, nrows=new_rows)
    new_rec.data[old_rows:] = new_entry
    return new_rec.data

def straylight_area(det1im, regfile, evf):
    '''
    Utility script to compute the area for a given region file. Uses the astropy.regions
    module to compute the mask. Allows for overlapping include/reject regions and multiple
    source regions
    
    Parameters
    ----------
    det1im : str
        Full path to the DET1 exposure image file.
        
    regfil : str
        Full path to the region file, assumed to be a ds9 region file
        in "IMAGE" coordinates
    evf : str
        Full path to the event file. Only used for pulling down
        the base exposure in the FITS header.

    Returns
    --------
    area: astropy-units float
        Illuminated detector area with units of cm2.

    '''
    
    from regions import Regions, PixCoord
    from astropy.io import fits
    from nustar_gen import info
    from astropy import units as u
    

    # Check to see that all files exist:
    assert os.path.isfile(det1im), \
        f'straylight_area: DET1 exposure map not found: {det1im}'
    assert os.path.isfile(regfile), \
        f'straylight_area: Region file not found: {regfile}'
    assert os.path.isfile(evf), \
        f'straylight_area: Event file not found: {evf}'


    reg = Regions.read(regfile)
    hdr = fits.getheader(evf)
    
    # exp = hdr['EXPOSURE']
    
    # Okay, so the DET1 exposure map is given in *clock* seconds spent in different
    # configurations, and not livetime. So we need to track the total duration
    # below via the "DURATION" keyword in the extensions.
    # Corrected after talking with KKM.
    exp = 0.
    
    ns = info.NuSTAR()
    det1_pixarea = (ns.pixel_um.to(u.cm))**2
    
    # Define pixel coordinates for use below    
    coords = np.array(np.meshgrid(range(360), range(360)))
    pix = PixCoord(coords[0, :], coords[1, :])


    # Loop over every exposure map HDU
    hdu = fits.open(det1im)
    exp_area = 0
    for ii, ihdu in enumerate(hdu):
        if ii == 0:
            continue    
        exp += ihdu.header['DURATION']
        expim = ihdu.data
        
        # Find the include regions first. Should be agnostic about what order you
        # put the regions in.
        
        # Loop over regions and find those that are include first:
        set=False

        for ri in reg:
            if (ri.meta['include']) is False:
                continue
            if set is False:
                # Generate a base mask:
                # Defaults to False
                all_mask = np.array(expim).astype(bool)
                all_mask[:] = False
                set=True
        
            # Turn on any pixels in the mask for the source
            rm = ri.contains(pix)
            all_mask = all_mask | rm

        # Mask out the exclusion regions:
        set = False
        for ri in reg:
            if (ri.meta['include']) is True:
                continue
            if set is False:
                # Defaults to False
                reg_mask = np.array(expim).astype(bool)
                reg_mask[:] = True
                set=True
            
            # Will return "True" for regions that are *not* in the exclusion region.
            rm = ri.contains(pix)
            reg_mask = reg_mask & rm
    
        # If you had an exclusion region, do an & since reg_mask will
        # have True for regions NOT in the exclusion region
        if set is True:
            # Invert reg_mask
            src_mask = all_mask & reg_mask
        else:
            src_mask = all_mask

        exp_area += expim[src_mask].sum()
        
    hdu.close()
    area = (exp_area * det1_pixarea) / exp
    return area
    


def make_straylight_arf(det1im, regfile, filt_file, mod, obs):
    '''
    Produces an ARF for a given observation. Currently uses counts-weighting to
        determine the contribution of the detabs parameters for each detector.
        
    Parameters
    ----------
    det1im : str
        Full path to the DET1 exposure image file.
        
    regfil : str
        Full path to the region file, assumed to be a ds9 region file in "IMAGE"
        coordinates
    
    filt_file : str
        Full path to the DET1 screened event file (i.e. output of
        nustar_gen.wrappers.extract_det1_events() ). This used to scale the ARF by
        the number of 3--10 keV counts on each detector.

    mod : str
        'A' or 'B'
        
    obs: nustar_gen.info.Observation
        Observation meta data


    Returns
    -------
    out_arf: string
        Full path to the new arf file
    
    '''   
    
    from astropy.io import fits
    import glob
    from nustar_gen import io
    from nustar_gen.utils import energy_to_chan
    import subprocess
    
        # Check to see that all files exist:
    assert os.path.isfile(det1im), \
        f'straylight_area: DET1 exposure map not found: {det1im}'
    assert os.path.isfile(regfile), \
        f'straylight_area: Region file not found: {regfile}'
    assert os.path.isfile(filt_file), \
        f'straylight_area: Event file not found: {filt_file}'

    outdir = obs.out_path
    
    out_arf = os.path.join(outdir, os.path.splitext(os.path.basename(filt_file))[0]+'.arf')

    # Get the effective are target:
    area = straylight_area(det1im, regfile, filt_file)
    
    # Use counts relative scaling to blend the ARFs:
    ev_hdu = fits.open(filt_file)
    evts = ev_hdu[1].data
    scale_evts = evts[ (evts['PI']>energy_to_chan(3.)) & (evts['PI'] < energy_to_chan(10.) )]
    dethist, detedges = np.histogram(scale_evts['DET_ID'], range=[0, 3],bins = 4)
    totcts = dethist.sum()
    scale = [ float(x) / totcts for x in dethist]
    ev_hdu.close()

    # Load the base ARF file:
    arf_file = io.find_arf_template()
    arf_hdu = fits.open(arf_file)
    arf = arf_hdu[1].data
    arf_header = arf_hdu[1].header
    # Trim header
    arf_hdu[1].header = arf_header[0:26]

    # Mock up an ARF that's flat with just the detector area
    arf['SPECRESP'] = [area.value for x in arf['SPECRESP']]


    cmdstring=f"quzcif nustar FPM{mod} DET0 - DETABS - - -"
    cmds = cmdstring.split()
    ret = subprocess.check_output(cmds)    
    detabs_file = (ret.split())[0]

    detabs_hdu = fits.open(detabs_file)

    for detind, isc in enumerate(scale):
        detabs = detabs_hdu[detind+1].data
        if detind == 0:
            arfabs = detabs['DETABS'] * isc
        else:
            arfabs += detabs['DETABS'] * isc
    arf['SPECRESP'] *= arfabs
    detabs_hdu.close()

    be_file = io.find_be_arf()
    be_arf = fits.getdata(be_file, 1)

    arf['SPECRESP'] *= be_arf['ATT']

    arf_hdu.writeto(out_arf, overwrite=True)
    arf_hdu.close()
    
    return out_arf
    
def validate_det1_region(regfile):
    """
    Code for making sure that region files that you're trying to use on the DET1
    analysis code are in "image" coordinates
        
    Parameters
    ----------
    regfile : str
        Full path to a ds9 region file
        


    Returns
    -------
    None
    
    """
    err=-1
    import regions
#    from regions.io.ds9.read import DS9Parser
    from regions import Regions
    assert os.path.isfile(regfile), f'{regfile} does not exist!'
    
#     with open(regfile) as fh: 
#         region_string = fh.read()
#     parser = DS9Parser(region_string)
#     assert parser.coordsys == 'image', \
#         f'Region coordinate system is {parser.coordsys}, not image!'

    reg = Regions.read(regfile)


    # Check and make sure this is a "pixel" region and not a "sky" region

    assert 'Pixel' in f'{type(reg[0])}', \
         f'Region coordinate system is not image coordinates for {regfile}\n'

    # Check to make sure tha the first region in the file is an "include" region
    for ri in reg:
        assert ri.meta['include'] == 1, \
            f'\n {regfile} has an exclusion region first! \n Put the source region first instead!'
        break

    return

def straylight_background(det1im_file='None', sky2det_file='None',
                          reg_file = 'None', det1_expo_file = 'None',
                          diag=False, src_rad = None):
    """
    Computes the average background rate in an input energy band. 
    
    Parameters
    ----------
    detim_file : str
        Path to the input DET1 image file. This should be the output of of 
        wrappers.make_det1_image
    
    sky2det_file : str
        Path to the sky2det file. This is generated by nusky2det and traces the
        source location as a function of time in DET1 and DET2 coordinates. Used for
        removing the primary source
    
    reg_file : str
        Region file for the image. Can have both source and exclusion regions. All regions
        listed will be treated as "exclusions" assuming that they contain all of the
        non-background counts in the image.
    
    det1_expo_file : str
        Path to the DET1 effective area file. Should be the DET1 output of
        wrappers.make_exposure_map 
        
    diag : bool
        Show diagnostic information and make some plots showing what is happening.
        Default is False
        
    src_rad : astropy.units
        Radius to exclude around the source region. Default is 2-arcmin.
    """
    
    from regions import Regions, PixCoord, CirclePixelRegion
    from astropy.io import fits
    from nustar_gen import info
    from astropy import units as u
    
    # Set up pixel coordinate grid:
    coords = np.array(np.meshgrid(range(360), range(360)))
    pix = PixCoord(coords[0, :], coords[1, :])
    


    assert os.path.isfile(sky2det_file) is True, \
        f'\n {sky2det_file} does not exist.'
    src = fits.getdata(sky2det_file)

    assert os.path.isfile(det1im_file) is True, \
        f'\n {det1im_file} does not exist.'
    im, im_hdr = fits.getdata(det1im_file, header=True)

    assert os.path.isfile(reg_file) is True, \
        f'\n {reg_file} does not exist.'
    reg = Regions.read(reg_file)


    assert os.path.isfile(det1_expo_file) is True, \
        f'\n {det1_expo_file} does not exist.'


    ns = info.NuSTAR()
    det1_pixarea = (ns.pixel_um.to(u.cm))**2
    exp_area = 0
    exp = 0
    if src_rad is None:
        radius = 2.0*u.arcmin
    else:
        radius = src_rad
    rad_pix = (radius / ns.pixel).cgs.value

    with fits.open(det1_expo_file) as hdu:

        # Hard work goes here
        start=True
        for ii, ihdu in enumerate(hdu):
            if ii == 0:
                continue
            exp += ihdu.header['DURATION']
            expim = ihdu.data
    
            if  start:
                all_mask = np.array(expim).astype(bool)
                all_mask[:] = True
                start=False
        
            # Find the source DET1X/DET1Y at this time:
            sx = np.interp(ihdu.header['REF_TIME'], src['TIME'], src['DET1X'])
            sy = np.interp(ihdu.header['REF_TIME'], src['TIME'], src['DET1Y'])
    
            # Construct a pixel region
            src_reg = CirclePixelRegion(PixCoord(x=sx, y=sy), radius=rad_pix)
    
            # Loop over regions and find those that are include first:
            set=False

            for ri in reg:
                if (ri.meta['include']) is False:
                    continue
                # Turn on any pixels in the mask for the source
                rm = ri.contains(pix)
                all_mask = all_mask & ~rm


            # Turn on any pixels in the mask for the source
            rm = src_reg.contains(pix)
            all_mask = all_mask & ~rm

    
        
        #     # Mask out the exclusion regions:
            set = False
            for ri in reg:
                if (ri.meta['include']) is True:
                    continue
                if set is False:
                    # Defaults to False
                    reg_mask = np.array(expim).astype(bool)
                    reg_mask[:] = True
                    set=True

        #         # Will return "True" for regions that are *not* in the exclusion region.
                rm = ri.contains(pix)
                reg_mask = reg_mask & rm

            if set is True:
                # Invert reg_mask
                src_mask = all_mask & reg_mask
            else:
                src_mask = all_mask
            

        exp_area = (expim*src_mask).sum() / expim.max()
        area = (exp_area * det1_pixarea)

    counts = (im*src_mask).sum() * u.ct * 1.0
    rate= counts / (area *im_hdr['EXPOSURE']*u.s)



    if diag is True:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from scipy import ndimage
        axs = plt.figure(figsize=(10, 10)).subplots(nrows=2, ncols=2)
        axs[0,0].imshow(src_mask, origin = 'lower')
        axs[0, 0].set_title('Mask')
        axs[0,1].imshow(expim*src_mask, origin='lower') 
        axs[0, 1].set_title('Exposure Map')
        axs[1, 0].imshow(im, origin = 'lower',
                        norm = LogNorm(vmin=0.1, vmax=im.max()*2))
        axs[1, 0].set_title('Counts Image')
        axs[1, 1].imshow(im*src_mask, origin='lower',
                        norm = LogNorm(vmin=0.1, vmax=(im).max()*2))
        axs[1, 1].set_title('Masked Counts Image')
        plt.show()
        print(f'straylight_background:')
        print(f'Background area: {area:8.2f}')
        print(f'Background rate: {rate:8.5f}')
        print('')
    return rate






