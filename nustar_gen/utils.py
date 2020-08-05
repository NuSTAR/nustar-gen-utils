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
    return en

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
            new_gti = gti_stub
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
    
    from regions import read_ds9, PixCoord
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


    reg = read_ds9(regfile)
    hdr = fits.getheader(evf)
    
    # exp = hdr['EXPOSURE']
    
    # Okay, so the DET1 exposure map is given in *clock* seconds spent in different
    # configurations, and not livetime. So we need to track the total duration
    # below via the "DURATION" keyword in the extensions.
    # Corrected after talking with KKM.
    exp = 0.
    
    ns = info.NuSTAR()
    det1_pixarea = (ns.pixel_um.to(u.cm))**2
    
    hdu = fits.open(det1im)
    exp_area = 0
    for ii, ihdu in enumerate(hdu):
        if ii == 0:
            continue    
        exp += ihdu.header['DURATION']
        expim = ihdu.data
        
         # Find the include regions first:
        # Loop over regions and find those that are include first:
        set=False

        for ri in reg:
            if (ri.meta['include']) is False:
                continue
            if set is False:
                all_reg=ri
                set=True
            else:
                # OR inclusive logic to add multiple regions if you want to
                all_reg = all_reg.union(ri)

        source_mask = all_reg.to_mask()
#        exp_area += source_mask.multiply(ihdu.data).sum()

        # Now loop over the excluded regions:
        set = False
        for ri in reg:
            if (ri.meta['include']) is True:
                continue
            if set is False:
                excl_reg = ri
                set = True
            else:
                excl_reg = excl_reg.union(ri)
        if set is True:
            merge_region = all_reg.intersection(excl_reg)
        else:
            merge_region = all_reg
        
        # Brute force a mask, because below was being weird
        coords = np.array(np.meshgrid(range(360), range(360)))
        pix = PixCoord(coords[0, :], coords[1, :])
        mask = merge_region.contains(pix)        
        exp_area += ihdu.data[mask].sum()
        
#            excl_area = excl_overlap.to_mask().multiply(ihdu.data).sum()
#            exp_area -= excl_area
        
    hdu.close()
    area = (exp_area * det1_pixarea) / exp
    return area
    


def make_straylight_arf(det1im, regfile, filt_file, mod, outpath=None):
    '''
    Produces an ARF for a given observation. Currently uses counts-weighting to
        determine the contribution of the detabs parameters for each detector.
        SHOULD PROBABLY BE MOVED TO WRAPPERS OR SOMETHING?
        
    Parameters
    ----------
    det1im : str
        Full path to the DET1 exposure image file.
        
    regfil : str
        Full path to the region file, assumed to be a ds9 region file in "IMAGE"
        coordinates
    
    filt_file : str
        Full path to the DET1 screened event file (i.e. output of
        nustar_gen.wrappers.extract_det1_events() ). This is only used for pulling out
        the exposure in straylight_area.

    mod : str
        'A' or 'B'

    Other parameters:
    -----------------
    
    outpath : str
        Output directory. If not specified, places this in the same location as evf

    Returns
    -------
    out_arf: string
        Full path to the new arf file
    
    '''   
    
    from astropy.io import fits
    import glob
    from nustar_gen import io
    
        # Check to see that all files exist:
    assert os.path.isfile(det1im), \
        f'straylight_area: DET1 exposure map not found: {det1im}'
    assert os.path.isfile(regfile), \
        f'straylight_area: Region file not found: {regfile}'
    assert os.path.isfile(filt_file), \
        f'straylight_area: Event file not found: {evf}'

    if outpath is None:
        outdir = os.path.dirname(evf)
    else:
        outdir=outpath
    
    out_arf = outdir+'/'+os.path.splitext(os.path.basename(filt_file))[0]+'.arf'

    # Get the effective are target:
    area = straylight_area(det1im, regfile, filt_file)
    
    # Use counts relative scaling to blend the ARFs:
    ev_hdu = fits.open(filt_file)
    evts = ev_hdu[1].data
    dethist, detedges = np.histogram(evts['DET_ID'], range=[0, 3],bins = 4)
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


    caldb = os.environ['CALDB']
    detabs_files = glob.glob(caldb+f'/**/nu{mod}detabs*', recursive=True)
    use = len(detabs_files)-1

    detabs_hdu = fits.open(detabs_files[use])

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

    arf['SPECRESP'] *= be_arf['SPECRESP']

    arf_hdu.writeto(out_arf, overwrite=True)
    arf_hdu.close()
    
    return out_arf
    
def validate_det1_region(regfile):
    """
    Code for making sure that region files that you're trying to use on the DET1
    analysis code are in "image" coordinates
    """
    err=-1
    from regions.io.ds9.read import DS9Parser
    assert os.path.isfile(regfile), f'{regilfe} does not exist!'
    
    with open(regfile) as fh: 
        region_string = fh.read()
    parser = DS9Parser(region_string)
    assert parser.coordsys == 'image', f'Region coordinate system is {parser.coordsys}, not image!'
