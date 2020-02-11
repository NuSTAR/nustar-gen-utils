import os
import warnings
import numpy as np



def energy_to_chan(keV):
    '''
    Convert keV to NuSTAR PI channels
    PI = (keV - 1.6) / 0.04
    
    chan = energy_to_chan(3.0)
    
    
     Example
    -------
    >>> from nustar_gen.utils import energy_to_chan
    >>> chan = energy_to_chan(10.)
    >>> np.isclose(chan, 210)
    True

    '''
    return int((keV - 1.6) / 0.04)


def make_usr_gti(input_gtis, outfile='usrgti.fits', **kwargs):
    '''
    Utility script to make a usrgti file. Uses the template distributed with the
    repo.
    
    Parameters
    ----------
    gtis :  dict
        Is assumed to have a 'tstart' and 'tend' key.
        Can be an array of entries.

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
    module to compute the mask.
    
    Parameters
    ----------
    det1im : str
        Full path to the DET1 exposure image file.
        
    regfil : str
        Full path to the region file, assumed to be a ds9 region file in "IMAGE"
            coordinates
    
    evf : str
        Full path to the event file. Only used for pulling down the base exposure in
            the FITS header.

    '''
    
    from regions import read_ds9
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
    exp = hdr['EXPOSURE']
    
    ns = info.NuSTAR()
    det1_pixarea = (ns.pixel_um.to(u.cm))**2
    
    hdu = fits.open(det1im)
    exp_area = 0
    for ii, ihdu in enumerate(hdu):
        if ii == 0:
            continue    
        expim = ihdu.data
        for ri in reg:
            mask = ri.to_mask()
            if ri.meta['include']:
                exp_area += mask.cutout(ihdu.data).sum()
            else:
                exp_area -= mask.cutout(ihdu.data).sum()
    hdu.close()
    area = (exp_area * det1_pixarea) / exp
    return area
    


def make_straylight_arf(det1im, regfile, filt_file, mod, outpath=None):
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
        nustar_gen.wrappers.extract_det1_events() ).

    mod : str
        'A' or 'B'

    Other parameters:
    -----------------
    
    outpath : str
        Output directory. If not specified, places this in the same location as evf

    Returns
    -------
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
