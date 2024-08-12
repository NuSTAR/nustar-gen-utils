import os
import warnings
import astropy.units as u

def find_be_arf():
    '''
    Find and return the ARF for the Be window (for straylight analysis only)
    
    Parameters
    ----------
    None
    
    Returns
    -------
    be_arf_file : str
        Full path to the Be window file in the CALDB
    
    Examples
    --------
    >>> be_file = find_be_arf()
    >>> print(os.path.isfile(be_file))
    True
    
    '''
    
    # Get the CALDB path:
    caldb=os.environ["CALDB"]    
    subdir='data/nustar/fpm/bcf/instrument'
    be_arf_file = os.path.join(caldb, subdir)
    be_file = 'nuCabsparBe20100101v001.fits'
    be_arf_file = os.path.join(be_arf_file, be_file)
    
    assert os.path.isfile(be_arf_file), \
        f'find_be_arf: Check your CALDB. Be ARF file not found {be_arf_file}'
    
    
#    curdir = os.path.dirname(__file__)
#    tempdir = os.path.join(curdir, 'templates')
#    be_arf_file = os.path.join(tempdir, 'be.arf')
    return be_arf_file
    
def find_arf_template():
    '''
    Find and return the template ARF file.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    arf_file : str
        Full path to the Be window file, which should be included when this repo is checked out
    
    Examples
    --------
    >>> arf_file = find_arf_template()
    >>> print(os.path.isfile(arf_file))
    True
    
    '''
    curdir = os.path.dirname(__file__)
    tempdir = os.path.join(curdir, 'templates')
    arf_file = os.path.join(tempdir, 'nu40101010002_srcA_sr.arf')
    return arf_file
    
def write_ds9(coord, rad=60*u.arcsec, outfile='ds9.reg'):
    '''
    Exports a ds9 region file for a circular extraction region
    based on an Astropy coordinates object.
    
    Parameters
    ----------
    coord : Astropy Coordinates object
        Coordinates for the source

    rad : Astropy unit, default value: 60*u.arcsec

    outfile : str, optional, default='ds9.reg'
        Full path to the output region file
    
    
    '''
    fk5 = coord.fk5
    with open(outfile, 'w') as f:
        f.write('# Region file format: DS9\n')
        f.write('fk5\n')
        
        regstring = f'circle({fk5.ra.deg}, {fk5.dec.deg}, {rad.to(u.arcsec).value}")\n'
        f.write(regstring)

    return    
    