import os
import warnings


def find_be_arf():
    '''
    Find and return the ARF for the Be window (for straylight analysis only)
    
    Parameters
    ----------
    None
    
    Returns
    -------
    Full path to the Be window file in the CALDB
    
    Examples
    --------
    >>> be_file = find_be_arf()
    >>> print(os.path.isfile(be_file))
    True
    
    '''
    
    # Get the CALDB path:
    caldb=os.environ["CALDB"]
    print(caldb)
    
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
    Full path to the Be window file, which should be included when this repo is
    checked out
    
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