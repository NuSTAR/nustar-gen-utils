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
	gtis : 	dict
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
		
	
	fits.writeto(outfile, gti, new_hdr, **kwargs)

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