import os
import glob
from skimage.feature import peak_local_max
from scipy.optimize import curve_fit
import scipy.ndimage as ndi
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
import numpy as np
import astropy.units as u
import warnings
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning, append=True)

def load_psf(rind, energy=4, filt_range = 5, show=False):
    '''
    Load the PSF from the CALDB and generate a radial profile
    '''
    caldb = os.environ['CALDB']
    psf_files = glob.glob(caldb+'/**/*2dpsfen*', recursive=True)
    psf_file = psf_files[energy]
    psf = fits.getdata(psf_file, 4)

    dr = rind[1] - rind[0]
    psf2 = ndi.maximum_filter(psf, size=5)
    psf_coord = peak_local_max(psf, min_distance=filt_range, threshold_abs=0.5*psf2.max())
    psf_profile_raw, psf_err = azimuthalAverage(psf, center= [psf_coord[0][1], psf_coord[0][0]])
    psf_rind = np.arange(len(psf_profile_raw))*dr

    psf_profile = np.interp(rind, psf_rind, psf_profile_raw)
    
    if show:
        fig = plt.figure(figsize=[12, 4])
        ax = fig.add_subplot(111)
        ax.imshow(psf, origin='lower', cmap=plt.cm.viridis)
        plt.show()

    
    return psf_profile
    


def profile_model(rind, psf_scale, bgd_offset, bgd_slope):
    '''
    Model for fitting the radial profile by scaling the PSF and then adding a background
    term that's assumed to be flat in cts / area
    
    '''
    return psf_scale * load_psf(rind, energy = 3) + bgd_offset + bgd_slope*rind
    
    
    

def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
             
    Returns - the radial profile in units of cts / area.
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr
    radial_err = np.sqrt(tbin) / nr

    return radial_prof, radial_err

def find_source(filename, show_image=False, filt_range =10, gauss_sig=3):
    '''
    Find a source in the image. Returns the location in pixel_coordiantes.
    
    filename should point to a file produed by nustar_gen.wrappers.make_image()
    
    '''
    
    hdu = fits.open(filename, uint=True)[0]
    wcs = WCS(hdu.header)

    im = hdu.data  
    im2 = ndi.gaussian_filter(im.astype(float), sigma=gauss_sig)
    im3 = ndi.maximum_filter(im2, size=filt_range)
    coordinates = peak_local_max(im2, min_distance=filt_range, threshold_abs=0.75*im2.max())
    
    
    if len(coordinates) == 0:
        print('No sources found')
    
    if show_image:
        fig = plt.figure(figsize=[12, 4])
        ax = fig.add_subplot(131, projection=wcs)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)





        ax.imshow(im, origin='lower', cmap=plt.cm.viridis)
        ax.set_title('Raw Image')

        imrange = [250, 750]
        
        ax2.imshow(im2, origin='lower', cmap=plt.cm.viridis)
        ax2.set_title('Smoothed Image')
        ax2.axis('off')
        
        ax3.imshow(im3, origin='lower', cmap=plt.cm.viridis)

        ax3.set_title('Filtered Image')
        ax2.axis('Off')
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        
        if len(coordinates) > 0:
            for coord in coordinates:
                ax.plot(coord[1], coord[0], 'rx')
                ax2.plot(coord[1], coord[0], 'rx')
                ax3.plot(coord[1], coord[0], 'rx')

        ax2.set_xbound(imrange)
        ax2.set_ybound(imrange)
        ax.set_xbound(imrange)
        ax.set_ybound(imrange)
        ax3.set_xbound(imrange)
        ax3.set_ybound(imrange)

        plt.show()
        
    # To switch back to [X, Y] nomenclature
    return coordinates



def make_radial_profile(filename, show_image=False, filt_range = 5, 
                       coordinates=False):
    '''
    
    Makes the radial profile.
    
    filename should show be the input image filename
    
    If `coordinates` is False, then does a maximum image filter and find the source.
    If True, then uses those coordinates instead and the filtered image is not used.
    
    The filtered image is NOT used for the radial profile itself, just for finding the
    source location.
    
    
    '''
    hdu = fits.open(filename, uint=True)[0]
    wcs = WCS(hdu.header)

    
    im = hdu.data
    
    
    im2 = ndi.maximum_filter(im, size=filt_range)
    if coordinates[0] is False:
        coordinates = find_source(filename, filt_range = filt_range)

    rad_profile, radial_err = azimuthalAverage(im, center= [coordinates[0][1], coordinates[0][0]])
    radius = np.arange(len(rad_profile))

    # Get the radius in arcsec
    dr = np.abs(wcs.pixel_scale_matrix[0,0])*u.deg.to(u.arcsec)
    rind = np.arange(len(rad_profile))*dr
    psf_profile = load_psf(rind, filt_range = filt_range)
    psf_scale = (rad_profile[0:5].sum()) / (psf_profile[0:5].sum())

    if show_image:
    
        fig = plt.figure(figsize=[16, 4])
        ax = fig.add_subplot(131, projection=wcs)
        ax2 = fig.add_subplot(132, projection=wcs)
        ax3 = fig.add_subplot(133)


 
        ax.imshow(im, origin='lower', cmap=plt.cm.viridis)
        ax2.imshow(im2, origin='lower', cmap=plt.cm.viridis)


        ax.plot(coordinates[0][1], coordinates[0][0], 'rx')
        ax2.plot(coordinates[0][1], coordinates[0][0], 'rx')

        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        ax2.set_xlabel('RA')
        ax2.set_ylabel('Dec')

        range = [-100, 100]
        ax.set_xbound(range + coordinates[0][1])
        ax.set_ybound(range+coordinates[0][0])
        ax.set_title('Raw Image')
        ax2.set_xbound(range + coordinates[0][1])
        ax2.set_ybound(range+coordinates[0][0])
        ax2.set_title('Filtered Image')


        ax3.set_xlabel('Radius (arcsec)')
        ax3.set_ylabel('Counts / arcsec**2')


        ax3.errorbar(rind, rad_profile, radial_err, fmt='r.')
        ax3.set_yscale('log')
        ax3.semilogy(rind, psf_profile*psf_scale, 'g+')
        ax3.set_xbound([0.9*rind.min(), 600])
        ax3.set_title('Radial Profile')
        plt.show()
        

    return rind, rad_profile, radial_err, psf_profile



def optimize_radius_snr(rind, rad_profile, radial_err, psf_profile,
                    bgd_rin=200, bgd_rout = 300, rlimit = 300,
                    source_ratio = 0.5, show=True):
    '''
    Assumes that you have already constructed the radial profile from the
    source location using make_radial_profile in a given energy band.
    
    This method fits the "source" and "background" components using the profile_model
    method above.
    
    Based on the relative strength of the source and the background this then computes
    the radius that optimizes the signal-to-noise in this band.

    Returns the "best" radius in arcseconds.

    If `show` is True, then generates diagnostic plots.
    

    '''

    psf_scale = (rad_profile[0:5].sum()) / (psf_profile[0:5].sum())
    bgd_scale = rad_profile[(rind > bgd_rin) & (rind<bgd_rout)].mean()

    mask = (rind<rlimit)
    popt, pcov = curve_fit(profile_model, rind[mask], rad_profile[mask], 
                           sigma=radial_err[mask], p0=[psf_scale, bgd_scale, 0])

    src = profile_model(rind, popt[0], 0, 0)
    bgd = profile_model(rind, 0, popt[1], popt[2])
    

    dr = rind[1]-rind[0] # in arcsec
    last_r = 0.
    rtest = 0
    src_cts = 0.
    bgd_cts = 0.
    old_snr = 0.
    all_snr = []

    # Dumb little loop to step out in radius and find the best radius.
    # Keeps all values for plotting regardless of whether you set show=True (to be fixed?)
    # I'm lazy...
    while rtest*dr < rlimit:
        rtest += 1
        area = np.pi*( rtest**2 - last_r**2) # area in pixels

        src_cts += src[rtest] * area
        bgd_cts += bgd[rtest] * area
        snr = src_cts / np.sqrt(src_cts + bgd_cts)
        if snr > old_snr:
            old_snr = snr
            best_radius = rtest
        last_r = rtest
        all_snr.append(snr)
    

    # If yes, make diagnostic plots
    if show:
        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(121)

        ax.errorbar(rind, rad_profile, radial_err, fmt='r.', label = 'Data')
        ax.plot(rind, profile_model(rind, *popt), label = 'Combined Model')
        ax.plot(rind, profile_model(rind, popt[0], 0,0), label = 'Source Only')
        ax.semilogy(rind, bgd, label = 'BGD Only')
        ax.set_ylim([0.5*bgd[mask].min(), 2*rad_profile.max()])
        ax.legend()
        ax.set_xlabel('Radius (arcsec)')
        ax.set_ylabel('Counts / arcsec**2')
        ax.set_xlim([0.9*rind.min(), rlimit])


        ax.set_title('Best Fit')
        ax.axvline(x=best_radius*dr)
        
        
        
        ax2 = fig.add_subplot(122)
        ax2.plot(range(len(all_snr))*dr, all_snr)
        ax2.axvline(x=best_radius*dr)
        ax2.set_xlabel('Radius (arcsec)')
        ax2.set_ylabel('Est Cumulative SNR')
        plt.show()

    return best_radius*dr




