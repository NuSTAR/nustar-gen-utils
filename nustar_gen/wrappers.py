import os, stat
import warnings

from nustar_gen.utils import energy_to_chan
from astropy import units as u



def make_lightcurve(infile, mod, src_reg,
    barycorr=False, time_bin=100*u.s, mode='01',
    bgd_reg='None', outpath='None', elow=3, ehigh=20):
    '''
    Generate a script to run nuproducts
    
    Parameters
    ----------
    
    infile: str
        Full path to the input event file.
    
    mod: str
        'A' or 'B'
        
    src_reg: str
        Full path to source region.
    
    
    Optional Parameters
    -------------------
    

    bgd_reg: str
        If not 'None', then must be the full path to the background region file

    barycorr: bool 
        Default is 'False'. If 'True', then queries the infile for the OBJ J2000
        coordinates and uses these for the barycenter correction.
        
    elow: float
        Low-eneryg bound. Default is 3 keV.
    
    ehigh: float
        High-energy bound. Default is 20 keV.
    
    outpath: str
        Optional. Default is to put the lightcurves in the same location as infile
        
    mode: str
        Optional. Used primarily if you're doing mode06 analysis and need to specify
        output names that are more complicated.

    '''

    from astropy.io.fits import getheader

    # Make sure environment is set up properly
    _check_environment()

    # Check to see that all files exist:
    assert os.path.isfile(infile), 'make_lightcurve: infile does not exist!'
    assert os.path.isfile(src_reg), 'make_lightcurve: src_reg does not exist!'

    if bgd_reg is not 'None':
        assert os.path.isfile(bgd_reg), 'make_lightcurve: bgd_reg does not exist!'
        bkgextract='yes'
    else:
        bkgextract='no'

    reg_base = os.path.basename(src_reg)
    reg_base = os.path.splitext(reg_base)[0]
    
    evdir = os.path.dirname(infile)
    
    seqid = os.path.basename(os.path.dirname(evdir))
    
    if outpath is 'None':
        outdir = evdir
    else:
        outdir = outpath
        
    time_bin = int((time_bin.to(u.s)).value)
    stemout = f'nu{seqid}{mod}{mode}_{reg_base}_{elow}to{ehigh}_{time_bin}s'
    lc_script = outdir+f'/runlc_{stemout}.sh'    
    
    
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    
    with open(lc_script, 'w') as f:
        f.write('nuproducts phafile=NONE bkgphafile=NONE imagefile=NONE ')
        f.write('runmkarf=no runmkrmf=no ')
        f.write(f'indir={evdir} outdir={outdir} instrument=FPM{mod} ')
        f.write(f'steminputs=nu{seqid} stemout={stemout} ')
        f.write(f'srcregionfile={src_reg} ')
        
        if bkgextract is 'no':
            f.write(f'bkgextract=no ')
        else:
            f.write(f'bkgextract=yes bkgregionfile={bgd_reg} ')
        f.write(f'binsize={time_bin} ')
        
        if barycorr:
            attorb=evdir+f'/nu{seqid}{mod}.attorb'
            hdr = getheader(infile)
            ra = hdr['RA_OBJ']
            dec = hdr['DEC_OBJ']
            f.write(f'barycorr=yes srcra_barycorr={ra} srcdec_barycorr={dec} ')
            f.write(f'orbitfile={attorb} ')
            
        f.write('clobber=yes')
        
    os.chmod(lc_script, stat.S_IRWXG+stat.S_IRWXU)
    return lc_script


def make_exposure_map(obs, mod, vign_energy = False,
    det_expo=False):
    '''
    Spawn an instance of nuexpomap to produce an exposure map.
    
    Parameters
    ----------
    
    obs: nustar_gen.info.Observation(), required
        A valid observation metadata.
                
    mod: str
        'A' or 'B'
    
    Other Parameters
    ----------------
    vign_energy: float, optional
        Energy where you want to apply the vignetting. Default is no vignetting.
    
    det_expo : boolean, optional, default=False
        Whether or not to retain the DET1 exposure map file

    '''
    import glob
    
    # Make sure environment is set up properly
    _check_environment()

   # Locate the mast file, attfile, which are what you need for inputs.
    
    auxdir = obs._auxdir    
    evdir = obs._evdir

    # Find the mast file. glob is necessary to handle .gz or .fits extensions:
    mastaspectfile = glob.glob(obs._evdir+'/nu'+obs.seqid+'*mast*')[0]
    # Find the attitude file:
    attfile = glob.glob(obs._auxdir+'/nu'+obs.seqid+'*att*')[0]
    # Find the det1reffile:
    det1reffile = glob.glob(obs._evdir+'/nu'+obs.seqid+mod+'*det1*')[0]
    
    # Only do this for A01, since that's all that matters
    evfile = obs.science_files[mod][0]
    assert '01' in evfile, f'make_exposure_map: Not an 01 event file: {evfile}'
    
    
    
    
    # Construct the nuexpomap call:
    expo_script = obs.out_path+'/nu'+obs.seqid+mod+'_runexpo.sh'
    expo = open(expo_script, 'w')
    
    cmd_string = 'nuexpomap '
    cmd_string += f'infile={evfile} '
    
    if vign_energy is not False:
        cmd_string+=f'vignflag=yes energy={vign_energy} '
    else:
        cmd_string += 'vignflag=no '
        
    cmd_string += f'mastaspectfile={mastaspectfile} '
    cmd_string += f'attfile={attfile} '
    cmd_string += f'det1reffile={det1reffile} '
    
    sky_expo_file = obs.out_path+'/nu'+obs.seqid+mod+'_sky_expo.fits'
    cmd_string += f'expomapfile={sky_expo_file} '

    if det_expo:
        det_expo_file = obs.out_path+'/nu'+obs.seqid+mod+'_det1_expo.fits'
        cmd_string += f'det1instrfile={det_expo_file} '

    cmd_string += 'clobber=yes '
    expo.write(cmd_string)
    expo.close()
    os.chmod(expo_script, stat.S_IRWXG+stat.S_IRWXU)
    
    return expo_script


def make_image(infile, elow = 3, ehigh = 20, clobber=True, outpath=False):
    '''
    Spawn an xselect instance that produces the image in the energy range.

    Parameters
    ----------
    infile: str
        Full path tot eh file that you want to process
    elow: float
        Low-energy band for the image
    ehigh: float
        High-energy band for the image
        
    Other Parameters
    ----------------
    
    clobber: boolean, optional, default=True
        Overwrite existing files?

    outpath: str, optional, default=os.path.dirname(infile)
        Set the destination for output. Defaults to same location as infile.
   
    Return
    -------
    outfile: str
        The full path to the output image.
   
    '''
    # Make sure environment is set up properly
    _check_environment()

    # Check if input file exists:
    try:
        with open(infile) as f:
            pass
    except IOError:
        raise IOError("make_image: File does not exist %s" % (infile))

    if not outpath:
        outdir=os.path.dirname(infile)
    else:
        outdir=outpath
    
    # Trime the filename:
    sname=os.path.basename(infile)
    if sname.endswith('.gz'):
        sname = os.path.splitext(sname)[0]
    sname = os.path.splitext(sname)[0]
    
    # Generate outfile name
    outfile = outdir + '/'+sname+f'_{elow}to{ehigh}keV.fits'
       
    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('make_image: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands(infile, outfile, elow, ehigh)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    return outfile


def extract_det1_events(infile, regfile, clobber=True, outpath=False):
    '''
    Spawn an xselect instance that produces a new event file screened using a det1 region
    file.
    
    Parameters
    ----------
    infile: str
        Full path tot eh file that you want to process
    regfile: str
        Full path to a ds9 region file (in physical coordinates) to be used to filter
        the events.
                
    Other Parameters
    ----------------
    
    clobber: boolean, optional, default=True
        Overwrite existing files?

    outpath: str, optional, default=os.path.dirname(infile)
        Set the destination for output. Defaults to same location as infile.
   
    Return
    -------
    outfile: str
        The full path to the output image.
   

    
    
    '''

    # Make sure environment is set up properly
    _check_environment()

    # Check if input file exists:
    try:
        with open(infile) as f:
            pass
    except IOError:
        raise IOError("extract_det1_events: File does not exist %s" % (infile))

    try:
        with open(regfile) as f:
            pass
    except IOError:
        raise IOError("extract_det1_events: File does not exist %s" % (regfile))
        
        
    if not outpath:
        outdir=os.path.dirname(infile)
    else:
        outdir=outpath
    
    # Trim the filename:
    sname=os.path.basename(infile)
    if sname.endswith('.gz'):
        sname = os.path.splitext(sname)[0]
    sname = os.path.splitext(sname)[0]

    rshort = os.path.basename(regfile)
    rname = os.path.splitext(rshort)[0]
    
    # Generate outfile name
    outfile = outdir + '/'+sname+f'_{rname}.evt'

    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('extract_det1_events: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands_det1_evts(infile, outfile, regfile)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    
    return outfile


def _make_xselect_commands_det1_evts(infile, outfile, regfile):
    '''
    Helper script to generate the xselect commands to make an image in a given NuSTAR range
    '''
    
    import glob
    for oldfile in glob.glob("session1*"):
        os.system(f"rm {oldfile}")
    
    xsel=open("xsel.xco","w")
    xsel.write("session1\n")
    xsel.write("read events \n")
    evdir=os.path.dirname(infile)
    xsel.write(f'{evdir} \n ' )
    evfile = os.path.basename(infile)
    xsel.write(f'{evfile} \n ')
    xsel.write('yes \n')
    xsel.write('set xyname\n')
    xsel.write('DET1X\n')
    xsel.write('DET1Y\n')
    xsel.write(f'filter region {regfile} \n')
    xsel.write("extract events\n")
    xsel.write("save events\n")
    xsel.write("%s \n" % outfile)
    xsel.write('n \n')
    xsel.write('exit\n')           
    xsel.write('n \n')
    xsel.close()
    return 'xsel.xco'




def make_det1_image(infile, elow = 3, ehigh = 20, clobber=True, outpath=False):
    '''
    Spawn an xselect instance that produces a DET1 image in the energy range.
    
    Parameters
    ----------
    infile: str
        Full path tot eh file that you want to process
    elow: float
        Low-energy band for the image
    ehigh: float
        High-energy band for the image
        
    Other Parameters
    ----------------
    
    clobber: boolean, optional, default=True
        Overwrite existing files?

    outpath: str, optional, default=os.path.dirname(infile)
        Set the destination for output. Defaults to same location as infile.
   
    Return
    -------
    outfile: str
        The full path to the output image.
    
    '''
    # Make sure environment is set up properly
    _check_environment()

    # Check if input file exists:
    try:
        with open(infile) as f:
            pass
    except IOError:
        raise IOError("make_image: File does not exist %s" % (infile))

    if not outpath:
        outdir=os.path.dirname(infile)
    else:
        outdir=outpath
    
    # Trime the filename:
    sname=os.path.basename(infile)
    if sname.endswith('.gz'):
        sname = os.path.splitext(sname)[0]
    sname = os.path.splitext(sname)[0]
    
    # Generate outfile name
    outfile = outdir + '/'+sname+f'_{elow}to{ehigh}keV_det1.fits'
       
    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('make_image: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands_det1(infile, outfile, elow, ehigh)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    return outfile


def _make_xselect_commands_det1(infile, outfile, elow, ehigh):
    '''
    Helper script to generate the xselect commands to make an image in a given NuSTAR range
    '''
    
    import glob
    for oldfile in glob.glob("session1*"):
        os.system(f"rm {oldfile}")
    
    xsel=open("xsel.xco","w")
    xsel.write("session1\n")
    xsel.write("read events \n")
    evdir=os.path.dirname(infile)
    xsel.write(f'{evdir} \n ' )
    evfile = os.path.basename(infile)
    xsel.write(f'{evfile} \n ')
    xsel.write('yes \n')
    xsel.write('set xyname\n')
    xsel.write('DET1X\n')
    xsel.write('DET1Y\n')
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    xsel.write('filter pha_cutoff {} {} \n'.format(pi_low, pi_high))
    xsel.write('set xybinsize 1\n')
    xsel.write("extract image\n")
    xsel.write("save image\n")
    xsel.write("%s \n" % outfile)
    xsel.write('exit\n')           
    xsel.write('n \n')
    xsel.close()
    return 'xsel.xco'

def _make_xselect_commands(infile, outfile, elow, ehigh):
    '''
    Helper script to generate the xselect commands to make an image in a given NuSTAR range
    '''
    
    xsel=open("xsel.xco","w")
    xsel.write("session1\n")
    xsel.write("read events \n")
    xsel.write('./ \n ' )
    xsel.write('%s \n ' % infile)
    xsel.write('yes \n')
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    xsel.write('filter pha_cutoff {} {} \n'.format(pi_low, pi_high))
    xsel.write('set xybinsize 1\n')
    xsel.write("extract image\n")
    xsel.write("save image\n")
    xsel.write("%s \n" % outfile)
    xsel.write('exit\n')           
    xsel.write('n \n')
    xsel.close()
    return 'xsel.xco'


def _check_environment():
    try:
        if ("CALDB" in os.environ) & ("HEADAS" in os.environ):
            pass
    except IOerror:
        raise IOError("Environment variables $CALDB and $HEADAS not set")

