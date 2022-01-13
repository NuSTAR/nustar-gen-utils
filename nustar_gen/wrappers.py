import os, stat
import warnings

from nustar_gen.utils import energy_to_chan, validate_det1_region
from astropy import units as u

def make_spectra(infile, mod, src_reg, usrgti=False,
    mode='01', bgd_reg='None', outpath='None', runmkarf='yes', extended='no'):
    '''
    Generate a script to run nuproducts to extract a source (and optionally
    a background) spectrum along with their response files.
    
    Always runs numkrmf and numkarf for now.
    
    Parameters
    ----------
    
    infile: str
        Full path to the input event file.
    
    mod: str
        'A' or 'B'
        
    src_reg: str
        Full path to source region.
    
    
    Other Parameters
    -------------------
    

    bgd_reg: str
        If not 'None', then must be the full path to the background region file

    barycorr: bool 
    
    outpath: str
        Optional. Default is to put the lightcurves in the same location as infile
        
    usrgti: str
        Optional. Full path to the GTI that you want to apply to the spectrum.
        
    mode: str
        Optional. Used primarily if you're doing mode06 analysis and need to specify
        output names that are more complicated.

    '''

    from astropy.io.fits import getheader

    # Make sure environment is set up properly
    _check_environment()

    # Check to see that all files exist:
    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)


    
    assert os.path.isfile(infile), 'make_spectra: infile does not exist!'
    assert os.path.isfile(src_reg), 'make_spectra: src_reg does not exist!'

    if bgd_reg != 'None':
        assert os.path.isfile(bgd_reg), 'make_spectra: bgd_reg does not exist!'
        bkgextract='yes'
    else:
        bkgextract='no'

    reg_base = os.path.basename(src_reg)
    reg_base = os.path.splitext(reg_base)[0]
    
    evdir = os.path.dirname(infile)
    
    seqid = os.path.basename(os.path.dirname(evdir))
    
    if outpath == 'None':
        outdir = evdir
    else:
        outdir = outpath
        try:
            os.makedirs(outpath)
        except FileExistsError:
        
    # directory already exists
            pass    
        
    stemout = f'nu{seqid}{mod}{mode}_{reg_base}'
    lc_script = os.path.join(outdir, f'runspec_{stemout}.sh')
    
        # Parse GTI file if necessary    
    if usrgti is not False:
        rshort = os.path.basename(usrgti)
        rname = os.path.splitext(rshort)[0]
        stemout += f'_{rname}'
   
    with open(lc_script, 'w') as f:
        f.write('nuproducts imagefile=NONE lcfile=NONE bkglcfile=NONE ')
        f.write(f'runmkarf={runmkarf} extended={extended} runmkrmf=yes ')
        f.write(f'indir={evdir} outdir={outdir} instrument=FPM{mod} ')
        f.write(f'steminputs=nu{seqid} stemout={stemout} ')
        f.write(f'srcregionfile={src_reg} ')
        
        if bkgextract == 'no':
            f.write(f'bkgextract=no ')
        else:
            f.write(f'bkgextract=yes bkgregionfile={bgd_reg} ')
        
        if usrgti is not False:
            f.write(f'usrgtifile={usrgti} ')
        f.write('clobber=yes')
        
    os.chmod(lc_script, stat.S_IRWXG+stat.S_IRWXU)
    return lc_script


def make_lightcurve(infile, mod, src_reg,
    barycorr=False, time_bin=100*u.s, mode='01',usrgti=False,
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
    
    
    Other Parameters
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
    
    time_bin: astropy unit
        Length of the time bin. Default is 100*u.s. Min is 1*u.s
    
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
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    
    assert os.path.isfile(infile), 'make_lightcurve: infile does not exist!'
    assert os.path.isfile(src_reg), 'make_lightcurve: src_reg does not exist!'

    assert time_bin.to(u.s) >= 1*u.s, f'make_lightcurve: {time_bin} too fast, use another tool.'

    if bgd_reg != 'None':
        assert os.path.isfile(bgd_reg), 'make_lightcurve: bgd_reg does not exist!'
        bkgextract='yes'
    else:
        bkgextract='no'

    reg_base = os.path.basename(src_reg)
    reg_base = os.path.splitext(reg_base)[0]
    
    evdir = os.path.dirname(infile)
    
    seqid = os.path.basename(os.path.dirname(evdir))
    
    if outpath == 'None':
        outdir = evdir
    else:
        outdir = outpath
        try:
            os.makedirs(outdir)
        except FileExistsError:
    # directory already exists
            pass    
    time_bin = (time_bin.to(u.s)).value
    
    
    
    stemout = f'nu{seqid}{mod}{mode}_{reg_base}_{elow}to{ehigh}_{time_bin:3.4}s'
    
    
            # Parse GTI file if necessary    
    if usrgti is not False:
        rshort = os.path.basename(usrgti)
        rname = os.path.splitext(rshort)[0]
        stemout += f'_{rname}'

    lc_script = os.path.join(outdir, f'runlc_{stemout}.sh')  
    
    
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    
    with open(lc_script, 'w') as f:
        f.write('nuproducts phafile=NONE bkgphafile=NONE imagefile=NONE ')
        f.write(f'runmkarf=no runmkrmf=no pilow={pi_low} pihigh={pi_high} ')
        f.write(f'indir={evdir} outdir={outdir} instrument=FPM{mod} ')
        f.write(f'steminputs=nu{seqid} stemout={stemout} ')
        f.write(f'srcregionfile={src_reg} ')
        
        if bkgextract == 'no':
            f.write(f'bkgextract=no ')
        else:
            f.write(f'bkgextract=yes bkgregionfile={bgd_reg} ')
        f.write(f'binsize={time_bin} ')
        
        
        if usrgti is not False:
            f.write(f'usrgtifile={usrgti} ')
        
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
    det_expo=False, evf=False):
    '''
    Create a script to run nuexpoomap. Returns the script name.
    
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
    
    evdir = obs.evdir

    # Find the mast file. glob is necessary to handle .gz or .fits extensions:
    mastaspectfile = glob.glob(evdir+'/nu'+obs.seqid+'*mast*')[0]
    # Find the attitude file:
    attfile = glob.glob(evdir+'/nu'+obs.seqid+'*att.fits')[0]
    # Find the det1reffile:
    det1reffile = glob.glob(evdir+'/nu'+obs.seqid+mod+'_det1.fits')[0]
    
    # Only do this for A01, since that's all that matters
    # Override this with evfile keyword:
    if evf is False:
        evfile = obs.science_files[mod][0]
        assert '01' in evfile, f'make_exposure_map: Not an 01 event file: {evfile}'
    else:
        evfile=evf
    
    
    
    # Construct the nuexpomap call:
    expo_script = os.path.join(obs.out_path, 'runexpo_'+obs.seqid+mod+'.sh')
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


def make_image(infile, elow = 3, ehigh = 20, clobber=True, outpath=False, usrgti=False):
    '''
    Spawn an xselect instance that produces the image in the energy range.

    Parameters
    ----------
    infile: str
        Full path to the file that you want to process
    elow: float
        Low-energy band for the image. Default is 3 keV.
    ehigh: float
        High-energy band for the image. Default is 20 keV.
        
    Other Parameters
    ----------------
    
    clobber: boolean, optional, default=True
        Overwrite existing files?

    outpath: str, optional, default=os.path.dirname(infile)
        Set the destination for output. Defaults to same location as infile.
   
    usrgti : str, optional, default = False
        Use a GTI file to time-fitler the data (see nustar_gen.utils.make_usr_gti)
        If False, do nothing.
   
    Return
    -------
    outfile: str
        The full path to the output image.
   
    '''
    # Make sure environment is set up properly
    _check_environment()

    # Check if input file exists:
    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    
    try:
        with open(infile) as f:
            pass
    except IOError:
        raise IOError("make_image: File does not exist %s" % (infile))

    if not outpath:
        outdir=os.path.dirname(infile)
    else:
        outdir=outpath
        try:
            os.makedirs(outdir)
        except FileExistsError:
    # directory already exists
            pass
    
    # Trim the filename:
    sname=os.path.basename(infile)
    if sname.endswith('.gz'):
        sname = os.path.splitext(sname)[0]
    sname = os.path.splitext(sname)[0]
    

    if usrgti is not False:
        rshort = os.path.basename(usrgti)
        rname = os.path.splitext(rshort)[0]
        sname += f'_{rname}'

    
    # Generate outfile name
    outfile = os.path.join(outdir, sname+f'_{elow}to{ehigh}keV.fits')
    
    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('make_image: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands(infile, outfile, elow, ehigh, usrgti=usrgti)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    return outfile


def extract_sky_events(infile, regfile, clobber=True, outpath=False):
    '''
    Spawn an xselect instance that produces a new event file screened using a sky ds9
    region file.
    
    Parameters
    ----------
    infile: str
        Full path to the event file that you want to process
    regfile: str
        Full path to a ds9 region file (in sky coordinates) to be used to filter
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
    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
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
    outfile = os.path.join(outdir, sname+f'_{rname}.evt')

    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('extract_sky_events: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands_sky_evts(infile, outfile, regfile)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    
    return outfile

def barycenter_events(obs, infile, mod='A', barycorr_pars=None):
    '''
    Run barycorr on an event file. 

    Parameters
    --------------------
    obs: nustar_gen.info.Observation
        An instance of the Observation class
    
    infile: str
        Full path to input file

    mod: str
        Module to use. 'A' or 'B'

    Other Parameters
    -------------------

    barycorr_pars: dict
        additional parameters to be passed to barycorr,
        in the format {'refframe': 'icrs', 'ra': 123.456}

    '''

    # Locate the attorb file:
    evdir = obs.evdir
    attorb = os.path.join(evdir, f'nu{obs.seqid}{mod}.attorb')

    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    
    
    # Trim the filename:
    if obs.out_path is False:
        outdir = os.path.dirname(infile)
        print(outdir)
    else:
        outdir = obs.out_path

    sname=os.path.basename(infile)
    sname=os.path.splitext(sname)[0]
    if barycorr_pars is None:
        barycorr_pars = {}
    if 'outfile' not in barycorr_pars:
        # Generate outfile name
        barycorr_pars['outfile'] = os.path.join(outdir, sname+f'_barycorr.fits')
    if 'ra' not in barycorr_pars:
        barycorr_pars['ra'] = obs.source_position.ra.deg
    if 'dec' not in barycorr_pars:
        barycorr_pars['dec'] = obs.source_position.dec.deg
        
    bary_sh = os.path.join(outdir, 'run_bary_'+sname+'.sh')

    with open(bary_sh, 'w') as f:
        f.write(f'barycorr infile={infile} clobber=yes ')
        f.write(f'orbitfiles={attorb} ')
        for par in barycorr_pars:
            f.write(f'{par}={barycorr_pars[par]} ')
    
    os.environ['HEADASNOQUERY'] = ""
    os.environ['HEADASPROMPT'] = "/dev/null"
    
    os.chmod(bary_sh, stat.S_IRWXG+stat.S_IRWXU)
    os.system(f'{bary_sh}')
 
    return barycorr_pars['outfile']

def apply_gti(infile, gtifile, clobber=True, outpath=False):
    '''
    Spawn an xselect instance that produces a new event file screened using GTI file
    
    Parameters
    ----------
    infile: str
        Full path to the event file that you want to process
    regfile: str
        Full path to a ds9 region file (in sky coordinates) to be used to filter
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
    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    
    try:
        with open(infile) as f:
            pass
    except IOError:
        raise IOError("apply_gti: File does not exist %s" % (infile))

    try:
        with open(gtifile) as f:
            pass
    except IOError:
        raise IOError("apply_gti: File does not exist %s" % (gtifile))
        
        
    if not outpath:
        outdir=os.path.dirname(infile)
    else:
        outdir=outpath
    
    # Trim the filename:
    sname=os.path.basename(infile)
    if sname.endswith('.gz'):
        sname = os.path.splitext(sname)[0]
    sname = os.path.splitext(sname)[0]

    rshort = os.path.basename(gtifile)
    rname = os.path.splitext(rshort)[0]
    
    # Generate outfile name
    outfile = os.path.join(outdir, sname+f'_{rname}.evt')

    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('apply_gti: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands_apply_gti(infile, outfile, gtifile)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    
    return outfile


def _make_xselect_commands_apply_gti(infile, outfile, gtifile):
    '''
    Helper script to generate the xselect commands to extract events from
    a given sky region.
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
    xsel.write(f'filter time \n')
    xsel.write('file \n')
    xsel.write(f'{gtifile}\n')
    xsel.write('extract events\n')
    xsel.write("save events\n")
    xsel.write("%s \n" % outfile)
    xsel.write('n \n')
    xsel.write('exit\n')           
    xsel.write('n \n')
    xsel.close()
    return 'xsel.xco'


def _make_xselect_commands_sky_evts(infile, outfile, regfile):
    '''
    Helper script to generate the xselect commands to extract events from
    a given sky region.
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
    xsel.write(f'filter region {regfile} \n')
    xsel.write("extract events\n")
    xsel.write("save events\n")
    xsel.write("%s \n" % outfile)
    xsel.write('n \n')
    xsel.write('exit\n')           
    xsel.write('n \n')
    xsel.close()
    return 'xsel.xco'


def _make_xselect_commands(infile, outfile, elow, ehigh, usrgti=False):
    '''
    Helper script to generate the xselect commands to make an image in a given NuSTAR range
    '''
    
    xsel=open("xsel.xco","w")
    xsel.write("session1\n")
    xsel.write("read events \n")
    evdir=os.path.dirname(infile)
    xsel.write(f'{evdir} \n ' )
    evfile = os.path.basename(infile)
    xsel.write(f'{evfile} \n ')
    xsel.write('yes \n')
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    if usrgti is not False:
        xsel.write(f'filter time \n')
        xsel.write('file \n')
        xsel.write(f'{usrgti}\n')
        xsel.write('extract events\n')
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


###
###
###
###
###
###

# From here down are DET1 methods

###
###
###
###
###
###


def make_det1_image(infile, elow = 3, ehigh = 20, clobber=True, outpath=False):
    '''
    Spawn an xselect instance that produces a DET1 image in the energy range.
    
    Parameters
    ----------
    infile: str
        Full path tot eh file that you want to process
    elow: float
        Low-energy band for the image. Default is 3 keV.
    ehigh: float
        High-energy band for the image. Default is 20 keV.
        
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
    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    
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
    outfile = os.path.join(outdir, sname+f'_{elow}to{ehigh}keV_det1.fits')
       
    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('make_image: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands_det1(infile, outfile, elow, ehigh)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    return outfile

def extract_det1_events(infile, regfile, elow=0., ehigh=165., clobber=True, outpath=False):
    '''
    Spawn an xselect instance that produces a new event file screened using a det1 region
    file.
    
    Parameters
    ----------
    infile: str
        Full path to the event file that you want to process

    regfile: str
        Full path to a ds9 region file (in physical coordinates) to be used to filter
        the events.

    elow: float
        Lowest energy to use (default is 0 keV)
    
    ehigh: float
        Highest energy to use (default is 165 keV)
        
                
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
    # Make sure region file is correctly formatted
    validate_det1_region(regfile)

    # Check if input file exists:
    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    
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
    outfile = os.path.join(outdir, sname+f'_{rname}.evt')

    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('extract_det1_events: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
        
    
        
#    xsel_file = _make_xselect_commands_det1_evts(infile, outfile, regfile)
    xsel_file = _make_xselect_commands_det1_evts_pifilt(infile, outfile,
                                                        regfile, elow, ehigh)

    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    
    return outfile

def make_det1_lightcurve(infile, mod, obs,
    time_bin=100*u.s, mode='01',
    elow=3, ehigh=20, stemout=False, gtifile=False):
    '''
    Generate a script to run nuproducts to make a lightcurve using the whole
    FoV and turning off all vignetting and PSF effects. Assumes that infile 
    has already been filtered using extract_det1_events().
    
    Parameters
    ----------
    
    infile: str
        Full path to the input event file. This should be pre-filtered by
        by extract_det1_events
    
    mod: str
        'A' or 'B'
        
    obs: nustar_gen.info.Observation
        Observation meta data
    
    
    Other Parameters
    -------------------
        
    elow: float, optional, default = 3 keV
        Low-energy bound
    
    ehigh: float, optional, default is 20 keV
        High-energy bound
        
    time_bin: astropy unit
        Length of the time bin. Default is 100*u.s
    
    mode: str, optional, default is '01'
        Optional. Used to specify stemout if you're doing mode06 analysis and want
        to specify output names that are more complicated.
        
    gtifile: str
        Path to a GTI file. If this is set, then this is passed to nuproducts.
    
    stemout: str, optional
        Use the specified stemout string when calling nuproducts. Otherwise
        uses the default value.

    '''

    from astropy.io.fits import getheader

    # Make sure environment is set up properly
    _check_environment()

    # Check to see that all files exist:
    
    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    
    assert os.path.isfile(infile), 'make_det1_lightcurve: infile does not exist!'

    assert time_bin.to(u.s) >= 1*u.s, f'make_det1_lightcurve: time_bin {time_bin} too fast, use another tool.'

    evdir = obs.evdir
    seqid = obs.seqid
    outdir = obs.out_path
    
    hdr = getheader(infile)
    ra = hdr['RA_OBJ']
    dec = hdr['DEC_OBJ']

    time_bin = int((time_bin.to(u.s)).value)
    if stemout is False:
        stemout = f'nu{seqid}{mod}{mode}_full_FoV_{elow}to{ehigh}_{time_bin}s'
    
    lc_script = os.path.join(outdir ,f'rundet1lc_{stemout}.sh')
    
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    
    with open(lc_script, 'w') as f:
        f.write('nuproducts phafile=NONE bkgphafile=NONE imagefile=NONE ')
        f.write(f'infile={infile} ')
        f.write('runmkarf=no runmkrmf=no ')
        f.write(f'indir={evdir} outdir={outdir} instrument=FPM{mod} ')
        f.write(f'steminputs=nu{seqid} stemout={stemout} ')
        f.write(f'pilow={pi_low} pihigh={pi_high} ')
        f.write(f'bkgextract=no ')
        f.write(f'binsize={time_bin} ')
        f.write(f'srcra={ra} srcdec={dec} srcregionfile=DEFAULT srcradius=299 ')
        
        # Turn off all of the time-dependent corrections for the pointing here
        f.write(f'lcpsfflag=no lcexpoflag=no lcvignflag=no ')
        
        if (gtifile != False):
            f.write(f'usrgtifile={gtifile} ')    
        f.write('clobber=yes')
        
    os.chmod(lc_script, stat.S_IRWXG+stat.S_IRWXU)

    return lc_script

def make_det1_spectra(infile, mod, obs,
    stemout=False, gtifile=False):
    '''
    Generate a script to run nuproducts to extract a source 
    spectrum along with the associated RMF.
    
    Assumes that infile has already been filtered using extract_det1_events().
 
    Always runs numkrmf, never runs numkarf. Never extract background.
    
    Parameters
    ----------
    
    infile: str
        Full path to the input event file. This should be pre-filtered by
        by extract_det1_events
    
    mod: str
        'A' or 'B'
        
    obs: nustar_gen.info.Observation
        Observation meta data
   
    Other Parameters
    -------------------
            
    
    stemout: str
        Optional. Use the specified stemout string when calling nuproducts. Otherwise
        uses the default value.
    
    gtifile: str
        Path to a GTI file. If this is set, then this is passed to nuproducts.
        **NOTE** As of now proper treatment of this being barycenter corrected (or not)
        is not supported. If you're doing pulsar analysis, please write your own version.
        
    
    '''

#    from astropy.io.fits import getheader
    from os.path import basename
    
    # Make sure environment is set up properly
    _check_environment()

    # Account for the case where you only get the filename (i.e. if you're
    # in the current working directory)
    infile = os.path.abspath(infile)
    

    # Check to see that all files exist:
    assert os.path.isfile(infile), f'make_det1_spectra: {infile} does not exist!'
#    assert os.path.isfile(src_reg), 'make_det1_spectra: src_reg does not exist!'


    bkgextract='no'
   
    evdir = obs.evdir
    seqid = obs.seqid
    outdir = obs.out_path
    # Construct the output file name:
    
    
    ra =obs.source_position.ra.deg
    dec = obs.source_position.dec.deg
    
    # Use the default stemout unless this is set
    if stemout is False:
        stemout = basename(infile).split('.')[0]
        
    lc_script = os.path.join(outdir, f'rundet1spec_{stemout}.sh')
    
   
    with open(lc_script, 'w') as f:
        f.write('nuproducts imagefile=NONE lcfile=NONE bkglcfile=NONE ')
        f.write(f'infile={infile} ')
        f.write('runmkarf=no runmkrmf=yes ')
        f.write(f'indir={evdir} outdir={outdir} instrument=FPM{mod} ')
        f.write(f'steminputs=nu{seqid} stemout={stemout} ')
        f.write(f'srcra={ra} srcdec={dec} srcregionfile=DEFAULT srcradius=299 ')
        if (gtifile != False):
            f.write(f'usrgtifile={gtifile} ')
        f.write(f'runbackscale=no ')
        f.write(f'bkgextract=no ')            
        f.write('clobber=yes')
        
    os.chmod(lc_script, stat.S_IRWXG+stat.S_IRWXU)
    return lc_script

def _make_xselect_commands_det1_evts(infile, outfile, regfile):
    '''
    Helper script to generate the xselect commands to extract events from
    a given region.
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

def _make_xselect_commands_det1_evts_pifilt(infile, outfile, regfile, elow, ehigh):
    '''
    Helper script to generate the xselect commands to extract events from
    a given region.
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
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    xsel.write('filter pha_cutoff {} {} \n'.format(pi_low, pi_high))
    xsel.write("extract events\n")
    xsel.write("save events\n")
    xsel.write("%s \n" % outfile)
    xsel.write('n \n')
    xsel.write('exit\n')           
    xsel.write('n \n')
    xsel.close()
    return 'xsel.xco'



def _make_xselect_commands_det1(infile, outfile, elow, ehigh):
    '''
    Helper script to generate the xselect commands to make an image in a
    given NuSTAR energy range
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
