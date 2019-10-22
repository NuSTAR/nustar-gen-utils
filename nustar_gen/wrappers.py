import os
import warnings

from nustar_gen.utils import energy_to_chan


def make_image(infile, elow = 3, ehigh = 20, clobber=True):
    '''
    Spawn an xselect instance that produces the image in the energy range.
    
    - Checks that input file exists
    - Checks if output file exists. If yes and clobber=False, just returns the output filename.
    - If outfile does not exist or file exists and clobber=True, then converts the input energy range
          to PI channels and then generates an Xselect command file.
    - Spawns an instance of Xselect running in the background.
    - Returns the outfile name
    
    
    '''
    # Make sure environment is set up properly
    _check_environment()

    # Check if input file exists:
    try:
        with open(infile) as f:
            pass
    except IOError:
        raise IOError("make_image: File does not exist %s" % (infile))
    
    
    # Generate outfile name
    outfile = (os.path.splitext(infile))[0]+'_{}to{}keV.fits'.format(int(elow), int(ehigh))
    if (os.path.exists(outfile)) & (~clobber):
        warnings.warn('make_image: %s exists, use clobber=True to regenerate' % (outfile))
    else:
        os.system("rm "+outfile)
    xsel_file = _make_xselect_commands(infile, outfile, elow, ehigh)
    os.system("xselect @"+xsel_file)
    os.system("rm -r -f "+xsel_file)
    
    return outfile


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

