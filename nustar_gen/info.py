from astropy.time import Time, TimeDelta
import astropy.units as u
import os
import glob


class NuSTAR():
    '''
    Class for holding constant attributes
    
    Parameters
    -----------
    
    mjdref: Astropy Time
        Contains the integer part of the MJDREF date (Jan 1, 2010) for
        the NuSTAR MET epoch.
    
    pixel: float
        Size of a sky pixel (2.54 arcsec)
    
    raw_pixel: float
        Contains the ideal pixel pitch in NuSTAR (604.8 microns)
    
    pixel_um: float
        Size of a DET1 pixel (raw_pixel / 5)
        
    '''
    
    def __init__(self):
        self.mjdref = Time(55197., format = 'mjd')
        self.raw_pixel = 604.8 * u.micron
        self.pixel_um = self.raw_pixel / 5.
        self.pixel = 2.54 * u.arcsec

    @classmethod
    def time_to_met(self, time):
        ''' 
        Convert a datetime object to a unitless NuSTAR MET second.
        '''
        dt = (time - self.mjdref).to(u.s).value
        return dt
    
    @classmethod
    def met_to_time(self, met):
        '''
        Assumes unitless MET seconds input. Need to catch this.
        '''
        this_time = TimeDelta(met*u.s) + self.mjdref
        return this_time




class Observation():
    '''
    Class for storing meta-data about a given NuSTAR observation.
    
    Parameters
    ----------
    path: str, optional, default './'
        The top-level working directory. All paths are assumed to be relative to this
        location.
    
    seqid: str
        The sequence id that you're going to be working with. Data are assumed to be
        obs.path+os.seqid+`/event_cl/`. Set via set_seqid.

    exposure: dict
        'A' and 'B' keys give a list of the exposure times for each event file.
                
    source_position: Astropy SkyCoord
        J2000 coordinates of the source based on
        the FITS header.

    
    '''
        
    def __init__(self, path='./', seqid=False):
        self.path=path
        self.modules = ['A', 'B']
        self._source_position=None
        self._seqid=False

        if seqid is not False:
            self.set_seqid = seqid
            
        self._out_path=False

    @property
    def seqid(self):
        '''
        Returns the current sequence ID. Data are assumed to be in obs.path+obs.method
        '''
        return self._seqid
        
    @seqid.setter
    def set_seqid(self, value):
        '''Set the sequence ID. Raise error if path+sequence doesn't exist. Finds
        clean event files and parses the input fits header to populate the other'''
        if os.path.isdir(self.path+value):
            self._seqid = value
        else:
            raise ValueError(f"Path does not exist! {self.path+value})")
        self._datapath=self.path+self._seqid
        
        # Set subdirectories
        self._hkdir=self._datapath+'/hk/'
        self._evdir=self._datapath+'/event_cl/'
        self._auxdir=self._datapath+'/auxil/'
        
        self._find_cleaned_files()
        
        self._parse_header()

    @property
    def exposure(self):
        '''
        Returns an dict (with 'A' and 'B' as keys) with lists of exposures for
        all event files.
        '''
        return self._exposure
        
    @property
    def source_position(self):
        '''
        Returns the current source RA/Dec from the FITS headers as Astropy
        SkyCoord object
        '''
        return self._source_position

    @property
    def science_files(self):
        '''
        Returns a list of science (01) event files
        '''
        
        if self._seqid is False:
            raise ValueError(f"Set sequence ID and path first!")
        
        science_files={}
        for mod in self._evtfiles:
            science_files[mod] = []
            for file in self._evtfiles[mod]:
                if ('01' in file) &  \
                    ( (file.endswith('gz') ) | \
                    (file.endswith('evt') ) ):
                    science_files[mod].extend([file])   
        return science_files

    @property
    def out_path(self):
        '''
        Returns the output path.
        '''
        return self._out_path

    @out_path.setter
    def set_outpath(self, value):
        '''
        Set the output path.
        '''
        
        if os.path.isdir(self.path+value):
            self._out_path = self.path+value
        else:
            raise ValueError(f"Output path does not exist! {self.path+value})")
        return self._out_path

    @classmethod
    def _find_cleaned_files(self):
        '''
        Uses self._evdir to find all of the cleaned event files.
        '''
        self._evtfiles = {}
        for mod in self.modules:
            self._evtfiles[mod] = sorted(glob.glob(self._evdir+f'nu*{mod}*cl.evt*'))
        return
        
    @classmethod
    def _parse_header(self):
        from astropy.io.fits import getheader
        from astropy.coordinates import SkyCoord
        
        self._exposure = {}
        
        for mod in self.modules:
            for evtfile in self._evtfiles[mod]:
                hdr = getheader(evtfile)
                if self._source_position is None:               
                    self._source_position = \
                            SkyCoord(hdr['RA_OBJ'], hdr['DEC_OBJ'], unit='deg')
                
                for ti in range(5):
                    keystr = f'{mod}'+f'{ti+1}'.zfill(2)

                    if keystr in evtfile:
                        self._exposure[keystr] = hdr['EXPOSURE']


        return
        

    @classmethod
    def exposure_report(self):
        '''
        Make pretty output of the exposure
        '''
        
        if self._seqid is False:
            raise ValueError(f"Set sequence ID and path first!")
        
        for mod in self.modules:
            for ti in range(5):
                keystr = f'{mod}'+f'{ti+1}'.zfill(2)
                if keystr in self._exposure:
                    print(f'Exposure for FPM{mod}, mode '+f'{ti+1}'.zfill(2)+f' is: {1e-3*self._exposure[keystr]:10.4} ks')

            print()

        return

    @classmethod
    def download_bgd_report(self):
        '''
        Wrappers to download the background report from the SOC:
        '''

        base_html = 'http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/'
        base_html += 'SAA_Filtering/nulyses_reports/'
        for mod in self.modules:
            pdf_html = base_html + f'{self.seqid}/nu{self.seqid}_SAA_Report_{mod}.pdf'
#            print(pdf_html)
            os.system(f'wget {pdf_html}')
            os.system(f'mv {os.path.basename(pdf_html)} {self.out_path}')
        return
    
    
    
    
    
    
    
    
    
    