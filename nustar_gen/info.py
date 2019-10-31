from astropy.time import Time, TimeDelta
import astropy.units as u
import os
import glob


class NuSTAR():
	'''
	Class for holding constant attributes
	
	'''
	
	def __init__(self):
		self.mjdref = Time(55197., format = 'mjd')


	def time_to_met(self, time):
		''' 
		Convert a datetime object to a unitless NuSTAR MET second.
		'''
		dt = (time - self.mjdref).to(u.s).value
		return dt
		
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
    path : string
       The working directory. Defaults to the current directory.

	Methods
	----------
	exposure_report:
		Finds all of the cleaned event files, parses the 

	
	Attributes
	-----------
	seqid: 'str'
	
	source_position: Astropy SkyCoord object
	
	
	'''
		
	def __init__(self, path='./'):
		self.path=path
		self.modules = ['A', 'B']
		self._source_position=None
		self._seqid=False
		self._out_path=False

	@property
	def seqid(self):
		'''
		Returns the current sequence ID
		'''
		return self._seqid
	
	@property
	def exposure(self):
		'''
		Returns the current sequence ID
		'''
		return self._exposure
		
	@property
	def source_position(self):
		'''
		Returns the current sequence ID
		'''
		return self._source_position

	
	@property
	def science_files(self):
		'''
		Returns the current sequence ID
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
		Returns the current sequence ID
		'''
		return self._out_path

	@out_path.setter
	def set_outpath(self, value):
		'''
		Returns the current sequence ID
		'''
		
		if os.path.isdir(self.path+value):
			self._out_path = self.path+value
		else:
			raise ValueError(f"Output path does not exist! {self.path+value})")
		return self._out_path


	@seqid.setter
	def set_seqid(self, value):
		'''Set the sequence ID. Raise error if path+sequence doesn't exist'''
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




	
	def _find_cleaned_files(self):
		'''
		Uses self._evdir to find all of the cleaned event files.
		'''
		self._evtfiles = {}
		for mod in self.modules:
			self._evtfiles[mod] = sorted(glob.glob(self._evdir+f'nu*{mod}*cl.evt*'))
	

	
	
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
					print(f'Exposure for FPM{mod} is: {1e-3*self._exposure[keystr]:10.4} ks')

			print()
		