from astropy.time import Time, TimeDelta
import astropy.units as u

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
	