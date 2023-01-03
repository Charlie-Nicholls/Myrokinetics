import os
from numpy import *
from ncdf2dict import ncdf2dict as readnc
from .plotting import Plotters

'''
GYROKINETIC SCAN ANALYSIS
'''

class myro_single(object):

	def __init__(self, filename = None, directory = "./"):
		self.directory = directory
		self.filename = filename
		self.run = {}
		self._open_file()
		self.run['gr'] = imag(self.run['omega'][-1,0,0])
		self.run['mf'] = real(self.run['omega'][-1,0,0])
	
	def __getitem__(self, key):
		return self.run[key]

	def keys(self):
		return self.run.keys()
				
	def _open_file(self, filename = None, directory = None):
		if directory:
			self.directory = directory
		if filename:
			self.filename = filename
		if self.filename is None:
			print("ERROR: filename not given")
			return
		data_in = {}
		if self.filename[-7:] == ".out.nc":
			try:
				data_in = readnc(os.path.join(self.directory,self.filename))
			except:
				print(f"Could not load file {os.path.join(self.directory,self.filename)}")
				return
		else:
			try:
				data_in = readnc(os.path.join(self.directory,self.filename) + ".out.nc")
			except:
				print(f"Could not load file {os.path.join(self.directory,self.filename)}.npz")
				return
		
		self.run = data_in
		
	def plot_omega(self):
		Plotters['Eigen_Single'](data = self.run, var = 0)
	
	def plot_phi(self):
		Plotters['Eigen_Single'](data = self.run, var = 1)
	
	def plot_apar(self):
		Plotters['Eigen_Single'](data = self.run, var = 2)
		
	def plot_phi2(self):
		Plotters['Eigen_Single'](data = self.run, var = 3)
	
	def _plot_diag(self, var = 0):
		Plotters['Eigen_Single'](data = self.run, var = var)
