import os
from numpy import imag, real
from .ncdf2dict import ncdf2dict as readnc
from .plotting import Plotters

'''
SINGLE RUN ANALYSIS
'''

class myro_single(object):

	def __init__(self, out_file = None, in_file = None, out_dict = None, in_nml = None, directory = "./"):
		self.directory = directory
		self.output_file = out_file
		self.input_file = in_file
		self.run = {}
		self.namelist = {}
		if out_dict:
			self.run = out_dict
		elif out_file:
			self.open_out()
		if in_nml:
			self.namelist = in_nml
		elif in_file:
			self.open_in()
		self.run['gr'] = imag(self.run['omega'][-1,0,0])
		self.run['mf'] = real(self.run['omega'][-1,0,0])
		
	def __getitem__(self, key):
		if key in self.run.keys():
			return self.run[key]
		elif key in self.namelist.keys():
			return self.namelist[key]
		else:
			for supkey in self.namelist.keys():
				if key in self.namelist[supkey].keys():
					return self.namelist[supkey][key]

	def keys(self):
		return self.run.keys()
				
	def open_out(self, filename = None, directory = None):
		if directory:
			self.directory = directory
		if filename:
			self.output_file = filename
		if self.output_file is None:
			print("ERROR: filename not given")
			return
		data_in = {}
		if self.output_file[-7:] == ".out.nc":
			try:
				data_in = readnc(os.path.join(self.directory,self.output_file))
			except:
				print(f"Could not load file {os.path.join(self.directory,self.output_file)}")
				return
		else:
			try:
				data_in = readnc(os.path.join(self.directory,self.output_file) + ".out.nc")
			except:
				print(f"Could not load file {os.path.join(self.directory,self.output_file)}.out.nc")
				return
		
		data_in['omega'] = data_in['omega'][:,0,0]
		data_in['phi'] = data_in['phi'][0,0,:]
		data_in['apar'] = data_in['apar'][0,0,:]
		data_in['bpar'] = data_in['bpar'][0,0,:]
		
		self.run = data_in
		
	def open_in(self, filename = None, directory = None):
		import f90nml
		if directory:
			self.directory = directory
		if filename:
			self.input_name = filename
		if self.input_file is None:
			print("ERROR: filename not given")
			return
		nml = f90nml.read(os.path.join(self.directory,self.input_file))
		self.namelist = nml
	
	@property
	def input(self):
		for key in self.namelist.keys():
			print(f"\n{key}:")
			for subkey in self.namelist[key].keys():
				print(f"\t{subkey} = {self.namelist[key][subkey]}")
				
	def write_input_file(self, filename = None):
		if self.namelist is None:
			print("ERROR: No Input Namelist Loaded")
		if filename is None and self.input_file is None:
			filename = input("Input File Name: ")
		elif filename is None:
			filename = self.input_file
		self.namelist.write(filename)
	
	def plot_omega(self):
		Plotters['Diag_Single'](data = self.run, var = 0)
	
	def plot_phi(self):
		Plotters['Diag_Single'](data = self.run, var = 1)
	
	def plot_apar(self):
		Plotters['Diag_Single'](data = self.run, var = 2)
		
	def plot_bpar(self):
		Plotters['Diag_Single'](data = self.run, var = 3)
		
	def plot_phi2(self):
		Plotters['Diag_Single'](data = self.run, var = 4)
	
	def _plot_diag(self, var = 0):
		Plotters['Diag_Single'](data = self.run, var = var)
		
	def plot_theta(self, var = 0):
		Plotters['Theta_Single'](run = self.run, var = var)
