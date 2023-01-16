import os
from numpy import *
from ncdf2dict import ncdf2dict as readnc
from .plotting import Plotters
from .myro_single import myro_single
from matplotlib.pyplot import * #TEMPORARY

'''
DIAGNOSTIC SET ANALYSIS
'''

class myro_set(object):

	def __init__(self, out_files = None, in_files = None, var_name = None, directory = "./"):
		self.directory = directory
		self.runs = {}
		if out_files or in_files:
			self.load_runs(out_files = out_files, in_files = in_files, directory = directory)
		
	def __getitem__(self, key):
		if key in self.runs.keys():
			return self.runs[key]

	def keys(self):
		return self.runs.keys()
		
	def _get_new_id(self):
		if self.runs == {}:
			return 0
		return list(self.runs.keys())[-1] + 1
				
	def load_runs(self, out_files = None, in_files = None, directory = None):
		if type(out_files) == str:
			out_files = [outfiles]
		if type(in_files) == str:
			in_files = [in_files]
			
		if out_files is None and in_files is None:
			print("ERROR: No output or input files passed")
			return
		elif in_files is None:
			in_files = [None] * len(out_files)
		elif out_files is None:
			out_files = [None] * len(in_files)
		elif len(out_files) != len(in_files):
			print("ERROR: Number of input and output files is not equal")
			return
			
		for out_file, in_file in zip(out_files, in_files):
			self.load_run(out_file = out_file, in_file = in_file, directory = directory)
		
		self._find_diff()
		
	def load_run(self, out_file = None, in_file = None, directory = None):
		if directory is None and self.directory is None:
			directory = "./"
		elif directory is None:
			directory = self.directory
		run_id = self._get_new_id()
		self.runs[run_id] = {'data': myro_single(out_file = out_file, in_file = in_file, directory = directory)}
	
	def _find_diff(self):
		nml = self.runs[0]['data'].namelist
		vs = set()
		for skey in nml:
			for key in nml[skey].keys():
				for i in self.runs.keys():
					if nml[skey][key] != self.runs[i]['data'].namelist[skey][key]:
						vs.add((skey,key))
		for i in self.runs.keys():
			for skey, key in vs:
				self.runs[i][key] = self.runs[i]['data'].namelist[skey][key]
		
	def plot_omega(self):
		self._plot_diag(var = 0)
	
	def plot_phi(self):
		self._plot_diag(var = 1)
	
	def plot_apar(self):
		self._plot_diag(var = 2)
		
	def plot_phi2(self):
		self._plot_diag(var = 3)
	
	def _plot_diag(self, var = 0):
		Plotters['Diag_Set'](runs = self.runs, var = var)
	
	def plot_set(self, var = None):
		if var is None:
			var = [x for x in self.runs[0].keys() if x != 'data'][0]
		if var not in self.runs[0].keys():
			print("ERROR: Invalid Key")
			return

		fig, ax = subplots(2,1)
		ax[1].plot([self.runs[i][var] for i in self.runs.keys()],[self.runs[i]['data'].run['gr'] for i in self.runs.keys()],'b.')
		ax[0].plot([self.runs[i][var] for i in self.runs.keys()],[self.runs[i]['data'].run['mf'] for i in self.runs.keys()],'r.')
		ax[1].set_xlabel(var)
		ax[1].set_ylabel("Growth Rate")
		ax[0].set_ylabel("Mode Frequency")
		show()
