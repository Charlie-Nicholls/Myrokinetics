import os
from numpy import load
from .plotting import Plotters
from .myro_single import myro_single

'''
DIAGNOSTIC SET ANALYSIS
'''

class myro_set_read(object):

	def __init__(self, filename = None, out_files = None, in_files = None,  directory = "./"):
		self.directory = directory
		self.runs = {}
		self.data = {}
		if filename:
			self.open_file(filename = filename, directory = directory)
		if out_files or in_files:
			self.load_runs(out_files = out_files, in_files = in_files, directory = directory)
	
	def __getitem__(self, key):
		if key in self.runs.keys():
			return self.runs[key]

	def keys(self):
		return self.runs.keys()
		
	def _get_new_id(self):
		if self.runs == {}:
			key = 0
		else:
			key = list(self.runs.keys())[-1] + 1
		self.runs[key] = {}
		return key
			
	def open_file(self, filename = None, directory = None):
		if directory:
			self.directory = directory
		if filename:
			self.filename = filename
		if self.filename is None:
			print("ERROR: filename not given")
			return
		data_in = {}
		if self.filename[-4:] == ".npz":
			try:
				data_in = load(os.path.join(self.directory,self.filename), allow_pickle=True)
			except:
				print(f"Could not load file {os.path.join(self.directory,self.filename)}")
				return
		else:
			try:
				data_in = load(os.path.join(self.directory,self.filename) + ".npz", allow_pickle=True)
			except:
				print(f"Could not load file {os.path.join(self.directory,self.filename)}.npz")
				return
		
		try:
			inputs = data_in['inputs'].item()
		except:
			print("ERROR: could not load Inputs")
			inputs = None
		try:
			info = data_in['run_info'].item()
		except:
			print("ERROR: could not load Run Info")
			info = None
		try:
			input_namelists = data_in['input_namelists']
		except:
			print("ERROR: could not load Input Files")
			input_namelists = None
		try:
			outputs = data_in['output_dicts']
		except:
			print("ERROR: could not load Outputs")
			outputs = None
		try:
			files = data_in['files'].item()
		except:
			print("ERROR: could not load eqbm files")
			files = None
			
		self.data = {'inputs': inputs, 'info': info, 'input_namelists': input_namelists, 'outputs': outputs, 'files': files}
		
		for p, psiN in enumerate(self.data['inputs']['psiNs']):
			for v, value in enumerate(self.data['inputs']['values']):
				for k, aky in enumerate(self.data['inputs']['aky_values']):
					run_id = self._get_new_id()
					if outputs[p][v][k] is not None and outputs[p][v][k] != {}:
						self.runs[run_id]['data'] = myro_single(out_dict = outputs[p][v][k], in_nml = input_namelists[p][v][k], directory = directory)
					else:
						self.runs[run_id]['data'] = None
					self.runs[run_id]['psiN'] = psiN
					self.runs[run_id]['aky'] = aky
					self.runs[run_id][self.data['inputs']['variable']] = value
			
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
		self.runs[run_id]['data'] = myro_single(out_file = out_file, in_file = in_file, directory = directory)
	
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
		
	def plot_omega(self, init = None):
		self._plot_diag(var = 0, init = init)
	
	def plot_phi(self, init = None):
		self._plot_diag(var = 1, init = init)
	
	def plot_apar(self, init = None):
		self._plot_diag(var = 2, init = init)
		
	def plot_phi2(self, init = None):
		self._plot_diag(var = 3, init = init)
	
	def _plot_diag(self, var = 0, init = None):
		Plotters['Diag_Set'](runs = self.runs, var = var, init = init)
	
	def plot_set(self, var = None, init = 0):
		Plotters['Set'](runs = self.runs, var = var, init = init)
