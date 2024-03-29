import os
from numpy import load, savez, nan, isfinite
from .plotting import Plotters
from .verify_runs import verify_scan
from .inputs import scan_inputs

'''
GYROKINETIC SCAN ANALYSIS
'''

class myro_read(object):

	def __init__(self, filename = None, directory = "./", verify = True):
		if directory == "./" or directory is None:
			directory = os.getcwd() 
		self.directory = directory
		self.filename = filename
		self.run = {}
		opened = self._open_file()
		self._gr_type = "Normalised"
		self.eqbm = None
		if opened and not self.verify and verify:		
			self._verify_run()
	
	def __call__(self, key, *indexes):
		if type(indexes[0]) == dict:
			run = indexes[0]
			ids = None
		elif type(indexes[0]) == int:
			run = None
			ids = list(indexes)
		else:
			run = None
			ids = list(indexes[0])
		
		key.lower()
		if key in ['gr','growth','gamma']:
			key = 'growth_rate'
		elif key in ['mf','mode','frequency']:
			key = 'mode_frequency'
		elif key in ['time']:
			key = 't'
		elif key in ['ql']:
			key = 'quasilinear'
		elif key in ['norm_growth_rate','normalised_gr','normalised_growth_rate']:
			key = 'norm_gr'
		elif key in ['abs_growth_rate','absolute_gr','absolute_growth_rate']:
			key = 'abs_gr'
		
		if key in ['quasilinear','norm_gr','abs_gr']:
			if f"_{key}_keys" not in self.data.keys():
				print(f"ERROR: {key} not calculated")
				return None
			
			if run is None:
				dim_order = [x for x in self.inputs.dim_order if x not in ['ky']] #Update to exclude theta0
				if len(ids) != len(dim_order):
					print(f"ERROR: ids must be length {len(dim_order)}")
				run = {}
				for idx, i in enumerate(ids):
					run[dim_order[idx]] = self.dimensions[dim_order[idx]].values[i]
			run_id = self.get_run_id(run, keys = f'_{key}_keys')
			if key == 'quasilinear':
				return self.data['quasilinear'][run_id]
			elif key == 'abs_gr':
				return self.data['gyro'][run_id]['growth_rate']
			elif key == 'norm_gr':
				return self.data['gyro'][run_id]['growth_rate']/self.data['gyro'][run_id]['ky']**2
		
		else:
			
			if run is None:
				if len(ids) != len(self.dimensions):
					print(f"ERROR: ids must be length {len(self.dimensions)}")
					return None
				run = {}
				for idx, i in enumerate(ids):
					run[self.inputs.dim_order[idx]] = self.dimensions[self.inputs.dim_order[idx]].values[i]
			run_id = self.get_run_id(run)
			if run_id is None:
				return None
			if key in self.data['gyro'][run_id].keys():
				return self.data['gyro'][run_id][key]
			elif 'group_key' in self.data['gyro'][run_id].keys() and key in self.data['group'][self.data['gyro'][run_id]['group_key']].keys():
				return self.data['group'][self.data['gyro'][run_id]['group_key']][key]
			elif key in ['data','all']:
				return self.data['gyro'][run_id]
			elif key in ['nt']:
				return self.verify.nts[run_id]
			else:
				print(f"ERROR: {key} not found")
			return None
		
	def __getitem__(self, key):
		key = key.lower()
		if key == "inputs":
			return self.inputs.inputs
			
		elif key in self.data.keys():
			return self.data[key]
		elif key in self.files.keys():
			return self.files[key]
		
		elif key in ["namelist_diffs", "nml_diffs", "namelist_diff", "nml_diff", "namelists", "nmls"]:
			return self.files['namelist_differences']
		elif key in ["ql", "quasi linear", "quasi_linear"]:
			return self.data['quasilinear']
			
		elif self.dimensions and key in self.dimensions:
			return self.dimensions.values
		
		elif self.verify and key in self.verify._all_keys():
			return self.verify[key]
			
		elif self.inputs:
			return self.inputs[key]
			#inputs also contains logic for if key not found
			
	def __len__(self):
		tot = 1
		for dim in self.dimensions.values():
			tot *= len(dim)
		return tot
	
	def get_run_from_id(self, run_id):
		run = {}
		for dim in self.dimensions:
			run[dim] = self.data['gyro'][run_id][dim]
		return run
	
	def get_run_id(self, run, keys = '_gyro_keys'):
		run_id = self.get_run_list(run, keys = keys)
		if run_id is None:
			return None
		if len(run_id) == 1:
			return run_id[0]
		else:
			return None
	
	def get_run_list(self, run, keys = '_gyro_keys'):
		if run == {}:
			return [x for x in self.data['gyro'].keys() if x != 'group']
		idlist = []
		for key, val in run.items():
			if key in self.data[keys] and val in self.data[keys][key]:
				idlist.append(self.data[keys][key][val])
		if len(idlist) == 0:
			return None
		elif len(idlist) == 1:
			return list(idlist[0])
		else:
			run_id = list(set.intersection(*map(set, idlist)))
			if len(run_id) > 0:
				return run_id
			else:
				return None
        		
	def print_inputs(self):
        	self.inputs.print_inputs()

	def print_info(self):
		for key, val in self.inputs['info'].items():
        		print(f"{key} = {val}")
		
	def _print_file(self, filetype = ''):
		if self.files is None:
			print("ERROR: No files found")
			return
		filetype = filetype.lower()
		if filetype == '':
			print("ERROR: filetype not given")
			return
		if filetype in ['template_file', 'template', 'gk_template']:
			print(self['template_file'])
			return
		elif filetype in ['eq','eq_file','geqdsk','geqdsk_file']:
			lines = self['eq_file']
		elif filetype in ['kin','kin_file','kinetics','kinetics_file','peqdsk','peqdsk_file']:
			lines = self['kin_file']
		else:
			print(f"ERROR: File Type {filetype} Not Found")
			return
		if lines is None:
			print("ERROR: File Not Found")
		for line in lines:
			print(line,end='')

	def template_file(self):
		self._print_file(filetype = 'template')

	def eq_file(self):
		self._print_file(filetype = 'eq')

	def kin_file(self):
		self._print_file(filetype = 'kin')
		
	def _write_file(self, filetype, filename = None, directory = None):
		if self.files is None:
			print("ERROR: No files found")
			return
		if directory is None:
			directory = self.directory
		filetype = filetype.lower()
		if filetype in ['template_file', 'template', 'gk_template']:
			lines = self['template_file']
			name = self.inputs['template_name']
		elif filetype in ['eq','eq_file','geqdsk','geqdsk_file']:
			lines = self['eq_file']
			name = self.inputs['eq_name']
		elif filetype in ['kin','kin_file','kinetics','kinetics_file','peqdsk','peqdsk_file']:
			lines = self['kin_file']
			name = self.inputs['kin_name']
		else:
			print(f"ERROR: File Type {filetype} Not Found")
			return
		if lines is None:
			print(f"ERROR: {filetype} lines not saved")
			return
		if filename:
			name = filename
		if not directory:
			directory = self.directory
		if os.path.exists(f"{directory}/{name}"):
			print(f"{os.path.join(directory,name)} Already Exists")
			return
		if type(lines) != list:
			lines.write(os.path.join(directory,name))
			print(f"Created {name} at {directory}")
			return
		f = open(os.path.join(directory,name),'w')
		for line in lines:
			f.write(line)
		f.close()
		print(f"Created {name} at {directory}")
		
	def write_template_file(self, filename = None, directory = None):
		self._write_file(filetype = 'template', filename = filename, directory = directory)
	
	def write_eq_file(self, filename = None, directory = None):
		self._write_file(filetype = 'eq', filename = filename, directory = directory)
	
	def write_kin_file(self, filename = None, directory = None):
		self._write_file(filetype = 'kin', filename = filename, directory = directory)
		
	def write_scan_input(self, filename = None, directory = "./", doPrint = True):
		self.inputs.write_scan_input(filename = filename, directory = directory, doPrint = doPrint)
		
	def _open_file(self, filename = None, directory = None):
		if directory:
			self.directory = directory
		if filename:
			self.filename = filename
		if self.filename is None:
			print("ERROR: filename not given")
			return

		if self.filename[-4:] == ".npz":
			try:
				data_in = load(os.path.join(self.directory,self.filename), allow_pickle=True)
			except:
				print(f"Could not load file {os.path.join(self.directory,self.filename)}")
				return False
		else:
			try:
				data_in = load(os.path.join(self.directory,self.filename) + ".npz", allow_pickle=True)
			except:
				print(f"Could not load file {os.path.join(self.directory,self.filename)}.npz")
				return False
		try:
			self.data = data_in['data'].item()
			possible_data = ['gyro','ideal','equilibrium','_gyro_keys','_ideal_keys','quasilinear']
			for key in [x for x in possible_data if x not in self.data.keys()]:
				self.data[key] = None
		except:
			print("ERROR: could not load Data")
			self.data = {}
		try:
			self.files = data_in['files'].item()
			possible_files = ['eq_file','kin_file','template_file','namelist_differences']
			for key in [x for x in possible_files if x not in self.files.keys()]:
				self.files[key] = {}
		except:
			print("ERROR: could not load Files")
			self.files = {}
		try:
			self.inputs = scan_inputs(input_dict = data_in['inputs'].item())
			self.dimensions = self.inputs.dimensions
			self.single_parameters = self.inputs.single_parameters
		except:
			print("ERROR: could not load Inputs")
			self.inputs = {}
			self.dimensions = {}
			self.single_parameters = {}
		try:
			self.verify = data_in['verify'].item()
		except:
			self.verify = None

		return True
	
	def save_file(self, filename = None, directory = "./"):
		if filename == None:
			filename = self['run_name']
		savez(f"{directory}/{filename}", inputs = self.inputs.inputs, data = self.data, files = self.files, verify = self.verify)
		
	def _verify_run(self):
		if not self['Gyro']:
			return
		self.verify = verify_scan(reader = self)
		self.data['gyro'] = self.verify.scan
		#self.calculate_gr()
	
	def get_all_runs(self, excludeDimensions = []):
		dim_order = [x for x in self.inputs.dim_order if x not in excludeDimensions]
		if len(dim_order) == 0:
		 	return [{}]
		def loop(n,variables={},runs=[]):
			dim = self.dimensions[dim_order[len(dim_order)-n]]
			for val in dim.values:
				variables[dim.name] = val
				if n>1:
					loop(n=n-1,variables=variables)
				else:
					runs.append(variables.copy())
			if n == len(dim_order):
				return runs
		return loop(n=len(dim_order))
	
	def calculate_ql(self):
		#reexclude theta0 when ql calculation is updated in line 322 & 327 & 330 & 56 & plotters
		from .quasilinear import QL
		from uuid import uuid4
		if 'ky' not in self.dimensions:
			print("ERROR: Requires ky dimension")
			return
		qls = {}
		ql_keys = {}
		for dim in [x for x in self.dimensions.values() if x.name not in ['ky']]:
			ql_keys[dim.name] = {}
			for val in dim.values:
				ql_keys[dim.name][val] = []
				
		for runs in self.get_all_runs(excludeDimensions = ['ky']):
			run_ids = self.get_run_list(runs)
			ql_key = str(uuid4())
			for dim_name, val in [(x, y) for x, y in self.data['gyro'][run_ids[0]].items() if (x in self.dimensions and x not in ['ky'])]:
				ql_keys[dim_name][val].append(ql_key)
			ql, [ql_norms,kys] = QL(run_ids,self.data['gyro'], returnlist = True)
			qls[ql_key] = ql
			for run_id in run_ids:
				if self.data['gyro'][run_id]['ky'] in kys:
					self.data['gyro'][run_id]['ql_norm'] = ql_norms[kys.index(self.data['gyro'][run_id]['ky'])]
				else:
					self.data['gyro'][run_id]['ql_norm'] = nan
			
		self.data['quasilinear'] = qls
		self.data['_quasilinear_keys'] = ql_keys
		
	def calculate_gr(self):
		abs_gr_keys = {}
		norm_gr_keys = {}
		for dim in [x for x in self.dimensions.values() if x.name not in ['ky','theta0']]:
			abs_gr_keys[dim.name] = {}
			norm_gr_keys[dim.name] = {}
			for val in dim.values:
				abs_gr_keys[dim.name][val] = []
				norm_gr_keys[dim.name][val] = []
		for runs in self.get_all_runs(excludeDimensions = ['ky','theta0']):
			run_ids = self.get_run_list(runs)
			abs_grs = []
			norm_grs = []
			for run_id in run_ids:
				abs_grs.append(self.data['gyro'][run_id]['growth_rate'])
				if 'ky' in self.data['gyro'][run_id]:
					ky = self.data['gyro'][run_id]['ky']
				elif 'ky' in self.single_parameters:
					ky = self.single_parameters['ky'].values[0]
				else:
					ky = nan
				if ky == 0:
					if self.data['gyro'][run_id]['growth_rate'] == 0:
						norm_gr = 0
					else:
						norm_gr = nan
				else:
					norm_gr = self.data['gyro'][run_id]['growth_rate']/(ky**2)
				norm_grs.append(norm_gr)
				self.data['gyro'][run_id]['growth_rate_norm'] = norm_gr
			
			if len([x for x in abs_grs if isfinite(x)]) == 0:
				if len(run_ids) == 0:
					abs_id = None
				else:
					abs_id = run_ids[0]
			else:
				abs_gr = max([x for x in abs_grs if isfinite(x)])
				abs_id = run_ids[abs_grs.index(abs_gr)]
			if len([x for x in norm_grs if isfinite(x)]) == 0:
				norm_id = run_ids[0]
			else:
				norm_gr = max([x for x in norm_grs if isfinite(x)])
				norm_id = run_ids[norm_grs.index(norm_gr)]

			for dim_name, val in [(x, y) for x, y in self.data['gyro'][abs_id].items() if (x in self.dimensions and x not in ['ky','theta0'])]:
				abs_gr_keys[dim_name][val].append(abs_id)
			for dim_name, val in [(x, y) for x, y in self.data['gyro'][norm_id].items() if (x in self.dimensions and x not in ['ky','theta0'])]:
				norm_gr_keys[dim_name][val].append(norm_id)
		
		self.data['_abs_gr_keys'] = abs_gr_keys
		self.data['_norm_gr_keys'] = norm_gr_keys
	
	'''
	def calculate_alpha(self):
		if not self.eqbm:
			self.load_equilibrium()
		from scipy.interpolate import InterpolatedUnivariateSpline
		qspline = InterpolatedUnivariateSpline(self.eqbm.eq_data['psiN'],self.eqbm.eq_data['qpsi'])
		alpha_axis = []
		alpha_values = []
		alpha_axis_ideal = []
		for rid, run in self.data['gyro'].items():
			self.data['gyro'][rid]['alpha'] = abs(run['beta_prime'])*self.eqbm.eq_data['rmaxis']*qspline(run['psin'])**2
		for 
			alpha_axis_ideal.append([abs(x)*self.eqbm.eq_data['rmaxis']*qspline(psiN)**2 for x in self['beta_prime_axis_ideal'][p]])
	'''
	
	def make_slider_axes(self, settings = {}, init = None):
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		return Plotters['Sliders'](reader = self, settings = settings)			
		
	def plot_aky(self, settings = {}, init = None):
		return self.plot_scan(init = init, aky = True, settings = settings, sliders = None)
		
	def plot_scan(self, settings = {}, init = None, aky = None, sliders = None):
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		if aky is not None:
			settings['aky'] = aky
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} Scan"
		return Plotters['Scan'](reader = self, settings = settings, sliders = sliders)
		
	def plot_kxky(self, settings = {}, init = None, sliders = None):
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} kxky"
		return Plotters['kxky'](reader = self, settings = settings, sliders = sliders)
	
	def plot_ql(self, settings = {}, init = None, sliders = None):
		if self['ql'] is None:
			self.calculate_ql()
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} QuasiLinear"
		return Plotters['2D'](reader = self, settings = settings, sliders = sliders)
	
	def plot_ideal(self, settings = {}, init = None):
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} Ideal Ballooning"
		return Plotters['Ideal'](reader = self, settings = settings)
	
	def plot_omega(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'omega', init = init, settings = settings, sliders = sliders)
	
	def plot_phi(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'phi', init = init, settings = settings, sliders = sliders)
	
	def plot_apar(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'apar', init = init, settings = settings, sliders = sliders)
		
	def plot_bpar(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'bpar', init = init, settings = settings, sliders = sliders)
		
	def plot_epar(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'epar', init = init, settings = settings, sliders = sliders)
		
	def plot_phi2(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'phi2', init = init, settings = settings, sliders = sliders)
	
	def plot_phi2_avg(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'phi2_avg', init = init, settings = settings, sliders = sliders)
		
	def plot_jacob(self, settings = {}, init = None, sliders = None):
		return self._plot_diag(var = 'jacob', init = init, settings = settings, sliders = sliders)
	
	def _plot_diag(self, settings = {}, init = None, var = None, sliders = None):
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		if var is not None:
			settings['var'] = var
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} {var}"
		return Plotters['Diag'](reader = self, settings = settings, sliders = sliders)
	
	def plot_box_phi(self, settings = {}, init = None):
		return self._plot_box_diag(var = 'phi', init = init, settings = settings)
	
	def plot_box_apar(self, settings = {}, init = None):
		return self._plot_box_diag(var = 'apar', init = init, settings = settings)
		
	def plot_box_bpar(self, settings = {}, init = None):
		return self._plot_box_diag(var = 'bpar', init = init, settings = settings)
		
	def plot_box_epar(self, settings = {}, init = None):
		return self._plot_box_diag(var = 'epar', init = init, settings = settings)
	
	def _plot_box_diag(self, settings = {}, init = None, var = None):
		if self.inputs['grid_option'] != 'box':
			print("ERROR: box plot only used for flux tube runs [inputs -> grid_option = box]")
			return
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		if var is not None:
			settings['var'] = var
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} {var}"
		return Plotters['Box_Diag'](reader = self, settings = settings)
	
	#def plot_epar_scan(self):
		#Plotters['Epar'](data = self.run['data'], inputs = self.inputs)
	
	def plot_theta(self, settings = {}, init = None, var = None, periods = None, polar = None):
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
					settings[f"slider_{i+1}"]['dimension_type'] = self.inputs.dim_order[i]
		if var is not None:
			settings['var'] = var
		if periods is not None:
			settings['periods'] = periods
		if polar is not None:
			settings['polar'] = polar
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} {var}"
		return Plotters['Theta'](reader = self, settings = settings)
		
	def plot_slice(self,  settings = {}, x_dim = None, y_dim = None, init = None, limit = None, sliders = None):
		if y_dim in ['quasilinear','ql_norm'] and self['ql'] is None:
			self.calculate_ql()
		if init is not None:
			init = list(init)
			for i, ini in enumerate(init):
				if i < len(self.inputs.dim_order):
					if f"slider_{i+1}" not in settings:
						settings[f"slider_{i+1}"] = {}
					settings[f"slider_{i+1}"]['id'] = ini
		if x_dim is not None:
			settings['x_axis_type'] = x_dim
		if y_dim is not None:
			settings['y_axis_type'] = y_dim
		if limit is not None:
			settings['limit'] = limit
		return Plotters['Slice'](reader = self, settings = settings, sliders = sliders)
	
	def plot_eq(self):
		if not self.eqbm:
			self.load_equilibrium()
		self.eqbm.plot_eq()
	
	def load_equilibrium(self, directory = None):
		from .equilibrium import equilibrium
		if directory is None:
			directory = self.directory
		if not os.path.exists(f"{os.path.join(directory,self.inputs['eq_name'])}"):
			self.write_eq_file(directory = directory)
		if not os.path.exists(f"{os.path.join(directory,self.inputs['kin_name'])}"):
			self.write_kin_file(directory = directory)
		if not os.path.exists(f"{os.path.join(directory,self.inputs['template_name'])}"):
			self.write_template_file(directory = directory)
			
		self.eqbm = self.equilibrium = equilibrium(inputs = self.inputs, directory = directory)
		self.eqbm.load_pyro(directory=directory)
	
	def write_gs2_input(self, indexes = None, run = None, filename = None, directory = None):
		if directory is None and self.directory is None:
			directory = "./"
		elif directory is None:
			directory = self.directory
		if run is None and indexes is None:
			print("ERROR: run and indexes can not both be None")
			return

		if filename is None:
			filename = f"itteration_{self.inputs['itteration']}.in"
		
		if run is None:
			if len(indexes) != len(self.inputs.dimensions):
				print(f"ERROR: indexes must be of length {len(self.dimensions)}, {[self.inputs.dim_order]}")
				return None
			run = {}
			for i, dim in zip(indexes,self.inputs.dimensions.values()):
				run[dim.name] = dim.values[i]
				
		run_id = self.get_run_id(run=run)
		
		if self.eqbm is None:
			self.load_equilibrium(directory = directory)
		if 'nml_diffs' in self.data['gyro'][run_id].keys():
			namelist_diff = self.data['gyro'][run_id]['nml_diffs']
		else:
			namelist_diff = {}
		
		nml = self.eqbm.get_gyro_input(run = run, indexes = indexes, namelist_diff = namelist_diff)
		
		if not nml:
			return
		
		nml.write(os.path.join(directory,filename), force=True)
		print(f"Created {filename} at {directory}")
	
	def merge_myro(self, myro):
		dims = set(self.dimensions.keys()).union(set(self.single_parameters.keys()))
		compare_dims = set(myro.dimensions.keys()).union(set(myro.single_parameters.keys()))
		if dims != compare_dims:
			print(f"ERROR: merging myro objects requires the same set of combined dimensions & single parameters; {dims} vs {compare_dims}")
			return
		from copy import deepcopy
		inputs_nml = deepcopy(self.inputs.inputs)
		del(inputs_nml['single_parameters'])
		new_dim_order = self.inputs.dim_order
		for dim in dims:
			if dim not in new_dim_order:
				new_dim_order.append(dim)
		for i, dim in enumerate(new_dim_order):
			in_dict = {'type': dim}
			opt = opt2 = None
			if dim in self.dimensions:	
				vals = set(self.dimensions[dim].values)
				opt = self.dimensions[dim].option
			elif dim in self.single_parameters:
				vals = set(self.single_parameters[dim].values)
				opt = self.single_parameters[dim].option
			if dim in myro.dimensions:	
				vals = list(vals.union(set(myro.dimensions[dim].values)))
				opt2 = myro.dimensions[dim].option
			elif dim in myro.single_parameters:
				vals = list(vals.union(set(myro.single_parameters[dim].values)))
				opt2 = myro.single_parameters[dim].option
			vals.sort()
			in_dict['values'] = vals
			if opt is not None and opt2 is not None and opt != opt2:
				print(f"Warning: option for {dim} differs between scans; {opt} vs {opt2}")
			in_dict['option'] = opt
			inputs_nml[f'dimension_{i}'] = in_dict
		new_inputs = scan_inputs(input_dict = inputs_nml)
		new_inputs.check_inputs()
		new_inputs.load_dimensions()
		for dim in new_inputs.dimensions:
			if dim not in self.data['_gyro_keys']:
				self.data['_gyro_keys'][dim] = {}
			for val in new_inputs.dimensions[dim].values:
				if val not in self.data['_gyro_keys'][dim]:
					self.data['_gyro_keys'][dim][val] = set()
		for dim in new_inputs.dimensions:
			if dim in self.single_parameters:
				for rid, run in self.data['gyro'].items():
					self.data['gyro'][rid][dim] = self.single_parameters[dim].values[0]
					self.data['_gyro_keys'][dim][self.single_parameters[dim].values[0]].add(rid)
			if dim in myro.single_parameters:
				for rid, run in myro.data['gyro'].items():
					myro.data['gyro'][rid][dim] = myro.single_parameters[dim].values[0]
		for rid, run in myro.data['gyro'].items():
			self.data['gyro'][rid] = run
			for dim in new_inputs.dimensions:
				self.data['_gyro_keys'][dim][run[dim]].add(rid)
		self.data['group'].update(myro.data['group'])
		self.inputs = new_inputs
		self.dimensions = self.inputs.dimensions
		self.single_parameters = self.inputs.single_parameters
		self._verify_run()
		
	'''
	def _return_mf_set(self, psi_id, ky_id, mf = None, mferr = None, mfmax = None, mfmin = None, smin_id = None, smax_id = None, bmin_id = None, bmax_id = None):
		if (mf is None or mferr is None) and (mfmax is None or mfmin is None):
			print("ERROR: insufficient mf specifications")
			return
		mfs = array(self['mode_frequencies_all'])[psi_id,:,:,ky_id].tolist()
		if type(mf) == list:
			mf = mfs[mf[0]][mf[1]]
		if mfmax is None:
			mfmax = mf + mferr
		if mfmin is None:
			mfmin = mf - mferr
		runs = set()
		for i in range(self['n_beta']):
			for j in range(self['n_shat']):
				if mfmin < mfs[i][j] and mfmax > mfs[i][j]:
					runs.add((psi_id,i,j,ky_id))
		return runs
	'''
			
	def _save_run_set(self, runs = None, filename = None):
		if runs is None:
			print("ERROR: runs not given")
			return
		if filename is None:
			print("ERROR: filename not given")
			return

		f = open(filename, 'w')
		for p,i,j,k in runs:
			f.write(f"{p}_{i}_{j}_{k}\n")
		f.close()
		
	def _load_run_set(self, filename = None):
		if filename is None:
			print("ERROR: filename not given")
			return
		
		runs = set()		
		with open(filename) as f:
			lines = f.readlines()
			for line in lines:
				p,i,j,k = [eval(x) for x in line.strip("\n").split("_")]
				runs.add((p,i,j,k))
		return runs
