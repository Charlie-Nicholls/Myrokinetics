import os
from numpy import load, savez, array, nan, full
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
		self.plots = {}
		opened = self._open_file()
		self._gr_type = "Normalised"
		self.eqbm = None
		#if opened and not self.verify and verify:		
		#	self._verify_run()
	
	def __call__(self, key, *indexes):
		if type(indexes[0]) == int:
			ids = list(indexes)
		else:
			ids = list(indexes[0])
	
		key.lower()
		if key in ['gr','growth','gamma']:
			key = 'growth_rate'
		if key in ['mf','mode','frequency']:
			key = 'mode_frequency'
		if key in ['time']:
			key = 't'
		
		if key in ['growth_rate','mode_frequency','omega','phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob','epar_norm']:
			if len(ids) != len(self.dimensions):
				print(f"ERROR: ids must be length {len(self.dimensions)}")
				return None
			run = {}
			for idx, i in enumerate(ids):
				run[self.inputs.dim_order[idx]] = self.dimensions[self.inputs.dim_order[idx]].values[i]
			run_id = self.get_run_id(run)
			return self.gyro_data[run_id][key]
		
		if key in ['ql','quasilinear']:
			if not self.data['quasilinear']:
				print("ERROR: quasilinear not calculated")
				return None
			
			dim_order = [x for x in self.inputs.dim_order if x not in ['ky','theta0']]
			if len(ids) != len(dim_order):
				print(f"ERROR: ids must be length {len(dim_order)}")
			run = {}
			for idx, i in enumerate(ids):
				run[dim_order[idx]] = self.dimensions[dim_order[idx]].values[i]
			run_id = self.get_run_id(run, keys = '_ql_keys')
			return self.data['quasilinear'][run_id]
			
		
	def __getitem__(self, key):
		key = key.lower()
		if key == "inputs":
			return self.inputs.inputs
			self.print_info()
			
		elif key in self.data.keys():
			return self.data[key]
		elif key in self.info.keys():
			return self.info[key]
		elif key in self.files.keys():
			return self.files[key]
		
		elif key in ["namelist_diffs", "nml_diffs", "namelist_diff", "nml_diff", "namelists", "nmls"]:
			return self.files['namelist_differences']
		elif key in ["ql", "quasi linear", "quasi_linear"]:
			return self.data['quasilinear']
			
		elif key in self.verify._all_keys():
			return self.verify[key]
			
		else:
			return self.inputs[key]
			#inputs also contains logic for if key not found
			
	def __len__(self):
		tot = 1
		for dim in self.dimensions.values():
			tot *= len(dim)
		return tot
	
	def get_run_id(self, run, keys = '_run_keys'):
		run_id = self.get_run_list(run, keys = keys)
		if len(run_id) == 1:
			return run_id[0]
		elif len(run_id) > 1:
			print("ERROR: multiple runs found")
			return run_id
		else:
			print("ERROR: run not found")
			return None
	
	def get_run_list(self, run, keys = '_run_keys'):
		idlist = []
		for key, val in run.items():
			idlist.append(self.data[keys][key][val])
		run_id = list(set.intersection(*map(set, idlist)))
		if len(run_id) > 0:
			return run_id
		else:
			print("ERROR: run not found")
			return None
        		
	def print_inputs(self):
        	self.inputs.print_inputs()

	def print_info(self):
		for key, val in self.run['info'].items():
        		print(f"{key} = {val}")
		
	def _print_file(self, filetype = ''):
		if self.run['files'] is None:
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
		
	def _write_file(self, filetype = '', filename = None, directory = None):
		if self.run['files'] is None:
			print("ERROR: No files found")
			return
		filetype = filetype.lower()
		if filetype == '':
			print("ERROR: filetype not given")
			return
		if filetype in ['template_file', 'template', 'gk_template']:
			lines = self['template_file']
			name = self['template_file_name']
		elif filetype in ['eq','eq_file','geqdsk','geqdsk_file']:
			lines = self['eq_file']
			name = self['eq_file_name']
		elif filetype in ['kin','kin_file','kinetics','kinetics_file','peqdsk','peqdsk_file']:
			lines = self['kin_file']
			name = self['kin_file_name']
		else:
			print(f"ERROR: File Type {filetype} Not Found")
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
		
	def write_input_file(self, filename = None, directory = "./", doPrint = True):
		self.inputs.write_input_file(filename = filename, directory = directory, doPrint = doPrint)
		
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
			self.info = data_in['run_info'].item()
			possible_info = ['run_name','run_uuid','data_path','input_file','eq_file_name','template_file_name','kin_file_name',
			'kinetics_type','run_data','_eq_file_path','_kin_file_path','_template_file_path']
			for key in [x for x in possible_info if x not in self.info.keys()]:
				self.info[key] = None
		except:
			print("ERROR: could not load Run Info")
			self.info = None
		try:
			self.gyro_data = data_in['gyro_data'].item()
		except:
			print("ERROR: could not load Gyro Data")
			self.gyro_data = None
		try:
			self.data = data_in['data'].item()
			possible_data = ['ideal_data','equillibrium','_run_keys','quasilinear','max_growth_rate']
			for key in [x for x in possible_data if x not in self.data.keys()]:
				self.data[key] = None
		except:
			print("ERROR: could not load Data")
			self.data = None
		try:
			self.files = data_in['files'].item()
			possible_files = ['eq_file','kin_file','template_file','namelist_differences']
			for key in [x for x in possible_files if x not in self.files.keys()]:
				self.files[key] = None
		except:
			print("ERROR: could not load Input Files")
			self.files = None
		try:
			self.inputs = scan_inputs(input_dict = data_in['inputs'].item())
			self.dimensions = self.inputs.dimensions
			if self.info:
				self.inputs.input_file = self.info['input_file']
		except:
			print("ERROR: could not load Inputs")
			self.inputs = None
			self.dimensions = None
		try:
			self.verify = data_in['verify'].item()
		except:
			self.verify = {}

		return True
	
	def save_file(self, filename = None, directory = "./"):
		if filename == None:
			filename = self['run_name']
		savez(f"{directory}/{filename}", inputs = self.inputs.inputs, gyro_data = self.gyro_data, run_info = self.info, files = selfiles, verify = self.verify)
		
	def _verify_run(self):
		if not self['Gyro']:
			return
		self.verify = verify_scan(scan = self.run)
		self.run = self.verify.scan
		self._convert_gr(gr_type = self._gr_type, doPrint = False)
	
	def get_all_runs(self, excludeDimensions = []):
		dim_order = [x for x in self.inputs.dim_order if x not in excludeDimensions]
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
		from .quasilinear import QL
		from uuid import uuid4
		if 'ky' not in self.dimensions:
			print("ERROR: Requires ky dimension")
			return
		qls = {}
		ql_keys = {}
		for dim in [x for x in self.dimensions.values() if x.name not in ['ky','theta0']]:
			ql_keys[dim.name] = {}
			for val in dim.values:
				ql_keys[dim.name][val] = []
				
		for runs in self.get_all_runs(excludeDimensions = ['ky','theta0']):
			run_ids = self.get_run_list(runs)
			ql_key = str(uuid4())
			for dim_name, val in [(x, y) for x, y in self.gyro_data[run_ids[0]].items() if (x in self.dimensions and x not in ['ky','theta0'])]:
				ql_keys[dim_name][val].append(ql_key)
			qls[ql_key] = QL(run_ids,self.gyro_data)
			
		self.data['quasilinear'] = qls
		self.data['_ql_keys'] = ql_keys
		
	def calculate_gr(self, gr_type = "Normalised", doPrint = True):
		GR = full((len(self['psiNs']),self['n_beta'],self['n_shat']),None).tolist()
		MF = full((len(self['psiNs']),self['n_beta'],self['n_shat']),None).tolist()
		KY = full((len(self['psiNs']),self['n_beta'],self['n_shat']),None).tolist()
		SYM = full((len(self['psiNs']),self['n_beta'],self['n_shat']),None).tolist()
		
		if gr_type == "Absolute":
			for psiN in range(len(self['psiNs'])):
				for i in range(self['n_beta']):
					for j in range(self['n_shat']):
						nonan = [x for x in self['growth_rates_all'][psiN][i][j] if (str(x) != 'nan')]
						if len(nonan) != 0:
							idx = array(self['growth_rates_all'][psiN][i][j]).tolist().index(max(nonan))
							GR[psiN][i][j] = self['growth_rates_all'][psiN][i][j][idx]
							MF[psiN][i][j] = self['mode_frequencies_all'][psiN][i][j][idx]
							KY[psiN][i][j] = self['aky_values'][idx]
							SYM[psiN][i][j] = self['parities_all'][psiN][i][j][idx]
						else:
							GR[psiN][i][j] = nan
							MF[psiN][i][j] = nan
							KY[psiN][i][j] = self['aky_values'][0]
							SYM[psiN][i][j] = 0
							
			self.run['data']['growth_rates'] = GR
			self.run['data']['mode_frequencies'] = MF
			self.run['data']['akys'] = KY
			self.run['data']['parities'] = SYM
			self._gr_type = "Absolute"
			if doPrint:
				print("Converted Growth Rates To Absolute")
		elif gr_type == "Normalised":
			for psiN in range(len(self['psiNs'])):
				for i in range(self['n_beta']):
					for j in range(self['n_shat']):
						grns = array(self['growth_rates_all'][psiN][i][j])/array(self['aky_values'])**2
						grns = grns.tolist()
						nonan = [x for x in grns if (str(x) != 'nan')]
						if len(nonan) != 0:
							idx = grns.index(max(nonan))
							GR[psiN][i][j] = grns[idx]
							MF[psiN][i][j] = self['mode_frequencies_all'][psiN][i][j][idx]
							KY[psiN][i][j] = self['aky_values'][idx]
							SYM[psiN][i][j] = self['parities_all'][psiN][i][j][idx]
						else:
							GR[psiN][i][j] = nan
							MF[psiN][i][j] = nan
							KY[psiN][i][j] = self['aky_values'][0]
							SYM[psiN][i][j] = 0
							
			self.run['data']['growth_rates'] = GR
			self.run['data']['mode_frequencies'] = MF
			self.run['data']['akys'] = KY
			self.run['data']['parities'] = SYM
			self._gr_type = "Normalised"
			if doPrint:
				print("Converted Growth Rates To Normalised")
		elif gr_type == "Unnormalised":
			for psiN in range(len(self['psiNs'])):
				for i in range(self['n_beta']):
					for j in range(self['n_shat']):
						grns = array(self['growth_rates_all'][psiN][i][j])/array(self['aky_values'])**2
						grns = grns.tolist()
						nonan = [x for x in grns if (str(x) != 'nan')]
						if len(nonan) != 0:
							idx = grns.index(max(nonan))
							GR[psiN][i][j] = self['growth_rates_all'][psiN][i][j][idx]
							MF[psiN][i][j] = self['mode_frequencies_all'][psiN][i][j][idx]
							KY[psiN][i][j] = self['aky_values'][idx]
							SYM[psiN][i][j] = self['parities_all'][psiN][i][j][idx]
						else:
							GR[psiN][i][j] = nan
							MF[psiN][i][j] = nan
							KY[psiN][i][j] = self['aky_values'][0]
							SYM[psiN][i][j] = 0
						
			self.run['data']['growth_rates'] = GR
			self.run['data']['mode_frequencies'] = MF
			self.run['data']['akys'] = KY
			self.run['data']['parities'] = SYM
			self._gr_type = "Unnormalised"
			if doPrint:
				print("Converted Growth Rates To Unnormalised")
	
	def calculate_alpha(self):
		if not self.eqbm:
			self.load_equillibrium()
		from scipy.interpolate import InterpolatedUnivariateSpline
		qspline = InterpolatedUnivariateSpline(self.eqbm.eq_data['psiN'],self.eqbm.eq_data['qpsi'])
		alpha_axis = []
		alpha_values = []
		alpha_axis_ideal = []
		for p, psiN in enumerate(self['psiNs']):
			alpha_axis.append([abs(x)*self.eqbm.eq_data['rmaxis']*qspline(psiN)**2 for x in self['beta_prime_axis'][p]])
			alpha_values.append(abs(self['beta_prime_values'][p])*self.eqbm.eq_data['rmaxis']*qspline(psiN)**2)
			alpha_axis_ideal.append([abs(x)*self.eqbm.eq_data['rmaxis']*qspline(psiN)**2 for x in self['beta_prime_axis_ideal'][p]])
		self.run['data']['alpha_axis'] = alpha_axis
		self.run['data']['alpha_values'] = alpha_values
		self.run['data']['alpha_axis_ideal'] = alpha_axis_ideal
		
	def plot_aky(self, init = None, settings = {}):
		self.plot_scan(init = init, aky = True, settings = settings)
		
	def plot_scan(self, init = None, aky = None, settings = {}):
		if init is not None and type(init) == int:
			settings['psi_id'] = init
		elif init is not None:
			settings['psi_id'] = int(init[0])
			if len(init) > 1:
				settings['ky_id'] = int(init[1])
		if aky is not None:
			settings['aky'] = aky
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} Scan"
		self.plots['scan'] = Plotters['Scan'](data = self.run['data'], inputs = self.inputs, verify = self.verify, settings = settings)
	
	def plot_ql(self, init = None, settings = {}):
		if self['ql'] is None:
			self.calculate_ql()
		if init is not None:
			settings['psi_id'] = int(init)
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} QuasiLinear"
		self.plots['ql'] = Plotters['QL'](data = self.run['data'], inputs = self.inputs, settings = settings)
	
	def plot_ideal(self, init = 0, settings = {}):
		if init is not None:
			settings['psi_id'] = int(init)
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} Ideal Ballooning"
		self.plots['ideal'] = Plotters['Ideal'](data = self.run['data'], inputs = self.inputs, settings = settings)
	
	def plot_omega(self, init = None, aky = None, settings = {}):
		self.plots['omega'] = self._plot_diag(var = 'omega', init = init, aky = aky, settings = settings)
	
	def plot_phi(self, init = None, aky = None, absolute = None, settings = {}):
		self.plots['phi'] = self._plot_diag(var = 'phi', init = init, aky = aky, absolute = absolute, settings = settings)
	
	def plot_apar(self, init = None, aky = None, absolute = None, settings = {}):
		self.plots['apar'] = self._plot_diag(var = 'apar', init = init, aky = aky, absolute = absolute, settings = settings)
		
	def plot_bpar(self, init = None, aky = None, absolute = None, settings = {}):
		self.plots['bpar'] = self._plot_diag(var = 'bpar', init = init, aky = aky, absolute = absolute, settings = settings)
		
	def plot_phi2(self, init = None, aky = None, settings = {}):
		self.plots['phi2'] = self._plot_diag(var = 'phi2', init = init, aky = aky, settings = settings)
	
	def _plot_diag(self, init = None, aky = None, var = None, absolute = None, settings = {}):
		if init is not None:
			settings['psi_id'] = init[0]
			settings['bp_id'] = init[1]
			settings['sh_id'] = init[2]
			settings['ky_id'] = init[3]
		if aky is not None:
			settings['aky'] = aky
		if var is not None:
			settings['var'] = var
		if 'title' not in settings:
			settings['suptitle'] = f"{self['run_name']} {var}"
		return Plotters['Diag'](data = self.run['data'], inputs = self.inputs, info = self.run['info'], verify = self.verify, settings = settings)
	
	def plot_epar(self):
		Plotters['Epar'](data = self.run['data'], inputs = self.inputs)
	
	def plot_theta(self, init = [0,0,0,0], aky = True, var = 0, n = 3, polar = False):
		Plotters['Theta'](data = self.run['data'], inputs = self.inputs, var = var, init = init, aky = aky, n = n, polar = polar)
		
	def plot_slice(self, init = [0,0], limit = None, settings = {}):
		if self['ql'] is None:
			self.calculate_ql()
		if init is not None:
			settings['psi_id'] = init[0]
			settings['sh_id'] = init[1]
		if limit is not None:
			settings['limit'] = limit
		self.plots['slice'] = Plotters['Slice'](data = self.run['data'], inputs = self.inputs, settings = settings)
	
	def plot_eq(self):
		if not self.eqbm:
			self.load_equillibrium()
		self.eqbm.plot_eq()
	
	def load_equillibrium(self, eq_file = None, kin_file = None, kinetics_type = None, template_file = None, directory = None):
		from .equillibrium import equillibrium
		if directory is None:
			directory = self.directory
		if eq_file is None:
			eq_file = self['eq_file_name']
			if not os.path.exists(f"{os.path.join(directory,eq_file)}"):
				self.write_eq_file(filename = eq_file, directory = directory)
		if kin_file is None:
			kin_file = self['kin_file_name']
			if not os.path.exists(f"{os.path.join(directory,kin_file)}"):
				self.write_kin_file(filename = kin_file, directory = directory)
		if kinetics_type is None:
			kinetics_type = self['kinetics_type']
		if template_file is None:
			template_file = self['template_file_name']
			if not os.path.exists(f"{os.path.join(directory,template_file)}"):
				self.write_template_file(filename = template_file, directory = directory)
		self.eqbm = self.equillibrium = equillibrium(eq_file = eq_file, kin_file = kin_file, kinetics_type = kinetics_type, template_file = template_file, directory = directory, inputs = self.inputs)
	
	def write_gs2_input(self, indexes = None, filename = None, eq_file = None, kin_file = None, template_file = None, directory = None):
		try:
			if len(indexes) != 4:
				print("ERROR: indexes must be of length 4, [psiN,beta_prime,shear,ky]")
				return
		except:
			print("ERROR: must pass list of indexes, [psiN,beta_prime,shear,ky]")
		if directory is None and self.directory is None:
			directory = "./"
		elif directory is None:
			directory = self.directory

		p,i,j,k = indexes
		if filename is None:
			filename = f"{p}_{i}_{j}_{k}.in"
		psiN = self['psiNs'][p]
		bp = -1*abs(self['beta_prime_axis'][p][i])
		sh = self['shear_axis'][p][j]
		aky = self['akyv'][k]
		
		if self.eqbm is None:
			self.load_equillibrium(eq_file = eq_file, kin_file = kin_file, directory = directory)
	
		namelist_diff = self['namelist_diffs'][p][i][j][k]
		nml = self.eqbm.get_gyro_input(psiN = psiN, bp = bp, sh = sh, aky = aky, namelist_diff = namelist_diff)
		
		nml.write(os.path.join(directory,filename), force=True)
		print(f"Created {filename} at {directory}")
		
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
