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
		if opened and not self.verify and verify:		
			self._verify_run()
		
	def __getitem__(self, key):
		if key not in ["Miller", "Ideal", "Gyro", "Viking", "Fixed_delt", "Epar","psiNs","eparN","eparN_all"]:
			key = key.lower()
		if key == "inputs":
			self.print_inputs()
		elif key in ["info","run_info","run info"]:
			self.print_info()
		elif key == "data":
			self.print_data()
			
		elif key in self.run['data'].keys():
			return self.run['data'][key]
		elif key in self.run['info'].keys():
			return self.run['info'][key]
		elif key in self.run['files'].keys():
			return self.run['files'][key]
		
		elif key in ["gr","growth"]:
			return self.run['data']['growth_rates']
		elif key in ["mf","mode"]:
			return self.run['data']['mode_frequencies']
		elif key in ["gra","gr_a","gr_all"]:
			return self.run['data']['growth_rates_all']
		elif key in ["mfa","mf_a","mf_all"]:
			return self.run['data']['mode_frequencies_all']
		elif key in ["namelist_diffs", "nml_diffs", "namelist_diff", "nml_diff", "namelists", "nmls"]:
			return self.run['files']['namelist_differences']
		elif key in ["ql", "quasi linear", "quasi_linear"]:
			return self.run['data']['quasilinear']
			
		elif key in self.verify._all_keys():
			return self.verify[key]
			
		else:
			return self.inputs[key]
			#inputs also contains logic for if key not found
			
	def __len__(self):
		return len(self['psiNs'])*self['n_beta']*self['n_shat']*len(self['aky_values'])
	
	def print_data(self):
        	for key, val in self.run['data'].items():
        		print(key)
        		
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
		data_in = {}
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
			info = data_in['run_info'].item()
			possible_info = ['run_name','run_uuid','data_path','input_file','eq_file_name','template_file_name','kin_file_name',
			'kinetics_type','run_data','_eq_file_path','_kin_file_path','_template_file_path']
			for key in [x for x in possible_info if x not in info.keys()]:
				info[key] = None
		except:
			print("ERROR: could not load Run Info")
			info = None
		try:
			data = data_in['data'].item()
			possible_data = ['beta_prime_values','shear_values','beta_prime_axis','shear_axis','beta_prime_axis_ideal','shear_axis_ideal',
			'growth_rates','mode_frequencies','growth_rates_all','mode_frequencies_all','parities',
			'parities_all','ideal_stabilities','eparN','eparN_all','akys','omega','phi','apar',
			'bpar','phi2','time','theta','quasilinear','alpha_axis','alpha_axis_ideal','alpha_values']
			for key in [x for x in possible_data if x not in data.keys()]:
				data[key] = None
		except:
			print("ERROR: could not load Data")
			data = None
		try:
			files = data_in['files'].item()
			possible_files = ['eq_file','kin_file','template_file','namelist_differences']
			for key in [x for x in possible_files if x not in files.keys()]:
				files[key] = None
		except:
			print("ERROR: could not load Input Files")
			files = None
		try:
			self.inputs = scan_inputs(input_dict = data_in['inputs'].item())
			if info:
				self.inputs.input_file = info['input_file']
		except:
			print("ERROR: could not load Inputs")
			self.inputs = None
		try:
			verify = data_in['verify'].item()
		except:
			verify = {}

		self.run = {'info': info, 'data': data, 'files': files}
		self.verify = verify
		return True
	
	def _save_file(self, filename = None, directory = "./"):
		if filename == None:
			filename = self['run_name']
		savez(f"{directory}/{filename}", inputs = self.inputs.inputs, data = self.run['data'], run_info = self.run['info'], files = self.run['files'], verify = self.verify)
		
	def _remove_akys(self, akys = None, ids = None):
		if akys is None and ids is None:
			print("ERROR: akys or ids must be given")
			return
		if akys:
			ids = []
			for aky in akys:
				ids.append(self['aky_values'].index(aky))
		grs = array(self['growth_rates_all'])
		for idx in ids:
			grs[:,:,:,idx] = nan
		self.run['data']['growth_rates_all'] = grs.tolist()
		self._convert_gr(gr_type = self._gr_type, doPrint = False)
		
	def _verify_run(self):
		if not self['Gyro']:
			return
		self.verify = verify_scan(scan = self.run)
		self.run = self.verify.scan
		self._convert_gr(gr_type = self._gr_type, doPrint = False)
	
	def calculate_ql(self):
		from .quasilinear import QL
		QLs = full((len(self['psiNs']),self['n_beta'],self['n_shat']),None).tolist()
		for p in range(len(self['psiNs'])):
			for i in range(self['n_beta']):
				for j in range(self['n_shat']):
					QLs[p][i][j] = QL((p,i,j),self.run['data'],self['akyv'])
		self.run['data']['quasilinear'] = QLs
		
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
		import f90nml
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
