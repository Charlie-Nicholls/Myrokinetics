import os
from numpy import *
from .plotting import Plotters
from .verify_runs import verify_scan

'''
GYROKINETIC SCAN ANALYSIS
'''

class myro(object):

	def __init__(self, filename = None, directory = "./"):
		if directory == "./":
			directory = os.getcwd() 
		self.directory = directory
		self.filename = filename
		self.run = {}
		self._open_file()
		self._gr_type = "Normalised"
		self.verify = {}		
		self._verify_run()
	
	def __getitem__(self, key):
		if key not in ["Miller", "Ideal", "Gyro", "Viking", "Fixed_delt", "Epar","psiNs","eparN","eparN_all"]:
			key = key.lower()
		if key == "inputs":
			self.inputs()
		elif key in ["info","run_info","run info"]:
			self.run_info()
		elif key == "data":
			self.data()
			
		elif key in self.run['data'].keys():
			return self.run['data'][key]
		elif key in self.run['inputs'].keys():
			return self.run['inputs'][key]
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
		elif key in ["akyv","aky_v",]:
			return self.run['inputs']['aky_values']
		elif key in ['bad_phi2', 'badphi2']:
			return self.verify['phi2']
		elif key in ['bad_nstep', 'badnstep']:
			return self.verify['nstep']
		elif key in ['bad_other', 'badother']:
			return self.verify['other']
		elif key in ['unconv', 'unconverged']:
			return self.verify['unconv']
		elif key in ['saveerrors','save_errors']:
			return self.verify['save_errors']
	@property	
	def data(self):
        	for key, val in self.run['data'].items():
        		print(key)
	@property
	def inputs(self):
        	for key, val in self.run['inputs'].items():
        		print(f"{key} = {val}")
	@property
	def info(self):
		for key, val in self.run['info'].items():
        		print(f"{key} = {val}")
	@property
	def keys(self):
		print("Data Keys:")
		for key in self.run['data'].keys():
			print(key)
		print("\nInput Keys:")
		for key in self.run['inputs'].keys():
			print(key)
		print("\nRun Info Keys:")
		for key in self.run['info'].keys():
			print(key)
	@property
	def help(self):
		#Needs Updating, split into help and _help
		print("Available Commands:\n\nkeys(): List available dictionary keys.\n\ninputs(): List inputs used for run.\n\ninfo(): List run information\n\data(): List available data keys.\n\ntemplate_file(): Print gyrokinetics template file.\n\nplot_scan(): Display S-Alpha plot for growth rate and mode frequency for each flux surface. Uses gyrokinetic data with options to overlay ideal_ball data and parities.\n\nplot_ideal():Display S-Alpha plot for ideal_ball stability only.\n\nplot_omega(): Display the complext frequency as a function of time for specific runs\n\nplot_phi(): Display the electric potential as a function of ballooning angle for specific runs\n\nplot_apar(): Display the parallel magnetic potential as a function of ballooning angle for specific runs.\n\nplot_phi2(): Display the electrostatic potential squared and averaged over theta, kx and ky, as a function of time for specific runs.\n\n_open_file(): Opens and reads data files. Arguments: filename - Name of data file, defaults to stored self.filename; directory - Location of data file, defaults to stored self.directory.")
		
	def _print_file(self, filetype = ''):
		if self.run['files'] is None:
			print("ERROR: No files found")
			return
		filetype = filetype.lower()
		if filetype == '':
			print("ERROR: filetype not given")
			return
		if filetype in ['template_file', 'template', 'gk_template']:
			lines = self.run['files']['template_file']
		elif filetype in ['eq','eq_file','geqdsk','geqdsk_file']:
			lines = self.run['files']['eq_file']
		elif filetype in ['kin','kin_file','kinetics','kinetics_file','peqdsk','peqdsk_file']:
			lines = self.run['files']['kin_file']
		else:
			print(f"ERROR: File Type {filetype} Not Found")
			return
		for line in lines:
			print(line,end='')
	@property	
	def template_file(self):
		self._print_file(filetype = 'template')
	@property
	def eq_file(self):
		self._print_file(filetype = 'eq')
	@property
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
			lines = self.run['files']['template_file']
			name = self.run['info']['template_file_name']
		elif filetype in ['eq','eq_file','geqdsk','geqdsk_file']:
			lines = self.run['files']['eq_file']
			name = self.run['info']['eq_file_name']
		elif filetype in ['kin','kin_file','kinetics','kinetics_file','peqdsk','peqdsk_file']:
			lines = self.run['files']['kin_file']
			name = self.run['info']['kin_file_name']
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
			run_info = None
		try:
			data = data_in['data'].item()
		except:
			print("ERROR: could not load Data")
			data = None
		try:
			files = data_in['files'].item()
		except:
			print("ERROR: could not load Input Files")
			files = None
		self.run = {'inputs': inputs, 'info': info, 'data': data, 'files': files}
	
	def _save_file(self, directory = "./", filename = None):
		if filename == None:
			filename = self.run['info']['run_name']
		savez(f"{directory}/{filename}", inputs = self.run['inputs'], data = self.run['data'], run_info = self.run['info'], files = self.run['files'])
		
	def _remove_akys(self, akys = None, ids = None):
		if akys is None and ids is None:
			print("ERROR: akys or ids must be given")
			return
		if akys:
			ids = []
			for aky in akys:
				ids.append(self.run['inputs']['aky_values'].index(aky))
		else:
			grs = array(self.run['data']['growth_rates_all'])
			for idx in ids:
				grs[:,:,:,idx] = nan
			self.run['data']['growth_rates_all'] = grs.tolist()
		self._convert_gr(gr_type = self._gr_type, doPrint = False)
		
	def _verify_run(self):
		self.verify = verify_scan(scan = self.run)
		self.run['data']['growth_rates_all'] = self.verify.new_data['gra']
		self.run['data']['mode_frequncies_all'] = self.verify.new_data['mfa']
		self._convert_gr(gr_type = self._gr_type, doPrint = False)
	
	def _reset_data(self):
		self.run['data'] = self.verify.scan['data']
		self._convert_gr(gr_type = "Normalised", doPrint = False)
	
	def _reset_all(self):
		self.run = self.verify.scan
		self._convert_gr(gr_type = "Normalised", doPrint = False)
		
	def _convert_gr(self, gr_type = None, doPrint = True):
		GR = full((len(self.run['inputs']['psiNs']),self.run['inputs']['n_beta'],self.run['inputs']['n_shat']),None).tolist()
		MF = full((len(self.run['inputs']['psiNs']),self.run['inputs']['n_beta'],self.run['inputs']['n_shat']),None).tolist()
		KY = full((len(self.run['inputs']['psiNs']),self.run['inputs']['n_beta'],self.run['inputs']['n_shat']),None).tolist()
		SYM = full((len(self.run['inputs']['psiNs']),self.run['inputs']['n_beta'],self.run['inputs']['n_shat']),None).tolist()
		
		if gr_type is None:
			if self._gr_type == "Unnormalised":
				gr_type = "Absolute"
			elif self._gr_type == "Absolute":
				gr_type = "Normalised"
			elif self._gr_type == "Normalised":
				gr_type = "Unnormalised"

		if gr_type == "Absolute":
			for psiN in range(len(self.run['inputs']['psiNs'])):
				for i in range(self.run['inputs']['n_beta']):
					for j in range(self.run['inputs']['n_shat']):
						idx = array(self.run['data']['growth_rates_all'][psiN][i][j]).tolist().index(max([x for x in self.run['data']['growth_rates_all'][psiN][i][j] if (str(x) != 'nan')]))
						GR[psiN][i][j] = self.run['data']['growth_rates_all'][psiN][i][j][idx]
						MF[psiN][i][j] = self.run['data']['mode_frequencies_all'][psiN][i][j][idx]
						KY[psiN][i][j] = self.run['inputs']['aky_values'][idx]
						SYM[psiN][i][j] = self.run['data']['parities_all'][psiN][i][j][idx]
			self.run['data']['growth_rates'] = GR
			self.run['data']['mode_frequencies'] = MF
			self.run['data']['akys'] = KY
			self.run['data']['parities'] = SYM
			self._gr_type = "Absolute"
			if doPrint:
				print("Converted Growth Rates To Absolute")
		elif gr_type == "Normalised":
			for psiN in range(len(self.run['inputs']['psiNs'])):
				for i in range(self.run['inputs']['n_beta']):
					for j in range(self.run['inputs']['n_shat']):
						grns = array(self.run['data']['growth_rates_all'][psiN][i][j])/array(self.run['inputs']['aky_values'])**2
						grns = grns.tolist()
						idx = grns.index(max([x for x in grns if (str(x) != 'nan')]))
						GR[psiN][i][j] = grns[idx]
						MF[psiN][i][j] = self.run['data']['mode_frequencies_all'][psiN][i][j][idx]
						KY[psiN][i][j] = self.run['inputs']['aky_values'][idx]
						SYM[psiN][i][j] = self.run['data']['parities_all'][psiN][i][j][idx]
			self.run['data']['growth_rates'] = GR
			self.run['data']['mode_frequencies'] = MF
			self.run['data']['akys'] = KY
			self.run['data']['parities'] = SYM
			self._gr_type = "Normalised"
			if doPrint:
				print("Converted Growth Rates To Normalised")
		elif gr_type == "Unnormalised":
			for psiN in range(len(self.run['inputs']['psiNs'])):
				for i in range(self.run['inputs']['n_beta']):
					for j in range(self.run['inputs']['n_shat']):
						grns = array(self.run['data']['growth_rates_all'][psiN][i][j])/array(self.run['inputs']['aky_values'])**2
						grns = grns.tolist()
						idx = grns.index(max([x for x in grns if (str(x) != 'nan')]))
						GR[psiN][i][j] = self.run['data']['growth_rates_all'][psiN][i][j][idx]
						MF[psiN][i][j] = self.run['data']['mode_frequencies_all'][psiN][i][j][idx]
						KY[psiN][i][j] = self.run['inputs']['aky_values'][idx]
						SYM[psiN][i][j] = self.run['data']['parities_all'][psiN][i][j][idx]
			self.run['data']['growth_rates'] = GR
			self.run['data']['mode_frequencies'] = MF
			self.run['data']['akys'] = KY
			self.run['data']['parities'] = SYM
			self._gr_type = "Unnormalised"
			if doPrint:
				print("Converted Growth Rates To Unnormalised")
	
	def plot_aky(self, init = [0,0]):
		Plotters['Scan'](scan = self.run, aky=True, init = init)
	
	def plot_scan(self, aky = False, init = [0,0,0,0]):
		Plotters['Scan'](scan = self.run, aky = aky, init = init)
	
	def plot_ideal(self):
		Plotters['Ideal'](scan = self.run)
	
	def plot_omega(self, aky = True, init = [0,0,0,0]):
		Plotters['Diag'](scan = self.run, var = 0, aky = aky, init = init, verify = self.verify)
	
	def plot_phi(self, aky = True, init = [0,0,0,0]):
		Plotters['Diag'](scan = self.run, var = 1, aky = aky, init = init, verify = self.verify)
	
	def plot_apar(self, aky = True, init = [0,0,0,0]):
		Plotters['Diag'](scan = self.run, var = 2, aky = aky, init = init, verify = self.verify)
		
	def plot_phi2(self, aky = True, init = [0,0,0,0]):
		Plotters['Diag'](scan = self.run, var = 3, aky = aky, init = init, verify = self.verify)
	
	def _plot_diag(self, var = 0, aky = True, init = [0,0,0,0]):
		Plotters['Diag'](scan = self.run, var = var, aky = aky, init = init, verify = self.verify)
	
	def plot_epar(self):
		Plotters['Epar'](scan = self.run)
	
	def write_gs2_input(self, indexes = None, eq_file = None, kin_file = None, template_file = None, filename = None, directory = None):
		from pathlib import Path
		from pyrokinetics import Pyro
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
		if directory == "./":
			directory = os.getcwd() 
		p,i,j,k = indexes
		if filename is None:
			filename = f"{indexes[0]}_{indexes[1]}_{indexes[2]}_{indexes[3]}.in"
		
		if eq_file is None:
			eq_file = f"{self.run['info']['run_name']}_eq"
			self.write_eq_file(filename = eq_file, directory = directory)
		if kin_file is None:
			kin_file = f"{self.run['info']['run_name']}_kin"
			self.write_kin_file(filename = kin_file, directory = directory)
		if template_file is None and self.run['files']['template_file'] is not None:	
			template_file = f"{self.run['info']['run_name']}_template"
			self.write_template_file(filename = template_file, directory = directory)

		if template_file is None:
			pyro = Pyro(
				eq_file=Path(directory) / eq_file,
			 	eq_type="GEQDSK",
			 	kinetics_file=Path(directory) / kin_file,
			 	kinetics_type=self.run['info']['kinetics_type'])
		else:
			pyro = Pyro(
				eq_file=Path(directory) / eq_file,
			 	eq_type="GEQDSK",
			 	kinetics_file=Path(directory) / kin_file,
			 	kinetics_type=self.run['info']['kinetics_type'],
			 	gk_file=Path(directory) / template_file)
	
		pyro.gk_code = "GS2"
		pyro.load_local_geometry(self.run['inputs']['psiNs'][p])
		pyro.load_local_species(self.run['inputs']['psiNs'][p])
		pyro.write_gk_file(os.path.join(directory,filename))
		nml = pyro._gk_input_record["GS2"].data
			
		if self.run['inputs']['Miller']:
			nml['theta_grid_eik_knobs']['iflux'] = 0
			nml['theta_grid_eik_knobs']['local_eq'] = True
		else:
			eq_dir = os.path.join(directory,eq_file)
			nml['theta_grid_eik_knobs']['eqfile'] = eq_dir
			nml['theta_grid_eik_knobs']['efit_eq'] =  True
			nml['theta_grid_eik_knobs']['iflux'] = 1
			nml['theta_grid_eik_knobs']['local_eq'] = False
		
		if self.run['inputs']['Epar']:
			nml['gs2_diagnostics_knobs']['write_final_epar'] = True
		else:
			nml['gs2_diagnostics_knobs']['write_final_epar'] = False
		
		bp = -self.run['data']['beta_prime_axis'][p][i]
		sh = self.run['data']['shear_axis'][p][j]
		
		nml['theta_grid_eik_knobs']['s_hat_input'] = sh
		nml['theta_grid_eik_knobs']['beta_prime_input'] = bp
		nml['kt_grids_single_parameters']['aky'] = self.run['inputs']['aky_values'][k]
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			mul = bp/(-2*(nml[spec]['tprim'] + nml[spec]['fprim'])*nml['parameters']['beta'])
			nml[spec]['tprim'] = nml[spec]['tprim']*mul
			nml[spec]['fprim'] = nml[spec]['fprim']*mul
		if self.run['inputs']['Fixed_delt'] is False:
			nml['knobs']['delt'] = 0.04/self.run['inputs']['aky_values'][k]
		
		nml.write(os.path.join(directory,filename), force=True)
		print(f"Created {filename} at {directory}")
