import os
import f90nml
from .inputs import scan_inputs
from .templates import template_dir, gs2_template, cgyro_template, tglf_template
from copy import deepcopy

class equilibrium(object):
	def __init__(self, inputs = None, directory = None):
		if directory == "./":
			directory = os.getcwd()
		self.path = directory
		if inputs is None:
			inputs = scan_inputs()
		self.load_inputs(inputs)
		self.eq_data = self.kin_data = self.pyro = self._eq_lines = self._kin_lines = self.beta_prime_profile = self.shear_profile = None
		self.surface_namelists = {}

	def load_geqdsk(self, eq_file = None, directory = None):
		from .geqdsk_reader import geqdsk
		if not self.inputs['eq_name'] and eq_file is None:
			print("ERROR: No GEQDSK file given")
			return
		elif eq_file is not None:
			self.inputs.inputs['files']['eq_name'] = eq_file
			self.inputs.check_inputs()
		
		if directory is None:
			if self.inputs['eq_path'] is None and self.path is None:
				self.inputs.inputs['files']['eq_path'] = os.getcwd()
			elif self.inputs['eq_path'] is None:
				self.inputs.inputs['files']['eq_path'] = self.path
			if self.inputs['eq_path'] == "./":
				self.inputs.inputs['files']['eq_path'] = os.getcwd()
			directory = self.inputs.inputs['files']['eq_path']
		if directory == "./":
			directory = os.getcwd()
		self.inputs.inputs['files']['eq_path'] = directory

		if self.surface_namelists:
			self.surface_namelists = {}
		if self.pyro:
			self.pyro = None
		
		with open(os.path.join(self.inputs['eq_path'],self.inputs['eq_name'])) as efile:
			self._eq_lines = efile.readlines()
			
		self.eq_data = geqdsk(filename = self.inputs['eq_name'], directory = self.inputs['eq_path'])
	
	def load_kinetics(self, kin_file = None, kinetics_type = None, directory = None):
		if self.inputs['kin_name'] is None and kin_file is None:
			print("ERROR: No Kinetics file given")
			return
		elif kin_file is not None:
			self.inputs.inputs['files']['kin_name'] = kin_file
			self.inputs.check_inputs()
		
		if directory is None:
			if self.inputs['kin_path'] is None and self.path is None:
				self.inputs.inputs['files']['kin_path'] = os.getcwd()
			elif self.inputs['kin_path'] is None:
				self.inputs.inputs['files']['kin_path'] = self.path
			if self.inputs['kin_path'] == "./":
				self.inputs.inputs['files']['kin_path'] = os.getcwd()
			directory = self.inputs.inputs['files']['kin_path']
		if directory == "./":
			directory = os.getcwd()
		self.inputs.inputs['files']['kin_path'] = directory

		
		if self.surface_namelists:
			self.surface_namelists = {}
		if self.pyro:
			self.pyro = None
		
		with open(os.path.join(self.inputs['kin_path'],self.inputs['kin_name'])) as kfile:
			self._kin_lines = kfile.readlines()
			
		if kinetics_type is not None:
			self.inputs.inputs['kinetics_type'] = kinetics_type
		if self.inputs['kinetics_type'] is None:
			print("Warning: kinetics_type Not Given, trying PEQDSK")
			self.inputs.inputs['kinetics_type'] = "PEQDSK"
		
		if self.inputs['kinetics_type'].upper() == "SCENE":
			import xarray as xr
			self.kin_data = xr.open_dataset(os.path.join(self.inputs['kin_path'],self.inputs['kin_name']))
		elif self.inputs['kinetics_type'].upper() in ["PEQDSK","PFILE"]:
			from .peqdsk_reader import peqdsk
			self.kin_data = peqdsk(filename = self.inputs['kin_name'], directory = self.inputs['kin_path'])
		else:
			print(f"ERROR: Kinetics type {self.inputs['kinetics_type']} not recognised. Currently supported: SCENE, PEQDSK/pFile")
	
	def load_pyro(self, template_file = None, directory = None):
		from pyrokinetics import Pyro
		from pathlib import Path
		if self.inputs['eq_name'] is None or self.inputs['kin_name'] is None:
			if self.inputs['eq_name'] is None:
				print("ERROR: No equilibrium loaded")
			if self.inputs['kin_name'] is None:
				print("ERROR: No kinetics file loaded")
			return
		
		if self.surface_namelists:
			self.surface_namelists = {}
		
		if template_file is not None:
			self.inputs.inputs['files']['template_name'] = template_file
			self.inputs.check_inputs()
		elif self.inputs['template_name'] is None:
			if self.inputs['gk_code'] == 'GS2':
				self.inputs.inputs['files']['template_name'] = gs2_template
			elif self.inputs['gk_code'] == 'CGYRO':
				self.inputs.inputs['files']['template_name'] = cgyro_template
			elif self.inputs['gk_code'] == 'TGLF':
				self.inputs.inputs['files']['template_name'] = tglf_template
			self.inputs.inputs['files']['template_path'] = template_dir
		
		if directory is None:
			if self.inputs['template_path'] is None and self.path is None:
				self.inputs.inputs['files']['template_path'] = os.getcwd()
			elif self.inputs['template_path'] is None:
				self.inputs.inputs['files']['template_path'] = self.path
			if self.inputs['template_path'] == "./":
				self.inputs.inputs['files']['template_path'] = os.getcwd()
			directory = self.inputs.inputs['files']['template_path']
		if directory == "./":
			directory = os.getcwd()
		self.inputs.inputs['files']['template_path'] = directory
		
		if not self.eq_data:
			self.load_geqdsk()
		if not self.kin_data:
			self.load_kinetics()
		
		if self.inputs['gk_code'] == 'GS2':
			self._template_lines = f90nml.read(os.path.join(self.inputs['template_path'],self.inputs['template_name']))
		elif self.inputs['gk_code'] in ['CGYRO','TGLF']:
			with open(os.path.join(self.inputs['template_path'],self.inputs['template_name']),'r') as f:
				self._template_lines = f.readlines()
		

		kin_type = 'pFile' if self.inputs['kinetics_type'].upper() == 'PEQDSK' else self.inputs['kinetics_type'].upper()
		self.pyro = Pyro(
			eq_file = Path(self.inputs['eq_path']) / self.inputs['eq_name'],
		 	eq_type = "GEQDSK",
		 	kinetics_file = Path(self.inputs['kin_path']) / self.inputs['kin_name'],
		 	kinetics_type = kin_type,
		 	gk_file = Path(self.inputs['template_path']) / self.inputs['template_name'],
			gk_code = self.inputs['gk_code']
		 	)
		#PYRO DOES NOT SEEM TO LOAD GK FILE PROPERLY, AT LEAST FOR CGYRO, UNSURE WHY
		
		
	def load_inputs(self, inputs):
		self.surface_namelists = {}
		if type(inputs) == str:
			self.inputs = inputs(input_file = inputs, directory = self.path)
		elif type(inputs) == scan_inputs:
			self.inputs = inputs
		else:
			print("ERROR: inputs must be of type string or scan_inputs")

	def AmendPEQDSK(self):
		if self.eq_data is None and self.inputs['eq_name'] is None:
			print("ERROR: No GEQDSK file provided")
			return
		elif self.eq_data is None:
			self.load_geqdsk()
		if self.kin_data is None and self.inputs['kin_name'] is None:
			print("ERROR: No Kinetics file provided")
			return
		elif self.kin_data is None:
			self.load_kinetics()
		
		if 'rhonorm' in self.kin_data.keys():
			return
		
		psi_n = self.kin_data['psinorm']
		rho = []
		for psiN in psi_n:
			if psiN < 1.0e-3:
				psi_n = psi_n[psi_n != psiN]	
			else:
				fs = self.eq_data.flux_surface(psiN = psiN)
				rho.append((max(fs['R']) - min(fs['R']))/2)
		rhonorm = rho/max(rho)
		
		f = open(self.inputs['kin_name'],'a')
		f.write(f"{len(rho)+1} psinorm rho rhonorm")
		f.write("\n 0.0000000   0.0000000   0.0000000")
		for i in range(len(rho)):
			f.write(f"\n {psi_n[i]:.7f}   {rho[i]:.7f}   {rhonorm[i]:.7f}")  
		f.close()
	
	def get_surface_input(self, psiN):
		if psiN in self.surface_namelists.keys():
			return deepcopy(self.surface_namelists[psiN])
		if self.pyro is None:
			self.load_pyro()

		if self.inputs['gk_code'] == 'GS2':
			return self._get_surface_input_gs2(psiN)
		elif self.inputs['gk_code'] == 'CGYRO':
			return self._get_surface_input_cgyro(psiN)
		elif self.inputs['gk_code'] == 'TGLF':
			return self._get_surface_input_tglf(psiN)
		else:
			print("ERROR: gk_code not found")
			return None

	def _get_surface_input_gs2(self, psiN):
		self.pyro.load_local(psi_n=psiN)
		self.pyro.update_gk_code()
		nml = deepcopy(self.pyro.gk_input.data)

		nml['theta_grid_parameters']['qinp'] = abs(nml['theta_grid_parameters']['qinp'])
		
		if 'parameters' in nml.keys() and 'beta' in nml['parameters'].keys():
			nml['knobs']['beta'] = nml['parameters']['beta']
			del(nml['parameters'])
		
		for dim in self.inputs.single_parameters.values():
			nml = dim.single_edit_nml(nml)
		beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
		
		beta = nml['knobs']['beta']
		bp_cal = sum([(nml[spec]['tprim'] + nml[spec]['fprim'])*nml[spec]['dens']*nml[spec]['temp'] for spec in [x for x in nml.keys() if 'species_parameters_' in x]])*beta*-1
		mul = beta_prim/bp_cal
		
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			nml[spec]['tprim'] = nml[spec]['tprim']*mul
			nml[spec]['fprim'] = nml[spec]['fprim']*mul
			
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			#Set bakdif to 0 for Electormagnetic Runs as a default
			nml[spec]['bakdif'] = 0
			#Force Uprim to 0, reccommendation by David
			nml[spec]['uprim'] = 0
			#TEMPORARY UNTIL PYRO FIXES BUG
			if 'bakdif' in nml[spec].keys():
				del(nml[spec]['bakdif'])
		#TEMPORARY UNTIL PYRO FIXES BUG
		if 'normalisations_knobs' in nml.keys():
			if 'qref' in nml['normalisations_knobs'].keys():
				del(nml['normalisations_knobs']['qref'])

		if 'ngauss' in nml['le_grids_knobs'] and 'npassing' in nml['le_grids_knobs']:
			del(nml['le_grids_knobs']['ngauss'])
		
		if self.inputs['force_zero_fs']:
			nml['dist_fn_knobs']['g_exb'] = 0
			
		nml['theta_grid_eik_knobs']['equal_arc'] = False
		nml['init_g_knobs']['ginit_option'] = 'random_sine'
		nml['dist_fn_knobs']['mach'] = 0
		
		if 'avail_cpu_time' not in nml['knobs'].keys():
				h, m, s = self.inputs['sbatch']['time'].split(':')
				nml['knobs']['avail_cpu_time'] = (int(h) * 3600) + (int(m) * 60) + int(s)
		
		if self.inputs['grid_option'] == 'single':
			nml['kt_grids_knobs']['grid_option'] = 'single'
			if 'kt_grids_single_parameters' not in nml:
				nml['kt_grids_single_parameters'] = {'aky': 0.1, 'theta0': 0}
			if 'kt_grids_box_parameters' in nml:
				del(nml['kt_grids_box_parameters'])
			if 'margin_cpu_time' not in nml['knobs'].keys():
				nml['knobs']['margin_cpu_time'] = 300
		elif self.inputs['grid_option'] == 'box':
			nml['kt_grids_knobs']['grid_option'] = 'box'
			nml['dist_fn_knobs']['boundary_option'] = 'linked'
			nml['dist_fn_knobs']['esv'] = True
			nml['fields_knobs']['field_option'] = 'local'
			if 'kt_grids_box_parameters' not in nml:
				nml['kt_grids_box_parameters'] = {'nx': 50, 'ny': 50, 'y0': -0.05, 'jtwist': 1}
			if 'kt_grids_single_parameters' in nml:
				del(nml['kt_grids_single_parameters'])
			nml['fields_knobs']['dump_response'] = True
			nml['fields_knobs']['response_dir'] = "response"
			nml['init_g_knobs']['restart_dir'] = "restart"
			nml['gs2_diagnostics_knobs']['save_for_restart'] = True
			nml['gs2_diagnostics_knobs']['nc_sync_freq'] = 1
			if nml['gs2_diagnostics_knobs']['nsave'] > 1000:
				nml['gs2_diagnostics_knobs']['nsave'] = 1
			if 'margin_cpu_time' not in nml['knobs'].keys():
				nml['knobs']['margin_cpu_time'] = 2400
			if 'nperiod' not in self.inputs.dimensions and 'nperiod' not in self.inputs.single_parameters:
					nml['theta_grid_parameters']['nperiod'] = 1
		else:
			print("ERROR: grid_option is invalid, valid: ['single','box']")
					
		if self.inputs['non_linear'] == True:
			if 'nonlinear_terms_knobs' not in nml.keys():
				nml['nonlinear_terms_knobs'] = {}
			nml['nonlinear_terms_knobs']['nonlinear_mode'] = 'on'
			if 'cfl' not in nml['nonlinear_terms_knobs'].keys():
				nml['nonlinear_terms_knobs']['cfl']	= 0.5
			if self.inputs['split_nonlinear'] == True:
				nml['nonlinear_terms_knobs']['split_nonlinear'] = True
				if 'split_nonlinear_terms_knobs' not in nml.keys():
					nml['split_nonlinear_terms_knobs'] = {'show_statistics': True}
			
		if self.inputs['Miller']:
			nml['theta_grid_eik_knobs']['local_eq'] = True
		else:
			nml['theta_grid_eik_knobs']['eqfile'] = os.path.join(self.inputs['data_path'],self.inputs['eq_name'])
			nml['theta_grid_eik_knobs']['efit_eq'] =  True
			nml['theta_grid_eik_knobs']['local_eq'] = False
		
		if self.inputs['Epar']:
			nml['gs2_diagnostics_knobs']['write_ascii'] = True
			nml['gs2_diagnostics_knobs']['write_final_epar'] = True
		else:
			nml['gs2_diagnostics_knobs']['write_ascii'] = False
			nml['gs2_diagnostics_knobs']['write_final_epar'] = False
		
		if 'ntheta_geometry' not in nml['theta_grid_eik_knobs'].keys():
			nml['theta_grid_eik_knobs']['ntheta_geometry'] = 4096

		if self.inputs['Ideal']:
			beta_mul = abs(self.inputs['beta_prime_max']/beta_prim)
			beta_div = abs(beta_prim/self.inputs['beta_prime_min'])

			nml['ballstab_knobs'] = {'n_shat': self.inputs['n_shat_ideal'], 'n_beta': self.inputs['n_beta_ideal'], 'shat_min': self.inputs['shat_min'], 'shat_max': self.inputs['shat_max'], 'beta_div': beta_div, 'beta_mul': beta_mul}

		nml['knobs']['wstar_units'] = False
		
		for dim_name, dim in self.inputs.single_parameters.items():
			nml = dim.single_edit_nml(nml)

		self.surface_namelists[psiN] = nml
		
		return deepcopy(self.surface_namelists[psiN])

	def _get_surface_input_cgyro(self, psiN):
		self.pyro.load_local(psi_n=psiN)
		self.pyro.update_gk_code()
		nml = deepcopy(self.pyro.gk_input.data)

		for dim in self.inputs.single_parameters.values():
			nml = dim.single_edit_nml(nml)
					
		if self.inputs['non_linear'] == True:
			nml['NONLINEAR_FLAG'] = 1
		else:
			nml['NONLINEAR_FLAG'] = 0
		
		nml['DELTA_T_METHOD'] = 1
		nml['EQUILIBRIUM_MODEL'] = 2
		
		nml['THETA_PLOT'] = nml['N_THETA'] #TEMPORARY UNTIL I FIND A BETTER SOLUTION TO ENSURING ALWAYS A FACTOR OF NTHETA AND/OR THE TEMPLATE LOADING ISSUE WITH PYRO
		
		#for dim in self.inputs.single_parameters.values():
		#	nml = dim.single_edit_nml(nml)

		self.surface_namelists[psiN] = nml
		
		return deepcopy(self.surface_namelists[psiN])
	
	def _get_surface_input_tglf(self, psiN):
		self.pyro.load_local(psi_n=psiN)
		self.pyro.update_gk_code()
		nml = deepcopy(self.pyro.gk_input.data)

		for dim in self.inputs.single_parameters.values():
			nml = dim.single_edit_nml(nml)
					
		self.surface_namelists[psiN] = nml
		
		return deepcopy(self.surface_namelists[psiN])

	def get_gyro_input(self, run = None, indexes = None, namelist_diff = {}):
		if self.inputs['gk_code'] == 'GS2':
			return self._get_gyro_input_gs2(run = run, indexes = indexes, namelist_diff = namelist_diff)
		if self.inputs['gk_code'] == 'CGYRO':
			return self._get_gyro_input_cgyro(run = run, indexes = indexes, namelist_diff = namelist_diff)
		if self.inputs['gk_code'] == 'TGLF':
			return self._get_gyro_input_tglf(run = run, indexes = indexes, namelist_diff = namelist_diff)
		else:
			print("ERROR: gk_code not found")
			return None

	def _get_gyro_input_gs2(self, run = None, indexes = None, namelist_diff = {}):
		if run is None and indexes is None:
			print("ERROR: Either indexes or run must be given")
			return None
		
		if run is None:
			if len(indexes) != len(self.inputs.dimensions):
				print(f"ERROR: indexes must be of length {len(self.dimensions)}, {[self.inputs.dim_order]}")
				return None
			run = {}
			for i, dim in zip(indexes,self.inputs.dimensions.values()):
				run[dim.name] = dim.values[i]
		
		if 'psin' in run:	
			psiN = run['psin']
		elif 'psin' in self.inputs.single_parameters:
			psiN = self.inputs.single_parameters['psin'].values[0]
		else:
			print("ERROR: psiN not defined")
			return None
			
		nml = self.get_surface_input(psiN)
		
		if self.inputs['knobs']['fixed_delt'] == False:
			ky = nml['kt_grids_single_parameters']['aky']
			delt = 0.04/ky
			if delt > 0.01:
				delt = 0.01
			nml['knobs']['delt'] = delt
		
		for dim_name, dim in self.inputs.single_parameters.items():
			nml = dim.single_edit_nml(nml)
		
		for dim_name, dim in self.inputs.dimensions.items():
			nml = dim.edit_nml(nml=nml,val=run[dim_name])
			
		for key in namelist_diff.keys():
			for skey in namelist_diff[key].keys():
				nml[key][skey] = namelist_diff[key][skey]
			
		return nml

	def _get_gyro_input_cgyro(self, run = None, indexes = None, namelist_diff = {}):
		if run is None and indexes is None:
			print("ERROR: Either indexes or run must be given")
			return None
		
		if run is None:
			if len(indexes) != len(self.inputs.dimensions):
				print(f"ERROR: indexes must be of length {len(self.dimensions)}, {[self.inputs.dim_order]}")
				return None
			run = {}
			for i, dim in zip(indexes,self.inputs.dimensions.values()):
				run[dim.name] = dim.values[i]
		
		if 'psin' in run:	
			psiN = run['psin']
		elif 'psin' in self.inputs.single_parameters:
			psiN = self.inputs.single_parameters['psin'].values[0]
		else:
			print("ERROR: psiN not defined")
			return None
			
		nml = self.get_surface_input(psiN)
		
		if self.inputs['knobs']['fixed_delt'] == False:
			ky = nml['kt_grids_single_parameters']['aky']
			delt = 0.04/ky
			if delt > 0.01:
				delt = 0.01
			nml['knobs']['delt'] = delt
		
		for dim_name, dim in self.inputs.dimensions.items():
			nml = dim.edit_nml(nml=nml,val=run[dim_name])
			
		for dim_name, dim in self.inputs.single_parameters.items():
			nml = dim.single_edit_nml(nml)
		
		for key in namelist_diff.keys():
			for skey in namelist_diff[key].keys():
				nml[key][skey] = namelist_diff[key][skey]
			
		return nml
	
	def _get_gyro_input_tglf(self, run = None, indexes = None, namelist_diff = {}):
		if run is None and indexes is None:
			print("ERROR: Either indexes or run must be given")
			return None
		
		if run is None:
			if len(indexes) != len(self.inputs.dimensions):
				print(f"ERROR: indexes must be of length {len(self.dimensions)}, {[self.inputs.dim_order]}")
				return None
			run = {}
			for i, dim in zip(indexes,self.inputs.dimensions.values()):
				run[dim.name] = dim.values[i]
		
		if 'psin' in run:	
			psiN = run['psin']
		elif 'psin' in self.inputs.single_parameters:
			psiN = self.inputs.single_parameters['psin'].values[0]
		else:
			print("ERROR: psiN not defined")
			return None
			
		nml = self.get_surface_input(psiN)
		
		return nml
	
	def write_nml(self, nml, directory = ".", filename = None):
		if self.inputs['gk_code'] == 'GS2':
			if filename is None:
				filename = "itteration_0.in"
			nml.write(f"{directory}/{filename}", force=True)
		elif self.inputs['gk_code'] == 'CGYRO':
			with open(f"{directory}/input.cgyro", "w") as f:
				for key, value in nml.items():
					f.write( f"{key} = {value}\n")
		elif self.inputs['gk_code'] == 'TGLF':
			with open(f"{directory}/input.tglf", "w") as f:
				for key, value in nml.items():
					f.write( f"{key} = {value}\n")		
		else:
			print(f"ERROR: gk_code {self.inputs['gk_code']} not found")

	def make_profiles(self):
		from scipy.interpolate import InterpolatedUnivariateSpline
		from numpy import linspace
		if not self.pyro:
			self.load_pyro()
		pyro = self.pyro
		psiNs = linspace(0.01,1,100)
		bp = []
		sh = []
		for psiN in psiNs:
			pyro.load_local_geometry(psi_n=psiN)
			bp.append(pyro.local_geometry['beta_prime'])
			sh.append(pyro.local_geometry['shat'])
		self.beta_prime_profile = InterpolatedUnivariateSpline(psiNs,bp)
		self.shear_profile = InterpolatedUnivariateSpline(psiNs,sh)
	
	def plot_eq(self):
		from matplotlib.pyplot import subplots, show, ion
		from numpy import linspace
		if not self.beta_prime_profile or not self.shear_profile:
			self.make_profiles()

		fig, ax = subplots(2,1)
		psiNs = linspace(0.01,1,100)
		bp = self.beta_prime_profile(psiNs)
		sh = self.shear_profile(psiNs)
		
		ion()
		ax[0].plot(psiNs, bp, 'b')
		ax[0].invert_yaxis()
		ax[1].plot(psiNs, sh, 'b')
		ax[0].set_xlabel("psiN")
		ax[1].set_xlabel("psiN")
		ax[0].set_ylabel("beta_prime")
		ax[1].set_ylabel("shear")
		show()
	
	def plot_kin(self):
		from matplotlib.pyplot import subplots, show, ion

		fig, ax = subplots(3,1)
		psiNs = self.kin_data['psinorm']
		ne = self.kin_data['ne']
		te = self.kin_data['te']
		ni = self.kin_data['ni']
		ti = self.kin_data['ti']
		
		ion()
		ax[2].plot(psiNs, ne, 'b')
		ax[1].plot(psiNs, te, 'b')
		ax[0].plot(psiNs, ne*te, 'b', label = "electron")
		ax[2].plot(psiNs, ni, 'r')
		ax[1].plot(psiNs, ti, 'r')
		ax[0].plot(psiNs, ni*ti, 'r',label = "ion")
		ax[2].set_xlabel("psiN")
		ax[1].set_xlabel("psiN")
		ax[0].set_xlabel("psiN")
		ax[2].set_ylabel("Density ($10^{20}m^{-3}$)")
		ax[1].set_ylabel("Temeperature (keV)")
		ax[0].set_ylabel("Pressure ($10^{20}keV m^{-3}$)")
		ax[0].legend()
		show()
