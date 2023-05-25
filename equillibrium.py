import os
import f90nml
from .templates import template_dir, gs2_template
from copy import deepcopy

class equillibrium(object):
	
	def __init__(self, eq_file = None, kin_file = None, kinetics_type = None, template_file = None, directory = None, inputs = {}):
		self.eq_name = eq_file
		self.kin_name = kin_file
		self.template_name = template_file
		self.kinetics_type = kinetics_type
		self.path = directory
		self.eq_data = self.kin_data = self.pyro = self._eq_path = self._kin_path = self._eq_lines = self._kin_lines = self.beta_prime_profile = self.shear_profile = None
		self.surface_namelists = {}
		self.inputs = {}
		self._inputs = {
		'shat_min': 0,
		'shat_max': 10,
		'beta_min': 0.01,
		'beta_max': 0.7,
		'shat_div': None,
		'shat_mul': None,
		'beta_div': None,
		'beta_mul': None,
		'n_shat_ideal': 100,
		'n_beta_ideal': 100,
		'Miller': False,
		'Ideal': True,
		'Fixed_delt': False,
		'Epar': False,
		}
		if eq_file:
			self.load_geqdsk()
		if kin_file:
			self.load_kinetics()
		if inputs:
			self.load_inputs(inputs)

	def load_geqdsk(self, eq_file = None, directory = None):
		from .geqdsk_reader import geqdsk
		
		if self.eq_name is None and eq_file is None:
			print("ERROR: No GEQDSK file given")
			return
		elif eq_file is None:
			eq_file = self.eq_name
		self.eq_name = eq_file
		
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		if directory == "./":
			directory = os.getcwd()
		self._eq_path = directory

		if self.surface_namelists:
			self.surface_namelists = {}
		if self.pyro:
			self.pyro = None
		
		with open(os.path.join(directory,self.eq_name)) as efile:
			self._eq_lines = efile.readlines()
			
		self.eq_data = geqdsk(filename = self.eq_name, directory = directory)
	
	def load_kinetics(self, kin_file = None, kinetics_type = None, directory = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		if directory == "./":
			directory = os.getcwd()
		self._kin_path = directory
		
		if self.kin_name is None and kin_file is None:
			print("ERROR: No Kinetics file given")
			return
		elif kin_file is None:
			kin_file = self.kin_name
		self.kin_name = kin_file
		
		if self.surface_namelists:
			self.surface_namelists = {}
		if self.pyro:
			self.pyro = None
		
		with open(os.path.join(directory,self.kin_name)) as kfile:
			self._kin_lines = kfile.readlines()
			
		if kinetics_type is not None:
			self.kinetics_type = kinetics_type
		if self.kinetics_type is None:
			print("Warning: kinetics_type Not Given, trying PEQDSK")
			self.kinetics_type = "PEQDSK"
		
		if self.kinetics_type.upper() == "SCENE":
			import xarray as xr
			self.kin_data = xr.open_dataset(os.path.join(self._kin_path,self.kin_name))
		elif self.kinetics_type.upper() == "PEQDSK":
			from .peqdsk_reader import peqdsk
			self.kin_data = peqdsk(filename = self.kin_name, directory = directory)
			if 'rhonorm' not in self.kin_data.keys():
				self.AmmendPEQDSK()
				self.kin_data = peqdsk(self.kin_name, directory)
				if 'rhonorm' not in self.kin_data.keys():
					print("ERROR: Could not load rho data from PEQDSK file")
					return			
		else:
			print(f"ERROR: Kinetics type {self.kinetics_type} not recognised. Currently supported: SCENE, PEQDSK")
	
	def load_pyro(self, template_file = None, directory = None):
		from pyrokinetics import Pyro
		from pathlib import Path
		if self.eq_name is None or self.kin_name is None:
			if self.eq_name is None:
				print("ERROR: No equillibrium loaded")
			if self.kin_name is None:
				print("ERROR: No kinetics file loaded")
			return
		
		if self.surface_namelists:
			self.surface_namelists = {}

		eq_file = Path(self._eq_path) / self.eq_name
		kin_file = Path(self._kin_path) / self.kin_name
		
		if template_file is not None:
				self.template_name = template_file
		if self.template_name is None:
			self.template_name = gs2_template
			directory = template_dir
		if directory is None and self.path is None:
			directory = os.getcwd()
		elif directory is None:
			directory = self.path
		if directory == "./":
			directory = os.getcwd()
		self._template_path = directory
		
		self._template_lines = f90nml.read(os.path.join(directory,self.template_name))

		self.pyro = Pyro(
			eq_file = eq_file,
		 	eq_type = "GEQDSK",
		 	kinetics_file = kin_file,
		 	kinetics_type = self.kinetics_type,
		 	gk_file = Path(directory) / self.template_name)
		self.pyro.gk_code = "GS2"
		
	def load_inputs(self, inputs):
		self.surface_namelists = {}
		for key, val in inputs.items():
			if key in self._inputs.keys():
				self.inputs[key] = val
	
	def edit_inputs(self, key = None, val = None):
		self.surface_namelists = {}
		if key is not None and key not in self._inputs.keys():
			print(f"ERROR: {key} is not a valid key. Valid keys: {self.inputs.keys()}")
			return
		elif key is not None and val is not None:
			self.inputs[key] = val
			return
		elif key is not None:
			keys = [key]
		else:
			keys = self._inputs.keys()
		for key in keys:
			if key in self.inputs.keys():
				if self.inputs[key] is None:
					val = input(f"Input value for {key} (Current value: None): ")
				else:
					val = input(f"Input value for {key} (Current value: {self.inputs[key]}): ")
			else:
				val = input(f"Input value for {key} (Current value: None): ")
			if val != "":
				self.inputs[key] = eval(val)

	def AmmendPEQDSK(self):
		if self.eq_data is None and self.eq_name is None:
			print("ERROR: No GEQDSK file provided")
			return
		elif self.eq_data is None:
			self.load_geqdsk()
		if self.kin_data is None and self.kin_name is None:
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
		
		f = open(self.kin_name,'a')
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
		if self.inputs.keys() != self._inputs.keys():
			idiff = [x for x in self._inputs.keys() if x not in self.inputs.keys()]
			if 'beta_min' not in idiff and 'beta_div' in idiff:
				idiff.remove('beta_div')
			if 'beta_max' not in idiff and 'beta_mul' in idiff:
				idiff.remove('beta_mul')
			if 'shear_min' not in idiff and 'shear_div' in idiff:
				idiff.remove('shear_div')
			if 'shear_max' not in idiff and 'shear_mul' in idiff:
				idiff.remove('shear_mul')
			if 'beta_min' in idiff and 'beta_div' not in idiff:
				idiff.remove('beta_min')
			if 'beta_max' in idiff and 'beta_mul' not in idiff:
				idiff.remove('beta_max')
			if 'shear_min' in idiff and 'shear_div' not in idiff:
				idiff.remove('shear_min')
			if 'shear_max' in idiff and 'shear_mul' not in idiff:
				idiff.remove('shear_max')
			
			if idiff:
				print("Warning: {idiff} input keys missing, using default")
				for key in idiff:
					self.inputs[key] = self._inputs[key]

		self.pyro.load_local(psi_n=psiN)
		self.pyro.update_gk_code()
		nml = deepcopy(self.pyro.gk_input.data)
		
		shear = nml['theta_grid_eik_knobs']['s_hat_input']
		beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
		beta =  nml['parameters']['beta']
		
		bp_cal = self.cal_bp(nml)

		mul = beta_prim/bp_cal
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			nml[spec]['tprim'] = nml[spec]['tprim']*mul
			nml[spec]['fprim'] = nml[spec]['fprim']*mul
		
		#Set bakdif to 0 for Electormagnetic Runs as a default
		nml['dist_fn_species_knobs_1']['bakdif'] = 0
		nml['dist_fn_species_knobs_2']['bakdif'] = 0
		nml['dist_fn_species_knobs_3']['bakdif'] = 0	
		
		if self.inputs['Miller']:
			nml['theta_grid_eik_knobs']['iflux'] = 0
			nml['theta_grid_eik_knobs']['local_eq'] = True
		else:
			nml['theta_grid_eik_knobs']['eqfile'] = os.path.join(self._eq_path,self.eq_name)
			nml['theta_grid_eik_knobs']['efit_eq'] =  True
			nml['theta_grid_eik_knobs']['iflux'] = 1
			nml['theta_grid_eik_knobs']['local_eq'] = False
			
		if self.inputs['Epar']:
			nml['gs2_diagnostics_knobs']['write_final_epar'] = True
		else:
			nml['gs2_diagnostics_knobs']['write_final_epar'] = False
		
		if 'ntheta_geometry' not in nml['theta_grid_eik_knobs'].keys():
			nml['theta_grid_eik_knobs']['ntheta_geometry'] = 4096

		if self.inputs['Ideal']:
			if self.inputs['beta_div'] is None:
				beta_div = abs(beta_prim/self.inputs['beta_min'])
			elif self.inputs['beta_min'] is None:
				beta_div = self.inputs['beta_div']
			else:
				beta_div = min(self.inputs['beta_div'],beta_prim/self.inputs['beta_min'])
				
			if self.inputs['beta_mul'] is None:
				beta_mul = abs(self.inputs['beta_max']/beta_prim)
			elif self.inputs['beta_max'] is None:
				beta_mul = self.inputs['beta_mul']
			else:
				beta_mul = min(self.inputs['beta_mul'],self.inputs['beta_max']/beta_prim)
				
			if self.inputs['shat_min'] is None:
				shat_min = shear/self.inputs['shat_div']
			elif self.inputs['shat_div'] is None:
				shat_min = self.inputs['shat_min']
			else:
				shat_min = max(self.inputs['shat_min'],shear/self.inputs['shat_div'])
				
			if self.inputs['shat_max'] is None:
				shat_max = self.inputs['shat_mul']*shear
			elif self.inputs['shat_mul'] is None:
				shat_max = self.inputs['shat_max']
			else:
				shat_max = min(self.inputs['shat_max'],self.inputs['shat_mul']*shear)
				
			nml['ballstab_knobs'] = {'n_shat': self.inputs['n_shat_ideal'], 'n_beta': self.inputs['n_beta_ideal'], 'shat_min': shat_min, 'shat_max': shat_max, 'beta_mul': beta_mul, 'beta_div': beta_div}
		try:
			if nml['kt_grids_knobs']['grid_option'] in ['single','default']:
				nml['knobs']['wstar_units'] = False
		except:
			nml['knobs']['wstar_units'] = False
		
		self.surface_namelists[psiN] = nml
		
		return deepcopy(self.surface_namelists[psiN])
				
	def get_gyro_input(self, psiN, bp, sh, aky, namelist_diff = {}):
		nml = self.get_surface_input(psiN)

		beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
		nml['theta_grid_eik_knobs']['s_hat_input'] = sh
		nml['theta_grid_eik_knobs']['beta_prime_input'] = -1*abs(bp)
		nml['kt_grids_single_parameters']['aky'] = aky
	
		bp_cal = self.cal_bp(nml)

		mul = bp/bp_cal
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			nml[spec]['tprim'] = nml[spec]['tprim']*mul
			nml[spec]['fprim'] = nml[spec]['fprim']*mul
					
		if self.inputs['Fixed_delt'] is False:
			delt = 0.04/aky
			if delt > 0.01:
				delt = 0.01
			nml['knobs']['delt'] = delt
		
		for key in namelist_diff.keys():
			for skey in namelist_diff[key].keys():
				nml[key][skey] = namelist_diff[key][skey]
			
		return nml
	
	def cal_bp(self, nml):
		bp_cal = sum((nml[spec]['tprim'] + nml[spec]['fprim'])*nml[spec]['dens']*nml[spec]['temp'] for spec in [x for x in nml.keys() if 'species_parameters_' in x])*nml['parameters']['beta']*-1
		return bp_cal
		
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
		from matplotlib.pyplot import subplots, show
		from numpy import linspace
		if not self.beta_prime_profile or not self.shear_profile:
			self.make_profiles()

		fig, ax = subplots(2,1)
		psiNs = linspace(0.01,1,100)
		bp = self.beta_prime_profile(psiNs)
		sh = self.shear_profile(psiNs)
		
		ax[0].plot(psiNs, bp, 'b')
		ax[0].invert_yaxis()
		ax[1].plot(psiNs, sh, 'b')
		ax[0].set_xlabel("psiN")
		ax[1].set_xlabel("psiN")
		ax[0].set_ylabel("beta_prime")
		ax[1].set_ylabel("shear")
		show()
