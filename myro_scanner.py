import os
from numpy import full, real, imag, nan, amax, array, isfinite, loadtxt, transpose, savez
from ncdf2dict import ncdf2dict as readnc
from .equillibrium import equillibrium
import f90nml

'''
GYROKINETIC SCAN PERFORMER
'''

class myro_scan(object):
	def __init__(self, input_file = None, eq_file = None, kin_file = None, template_file = None, kinetics_type = "PEQDSK", directory = "./", run_name = None):
		self._create_empty_inputs()
		self.template_name = template_file
		self.input_name = input_file
		self.run_name = run_name
		self._template_path = self.info = self.pyro = self._template_lines = self.dat = self.file_lines = self.verify = self.namelist_diffs = None

		if directory == "./":
			directory = os.getcwd() 
		self.path = directory
		
		self.eqbm = self.equillibrium = equillibrium(eq_file = eq_file, kin_file = kin_file, kinetics_type = kinetics_type, directory = directory)
		
		if input_file is not None:
			self.load_inputs()
	
	def __getitem__(self, key):
		if key == "inputs":
			self.inputs()
		else:
        		return self.inputs[key]

	@property
	def print_inputs(self):
        	for key, val in self.inputs.items():
        		print(f"{key} = {val}")
        	
	def keys(self):
		return self.inputs.keys()
		
	def _create_empty_inputs(self):
		self.inputs = {
		'shat_min': None,
		'shat_max': None,
		'beta_min': None,
		'beta_max': None,
		'shat_div': None,
		'shat_mul': None,
		'beta_div': None,
		'beta_mul': None,
		'n_shat': None,
		'n_beta': None,
		'n_shat_ideal': None,
		'n_beta_ideal': None,
		'psiNs': None,
		'aky_values': None,
		'Miller': True,
		'Ideal': True,
		'Gyro': True,
		'Viking': False,
		'Fixed_delt': False,
		'Epar': False,
		'grad_type': 0
		}

	def create_inputs(self, key = None, val = None):
		self.edit_inputs(key = key, val = val)
	
	def edit_inputs(self, key = None, val = None):
		if key is not None and key not in self.inputs.keys():
			print(f"ERROR: {key} is not a valid key. Valid keys: {self.inputs.keys()}")
			return
		elif key is not None and val is not None:
			self.inputs[key] = val
			print(f"{key} set to {val}")
			return
		elif key is not None:
			keys = [key]
		else:
			keys = self.inputs.keys()
		for key in keys:
			if self.inputs[key] is None:
				val = input(f"Input value for {key} (Current value: None): ")
			else:
				val = input(f"Input value for {key} (Current value: {self.inputs[key]}): ")
			if key == 'psiNs' or key == 'aky_values':
				if val[0] != "[":
					val = "[" + val
				if val[-1] != "]":
					val = val+ "]"
			if val not in ["","[]"]:
				self.inputs[key] = eval(val)
	
	def load_geqdsk(self, eq_file = None, directory = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
			
		if eq_file is None:
			print("ERROR: eq_file Not Given")
			return
		
		self.eqbm.load_geqdsk(eq_file = eq_file, directory = directory)
	
	def load_kinetics(self, kin_file = None, kinetics_type = None, directory = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path

		if kin_file is None:
			print("ERROR: kin_file Not Given")
			return
		
		if kinetics_type is None:
			print("ERROR: kinetics_type Not Given, trying PEQDSK")
			kinetics_type = "PEQDSK"
			
		self.eqbm.load_kinetics(self, kin_file = kin_file, kinetics_type = kinetics_type, directory = directory)
	
	def load_pyro(self, template_file = None, directory = None):
		if template_file is not None:
				self.template_name = template_file
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		if directory == "./":
			directory = os.getcwd()
		self._template_path = directory
		
		if self.template_name:
			self._template_lines = f90nml.read(os.path.join(directory,self.template_name))
		
		self.pyro = self.eqbm.load_pyro(template_file = self.template_name, directory = directory)
		
		
	def load_inputs(self, filename = None, directory = "./"):
		if self.input_name is None and filename is None:
			filename = input("Input File Name: ")
			if "." not in filename:
				filename = filename + ".in"
		elif self.input_name is None:
			if "." not in filename:
				filename = filename + ".in"
			self.input_name = filename

		if self.run_name is None:
			self.run_name = self.input_name.split(".")[0]

		with open(os.path.join(self.path,self.input_name)) as in_file:
			lines = in_file.readlines()
			for line in lines:
				key = line.split(" = ")[0].strip("\t\n ")
				val = line.split(" = ")[1].strip("\t\n ")
				if key == 'akys':
					key = 'aky_values'
				if key == 'psiNs' or key == 'aky_values':
					if val[0] != "[":
						val = "[" + val
					if val[-1] != "]":
						val = val + "]"
				if key not in self.inputs.keys():
					print(f"ERROR: key '{key}' not found")
				
				else:
					self.inputs[key] = eval(val)
			
	def write_scan_input(self, filename = None, directory = "./"):
		if directory is None and self.directory is None:
			directory = "./"
		elif directory is None:
			directory = self.directory
			
		if self.input_name is None and filename is None:
			filename = input("Input File Name: ")
		elif self.input_name is None:
			self.input_name = filename
		if "." not in self.input_name:
			self.input_name = self.input_name + ".in"
		
		with open(os.path.join(directory,self.input_name),'w') as in_file:
			for key in self.inputs.keys():
				in_file.write(f"{key} = {self.inputs[key]}\n")
	
	def _make_fs_in(self, run_path = None, psiN = None):
		if self.pyro is None:
			self.load_pyro()
		if run_path is None:
			run_path = self.path
		if psiN is None:
			print("ERROR: please speicify psiN")
			return
			
		self.pyro.load_local_geometry(psi_n=psiN)
		self.pyro.load_local_species(psi_n=psiN)
		file_name=f"{run_path}/{psiN}.in"
		self.pyro.write_gk_file(file_name)
		nml = self.pyro._gk_input_record["GS2"].data
		
		shear = nml['theta_grid_eik_knobs']['s_hat_input']
		beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
		tprim =  nml['species_parameters_1']['tprim']
		fprim =  nml['species_parameters_1']['fprim']
		beta =  nml['parameters']['beta']
		#Set bakdif to 0 for Electormagnetic Runs as a default
		nml['dist_fn_species_knobs_1']['bakdif'] = 0
		nml['dist_fn_species_knobs_2']['bakdif'] = 0
		nml['dist_fn_species_knobs_3']['bakdif'] = 0
		
		if self.inputs['Miller']:
			nml['theta_grid_eik_knobs']['iflux'] = 0
			nml['theta_grid_eik_knobs']['local_eq'] = True
		else:
			nml['theta_grid_eik_knobs']['eqfile'] = os.path.join(self.path,self.eqbm.eq_name)
			nml['theta_grid_eik_knobs']['efit_eq'] =  True
			nml['theta_grid_eik_knobs']['iflux'] = 1
			nml['theta_grid_eik_knobs']['local_eq'] = False
			
		if self.inputs['Epar']:
			nml['gs2_diagnostics_knobs']['write_final_epar'] = True
		else:
			nml['gs2_diagnostics_knobs']['write_final_epar'] = False

		if self.inputs['Ideal']:
			if self.inputs['beta_div'] is None:
				beta_div = beta_prim/self.inputs['beta_min']
			elif self.inputs['beta_min'] is None:
				beta_div = self.inputs['beta_div']
			else:
				beta_div = min(self.inputs['beta_div'],beta_prim/self.inputs['beta_min'])
				
			if self.inputs['beta_mul'] is None:
				beta_mul = self.inputs['beta_max']/beta_prim
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
		nml.write(file_name, force=True)
		return shear, beta_prim, tprim, fprim, beta, nml
	
	def run_scan(self, gyro = None, ideal = None, directory = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		
		if gyro is None:
			gyro = self.inputs['Gyro']
		elif type(gyro) == bool:	
			self.inputs['Gyro'] = gyro
		else:
			print("ERROR: gyro must be boolean")
			return
		if ideal is None:
			ideal = self.inputs['Ideal']
		elif type(ideal) == bool:
			self.inputs['Ideal'] = ideal
		else:
			print("ERROR: ideal must be boolean")
			return
		
		if not self.info:
			self._create_run_info()
		else:
			self.info['itteration'] += 1
		run_path = self.info['data_path']
		
		try:
			os.mkdir(run_path)
		except:
			pass	
		
		if not self._check_setup(ideal = ideal, gyro = gyro):
			return
		if self.inputs['Gyro']:
			self._run_gyro(directory = run_path)
		if self.inputs['Ideal']:
			self._run_ideal(directory = run_path)
		if not self.inputs['Viking']:
			self.save_out(directory = run_path)
	
	def _check_setup(self, ideal = None, gyro = None):
		if gyro is None:
			gyro = self.inputs['Gyro']
		elif type(gyro) == bool:	
			self.inputs['Gyro'] = gyro
		else:
			print("ERROR: gyro must be boolean")
			return False
		if ideal is None:
			ideal = self.inputs['Ideal']
		elif type(ideal) == bool:
			self.inputs['Ideal'] = ideal
		else:
			print("ERROR: ideal must be boolean")
			return False
			
		if self.info is None:
			self._create_run_info()
			
		if not self.eqbm.eq_name or not self.eqbm.kin_name:
			if not self.eqbm.eq_name:
				print("ERROR: No eq_file loaded")
			if not self.eqbm.eq_name:
				print("ERROR: No kin_file loaded")
			return False
			
		if not self.pyro:
			self.load_pyro()
			
		empty_elements = []
		if type(self.inputs['psiNs']) in [int,float]:
			self.inputs['psiNs'] = [self.inputs['psiNs']]
		elif self.inputs['psiNs'] is None:
			empty_elements.append('psiNs')
		else:
			self.inputs['psiNs'].sort()
		
		if gyro:
			if type(self.inputs['aky_values']) in [int,float]:
				self.inputs['aky_values'] = [self.inputs['aky_values']]
			elif self.inputs['aky_values'] is None:
				empty_elements.append('aky_values')
			else:
				self.inputs['aky_values'].sort()

		if self.inputs['beta_min'] is None and self.inputs['beta_div'] is None:
			empty_elements.append('beta_min/beta_div')
		if self.inputs['beta_max'] is None and self.inputs['beta_mul'] is None:
			empty_elements.append('beta_max/beta_mul')
		if self.inputs['shat_min'] is None and self.inputs['shat_div'] is None:
			empty_elements.append('shat_min/shat_div')
		if self.inputs['shat_max'] is None and self.inputs['shat_mul'] is None:
			empty_elements.append('shat_max/shat_mul')
		
		if self.inputs['n_shat_ideal'] is None and self.inputs['n_shat'] is None:
			if ideal:
				empty_elements.append('n_shat_ideal')
			if gyro:
				empty_elements.append('n_shat')
		elif gyro and not self.inputs['n_shat']:
			empty_elements.append('n_shat')
		elif ideal and not self.inputs['n_shat_ideal']:
			self.inputs['n_shat_ideal'] = self.inputs['n_shat']
			print("n_shat_ideal is empty, setting equal to n_shat")
		
		if self.inputs['n_beta_ideal'] is None and self.inputs['n_beta'] is None:
			if ideal:
				empty_elements.append('n_beta_ideal')
			if gyro:
				empty_elements.append('n_beta')
		elif gyro and not self.inputs['n_beta']:
			empty_elements.append('n_beta')
		elif ideal and not self.inputs['n_beta_ideal']:
			self.inputs['n_beta_ideal'] = self.inputs['n_beta']
			print("n_beta_ideal is empty, setting equal to n_beta")
		
		if empty_elements:
			print(f"ERROR: the following inputs are empty: {empty_elements}")
			return False
		
		if gyro and not self.namelist_diffs:
			self.namelist_diffs = self.namelist_diffs = [[[[{} for _ in range(len(self.inputs['aky_values']))] for _ in range(self.inputs['n_shat'])] for _ in range(self.inputs['n_beta'])] for _ in range(len(self.inputs['psiNs']))]
		
		if not os.path.exists(self.info['data_path']):
				os.mkdir(self.info['data_path'])

		if self.template_name is not None:
			os.system(f"cp \"{self._template_path}/{self.template_name}\" \"{self.info['data_path']}/{self.template_name}\"")
		os.system(f"cp \"{self.eqbm._kin_path}/{self.eqbm.kin_name}\" \"{self.info['data_path']}/{self.eqbm.kin_name}\"")
		os.system(f"cp \"{self.eqbm._eq_path}/{self.eqbm.eq_name}\" \"{self.info['data_path']}/{self.eqbm.eq_name}\"")
		if self.input_name is not None:		
			os.system(f"cp \"{self.path}/{self.input_name}\" \"{self.info['data_path']}/{self.input_name}\"")
		
		return True
	
	def _run_ideal(self, directory = None, checkSetup = True):
		self.inputs['Ideal'] = True
		if checkSetup:
			if not self._check_setup():
				return
		if directory is None:
			directory = self.info['data_path']
		check = self.check_complete(directory = directory, doPrint = False, gyro = False, ideal = True)
		if check['ideal_complete']:
			print(f"{len(check['ideal_complete'])} Existing Ideal Runs Detected")
		for psiN in check['ideal_incomplete']:
			run_path = os.path.join(directory, str(psiN))
			try:
				os.mkdir(run_path)
			except:
				pass
			file_name = os.path.join(run_path, f"{psiN}.in")
			if not os.path.isfile(file_name):
				self._make_fs_in(run_path=run_path,psiN=psiN)
			
			if not self.inputs['Viking']:
				os.system(f"ideal_ball \"{run_path}/{psiN}.in\"")
			else:
				
				jobfile = open(f"{run_path}/{psiN}.job",'w')
				jobfile.write(f"#!/bin/bash\n#SBATCH --time=01:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n#SBATCH --output={psiN}.slurm\n\nmodule purge\nmodule load tools/git\nmodule load compiler/ifort\nmodule load mpi/impi\nmodule load numlib/FFTW\nmodule load data/netCDF/4.6.1-intel-2018b\nmodule load data/netCDF-Fortran/4.4.4-intel-2018b\nmodule load numlib/imkl/2018.3.222-iimpi-2018b\nmodule load lang/Python/3.7.0-intel-2018b\nexport GK_SYSTEM=viking\nexport MAKEFLAGS=-IMakefiles\nexport PATH=$PATH:$HOME/gs2/bin\n\nwhich gs2\n\ngs2 --build-config\n\nideal_ball \"{run_path}/{psiN}.in\"")
				jobfile.close()
				os.chdir(f"{run_path}")
				os.system(f"sbatch \"{run_path}/{psiN}.job\"")
				os.chdir(f"{self.path}")

	def _run_gyro(self, directory = None, checkSetup = True, specificRuns = None):
		self.inputs['Gyro'] = True
		if checkSetup:
			if not self._check_setup():
				return
		if directory is None:
			directory = self.info['data_path']
		if not specificRuns:
			check = self.check_complete(directory = directory, doPrint = False, gyro = True, ideal = False)
			if check['gyro_complete']:
				print(f"{len(check['gyro_complete'])} Existing Gyro Runs Detected")
			runs = check['gyro_incomplete']
		else:
			runs = specificRuns

		if len(runs) > 10000:
			group_runs = True
		else:
			group_runs = False
		
		for p, psiN in enumerate(self.inputs['psiNs']):
			run_path = os.path.join(directory, str(psiN))
			try:
				os.mkdir(run_path)
			except:
				pass
				
			shear, beta_prim, tprim, fprim, beta, nml = self._make_fs_in(run_path=run_path, psiN=psiN)
			
			if self.inputs['beta_min'] is None:
				beta_min = beta_prim/self.inputs['beta_div']
			elif self.inputs['beta_div'] is None:
				beta_min = self.inputs['beta_min']
			else:
				beta_min = max(self.inputs['beta_min'],beta_prim/self.inputs['beta_div'])
				
			if self.inputs['beta_max'] is None:
				beta_max = self.inputs['beta_mul']*beta_prim
			elif self.inputs['beta_mul'] is None:
				beta_max = self.inputs['beta_max']
			else:
				beta_max = min(self.inputs['beta_max'],self.inputs['beta_mul']*beta_prim)
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
			
			
			for i in range(self.inputs['n_beta']):
				for j in range(self.inputs['n_shat']):
					group_kys = []
					fol = f"{i}_{j}"
					sub_path = os.path.join(run_path,fol)
					try:
						os.mkdir(sub_path)
					except:
						pass
					bp = (beta_max - beta_min)*i/(self.inputs['n_beta']-1) + beta_min
					sh = (shat_max - shat_min)*j/(self.inputs['n_shat']-1) + shat_min
					if sh < 1e-4:
						sh = 1e-4

					for k, aky in enumerate(self.inputs['aky_values']):
						if (p,i,j,k) in runs:
							if os.path.exists(f"{sub_path}/{p}_{fol}_{k}_old.out.nc"):
								os.remove(f"{sub_path}/{p}_{fol}_{k}_old.out.nc")
							if os.path.exists(f"{sub_path}/{p}_{fol}_{k}.out.nc"):
								os.rename(f"{sub_path}/{p}_{fol}_{k}.out.nc",f"{sub_path}/{p}_{fol}_{k}_old.out.nc")
							if os.path.exists(f"{sub_path}/{p}_{fol}_{k}.slurm"):
								os.rename(f"{sub_path}/{p}_{fol}_{k}.slurm",f"{sub_path}/{p}_{fol}_{k}_old.slurm")
							if os.path.exists(f"{sub_path}/{p}_{fol}.slurm"):
								os.rename(f"{sub_path}/{p}_{fol}.slurm",f"{sub_path}/{p}_{fol}_old.slurm")
							subnml = nml
							subnml['theta_grid_eik_knobs']['s_hat_input'] = sh
							subnml['theta_grid_eik_knobs']['beta_prime_input'] = bp
							subnml['kt_grids_single_parameters']['aky'] = aky
							for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
								if self.inputs['grad_type'] == 2:
									mul = (bp/(beta*-2) - nml[spec]['tprim'])/nml[spec]['fprim']
									subnml[spec]['fprim'] = nml[spec]['fprim']*mul
								elif self.inputs['grad_type'] == 1:
									mul = (bp/(beta*-2) - nml[spec]['tprim'])/nml[spec]['fprim']
									subnml[spec]['tprim'] = nml[spec]['fprim']*mul
								else:
									mul = bp/(-2*(nml[spec]['tprim'] + nml[spec]['fprim'])*beta)
									subnml[spec]['tprim'] = nml[spec]['tprim']*mul
									subnml[spec]['fprim'] = nml[spec]['fprim']*mul
							if self.inputs['Fixed_delt'] is False:
								delt = 0.04/aky
								if delt > 0.1:
									delt = 0.1
								subnml['knobs']['delt'] = delt
							
							if self.namelist_diffs[p][i][j][k]:
								for key in self.namelist_diffs[p][i][j][k].keys():
									for skey in self.namelist_diffs[p][i][j][k][key].keys():
										subnml[key][skey] = self.namelist_diffs[p][i][j][k][key][skey]
							
							subnml.write(f"{sub_path}/{p}_{fol}_{k}.in", force=True)
							
							if not self.inputs['Viking']:
								os.system(f"mpirun -np 8 gs2 \"{sub_path}/{p}_{fol}_{k}.in\"")
							elif not group_runs:
								
								jobfile = open(f"{sub_path}/{p}_{fol}_{k}.job",'w')
								jobfile.write(f"#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n#SBATCH --output={p}_{fol}_{k}.slurm\n\nmodule purge\nmodule load tools/git\nmodule load compiler/ifort\nmodule load mpi/impi\nmodule load numlib/FFTW\nmodule load data/netCDF/4.6.1-intel-2018b\nmodule load data/netCDF-Fortran/4.4.4-intel-2018b\nmodule load numlib/imkl/2018.3.222-iimpi-2018b\nmodule load lang/Python/3.7.0-intel-2018b\nexport GK_SYSTEM=viking\nexport MAKEFLAGS=-IMakefiles\nexport PATH=$PATH:$HOME/gs2/bin\n\necho \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"\necho \"Input: {sub_path}/{p}_{fol}_{k}.in\"\necho \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"\n\nwhich gs2\n\ngs2 --build-config\n\ngs2 \"{sub_path}/{p}_{fol}_{k}.in\"")
								
								jobfile.close()
								os.chdir(f"{sub_path}")
								os.system(f"sbatch \"{sub_path}/{p}_{fol}_{k}.job\"")
								os.chdir(f"{self.path}")
							
							else:
								group_kys.append(k)
							
						else:
							pass
				
					if self.inputs['Viking'] and group_runs and group_kys:
						hours = 4*len(self.inputs['aky_values'])
						jobfile = open(f"{sub_path}/{p}_{fol}.job",'w')
						jobfile.write(f"#!/bin/bash\n#SBATCH --time={hours}:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n#SBATCH --output={p}_{fol}.slurm\n\nmodule purge\nmodule load tools/git\nmodule load compiler/ifort\nmodule load mpi/impi\nmodule load numlib/FFTW\nmodule load data/netCDF/4.6.1-intel-2018b\nmodule load data/netCDF-Fortran/4.4.4-intel-2018b\nmodule load numlib/imkl/2018.3.222-iimpi-2018b\nmodule load lang/Python/3.7.0-intel-2018b\nexport GK_SYSTEM=viking\nexport MAKEFLAGS=-IMakefiles\nexport PATH=$PATH:$HOME/gs2/bin\n\nwhich gs2\n\ngs2 --build-config\n\n")

						for k in group_kys:
							jobfile.write(f"echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"\necho \"Input: {sub_path}/{p}_{fol}_{k}.in\"\necho \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"\n\ngs2 \"{sub_path}/{p}_{fol}_{k}.in\"\n")
						jobfile.close()
						os.chdir(f"{sub_path}")
						os.system(f"sbatch \"{sub_path}/{fol}.job\"")
						os.chdir(f"{self.path}")
	def _create_run_info(self):
		try:
			import uuid
			ID = uuid.uuid4()
		except:
			print("ERROR: unable to import uuid module, setting ID to None")
			ID = None
			
		path = self.path
		if path[:2] == "./":
			path = os.getcwd() + path[1:]
			
		if self.run_name is None:
			run_path = os.path.join(path, str(ID))
		else:	
			run_path = os.path.join(path, self.run_name)
		try:
			from datetime import datetime as dt
			date = dt.now().strftime("%d/%m/%Y, %H:%M:%S")
		except:
			print("ERROR: unable to import datetime module, setting run date to None")
			date = None
		itt = 0
		self.info = {'run_name': self.run_name, 'run_uuid': str(ID), 'data_path': run_path, 'input_file': self.input_name, 'eq_file_name': self.eqbm.eq_name, 'template_file_name': self.template_name, 'kin_file_name': self.eqbm.kin_name, 'kinetics_type': self.eqbm.kinetics_type, 'run_data': date, '_eq_file_path': self.eqbm._eq_path, '_kin_file_path': self.eqbm._kin_path, '_template_file_path': self._template_path, 'itteration': 0}
	
	def rerun_cancelled(self, directory = None, checkSetup = True):
		if self.info is None:
			self._create_run_info()
		if directory is None:
			directory = self.info['data_path']
		cancelled = self.check_cancelled(directory = directory, doPrint = False)
		if len(cancelled) == 0:
			print("No Cancelled Runs")
			return
		for p, i, j, k in cancelled:
			if os.path.exists(f"{directory}/{self.inputs['psiNs'][p]}/{i}_{j}/{p}_{i}_{j}_{k}.out.nc"):
				os.remove(f"{directory}/{self.inputs['psiNs'][p]}/{i}_{j}/{p}_{i}_{j}_{k}.out.nc")
		self._run_gyro(directory = directory, specificRuns = cancelled)
		
	def check_cancelled(self, directory = None, doPrint = True):
		if self.info is None:
			self._create_run_info()		
		if directory is None:
			directory = self.info['data_path']
		if not self.inputs['Gyro']:
			print("Only Used For Gyro Runs")
			return
		
		cancelled = set()
		os.system(f"grep --include=\*slurm -rnwl {directory} -e \"CANCELLED\" > {directory}/grep.out")
		grep = open(f"{directory}/grep.out")
		lines = grep.readlines()
		grep.close()
		if len(lines) > 0:
			for line in lines:
				if "_old.slurm" not in line:
					ids = line.split("/")[-1].strip("\n").strip(".slurm").split("_")
					if len(ids) == 4:
						p, i, j, k = [eval(x) for x in ids]
						cancelled.add((p,i,j,k))
					elif len(ids) == 3:
						f = open(line.strip("\n"))
						lins = f.readlines()
						f.close()
						for l in lins:
							if ".in" in l:
								inp = l.split("/")[-1].split(".")[0]
						p, i, j, k = [eval(x) for x in inp.split("_")]
						for ki in range(k, len(self.inputs['aky_values'])):
							cancelled.add((p,i,j,ki))
		if doPrint:		
			print(f"{len(cancelled)} Cancelled Runs")
			return
		else:
			return cancelled
	
	def check_complete(self, directory = None, doPrint = True, ideal = None, gyro = None):
		if self.info is None:
			self._create_run_info()		
		if directory is None:
			directory = self.info['data_path']
			
		if gyro is None:
			gyro = self.inputs['Gyro']
		elif type(gyro) != bool:	
			print("ERROR: gyro must be boolean")
			return
		if ideal is None:
			ideal = self.inputs['Ideal']
		elif type(ideal) != bool:
			print("ERROR: ideal must be boolean")
			return
		
		unfinished_gyro = set()
		finished_gyro = set()
		if gyro:
			for p, psiN in enumerate(self.inputs['psiNs']):
				for i in range(self.inputs['n_beta']):
					for j in range(self.inputs['n_shat']):
						for k, aky in enumerate(self.inputs['aky_values']):
							if os.path.exists(f"{directory}/{psiN}/{i}_{j}/{p}_{i}_{j}_{k}.out.nc"):
								finished_gyro.add((p,i,j,k))
							else:
								unfinished_gyro.add((p,i,j,k))

		unfinished_ideal = set()
		finished_ideal = set()
		if ideal:
			for p, psiN in enumerate(self.inputs['psiNs']):
				if os.path.exists(f"{directory}/{psiN}/{psiN}.ballstab_2d"):
					finished_ideal.add(psiN)
				else:
					unfinished_ideal.add(psiN)
		
		if doPrint:
			print(f"Gyro Runs Complete: {len(finished_gyro)} | Incomplete : {len(unfinished_gyro)}")
			print(f"Ideal Runs Complete: {len(finished_ideal)} | Incomplete : {len(unfinished_ideal)}")
			return
		else:
			return {'gyro_complete': finished_gyro, 'gyro_incomplete': unfinished_gyro, 'ideal_complete': finished_ideal, 'ideal_incomplete': unfinished_ideal}
	
	def _save_obj(self, filename = None, directory = None):
		if filename is None:
			filename = "scan.obj"
		if directory is None:
			directory = self.path
		import pickle
		temp = self.pyro
		self.pyro = None
		with open(filename,'wb') as obj:
			pickle.dump(self,obj)
		self.pyro = temp
	
	def quick_save(self, filename = None, directory = None, VikingSave = False):
		self.save_out(filename = filename, directory = directory, VikingSave = VikingSave, QuickSave = True)
	def save_out(self, filename = None, directory = None, VikingSave = False, QuickSave = False):
		if filename is None and self.run_name is None:
			filename = input("Output File Name: ")
			filename = filename.split(".")[0]
		elif filename is None:
			filename = self.run_name
			
		if self.info is None:
			self._create_run_info()
		if directory is None:
			directory = self.info['data_path']
		
		if not self.inputs['Gyro'] and not self.inputs['Ideal']:
			print("Error: Both Gyro and Ideal are False")
			return
		
		if self.inputs['Viking'] and not VikingSave:
			os.chdir(f"{directory}")
			self._save_obj(filename = "scan.obj", directory = directory)
			job = open(f"save_out.job",'w')
			job.write(f"#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n#SBATCH --mem=10gb\n#SBATCH --output=save_out.slurm\n\nmodule load lang/Python/3.7.0-intel-2018b\nmodule swap lang/Python lang/Python/3.10.4-GCCcore-11.3.0\n\nsource $HOME/pyroenv2/bin/activate\n\npython {directory}/save_out.py")
			job.close()
			pyth = open(f"save_out.py",'w')
			#pyth.write(f"from Myrokinetics import myro_scan\nimport pickle\nscan = open(\"{directory}/scan.obj\",\'rb\')\nrun = pickle.load(scan)\nscan.close()\nrun.save_out(filename = \"{filename}\", directory = \"{directory}\", VikingSave = True, QuickSave = {QuickSave})")
			pyth.write(f"from Myrokinetics import myro_scan\nimport pickle\nwith open(\"{directory}/scan.obj\",\'rb\') as obj\n\trun = pickle.load(obj)\nrun.save_out(filename = \"{filename}\", directory = \"{directory}\",VikingSave = True,QuickSave = {QuickSave})")
			pyth.close()
			os.system(f"sbatch \"save_out.job\"")
			os.chdir(f"{self.path}")
			return
			
		if not self._check_setup():
			return
			
		psiNs = self.inputs['psiNs']
		beta_prime_values = full((len(psiNs)),None)
		shear_values = full((len(psiNs)),None)
		
		if self.inputs['Gyro']:
			beta_prime_axis = full((len(psiNs),self.inputs['n_beta']),None).tolist()
			shear_axis = full((len(psiNs),self.inputs['n_beta']),None).tolist()
			gr = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat']),None).tolist()
			mf = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat']),None).tolist()
			sym = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat']),None).tolist()
			akys = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat']),None).tolist()
			grs = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
			mfs = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
			syms = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
			
			if self.inputs['Epar']:
				eparN = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat']),None).tolist()
				eparNs = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
			else:
				eparN = None
				eparNs = None
			
			if not QuickSave:
				omega = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				phi = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				apar = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				bpar = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				phi2 = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				time = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				theta = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
			else:
				omega = phi = apar = bpar = phi2 = time = theta = None
		
		else:
				gr = mf = grs = mfs = sym = syms = beta_prime_axis = shear_axis = akys = self.inputs['n_shat'] = self.inputs['n_beta'] = self.inputs['aky_values'] = eparN = eparNs = omega = phi = apar = bpar = phi2 = time = theta = None
		
		if self.inputs['Ideal']:
			beta_prime_axis_ideal = full((len(psiNs),self.inputs['n_beta_ideal']),None).tolist()
			shear_axis_ideal = full((len(psiNs),self.inputs['n_shat_ideal']),None).tolist()
			stabilities = full((len(psiNs),self.inputs['n_beta_ideal'],self.inputs['n_shat_ideal']),None).tolist()
		else:
			beta_prime_axis_ideal = shear_axis_ideal = self.inputs['n_shat_ideal'] = self.inputs['n_beta_ideal'] = stabilities = None
		
		for p, psiN in enumerate(psiNs):
			run_path = os.path.join(directory,str(psiN))
			nml = f90nml.read(f"{run_path}/{psiN}.in")
			shear = nml['theta_grid_eik_knobs']['s_hat_input']
			beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
			shear_values[p] = shear
			beta_prime_values[p] = beta_prim
			
			if self.inputs['Gyro']:
				for i in range(self.inputs['n_beta']):
					per = 100*(p*self.inputs['n_beta']+i)/(self.inputs['n_beta']*len(psiNs))
					print(f"Saving {per:.2f}%", end='\r')
					for j in range(self.inputs['n_shat']):
						fol = str(i) + "_" + str(j)
						gr_aky = []
						gr_list = []
						mf_list = []
						sym_list = []
						epars = []
						for k, ky in enumerate(self.inputs['aky_values']):
							try:
							
								if self.inputs['Epar'] and not QuickSave:
									data = readnc(f"{run_path}/{fol}/{p}_{fol}_{k}.out.nc",only=['omega','phi','bpar','apar','phi2','t','theta'])
								elif self.inputs['Epar']:
									data = readnc(f"{run_path}/{fol}/{p}_{fol}_{k}.out.nc",only=['omega','phi','bpar'])
								elif not QuickSave:
									data = readnc(f"{run_path}/{fol}/{p}_{fol}_{k}.out.nc",only=['omega','phi','apar','bpar','phi2','t','theta'])
								else:
									data = readnc(f"{run_path}/{fol}/{p}_{fol}_{k}.out.nc",only=['omega','phi'])	
								
								om = data['omega'][-1,0,0]
								gr_list.append(imag(om))
								mf_list.append(real(om))
								gr_aky.append(imag(om)/(ky**2))
								
								symsum = sum(abs(data['phi'][0,0,:] + data['phi'][0,0,::-1]))/sum(abs(data['phi'][0,0,:]))
								if  symsum > 1.9:
									sym_list.append(1)
								elif symsum < 1:
									sym_list.append(-1)
								else:
									sym_list.append(0)
								
								if not QuickSave:
									try:
										omega[p][i][j][k] = data['omega'][:,0,0].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: omega")
									try:
										phi[p][i][j][k] = data['phi'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: phi")
									try:
										apar[p][i][j][k] = data['apar'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: apar")
									try:
										bpar[p][i][j][k] = data['bpar'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: bpar")
									try:
										phi2[p][i][j][k] = data['phi2'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: phi2")
									try:
										time[p][i][j][k] = data['t'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: time")
									try:
										theta[p][i][j][k] = data['theta'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: theta")
								
								if self.inputs['Epar']:
									epar_path = f"{run_path}/{fol}/{p}_{fol}_{k}.epar"
									bpar = data['bpar'][0,0,:]
									try:
										epar_data = loadtxt(epar_path)
										epar = []
										for k in range(len(epar_data[:,3])):
											epar.append(complex(epar_data[k,3],epar_data[k,4]))
										epar = array(epar)
										epar_norm = max(abs(epar))/max(abs(bpar))
										
										epars.append(epar_norm)
									except:
										print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: epar")
								
							except Exception as e:
								print(f"Save Error {psiN}/{fol}/{p}_{fol}_{k}: {e}")
								gr_aky.append(nan)
								gr_list.append(nan)
								mf_list.append(nan)
								sym_list.append(nan)
								epars.append(nan)
						try:
							aky_idx = gr_aky.index(amax(array(gr_aky)[isfinite(gr_aky)]))

							grs[p][i][j] = gr_list
							mfs[p][i][j] = mf_list
							syms[p][i][j] = sym_list
							akys[p][i][j] = self.inputs['aky_values'][aky_idx]
							gr[p][i][j] = gr_list[aky_idx]
							mf[p][i][j] = mf_list[aky_idx]
							sym[p][i][j] = sym_list[aky_idx]
							if self.inputs['Epar']:
								eparNs[p][i][j] = epars
								eparN[p][i][j] = epars[aky_idx]
						except:
							grs[p][i][j] = gr_list
							mfs[p][i][j] = gr_list
							syms[p][i][j] = None
							akys[p][i][j] = None
							gr[p][i][j] = nan
							mf[p][i][j] = nan
							sym[p][i][j] = None
							if self.inputs['Epar']:
								eparNs[p][i][j] = None
								eparN[p][i][j] = None
							
				if self.inputs['beta_min'] is None:
					beta_min = beta_prim/self.inputs['beta_div']
				else:
					beta_min = self.inputs['beta_min']
				if self.inputs['beta_max'] is None:
					beta_max = beta_prim*self.inputs['beta_mul']
				else:
					beta_max = self.inputs['beta_max']
				for i in range(self.inputs['n_beta']):
					beta_prime_axis[p][i] = abs((beta_max - beta_min)*i/(self.inputs['n_beta']-1) + beta_min)
					
				if self.inputs['shat_min'] is None:
					shat_min = shear/self.inputs['shat_div']
				else:
					shat_min = self.inputs['shat_min']
				if self.inputs['shat_max'] is None:
					shat_max = shear*self.inputs['shat_mul']
				else:
					shat_max = self.inputs['shat_max']
				for i in range(self.inputs['n_shat']):
					sh = (shat_max - shat_min)*i/(self.inputs['n_shat']-1) + shat_min
					if sh == 0:
						sh = 1e-4
					shear_axis[p][i] = sh

			if self.inputs['Ideal']:
				shear = loadtxt(f"{run_path}/{psiN}.ballstab_shat")
				bp = loadtxt(f"{run_path}/{psiN}.ballstab_bp")
				stab = loadtxt(f"{run_path}/{psiN}.ballstab_2d")
				
				bp = [abs(x) for x in bp]
				beta_prime_axis_ideal[p] = bp
				shear_axis_ideal[p] = shear
				stabilities[p] = transpose(stab)
			
		self.dat = {'beta_prime_values':beta_prime_values,
		'shear_values': shear_values,
		'beta_prime_axis': beta_prime_axis,
		'shear_axis': shear_axis,
		'beta_prime_axis_ideal': beta_prime_axis_ideal,
		'shear_axis_ideal': shear_axis_ideal,
		'growth_rates': gr,
		'mode_frequencies': mf,
		'growth_rates_all': grs,
		'mode_frequencies_all': mfs,
		'parities': sym,
		'parities_all': syms,
		'ideal_stabilities': stabilities,
		'eparN': eparN,
		'eparN_all': eparNs,
		'akys': akys,
		'omega': omega,
		'phi': phi,
		'apar': apar,
		'bpar': bpar,
		'phi2': phi2,
		'time': time,
		'theta': theta
		}
		self.file_lines = {'eq_file': self.eqbm._eq_lines, 'kin_file': self.eqbm._kin_lines, 'template_file': self._template_lines, 'namelist_differences': self.namelist_diffs}
		savez(f"{self.path}/{filename}", inputs = self.inputs, data = self.dat, run_info = self.info, files = self.file_lines)
		
	def rerun_errors(self, save = None, specificRuns = None, directory = None):
		if specificRuns is None:
			if self.dat:
				scan = {'inputs': self.inputs, 'data': self.dat, 'info': self.info, 'files': self.file_lines}
				self.verify = verify_scan(scan = scan)
			if not self.dat:
				if not save:
					print("ERROR: No Run Data or Save File")
					return
				from .myro_reader import myro_read
				myro = myro_read(filename = save, directory = directory)
				self.verify = myro.verify
			
			specificRuns = self.verify.runs_with_errors
			
		if not self.pyro:
			self.load_pyro()
			
		
		if not self.namelist_diffs:
			self.namelist_diffs = [[[[{} for _ in range(len(self.inputs['aky_values']))] for _ in range(self.inputs['n_shat'])] for _ in range(self.inputs['n_beta'])] for _ in range(len(self.inputs['psiNs']))]
		for [p,i,j,k] in specificRuns:
			if 'knobs' not in self.namelist_diffs[p][i][j][k].keys():
				self.namelist_diffs[p][i][j][k]['knobs'] = {}
			if 'theta_grid_parameters' not in self.namelist_diffs[p][i][j][k].keys():
				self.namelist_diffs[p][i][j][k]['theta_grid_parameters'] = {}
			self.namelist_diffs[p][i][j][k]['knobs']['nstep'] = 2*self._template_lines['knobs']['nstep']
			self.namelist_diffs[p][i][j][k]['theta_grid_parameters']['ntheta'] = 2*self._template_lines['theta_grid_parameters']['ntheta']
			self.namelist_diffs[p][i][j][k]['theta_grid_parameters']['nperiod'] = 2*self._template_lines['theta_grid_parameters']['nperiod']
			if self.inputs['Fixed_delt'] is False:
				delt = 0.004/self.inputs['aky_values'][k]
				if delt > 0.01:
					delt = 0.01
				self.namelist_diffs[p][i][j][k]['knobs']['delt'] = delt
			else:
				self.namelist_diffs[p][i][j][k]['knobs']['delt'] = self._template_lines['knobs']['delt']/10
		self.info['itteration'] += 1
		self._run_gyro(specificRuns = specificRuns)
	
	def rerun(self, specificRuns = None, nml = None, directory = None):
		if specificRuns is None:
			print("ERROR: specificRuns not given")
			return
		if nml is None:
			print("ERROR: nml not given")
			return
		if not self.namelist_diffs:
			self.namelist_diffs = [[[[{} for _ in range(len(self.inputs['aky_values']))] for _ in range(self.inputs['n_shat'])] for _ in range(self.inputs['n_beta'])] for _ in range(len(self.inputs['psiNs']))]
			
		for p,i,j,k in specificRuns:
			self.namelist_diffs[p][i][j][k] = nml
		self.info['itteration'] += 1
		self._run_gyro(specificRuns = specificRuns, directory = directory)
	
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
		
