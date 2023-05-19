import os
from numpy import full, real, imag, nan, amax, array, isfinite, loadtxt, transpose, savez
from .ncdf2dict import ncdf2dict as readnc
from .equillibrium import equillibrium
from .templates import modules, save_modules
import f90nml

'''
GYROKINETIC SCAN PERFORMER
'''

class myro_scan(object):
	def __init__(self, input_file = None, eq_file = None, kin_file = None, template_file = None, kinetics_type = "PEQDSK", directory = "./", run_name = None):
		self._create_empty_inputs()
		self.input_name = input_file
		self.template_name = template_file
		self.run_name = run_name
		self.namelist_diffs = None
		self._template_path = self.info = self.dat = self.file_lines = self.verify = self.namelist_diffs = None
		self.jobs = []

		if directory == "./":
			directory = os.getcwd() 
		self.path = directory

		self.eqbm = self.equillibrium = equillibrium(eq_file = eq_file, kin_file = kin_file, kinetics_type = kinetics_type, template_file = template_file, directory = directory)
		
		if input_file is not None:
			self.load_inputs()
	
	def __getitem__(self, key):
		if key == "inputs":
			self.inputs()
		else:
        		return self.inputs[key]

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
		'Miller': False,
		'Ideal': True,
		'Gyro': True,
		'Viking': False,
		'Fixed_delt': False,
		'Epar': False,
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
			if (key == 'psiNs' or key == 'aky_values') and val != "":
				if val[0] != "[":
					val = "[" + val
				if val[-1] != "]":
					val = val+ "]"
			if val != "":
				self.inputs[key] = eval(val)

		self.eqbm.load_inputs(self.inputs)
	
	def load_geqdsk(self, eq_file, directory = None):
		self.eqbm.load_geqdsk(eq_file = eq_file, directory = directory)
	
	def load_kinetics(self, kin_file, kinetics_type = None, directory = None):
		self.eqbm.load_kinetics(self, kin_file = kin_file, kinetics_type = kinetics_type, directory = directory)
	
	def load_pyro(self, template_file = None, directory = None):
		self.eqbm.load_pyro(template_file = template_file, directory = directory)
		
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

		self.eqbm.load_inputs(self.inputs)
			
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
	
	def run_scan(self, gyro = None, ideal = None, directory = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		
		if gyro is None:
			gyro = self['Gyro']
		elif type(gyro) == bool:	
			self.inputs['Gyro'] = gyro
		else:
			print("ERROR: gyro must be boolean")
			return
		if ideal is None:
			ideal = self['Ideal']
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
		if self['Ideal']:
			self._make_ideal_files(directory = run_path)
		if self['Gyro']:
			self._make_gyro_files(directory = run_path)
		self._run_jobs()
		if not self['Viking']:
			self.save_out(directory = run_path)
	
	def _check_setup(self, ideal = None, gyro = None):
		if gyro is None:
			gyro = self['Gyro']
		elif type(gyro) == bool:	
			self.inputs['Gyro'] = gyro
		else:
			print("ERROR: gyro must be boolean")
			return False
		if ideal is None:
			ideal = self['Ideal']
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
			
		empty_elements = []
		if type(self['psiNs']) in [int,float]:
			self.inputs['psiNs'] = [self['psiNs']]
		elif self['psiNs'] is None:
			empty_elements.append('psiNs')
		else:
			self.inputs['psiNs'].sort()
		
		if gyro:
			if type(self['aky_values']) in [int,float]:
				self.inputs['aky_values'] = [self['aky_values']]
			elif self['aky_values'] is None:
				empty_elements.append('aky_values')
			else:
				self.inputs['aky_values'].sort()

		if self['beta_min'] is None and self['beta_div'] is None:
			empty_elements.append('beta_min/beta_div')
		if self['beta_max'] is None and self['beta_mul'] is None:
			empty_elements.append('beta_max/beta_mul')
		if self['shat_min'] is None and self['shat_div'] is None:
			empty_elements.append('shat_min/shat_div')
		if self['shat_max'] is None and self['shat_mul'] is None:
			empty_elements.append('shat_max/shat_mul')
		
		if self['n_shat_ideal'] is None and self['n_shat'] is None:
			if ideal:
				empty_elements.append('n_shat_ideal')
			if gyro:
				empty_elements.append('n_shat')
		elif gyro and not self['n_shat']:
			empty_elements.append('n_shat')
		elif ideal and not self['n_shat_ideal']:
			self.inputs['n_shat_ideal'] = self['n_shat']
			print("n_shat_ideal is empty, setting equal to n_shat")
		
		if self['n_beta_ideal'] is None and self['n_beta'] is None:
			if ideal:
				empty_elements.append('n_beta_ideal')
			if gyro:
				empty_elements.append('n_beta')
		elif gyro and not self['n_beta']:
			empty_elements.append('n_beta')
		elif ideal and not self['n_beta_ideal']:
			self.inputs['n_beta_ideal'] = self['n_beta']
			print("n_beta_ideal is empty, setting equal to n_beta")
		
		if empty_elements:
			print(f"ERROR: the following inputs are empty: {empty_elements}")
			return False
		
		if gyro and self.namelist_diffs is None:
			self.namelist_diffs = self.namelist_diffs = [[[[{} for _ in range(len(self['aky_values']))] for _ in range(self['n_shat'])] for _ in range(self['n_beta'])] for _ in range(len(self['psiNs']))]
		
		if not os.path.exists(self.info['data_path']):
				os.mkdir(self.info['data_path'])

		if not self.eqbm.pyro:
			self.eqbm.load_pyro()

		if self.template_name is not None:
			os.system(f"cp \"{self.eqbm._template_path}/{self.template_name}\" \"{self.info['data_path']}/{self.template_name}\"")
		os.system(f"cp \"{self.eqbm._kin_path}/{self.eqbm.kin_name}\" \"{self.info['data_path']}/{self.eqbm.kin_name}\"")
		os.system(f"cp \"{self.eqbm._eq_path}/{self.eqbm.eq_name}\" \"{self.info['data_path']}/{self.eqbm.eq_name}\"")
		if self.input_name is not None:		
			os.system(f"cp \"{self.path}/{self.input_name}\" \"{self.info['data_path']}/{self.input_name}\"")
		
		return True
	
	def _run_jobs(self):
		for job in self.jobs:
			os.system(job)
		self.jobs = []
	
	def _make_ideal_files(self, directory = None, checkSetup = True):
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
			nml = self.eqbm.get_surface_input(psiN = psiN)
			nml.write(f"{run_path}/{psiN}.in", force=True)
			if not self['Viking']:
				self.jobs.append(f"ideal_ball \"{run_path}/{psiN}.in\"")
			else:
				jobfile = open(f"{run_path}/{psiN}.job",'w')
				jobfile.write(f"""#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --job-name={self.info['run_name']}
#SBATCH --ntasks=1
#SBATCH --output={run_path}/{psiN}.slurm

{modules}

which gs2
gs2 --build-config

ideal_ball \"{run_path}/{psiN}.in\"""")
				jobfile.close()
				self.jobs.append(f"sbatch \"{run_path}/{psiN}.job\"")
				
	def _make_gyro_files(self, directory = None, checkSetup = True, specificRuns = None):
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

		if len(runs) > 5000:
			group_runs = True
		else:
			group_runs = False
		
		for p, psiN in enumerate(self['psiNs']):
			run_path = os.path.join(directory, str(psiN))
			try:
				os.mkdir(run_path)
			except:
				pass
			
			nml = self.eqbm.get_surface_input(psiN = psiN)
			shear = nml['theta_grid_eik_knobs']['s_hat_input']
			beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
			beta =  nml['parameters']['beta']
			
			if self['beta_min'] is None:
				beta_min = beta_prim/self['beta_div']
			elif self['beta_div'] is None:
				beta_min = -1*abs(self['beta_min'])
			else:
				beta_min = max(self['beta_min'],beta_prim/self['beta_div'])
				
			if self['beta_max'] is None:
				beta_max = self['beta_mul']*beta_prim
			elif self['beta_mul'] is None:
				beta_max = -1*abs(self['beta_max'])
			else:
				beta_max = min(self['beta_max'],self['beta_mul']*beta_prim)
			if self['shat_min'] is None:
				shat_min = shear/self['shat_div']
			elif self['shat_div'] is None:
				shat_min = self['shat_min']
			else:
				shat_min = max(self['shat_min'],shear/self['shat_div'])
				
			if self['shat_max'] is None:
				shat_max = self['shat_mul']*shear
			elif self['shat_mul'] is None:
				shat_max = self['shat_max']
			else:
				shat_max = min(self['shat_max'],self['shat_mul']*shear)
			
			for i in range(self['n_beta']):
				for j in range(self['n_shat']):
					group_kys = []
					fol = f"{i}_{j}"
					sub_path = os.path.join(run_path,fol)
					try:
						os.mkdir(sub_path)
					except:
						pass
					bp = (beta_max - beta_min)*i/(self['n_beta']-1) + beta_min
					sh = (shat_max - shat_min)*j/(self['n_shat']-1) + shat_min
					if sh < 1e-4:
						sh = 1e-4

					for k, aky in enumerate(self['aky_values']):
						if (p,i,j,k) in runs:
							if os.path.exists(f"{sub_path}/{fol}_{k}_old.out.nc"):
								os.remove(f"{sub_path}/{fol}_{k}_old.out.nc")
							if os.path.exists(f"{sub_path}/{fol}_{k}.out.nc"):
								os.rename(f"{sub_path}/{fol}_{k}.out.nc",f"{sub_path}/{fol}_{k}_old.out.nc")
							if os.path.exists(f"{sub_path}/{fol}_{k}.slurm"):
								os.rename(f"{sub_path}/{fol}_{k}.slurm",f"{sub_path}/{fol}_{k}_old.slurm")
							if os.path.exists(f"{sub_path}/{fol}.slurm"):
								os.rename(f"{sub_path}/{fol}.slurm",f"{sub_path}/_{fol}_old.slurm")
							
							subnml = self.eqbm.get_gyro_input(psiN = psiN, bp = bp, sh = sh, aky = aky, namelist_diff = self.namelist_diffs[p][i][j][k])
							subnml.write(f"{sub_path}/{fol}_{k}.in", force=True)
							
							if not self['Viking']:
								self.jobs.append(f"mpirun -np 8 gs2 \"{sub_path}/{fol}_{k}.in\"")
							elif not group_runs:
								jobfile = open(f"{sub_path}/{fol}_{k}.job",'w')
								jobfile.write(f"""#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --job-name={self.info['run_name']}
#SBATCH --ntasks=1
#SBATCH --output={sub_path}/{fol}_{k}.slurm

{modules}

which gs2
gs2 --build-config

echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"
echo \"Input: {sub_path}/{fol}_{k}.in\"
gs2 \"{sub_path}/{fol}_{k}.in\"
echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"""")
								jobfile.close()
								self.jobs.append(f"sbatch \"{sub_path}/{fol}_{k}.job\"")
							else:
								group_kys.append(k)
							
						else:
							pass
				
					if self['Viking'] and group_runs and group_kys:
						hours = 2*len(self['aky_values'])
						if hours > 36:
							hours = 36
						jobfile = open(f"{sub_path}/{fol}.job",'w')
						jobfile.write(f"""#!/bin/bash
#SBATCH --time={hours}:00:00
#SBATCH --job-name={self.info['run_name']}
#SBATCH --ntasks=1
#SBATCH --output={sub_path}/{fol}.slurm

{modules}

which gs2
gs2 --build-config""")
						for k in group_kys:
							jobfile.write(f"""
echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"
echo \"Input: {sub_path}/{fol}_{k}.in\"
gs2 \"{sub_path}/{fol}_{k}.in\"
echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"""")
						jobfile.close()
						self.jobs.append(f"sbatch \"{sub_path}/{fol}.job\"")

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
			if os.path.exists(f"{directory}/{self['psiNs'][p]}/{i}_{j}/{i}_{j}_{k}.out.nc"):
				os.remove(f"{directory}/{self['psiNs'][p]}/{i}_{j}/{i}_{j}_{k}.out.nc")
		self._make_gyro_files(directory = directory, specificRuns = cancelled)
		self._run_jobs()
		
	def check_cancelled(self, directory = None, doPrint = True):
		if self.info is None:
			self._create_run_info()		
		if directory is None:
			directory = self.info['data_path']
		if not self['Gyro']:
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
					if len(ids) != 1:
						psiN = eval(line.split("/")[-3])
						p = self['psiNs'].index(psiN)
						if len(ids) == 3:
							i, j, k = [eval(x) for x in ids]
							cancelled.add((p,i,j,k))
						elif len(ids) == 2:
							f = open(line.strip("\n"))
							lins = f.readlines()
							f.close()
							for l in lins:
								if ".in" in l:
									inp = l.split("/")[-1].split(".")[0]
							i, j, k = [eval(x) for x in inp.split("_")]
							for ki in range(k, len(self['aky_values'])):
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
			gyro = self['Gyro']
		elif type(gyro) != bool:	
			print("ERROR: gyro must be boolean")
			return
		if ideal is None:
			ideal = self['Ideal']
		elif type(ideal) != bool:
			print("ERROR: ideal must be boolean")
			return
		
		unfinished_gyro = set()
		finished_gyro = set()
		if gyro:
			for p, psiN in enumerate(self['psiNs']):
				for i in range(self['n_beta']):
					for j in range(self['n_shat']):
						for k, aky in enumerate(self['aky_values']):
							if os.path.exists(f"{directory}/{psiN}/{i}_{j}/{i}_{j}_{k}.out.nc"):
								finished_gyro.add((p,i,j,k))
							else:
								unfinished_gyro.add((p,i,j,k))

		unfinished_ideal = set()
		finished_ideal = set()
		if ideal:
			for p, psiN in enumerate(self['psiNs']):
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
		temp = self.eqbm.pyro
		self.eqbm.pyro = None
		with open(filename,'wb') as obj:
			pickle.dump(self,obj)
		self.eqbm.pyro = temp

	def _save_for_save(self, filename = None, directory = None):
		if filename is None:
			filename = "save_info"
		if directory is None:
			directory = self.path
		savez(f"{directory}/{filename}", name_diffs = self.namelist_diffs, info = self.info)
	
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
		
		if not self['Gyro'] and not self['Ideal']:
			print("Error: Both Gyro and Ideal are False")
			return
		
		if self['Viking'] and not VikingSave:
			os.chdir(f"{directory}")
			self._save_for_save(filename = "save_info", directory = directory)
			job = open(f"save_out.job",'w')
			job.write(f"""#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --job-name={self.info['run_name']}
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --output=save_out.slurm

{save_modules}

python {directory}/save_out.py""")
			job.close()
			pyth = open(f"save_out.py",'w')
			pyth.write(f"""from Myrokinetics import myro_scan
from numpy import load
with load(\"{directory}/save_info.npz\",allow_pickle = True) as obj:
	nd = obj['name_diffs']
	info = obj['info'].item()
	run = myro_scan(eq_file = \"{self.eqbm.eq_name}\", kin_file = \"{self.eqbm.kin_name}\", input_file = \"{self.input_name}\", kinetics_type = \"{self.eqbm.kinetics_type}\", template_file = \"{self.template_name}\", directory = \"{self.path}\", run_name = \"{self.run_name}\")
	run.info = info
	run.namelist_diffs = nd
	run.save_out(filename = \"{filename}\", directory = \"{directory}\",VikingSave = True,QuickSave = {QuickSave})""")
			pyth.close()
			os.system(f"sbatch \"save_out.job\"")
			os.chdir(f"{self.path}")
			return
			
		if not self._check_setup():
			return
			
		psiNs = self['psiNs']
		beta_prime_values = full((len(psiNs)),None)
		shear_values = full((len(psiNs)),None)
		
		if self['Gyro']:
			beta_prime_axis = full((len(psiNs),self['n_beta']),None).tolist()
			shear_axis = full((len(psiNs),self['n_shat']),None).tolist()
			gr = full((len(psiNs),self['n_beta'],self['n_shat']),None).tolist()
			mf = full((len(psiNs),self['n_beta'],self['n_shat']),None).tolist()
			sym = full((len(psiNs),self['n_beta'],self['n_shat']),None).tolist()
			akys = full((len(psiNs),self['n_beta'],self['n_shat']),None).tolist()
			grs = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
			mfs = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
			syms = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
			
			if self['Epar']:
				eparN = full((len(psiNs),self['n_beta'],self['n_shat']),None).tolist()
				eparNs = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
			else:
				eparN = None
				eparNs = None
			
			if not QuickSave:
				omega = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				phi = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				apar = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				bpar = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				phi2 = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				time = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				theta = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				jacob = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
				gds2 = full((len(psiNs),self['n_beta'],self['n_shat'],len(self['aky_values'])),None).tolist()
			else:
				omega = phi = apar = bpar = phi2 = time = theta = None
		
		else:
				gr = mf = grs = mfs = sym = syms = beta_prime_axis = shear_axis = akys = self['n_shat'] = self['n_beta'] = self['aky_values'] = eparN = eparNs = omega = phi = apar = bpar = phi2 = time = theta = None
		
		if self['Ideal']:
			beta_prime_axis_ideal = full((len(psiNs),self['n_beta_ideal']),None).tolist()
			shear_axis_ideal = full((len(psiNs),self['n_shat_ideal']),None).tolist()
			stabilities = full((len(psiNs),self['n_beta_ideal'],self['n_shat_ideal']),None).tolist()
		else:
			beta_prime_axis_ideal = shear_axis_ideal = self['n_shat_ideal'] = self['n_beta_ideal'] = stabilities = None
		
		for p, psiN in enumerate(psiNs):
			run_path = os.path.join(directory,str(psiN))
			nml = f90nml.read(f"{run_path}/{psiN}.in")
			shear = nml['theta_grid_eik_knobs']['s_hat_input']
			beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
			shear_values[p] = shear
			beta_prime_values[p] = beta_prim
			
			if self['Gyro']:
				for i in range(self['n_beta']):
					per = 100*(p*self['n_beta']+i)/(self['n_beta']*len(psiNs))
					print(f"Saving {per:.2f}%", end='\r')
					for j in range(self['n_shat']):
						fol = str(i) + "_" + str(j)
						gr_aky = []
						gr_list = []
						mf_list = []
						sym_list = []
						epars = []
						for k, ky in enumerate(self['aky_values']):
							try:
							
								if self['Epar'] and not QuickSave:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob'])
								elif self['Epar']:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi','bpar'])
								elif not QuickSave:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi','apar','bpar','phi2','t','theta', 'gds2', 'jacob'])
								else:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi'])	
								
								om = data['omega'][-1,0,0]
								if type(om) != complex:
									om = data['omega'][-2,0,0]
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
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: omega")
									try:
										phi[p][i][j][k] = data['phi'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: phi")
									try:
										apar[p][i][j][k] = data['apar'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: apar")
									try:
										bpar[p][i][j][k] = data['bpar'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: bpar")
									try:
										phi2[p][i][j][k] = data['phi2'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: phi2")
									try:
										time[p][i][j][k] = data['t'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: time")
									try:
										theta[p][i][j][k] = data['theta'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: theta")
									try:
										jacob[p][i][j][k] = data['jacob'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: jacob")
									try:
										gds2[p][i][j][k] = data['gds2'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: gds2")
								
								if self['Epar']:
									epar_path = f"{run_path}/{fol}/{fol}_{k}.epar"
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
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: epar")
								
							except Exception as e:
								print(f"Save Error {psiN}/{fol}/{fol}_{k}: {e}")
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
							akys[p][i][j] = self['aky_values'][aky_idx]
							gr[p][i][j] = gr_list[aky_idx]
							mf[p][i][j] = mf_list[aky_idx]
							sym[p][i][j] = sym_list[aky_idx]
							if self['Epar']:
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
							if self['Epar']:
								eparNs[p][i][j] = None
								eparN[p][i][j] = None
							
				if self['beta_min'] is None:
					beta_min = beta_prim/self['beta_div']
				else:
					beta_min = self['beta_min']
				if self['beta_max'] is None:
					beta_max = beta_prim*self['beta_mul']
				else:
					beta_max = self['beta_max']
				for i in range(self['n_beta']):
					beta_prime_axis[p][i] = abs((beta_max - beta_min)*i/(self['n_beta']-1) + beta_min)
					
				if self['shat_min'] is None:
					shat_min = shear/self['shat_div']
				else:
					shat_min = self['shat_min']
				if self['shat_max'] is None:
					shat_max = shear*self['shat_mul']
				else:
					shat_max = self['shat_max']
				for i in range(self['n_shat']):
					sh = (shat_max - shat_min)*i/(self['n_shat']-1) + shat_min
					if sh == 0:
						sh = 1e-4
					shear_axis[p][i] = sh

			if self['Ideal']:
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
		'theta': theta,
		'jacob': jacob,
		'gds2': gds2
		}
		self.file_lines = {'eq_file': self.eqbm._eq_lines, 'kin_file': self.eqbm._kin_lines, 'template_file': self.eqbm._template_lines, 'namelist_differences': self.namelist_diffs}
		savez(f"{self.path}/{filename}", inputs = self.inputs, data = self.dat, run_info = self.info, files = self.file_lines)
		
	def rerun_errors(self, save = None, runs = None, directory = None):
		if runs is None:
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
			
			runs = self.verify.runs_with_errors

		if not self.namelist_diffs:
			self.namelist_diffs = [[[[{} for _ in range(len(self['aky_values']))] for _ in range(self['n_shat'])] for _ in range(self['n_beta'])] for _ in range(len(self['psiNs']))]
		for [p,i,j,k] in runs:
			if 'knobs' not in self.namelist_diffs[p][i][j][k].keys():
				self.namelist_diffs[p][i][j][k]['knobs'] = {}
			if 'theta_grid_parameters' not in self.namelist_diffs[p][i][j][k].keys():
				self.namelist_diffs[p][i][j][k]['theta_grid_parameters'] = {}
			self.namelist_diffs[p][i][j][k]['knobs']['nstep'] = 2*self.eqbm._template_lines['knobs']['nstep']
			self.namelist_diffs[p][i][j][k]['theta_grid_parameters']['ntheta'] = 2*self.eqbm._template_lines['theta_grid_parameters']['ntheta']
			self.namelist_diffs[p][i][j][k]['theta_grid_parameters']['nperiod'] = 2*self.eqbm._template_lines['theta_grid_parameters']['nperiod']
			if self['Fixed_delt'] is False:
				delt = 0.004/self['aky_values'][k]
				if delt > 0.01:
					delt = 0.01
				self.namelist_diffs[p][i][j][k]['knobs']['delt'] = delt
			else:
				self.namelist_diffs[p][i][j][k]['knobs']['delt'] = self.eqbm._template_lines['knobs']['delt']/10
		self.info['itteration'] += 1
		self._make_gyro_files(specificRuns = runs)
		self._run_jobs
	
	def rerun(self, runs = None, nml = None, directory = None):
		if specificRuns is None:
			print("ERROR: specificRuns not given")
			return
		if nml is None:
			print("ERROR: nml not given")
			return
		if not self.namelist_diffs:
			self.namelist_diffs = [[[[{} for _ in range(len(self['aky_values']))] for _ in range(self['n_shat'])] for _ in range(self['n_beta'])] for _ in range(len(self['psiNs']))]
		
		if type(nml) == str:
			nml = f90nml.read(nml)
		for p,i,j,k in runs:
			self.namelist_diffs[p][i][j][k] = nml
		if self.info:
			self.info['itteration'] += 1
		self._make_gyro_files(specificRuns = runs, directory = directory)
		self._run_jobs()
	
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
	
	def load_info(self, directory = None, filename = "save_info.npz"):
		from numpy import load
		if directory is None:
			directory = os.path.join(self.path, self.run_name)
		with load(f"{directory}/{filename}",allow_pickle = True) as obj:
			nd = obj['name_diffs']
			info = obj['info'].item()
			self.info = info
			self.namelist_diffs = nd
