import os
from numpy import full, real, imag, array, loadtxt, transpose, savez
from .ncdf2dict import ncdf2dict as readnc
from .equilibrium import equilibrium
from .templates import systems
from .inputs import scan_inputs
import f90nml
import glob
from uuid import uuid4
from copy import deepcopy

'''
GYROKINETIC SCAN PERFORMER
'''

class myro_scan(object):
	def __init__(self, input_file = None, directory = "./"):
		if directory == "./":
			directory = os.getcwd()
		self.path = directory
		self.dat = self.file_lines = self.verify = self.dimensions = self.namelist_diffs = self.eqbm =  None
		self._input_files = set()
		self._ideal_input_files = set()
		self.load_inputs(input_file = input_file, directory = directory)
		self.eqbm = self.equilibrium = equilibrium(inputs = self.inputs, directory = directory)
	
	def __getitem__(self, key):
		if key == "inputs":
			self.inputs()
		else:
        		return self.inputs[key]
        
	def __len__(self):
		if self.dimensions is None:
			return None
		tot = 1
		for dim in self.dimensions.values():
			tot *= len(dim)
		return tot

	def print_inputs(self):
        	self.inputs.print_inputs()
        	
	def keys(self):
		return self.inputs.keys()
	
	def load_geqdsk(self, eq_file, directory = None):
		if directory:
			self.inputs.inputs['files']['eq_path'] = directory
		self.eqbm.load_geqdsk(eq_file = eq_file, directory = directory)
	
	def load_kinetics(self, kin_file, kinetics_type = None, directory = None):
		if directory:
			self.inputs.inputs['files']['kin_path'] = directory
		self.eqbm.load_kinetics(self, kin_file = kin_file, kinetics_type = kinetics_type, directory = directory)
	
	def load_pyro(self, template_file = None, directory = None):
		if directory:
			self.inputs.inputs['files']['template_path'] = directory
		self.eqbm.load_pyro(template_file = template_file, directory = directory)
		
	def load_inputs(self, input_file = None, directory = None):
		if input_file is None:
			self.inputs = None
			return
		self.inputs = scan_inputs(input_file = input_file, directory = directory)
		self.dimensions = self.inputs.dimensions
		self.single_parameters = self.inputs.single_parameters
		if self.eqbm:
			self.eqbm.load_inputs(self.inputs)
			
	def write_scan_input(self, filename = None, directory = "./", doPrint = True):
		self.inputs.write_scan_input(filename = filename, directory = directory, doPrint = doPrint)
	
	def run_scan(self, n_jobs = None, n_par = None, n_sim = None, gyro = None, ideal = None, directory = None, group_runs = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		
		if not self.check_setup(ideal = ideal, gyro = gyro):
			return
			
		run_path = self.inputs['data_path']
		
		if not os.path.exists(run_path):
			os.mkdir(run_path)
		
		if self['gyro']:
			self.make_gyro_files(directory = run_path, group_runs = group_runs)
			self.run_jobs(n_jobs = n_jobs, n_par = n_par, n_sim = n_sim)
		if self['ideal']:
			self.make_ideal_files(directory = run_path)
			if self['gyro']:
				print("Gyro scan currently running, use run_ideal_jobs when completed to run ideal scan")
			else:
				self.run_ideal_jobs(n_jobs = n_jobs, n_par = n_par, n_sim = n_sim)
	
	def check_setup(self, ideal = None, gyro = None):
		if not self.inputs.check_scan():
			return False
		
		if gyro is None:
			gyro = self['gyro']
		elif type(gyro) == bool:	
			self.inputs['gyro'] = gyro
		else:
			print("ERROR: gyro must be boolean")
			return False
		if ideal is None:
			ideal = self['ideal']
		elif type(ideal) == bool:
			self.inputs['ideal'] = ideal
		else:
			print("ERROR: ideal must be boolean")
			return False
			
		self.inputs.create_run_info()

		if not self.inputs['eq_name'] or not self.inputs['kin_name']:
			if not self.inputs['eq_name']:
				print("ERROR: No eq_file loaded")
			if not self.inputs['kin_name']:
				print("ERROR: No kin_file loaded")
			return False
		
		if gyro and self.namelist_diffs is None:
			self.namelist_diffs = {}
		
		if not os.path.exists(self.inputs['data_path']):
			os.mkdir(self.inputs['data_path'])

		if not self.eqbm.pyro:
			self.load_pyro()
		
		if self.inputs['data_path'] != self.inputs['template_path']:
			os.system(f"cp \"{self.inputs['template_path']}/{self.inputs['template_name']}\" \"{self.inputs['data_path']}/{self.inputs['template_name']}\"")
		if self.inputs['data_path'] != self.inputs['kin_path']:
			os.system(f"cp \"{self.inputs['kin_path']}/{self.inputs['kin_name']}\" \"{self.inputs['data_path']}/{self.inputs['kin_name']}\"")
		if self.inputs['data_path'] != self.inputs['eq_path']:
			os.system(f"cp \"{self.inputs['eq_path']}/{self.inputs['eq_name']}\" \"{self.inputs['data_path']}/{self.inputs['eq_name']}\"")
		if self.inputs['input_name']:
			input_name = self.inputs['input_name']
		elif self.inputs['info']['run_name']:
			input_name = f"{self.inputs['info']['run_name']}.in"
		else:
			input_name = "myro.in"
		input_path = self.inputs['input_path'] if self.inputs['input_path'] else self.path
		self.inputs.inputs['files']['input_name'] = input_name
		self.inputs.inputs['files']['input_path'] = input_path
		self.write_scan_input(filename = input_name, directory = input_path, doPrint=False)
		if input_path != self.inputs['data_path']:
			self.write_scan_input(filename = input_name, directory = self.inputs['data_path'], doPrint=False)

		return True
	
	def clear_jobs(self):
		self._input_files = set()
	
	def run_jobs(self, n_jobs = None, n_par = None, n_sim = None):
		if self['system'] in ['viking','archer2']:
			cwd = os.getcwd()
			compile_modules = systems[self['system']]['modules']
			sbatch = "#!/bin/bash"
			for key, val in self.inputs['sbatch'].items():
				if key == 'output' and '/' not in val:
					val = f"{self.inputs['data_path']}/submit_files/{val}"
				sbatch = sbatch + f"\n#SBATCH --{key}={val}"
	
		if self['system'] == 'ypi_server':
			if n_jobs is None or n_jobs > len(self._input_files):
				n_jobs = len(self._input_files)
			while n_jobs > 0:
				for input_file in self._input_files:
					os.system(f"mpirun -np 8 gs2 \"{input_file}\"")
					self._input_files.remove(input_file)
					n_jobs -= 1
		elif self['system'] == 'viking':
			if n_par is None:
				n_par = 1
			if n_sim is None:
				n_sim = n_par
			os.makedirs(f"{self.inputs['data_path']}/submit_files/",exist_ok=True)
			input_lists = {}
			for n in range(n_par):
				input_lists[n] = []
			if n_jobs == None or n_jobs*n_par > len(self._ideal_input_files):
				total_jobs = len(self._ideal_input_files)
			else:
				total_jobs = n_jobs*n_par
			input_list = list(self._ideal_input_files)
			for i in range(total_jobs):
				input_lists[i%n_par].append(input_list[i])
				self._ideal_input_files.remove(input_list[i])
			for n in range(n_par):
				sbatch_n = sbatch.replace(f"{self.inputs['sbatch']['output']}",f"{self.inputs['sbatch']['output']}_ideal_{n}")
				sbatch_n = sbatch_n.replace(f"{self.inputs['sbatch']['error']}",f"{self.inputs['sbatch']['error']}_ideal_{n}")
				filename = f"gyro_{n}"
				pyth = open(f"{self.inputs['data_path']}/submit_files/{filename}.py",'w')
				pyth.write(f"""import os, sys
				
input_files = {input_lists[n]}

if __name__ == '__main__':
	slurm_id = int(sys.argv[1])
	input_file = input_files[slurm_id]
	os.system(f"echo \\\"Input: {{input_file}}\\\"")
	os.system(f"srun --ntasks={self.inputs['sbatch']['cpus-per-task']} \\\"{{input_file}}\\\"")
	if os.path.exists(f\"{{input_file[:-3]}}.out.nc\"):
		os.system(f"touch \\\"{{input_file[:-3]}}.fin\\\"")""")
				pyth.close()
				jobfile = open(f"{self.inputs['data_path']}/submit_files/{filename}.job",'w')
				jobfile.write(f"""{sbatch_n}
#SBATCH --array=0-{len(input_lists[n])}

{compile_modules}

which gs2
gs2 --build-config

python {self.inputs['data_path']}/submit_files/{filename}.py &

wait""")
				if n_par > n_sim and n + n_sim < n_par:
					jobfile.write(f"\nsbatch {self.inputs['data_path']}/submit_files/gyro_{n+n_sim}.job")
				jobfile.close()
			for n in range(n_sim):
				os.system(f"sbatch \"{self.inputs['data_path']}/submit_files/gyro_{n}.job\"")	
		if self['system'] == 'archer2':
			if n_par is None:
				n_par = 1
			if n_sim is None:
				n_sim = n_par if n_par < 8 else 8
			if n_sim > 8:
				print("Archer supports a maximum of n_sim = 8")
				n_sim = 8
			os.makedirs(f"{self.inputs['data_path']}/submit_files/",exist_ok=True)
			input_lists = {}
			for n in range(n_par):
				input_lists[n] = []
			if n_jobs == None or n_jobs*n_par > len(self._input_files):
				total_jobs = len(self._input_files)
			else:
				total_jobs = n_jobs*n_par
			input_list = list(self._input_files)
			for i in range(total_jobs):
				input_lists[i%n_par].append(input_list[i])
				self._input_files.remove(input_list[i])
			for n in range(n_par):
				sbatch_n = sbatch.replace(f"{self.inputs['sbatch']['output']}",f"{self.inputs['sbatch']['output']}_{n}")
				filename = f"gyro_{n}"
				pyth = open(f"{self.inputs['data_path']}/submit_files/{filename}.py",'w')
				pyth.write(f"""import os
from joblib import Parallel, delayed
from time import sleep

input_files = {input_lists[n]}

def start_run(run):
	os.system(f"echo \\\"Input: {{run}}\\\"")
	os.system(f"srun --nodes={self.inputs['sbatch']['nodes']} --ntasks={self.inputs['sbatch']['ntasks-per-node']} gs2 \\\"{{run}}\\\"")
	if os.path.exists(f\"{{run[:-3]}}.out.nc\"):
		os.system(f"touch \\\"{{run[:-3]}}.fin\\\"")
	else:
		sleep(60)
		start_run(run)

Parallel(n_jobs={self.inputs['sbatch']['nodes']})(delayed(start_run)(run) for run in input_files)""")
				pyth.close()
				jobfile = open(f"{self.inputs['data_path']}/submit_files/{filename}.job",'w')
				jobfile.write(f"""{sbatch_n}

{compile_modules}

which gs2
gs2 --build-config

python {self.inputs['data_path']}/submit_files/{filename}.py &

wait""")
				if n_par > n_sim and n + n_sim < n_par:
					jobfile.write(f"\nsbatch {self.inputs['data_path']}/submit_files/gyro_{n+n_sim}.job")
				jobfile.close()
			for n in range(n_sim):
				os.system(f"sbatch \"{self.inputs['data_path']}/submit_files/gyro_{n}.job\"")
	
	def run_ideal_jobs(self, n_jobs = None, n_par = None, n_sim = None):
		if self['system'] in ['viking','archer2']:
			cwd = os.getcwd()
			compile_modules = systems[self['system']]['modules']
			sbatch = "#!/bin/bash"
			for key, val in self.inputs['sbatch'].items():
				if key == 'output' and '/' not in val:
					val = f"{self.inputs['data_path']}/submit_files/{val}"
				sbatch = sbatch + f"\n#SBATCH --{key}={val}"
	
		if self['system'] == 'ypi_server':
			if n_jobs is None or n_jobs > len(self._input_files):
				n_jobs = len(self._input_files)
			while n_jobs > 0:
				for input_file in self._ideal_input_files:
					os.system(f"ideal_ball \"{input_file}\"")
					self._ideal_input_files.remove(input_file)
					n_jobs -= 1
		elif self['system'] == 'viking':
			if n_par is None:
				n_par = 1
			if n_sim is None:
				n_sim = n_par
			os.makedirs(f"{self.inputs['data_path']}/submit_files/",exist_ok=True)
			input_lists = {}
			for n in range(n_par):
				input_lists[n] = []
			if n_jobs == None or n_jobs*n_par > len(self._ideal_input_files):
				total_jobs = len(self._ideal_input_files)
			else:
				total_jobs = n_jobs*n_par
			input_list = list(self._ideal_input_files)
			for i in range(total_jobs):
				input_lists[i%n_par].append(input_list[i])
				self._ideal_input_files.remove(input_list[i])
			for n in range(n_par):
				sbatch_n = sbatch.replace(f"{self.inputs['sbatch']['output']}",f"{self.inputs['sbatch']['output']}_ideal_{n}")
				sbatch_n = sbatch_n.replace(f"{self.inputs['sbatch']['error']}",f"{self.inputs['sbatch']['error']}_ideal_{n}")
				sbatch_n = sbatch_n.replace(f"--cpus-per-task={self.inputs['sbatch']['cpus-per-task']}",f"--cpus-per-task=1")
				filename = f"ideal_{n}"
				pyth = open(f"{self.inputs['data_path']}/submit_files/{filename}.py",'w')
				pyth.write(f"""import os, sys
				
input_files = {input_lists[n]}

if __name__ == '__main__':
	slurm_id = int(sys.argv[1])
	input_file = input_files[slurm_id]
	os.system(f"echo \\\"Ideal Input: {{input_file}}\\\"")
	os.system(f"srun ideal_ball \\\"{{input_file}}\\\"")
	if os.path.exists(f\"{{input_file[:-3]}}.ballstab2d\"):
		os.system(f"touch \\\"{{input_file[:-3]}}.fin\\\"")""")
				pyth.close()
				jobfile = open(f"{self.inputs['data_path']}/submit_files/{filename}.job",'w')
				jobfile.write(f"""{sbatch_n}
#SBATCH --array=0-{len(input_lists[n])}

{compile_modules}

which gs2
gs2 --build-config

python {self.inputs['data_path']}/submit_files/{filename}.py $SLURM_ARRAY_TASK_ID""")
				if n_par > n_sim and n + n_sim < n_par:
					jobfile.write(f"\nsbatch {self.inputs['data_path']}/submit_files/ideal_{n+n_sim}.job")
				jobfile.close()
			for n in range(n_sim):
				os.system(f"sbatch \"{self.inputs['data_path']}/submit_files/ideal_{n}.job\"")	
		if self['system'] == 'archer2':
			if n_par is None:
				n_par = 1
			if n_sim is None:
				n_sim = n_par if n_par < 8 else 8
			if n_sim > 8:
				print("Archer supports a maximum of n_sim = 8")
				n_sim = 8
			os.makedirs(f"{self.inputs['data_path']}/submit_files/",exist_ok=True)
			input_lists = {}
			for n in range(n_par):
				input_lists[n] = []
			if n_jobs == None or n_jobs*n_par > len(self._ideal_input_files):
				total_jobs = len(self._ideal_input_files)
			else:
				total_jobs = n_jobs*n_par
			input_list = list(self._ideal_input_files)
			for i in range(total_jobs):
				input_lists[i%n_par].append(input_list[i])
				self._ideal_input_files.remove(input_list[i])
			for n in range(n_par):
				sbatch_n = sbatch.replace(f"{self.inputs['sbatch']['output']}",f"{self.inputs['sbatch']['output']}_ideal_{n}")
				sbatch_n = sbatch_n.replace(f"#SBATCH --nodes = {self.inputs['sbatch']['nodes']}",f"#SBATCH --nodes = 1")
				filename = f"ideal_{n}"
				pyth = open(f"{self.inputs['data_path']}/submit_files/{filename}.py",'w')
				pyth.write(f"""import os
from joblib import Parallel, delayed
from time import sleep

input_files = {input_lists[n]}

def start_run(run):
	os.system(f"echo \\\"Ideal Input: {{run}}\\\"")
	os.system(f"srun --nodes=1 --ntasks=1 ideal_ball \\\"{{run}}\\\"")
	if os.path.exists(f\"{{run[:-3]}}.ballstab_2d\"):
		os.system(f"touch \\\"{{run[:-3]}}.fin\\\"")
	else:
		sleep(60)
		start_run(run)

Parallel(n_jobs={self.inputs['sbatch']['ntasks-per-node']})(delayed(start_run)(run) for run in input_files)""")
				pyth.close()
				jobfile = open(f"{self.inputs['data_path']}/submit_files/{filename}.job",'w')
				jobfile.write(f"""{sbatch_n}

{compile_modules}

which gs2
gs2 --build-config

python {self.inputs['data_path']}/submit_files/{filename}.py &

wait""")
				if n_par > n_sim and n + n_sim < n_par:
					jobfile.write(f"\nsbatch {self.inputs['data_path']}/submit_files/ideal_{n+n_sim}.job")
				jobfile.close()
			for n in range(n_sim):
				os.system(f"sbatch \"{self.inputs['data_path']}/submit_files/ideal_{n}.job\"")
	
	def make_ideal_files(self, directory = None, specificRuns = None, checkSetup = True):
		if checkSetup:
			if not self.check_setup():
				return
		if directory is None:
			directory = self.inputs['data_path']
		if specificRuns:
			runs = specificRuns
		else:
			check = self.check_complete(directory = directory, doPrint = False, gyro = False, ideal = True)
			if check['ideal_complete']:
				print(f"{len(check['ideal_complete'])} Existing Ideal Runs Detected")
			runs = check['ideal_incomplete']
			
		for run in runs:
			sub_dir = self.get_ideal_run_directory(run)
			os.makedirs(sub_dir,exist_ok=True)
			
			existing_inputs = [] 
			for f in glob.glob(r'itteration_*.in'):
				existing_inputs.append([x for x in f if x.isdigit()])
			itt = max([eval("".join(x)) for x in existing_inputs],default=-1) + 1
			filename = f"itteration_{itt}"
			
			nml = self.eqbm.get_surface_input(psiN = run['psin'])
			nml['ballstab_knobs']['theta0'] = run['theta0']
			nml.write(f"{sub_dir}/{filename}.in", force=True)
			self._ideal_input_files.add(f"{sub_dir}/{filename}.in")
	
	def get_all_runs(self):
		def loop(n,variables={},runs=[]):
			if n == 0:
				return [{}]
			dim = self.dimensions[self.inputs.dim_order[len(self.dimensions)-n]]
			for val in dim.values:
				variables[dim.name] = val
				if n>1:
					loop(n=n-1,variables=variables)
				else:
					runs.append(variables.copy())
			if n == len(self.dimensions):
				return runs
		return loop(n=len(self.dimensions))
	
	def get_all_ideal_runs(self):
		runs = []
		if 'theta0' in self.dimensions:
			theta0s = self.dimensions['theta0'].values
		elif 'theta0' in self.single_parameters:
			theta0s = self.single_parameters['theta0'].values
		else:
			theta0s = [0]
		
		if 'psin' in self.dimensions:
			psins = self.dimensions['psin'].values
		else:
			psins = self.single_parameters['psin'].values
		
		for psiN in psins:
			for theta0 in theta0s:
				runs.append({'psin': psiN, 'theta0': theta0})
		return runs
				
	def make_gyro_files(self, directory = None, checkSetup = True, specificRuns = None, group_runs = None):
		if checkSetup:
			if not self.check_setup():
				return
		if directory is None:
			directory = self.inputs['data_path']
		if not specificRuns:
			check = self.check_complete(directory = directory, doPrint = False, gyro = True, ideal = False)
			if check['gyro_complete']:
				print(f"{len(check['gyro_complete'])} Existing Gyro Runs Detected")
			runs = check['gyro_incomplete']
		else:
			runs = specificRuns
		
		
			
		for run in runs:
			sub_dir = self.get_run_directory(run)
			os.makedirs(sub_dir,exist_ok=True)
			existing_inputs = [] 
			for f in glob.glob(r'itteration_*.in'):
				existing_inputs.append([x for x in f if x.isdigit()])
			itt = max([eval("".join(x)) for x in existing_inputs],default=-1)
			if itt < self.inputs['itteration']:
				filename = f"itteration_{self.inputs['itteration']}"
				subnml = self.eqbm.get_gyro_input(run = run)
				subnml.write(f"{sub_dir}/{filename}.in", force=True)
			else:
				filename = f"itteration_{itt}"
				
			self._input_files.add(f"{sub_dir}/{filename}.in")
	
	def get_run_directory(self, run):
		sub_dir = f"{self.inputs['data_path']}/gyro_files/" + "/".join([f"{name} = {run[name]:.4g}" for name in self.inputs.dim_order])
		return sub_dir
	
	def get_ideal_run_directory(self, run):
		if 'psin' not in run and 'psin' not in self.single_parameters:
			print("ERROR: psin not given")
			return None
		elif 'psin' not in run and 'psin' in self.single_parameters:
			run['psin'] = self.single_parameters['psin'].values[0]
		if 'theta0' not in run and 'theta0' not in self.single_parameters and 'theta0' not in self.dimensions:
			run['theta0'] = 0
		elif 'theta0' not in run and 'theta0' in self.single_parameters:
			run['theta0'] = self.single_parameters['theta0'].values[0]
		elif 'theta0' not in run and 'theta0' in self.dimensions:
			print("ERROR: theta0 not given")
			return None
		
		sub_dir = f"{self.inputs['data_path']}/ideal_files/" + "/".join([f"{name} = {run[name]:.4g}" for name in ['psin','theta0']])
		return sub_dir
	
	def update_itteration(self):
		self.inputs['info']['itteration'] = self.inputs['itteration'] + 1
		print(f"Updated to itteration {self.inputs['itteration']}")
	
	def create_run_info(self):
		self.inputs.create_run_info()
	
	def check_complete(self, directory = None, doPrint = True, ideal = None, gyro = None):
		if self.inputs['data_path'] is None:
			self.inputs.create_run_info()
		if directory is None:
			directory = self.inputs['data_path']
			
		if gyro is None:
			gyro = self['gyro']
		elif type(gyro) != bool:	
			print("ERROR: gyro must be boolean")
			return
		if ideal is None:
			ideal = self['ideal']
		elif type(ideal) != bool:
			print("ERROR: ideal must be boolean")
			return
		
		unfinished_gyro = []
		finished_gyro = []
		if gyro:
			for run in self.get_all_runs():
				sub_dir = self.get_run_directory(run)
				if self['system'] != 'archer2' and os.path.exists(f"{sub_dir}/itteration_0.out.nc"):
					finished_gyro.append(run)
				elif self['system'] == 'archer2' and os.path.exists(f"{sub_dir}/itteration_0.fin"):
					finished_gyro.append(run)
				else:
					unfinished_gyro.append(run)

		unfinished_ideal = []
		finished_ideal = []
		if ideal:
			for run in self.get_all_ideal_runs():
				sub_dir = self.get_ideal_run_directory(run)
				if os.path.exists(f"{sub_dir}/itteration_0.fin"):
					finished_ideal.append(run)
				else:
					unfinished_ideal.append(run)
		
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

	def _save_nml_diff(self, filename = None, directory = None):
		if filename is None:
			filename = "nml_diffs"
		if directory is None:
			directory = self.inputs['data_path']
		savez(f"{directory}/{filename}", name_diffs = self.namelist_diffs)
	
	def quick_save(self, filename = None, directory = None, SlurmSave = False):
		self.save_out(filename = filename, directory = directory, SlurmSave = SlurmSave, QuickSave = True)
	def save_out(self, filename = None, directory = None, SlurmSave = False, QuickSave = False):
		if filename is None and self.inputs['run_name'] is None:
			filename = input("Output File Name: ")
			filename = filename.split(".")[0]
		elif filename is None:
			filename = self.inputs['run_name']
			
		if self.inputs['data_path'] is None:
			self.inputs.create_run_info()
		if directory is None:
			directory = self.path
		
		if not self['gyro'] and not self['ideal']:
			print("Error: Both Gyro and Ideal are False")
			return
		
		if self['system'] in ['viking','archer2'] and not SlurmSave:
			save_modules = systems[self['system']]['save_modules']
			self._save_nml_diff()
			sbatch = "#!/bin/bash"
			for key, val in self.inputs['sbatch_save'].items():
				if key == 'output' and '/' not in val:
					val = f"{self.inputs['data_path']}/submit_files/{val}"
				sbatch = sbatch + f"\n#SBATCH --{key}={val}"
			job = open(f"{self.inputs['data_path']}/submit_files/save_out.job",'w')
			job.write(f"""{sbatch}

{save_modules}

python {self.inputs['data_path']}/submit_files/save_out.py""")
			job.close()
			pyth = open(f"{self.inputs['data_path']}/submit_files/save_out.py",'w')
			pyth.write(f"""from Myrokinetics import myro_scan
from numpy import load
with load(\"{self.inputs['data_path']}/nml_diffs.npz\",allow_pickle = True) as obj:
	nd = obj['name_diffs']
	run = myro_scan(input_file = \"{self.inputs.input_name}\", directory = \"{self.inputs['files']['input_path']}\")
	run.namelist_diffs = nd
	run.save_out(filename = \"{filename}\", directory = \"{directory}\",SlurmSave = True,QuickSave = {QuickSave})""")
			pyth.close()
			os.system(f"sbatch \"{self.inputs['data_path']}/submit_files/save_out.job\"")
			return
			
		if not self.check_setup():
			return
			
		
		psi_itt = self.single_parameters['psin'].values if 'psin' in self.single_parameters else self.dimensions['psin'].values
		equilibrium = {}
		for psiN in psi_itt:
			equilibrium[psiN] = {}
			nml = self.eqbm.get_surface_input(psiN)
			equilibrium[psiN]['shear'] = nml['theta_grid_eik_knobs']['s_hat_input']
			equilibrium[psiN]['beta_prime'] = nml['theta_grid_eik_knobs']['beta_prime_input']
		
		if self['gyro']:
			gyro_data = {}
			gyro_data['group'] = {}
			only = set({'omega','kx','ky'})
			if not QuickSave:
				only = only | set({'phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob','ql_metric_by_mode', 'phi2_by_mode'})
			#if self.inputs['epar']:
				#only = only | set({'epar'}) NOT CURRENTLY WORKING
			data_keys = ['growth_rate','mode_frequency','omega','phi','bpar','apar','epar','phi2','parity','ql_metric']
			group_keys = ['phi2_avg','t','theta', 'gds2', 'jacob']
			gyro_keys = {}
			for dim in self.dimensions.values():
				gyro_keys[dim.name] = {}
				for val in dim.values:
					gyro_keys[dim.name][val] = set()
			if self.inputs['grid_option'] == 'box':
				kxs = set()
				kys = set()
				gyro_keys['ky'] = {}
				gyro_keys['kx'] = {}
			
			runs = self.get_all_runs()
			for run in runs:
				sub_dir = self.get_run_directory(run)
				try:
					existing_inputs = [] 
					for f in glob.glob(r'itteration_*.in'):
			         			existing_inputs.append([x for x in f if x.isdigit()])
					itt = max([eval("".join(x)) for x in existing_inputs],default=0)
					run_data = readnc(f"{sub_dir}/itteration_{itt}.out.nc",only=only)	
					group_key = run_data['attributes']['id']
					gyro_data['group'][group_key] = {}
					for key in group_keys:
						gyro_data['group'][group_key][key] = None
					for xi, kx in enumerate(run_data['kx']):
						kxs.add(kx)
						for yi, ky in enumerate(run_data['ky']):
							kys.add(ky)
							run_key = str(uuid4())
							gyro_data[run_key] = deepcopy(run)
							for key in run:
								gyro_keys[key][run[key]].add(run_key)
							gyro_data[run_key]['group_key'] = group_key
							if self.inputs['grid_option'] == 'box':
								if ky not in gyro_keys['ky']:
									gyro_keys['ky'][ky] = set()
								if kx not in gyro_keys['kx']:
									gyro_keys['kx'][kx] = set()
								gyro_keys['ky'][ky].add(run_key)
								gyro_keys['kx'][kx].add(run_key)
							if 'kx' not in gyro_data[run_key]:
								gyro_data[run_key]['kx'] = kx
							if 'ky' not in gyro_data[run_key]:
								gyro_data[run_key]['ky'] = ky
							#gyro_data['nml_diffs'] = self.namelist_diffs[?]
							for key in data_keys:
								gyro_data[run_key][key] = None
							for key in only:
								try:
									key_data = run_data[key]
									
									if key == 'omega':
										om = key_data[-1,yi,xi]
										if type(om) != complex:
											om = key_data[-2,yi,xi]
										gyro_data[run_key]['growth_rate'] = imag(om)
										gyro_data[run_key]['mode_frequency'] = real(om)
										gyro_data[run_key]['omega'] = key_data[:,yi,xi].tolist()
									elif key in ['phi','apar','bpar']:
										gyro_data[run_key][key] = key_data[yi,xi,:].tolist()
										if key == 'phi':
											try:
												symsum = sum(abs(key_data[yi,xi,:] + key_data[yi,xi,::-1]))/sum(abs(key_data[yi,xi,:]))
											except:
												symsum = 1
											if  symsum > 1.5:
												gyro_data[run_key]['parity'] = 1
											elif symsum < 0.5:
												gyro_data[run_key]['parity'] = -1
											else:
												gyro_data[run_key]['parity'] = 0
									elif key in ['t','theta', 'gds2', 'jacob']:
										gyro_data['group'][group_key][key] = key_data.tolist()
									elif key in ['phi2']:
										gyro_data['group'][group_key]['phi2_avg'] = key_data.tolist()
									elif key in ['ql_metric_by_mode']:
										gyro_data[run_key]['ql_metric'] = key_data[-1,yi,xi]
									elif key in ['phi2_metric_by_mode']:
										gyro_data[run_key]['phi2'] = key_data[:,yi,xi]
									elif key in ['epar']:
										epar_path = f"{sub_dir}/itteration_{itt}.epar"
								
										bpar = key_data['bpar'][yi,xi,:]
										epar_data = loadtxt(epar_path)
										epar = []
										for l in range(len(epar_data[:,3])):
											epar.append(complex(epar_data[l,3],epar_data[l,4]))
										epar = array(epar)
										gyro_data[run_key]['epar'] = epar
								except Exception as e:
									print(f"Save Error in {sub_dir}/itteration_{itt}: {e}")
									if key == 'omega':
										gyro_data[run_key]['growth_rate'] = nan
										gyro_data[run_key]['mode_frequency'] = nan
										
				except Exception as e:
					print(f"Save Error {sub_dir}/itteration_{itt}: {e}")
			if self.inputs['grid_option'] == 'box':
				existing_dim_keys = []
				for key in [x for x in self.inputs.inputs.keys() if 'dimension_' in x]:
	         			existing_dim_keys.append([x for x in key if x.isdigit()])
				dim_n = max([eval("".join(x)) for x in existing_dim_keys],default=1) + 1
				kxs = list(kxs)
				kxs.sort()
				self.inputs.inputs[f'dimension_{dim_n}'] = {'type': 'kx', 'values': kxs, 'min': min(kxs), 'max': max(kxs), 'num': len(kxs), 'option': None}
				kys = list(kys)
				kys.sort()
				self.inputs.inputs[f'dimension_{dim_n+1}'] = {'type': 'ky', 'values': kys, 'min': min(kys), 'max': max(kys), 'num': len(kys), 'option': None}
				self.inputs.load_dimensions()
		else:
			gyro_data = None
			gyro_keys = None

		if self['ideal']:
			ideal_keys = {}
			if 'theta0' in self.single_parameters:
				theta0_itt = self.single_parameters['theta0'].values  
			if 'theta0' in self.dimensions:
				theta0_itt = self.dimensions['theta0'].values
			else:
				theta0_itt = [0]
			
			ideal_keys['psin'] = {}
			ideal_keys['theta0'] = {}
			for val in psi_itt:
				ideal_keys['psin'][val] = set()
			for val in theta0_itt:
				ideal_keys['theta0'][val] = set()

			ideal_data = {}
			for run in self.get_all_ideal_runs():
				run_id = str(uuid4())
				for key in run:
					ideal_keys[key][run[key]].add(run_id)
				ideal_data[run_id] = {}
				try:
					sub_dir = get_ideal_run_directory(run)
					existing_inputs = [] 
					for f in glob.glob(r'itteration_*.in'):
						existing_inputs.append([x for x in f if x.isdigit()])
					itt = max([eval("".join(x)) for x in existing_inputs],default=0)

					shear = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_shat")
					bp = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_bp")
					stab = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_2d")
					
					ideal_data[run_id]['beta_prime'] = [abs(x) for x in bp]
					ideal_data[run_id]['shear'] = shear.tolist()
					ideal_data[run_id]['stabilities'] = transpose(stab).tolist()
				except:
					ideal_data[run_id]['beta_prime'] = None
					ideal_data[run_id]['shear'] = None
					ideal_data[run_id]['stabilities'] = None
					print(f"Save Error for ideal run: {run}")
		else:
			ideal_data = None
			ideal_keys = None
		
		data = {'gyro': gyro_data,
			'ideal': ideal_data,
			'equilibrium': equilibrium,
			'_gyro_keys': gyro_keys,
			'_ideal_keys': ideal_keys,
			}
		
		self.file_lines = {'eq_file': self.eqbm._eq_lines, 'kin_file': self.eqbm._kin_lines, 'template_file': self.eqbm._template_lines}
		savez(f"{directory}/{filename}", inputs = self.inputs.inputs, data = data, files = self.file_lines)
	
	'''
	def rerun(self, runs = None, nml = None, directory = None, group_runs = None):
		if runs is None:
			print("ERROR: runs not given")
			return
		if nml is None:
			print("ERROR: nml not given, if you wish to rerun with no changes please use nml = {}")
			return
			
		self.check_setup()
		
		if type(nml) == str:
			nml = f90nml.read(nml)
		for p,i,j,k,t in runs:
			self.namelist_diffs[p][i][j][k][t] = nml
		self.inputs.inputs['itteration'] += 1
		self.make_gyro_files(specificRuns = runs, directory = directory, group_runs = group_runs)
		self.run_jobs()
		
	def load_run_set(self, filename = None):
		if filename is None:
			print("ERROR: filename not given")
			return
		
		runs = set()		
		with open(filename) as f:
			lines = f.readlines()
			for line in lines:
				p,i,j,k,t = [eval(x) for x in line.strip("\n").split("_")]
				runs.add((p,i,j,k,t))
		return runs
	'''
