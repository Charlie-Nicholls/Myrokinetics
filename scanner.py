import os
from numpy import real, imag, array, loadtxt, transpose, savez, nan, ceil, trapz
from .ncdf2dict import ncdf2dict as readnc
from .equilibrium import equilibrium
from .codes import codes
from .systems import systems
from .inputs import scan_inputs
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
		self._jobs = set()
		self._ideal_input_files = set()
		self._ideal_jobs = set()
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
		if self.inputs['scan_format'] == 'grid':
			tot = 1
			for dim in self.dimensions.values():
				tot *= len(dim)
		elif self['scan_format'] == 'point':
			return len(self.get_all_runs())
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
	
	def write_scan_input_copy(self, filename = None, directory = "./", doPrint = True):
		self.inputs.write_scan_input_copy(filename = filename, directory = directory, doPrint = doPrint)
	
	def run_scan(self, n_jobs = None, n_par = None, n_sim = None, gyro = None, ideal = None, directory = None, specificRuns = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		
		if gyro is None:
			gyro = self['gyro']
		if type(gyro) != bool:	
			print("ERROR: gyro must be boolean")
			return
		if ideal is None:
			ideal = self['ideal']
		if type(ideal) != bool:
			print("ERROR: ideal must be boolean")
			return
		
		if not self.check_setup():
			return
			
		run_path = self.inputs['data_path']
		
		if not os.path.exists(run_path):
			os.mkdir(run_path)
		
		if gyro:
			self.make_gyro_files(directory = run_path, specificRuns = specificRuns)
			self.make_job_files(n_jobs = n_jobs, n_par = n_par, n_sim = n_sim)
			self.run_jobs()
		if ideal:
			self.make_ideal_files(directory = run_path)
			if gyro:
				self.make_ideal_job_files()
			else:
				self.make_ideal_job_files(n_jobs = n_jobs, n_par = n_par, n_sim = n_sim)
			self.run_ideal_jobs()
	
	def check_setup(self, gyro = None, ideal = None):
		if not self.inputs.check_scan():
			return False
		
		if gyro is None:
			gyro = self['gyro']
		if ideal is None:
			ideal = self['ideal']
		
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
		self._ideal_input_files = set()

	def make_job_files(self, n_jobs = None, n_par = None, n_sim = None):
		systems[self.inputs['system']].make_job_files(self, n_jobs=n_jobs, n_par=n_par, n_sim=n_sim)
		
	def make_ideal_job_files(self, n_jobs = None, n_par = None, n_sim = None):
		systems[self.inputs['system']].make_ideal_job_files(self, n_jobs=n_jobs, n_par=n_par, n_sim=n_sim)
	
	def run_jobs(self):
		cwd = os.getcwd()
		os.chdir(f"{self.inputs['data_path']}/submit_files")
		for job in self._jobs:
			os.system(f"sbatch {job}")
		os.chdir(cwd)
		self._jobs = set()
	
	def run_ideal_jobs(self):
		cwd = os.getcwd()
		os.chdir(f"{self.inputs['data_path']}/submit_files")
		for job in self._ideal_jobs:
			os.system(f"sbatch {job}")
		os.chdir(cwd)
		self._ideal_jobs = set()
	
	def restart_run(self, run = {}, itt = None):
		import f90nml
		if self.inputs['grid_option'] != 'box':
			print("ERROR: restart_run only supported for grid_option = True")
			return
		if run not in self.get_all_runs(excludeDimensions = ['kx','ky']):
			print("ERROR: run not found")
			return
		file_dir = self.get_run_directory(run)
		if itt is None:
			itt = self['itteration']
			if not os.path.exists(f"{file_dir}/itteration_{itt}.out.nc"):
				print(f"ERROR: itteration {itt} not found, please specify itt")
				return
		nml = f90nml.read(f"{file_dir}/itteration_{itt}.in")
		nml['knobs']['delt_option'] = 'check_restart'
		h, m, s = self.inputs['sbatch']['time'].split(':')
		time = (int(h) * 3600) + (int(m) * 60) + int(s)
		nml['knobs']['avail_cpu_time'] = time
		nml['knobs']['margin_cpu_time'] = time // 20
		nml['knobs']['delt_option'] = 'check_restart'
		nml['init_g_knobs']['ginit_option'] = 'many'
		nml['gs2_diagnostics_knobs']['append_old'] = True
		nml.write(f"{file_dir}/itteration_{itt}.in",force=True)
		self._input_files.add(f"{file_dir}/itteration_{itt}.in")
		self.make_job_files()
		self.run_jobs()
	
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
		if self.inputs['scan_format'] == 'grid':
			return loop(n=len(dim_order))
		elif self.inputs['scan_format'] == 'point':
			return self.load_run_set(os.path.join(self.inputs['point_path'],self.inputs['point_name']))
	
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
				
	def make_gyro_files(self, directory = None, checkSetup = True, specificRuns = None):
		self._input_files = set()
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
			codes[self.inputs['gk_code']].make_gyro_file(eqbm=self.eqbm,run=run,sub_dir=sub_dir)
	
	def make_ideal_files(self, directory = None, specificRuns = None, checkSetup = True):
		self._ideal_input_files = set()
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
	
	def get_run_directory(self, run):
		dims = [x for x in self.inputs.dim_order if x not in ['kx','ky']] if self.inputs['nonlinear'] == True else self.inputs.dim_order
		dir_list = [f"{name}={run[name]:.4g}" for name in dims]
		dir_list.insert(0,f"{self.inputs['data_path']}/gyro_files")
		sub_dir = "/".join(dir_list)
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
		
		sub_dir = f"{self.inputs['data_path']}/ideal_files/" + "/".join([f"{name}={run[name]:.4g}" for name in ['psin','theta0']])
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
		all_runs = self.get_all_runs(excludeDimensions=['kx','ky']) if self.inputs['nonlinear'] else self.get_all_runs()
		if gyro:
			for run in all_runs:
				sub_dir = self.get_run_directory(run)
				if os.path.exists(f"{sub_dir}/run.fin"):
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
	
	def save_out(self, filename = None, directory = None, specificRuns = None, SlurmSave = False, QuickSave = False):
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
		
		if not scanner.check_setup():
			return
		
		if systems[scanner.inputs['system']].requires_slurm_save and not SlurmSave:
			save_modules = systems[scanner['system']].save_modules
			scanner._save_nml_diff()
			sbatch = "#!/bin/bash"
			for key, val in scanner.inputs['sbatch_save'].items():
				if key == 'output' and '/' not in val:
					val = f"{scanner.inputs['data_path']}/submit_files/{val}"
				sbatch = sbatch + f"\n#SBATCH --{key}={val}"
			job = open(f"{scanner.inputs['data_path']}/submit_files/save_out.job",'w')
			job.write(f"""{sbatch}

{save_modules}

python {scanner.inputs['data_path']}/submit_files/save_out.py""")
			job.close()
			pyth = open(f"{scanner.inputs['data_path']}/submit_files/save_out.py",'w')
			pyth.write(f"""from Myrokinetics import myro_scan
from numpy import load
specificRuns = {specificRuns}
with load(\"{scanner.inputs['data_path']}/nml_diffs.npz\",allow_pickle = True) as obj:
	nd = obj['name_diffs']
	run = myro_scan(input_file = \"{scanner.inputs.input_name}\", directory = \"{scanner.inputs['files']['input_path']}\")
	run.namelist_diffs = nd
	run.save_out(filename = \"{filename}\", directory = \"{directory}\", specificRuns = specificRuns, SlurmSave = True, QuickSave = {QuickSave})""")
			pyth.close()
			os.system(f"sbatch \"{scanner.inputs['data_path']}/submit_files/save_out.job\"")
			return
		
		data = codes[self.inputs['gk_code']].save_out(scanner, filename = filename, directory = directory, specificRuns = specificRuns, QuickSave = QuickSave)
		
		self.file_lines = {'eq_file': self.eqbm._eq_lines, 'kin_file': self.eqbm._kin_lines, 'template_file': self.eqbm._template_lines}
		savez(f"{directory}/{filename}", inputs = self.inputs.inputs, data = data, files = self.file_lines)
		return
	
	def print_run_input(self, run = {}, itt = None):
		import f90nml
		if run not in self.get_all_runs():
			print("ERROR: run not found")
			return
		file_dir = self.get_run_directory(run)
		if self.inputs['gk_code'] == 'GS2':
			if itt is None:
				itt = self['itteration']
				filepath = f"{file_dir}/itteration_{itt}.in"
				if not os.path.exists(filepath):
					print(f"ERROR: input \"{filepath}\" not found, please specify itt and ensure make files has been run")
					return
			nml = f90nml.read(filepath)
			print(nml)
		elif self.inputs['gk_code'] == 'CGYRO':
				with open(f"{file_dir}/input.cgyro") as f:
					lines = f.readlines()
					for line in lines:
						print(line)
		
	
	def print_submit_file(self, n = 0):
		filepath = f"{self['data_path']}/submit_files/"
		if self.inputs['nonlinear'] == True:
			filepath += "submit.job"
		else:
			filepath += f"gyro_{n}/gyro_{n}.job"
		if os.path.exists(filepath):
			sfile = open(filepath)
		else:
			print(f"ERROR: submit file \"{filepath}\" not found")
			return
		slines = sfile.readlines()
		sfile.close()
		for line in slines:
			print(line, end='')

	def print_ideal_submit_file(self, n = 0):
		filepath = f"{self['data_path']}/submit_files/ideal_{n}/ideal_{n}.job"
		if os.path.exists(filepath):
			sfile = open(filepath)
		else:
			print(f"ERROR: submit file \"{filepath}\" not found")
			return
		slines = sfile.readlines()
		sfile.close()
		for line in slines:
			print(line)
	
	def print_slurm(self, n = 0, i = 1):
		filepath = f"{self['data_path']}/submit_files/"
		if self.inputs['nonlinear']:
			filepath += f"{self.inputs['sbatch']['output']}"
		elif self.inputs['system'] == 'viking':
			filepath += f"gyro_{n}/{self.inputs['sbatch']['output']}_{i}"
		elif self.inputs['system'] == 'archer2':
			filepath += f"gyro_{n}/{self.inputs['sbatch']['output']}_{n}"
		if os.path.exists(filepath):
			sfile = open(filepath)
		else:
			print(f"ERROR: slurm file \"{filepath}\" not found")
			return
		slines = sfile.readlines()
		sfile.close()
		for line in slines:
			print(line, end='')
		
	def print_ideal_slurm(self, n = 0, i = 1):
		filepath = f"{self['data_path']}/submit_files/"
		if self.inputs['nonlinear']:
			filepath += f"{self.inputs['sbatch']['output']}_{n}"
		else:
			filepath += f"ideal_{n}/{self.inputs['sbatch']['output']}_{i}"
		if os.path.exists(filepath):
			sfile = open(filepath)
		else:
			print(f"ERROR: slurm file \"{filepath}\" not found")
			return
		slines = sfile.readlines()
		sfile.close()
		for line in slines:
			print(line, end='')
	
	def run_ingen(self, run = {}, itt = None):
		if run not in self.get_all_runs():
			print("ERROR: run not found")
			return
		file_dir = self.get_run_directory(run)
		if self.inputs['gk_code'] == 'GS2':
			if itt is None:
				itt = self['itteration']
				filepath = f"{file_dir}/itteration_{itt}.in"
				if not os.path.exists(filepath):
					print(f"ERROR: input \"{filepath}\" not found")
					return
			os.system(filepath)
			self.print_ingen(run = run, itt = itt)
		elif self.inputs['gk_code'] == 'CGYRO' and self.inputs['system'] == 'archer2':
			f = open(f"{file_dir}/ingen.job",'w')
			f.write(f"""#!/bin/bash
#SBATCH --time=00:00:30
#SBATCH --job-name=ingen
#SBATCH --nodes=1
#SBATCH --output=ingen.slurm
#SBATCH --account={self.inputs['sbatch_save']['account']}
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --distribution=block:block
#SBATCH --hint=nomultithread

{systems[self['system']]['modules']['CGYRO']}

cgyro -i "./" >& ingen.out
""")
			f.close()
			cwd = os.getcwd()
			os.chdir(file_dir)
			os.system(f"sbatch {file_dir}/ingen.job")
			os.chdir(cwd)
	
	def print_ingen(self, run = {}, itt = None):
		if run not in self.get_all_runs():
			print("ERROR: run not found")
			return
		file_dir = self.get_run_directory(run)
		if self.inputs['gk_code'] == 'GS2':
			if itt is None:
				itt = self['itteration']
			filepath = f"{file_dir}/itteration_{itt}.report"
		elif self.inputs['gk_code'] == 'CGYRO':
			filepath = f"{file_dir}/ingen.out"
		if os.path.exists(filepath):
			sfile = open(filepath)
		else:
			print(f"ERROR: report \"{filepath}\" not found, please specify itt and ensure ingen has been run")
			return
		slines = sfile.readlines()
		sfile.close()
		for line in slines:
			print(line, end='')
			
	def load_run_out(self, run = {}, itt = None):
		if run not in self.get_all_runs():
			print("ERROR: run not found")
			return
		file_dir = self.get_run_directory(run)
		if itt is None:
			itt = self['itteration']
		filepath = f"{file_dir}/itteration_{itt}.out.nc"
		if os.path.exists(filepath):
			sfile = readnc(filepath)
		else:
			print(f"ERROR: report \"{filepath}\" not found, please specify itt and ensure simulation has run")
			return
		return sfile

	def load_run_set(self, filename = None):
		if filename is None:
			print("ERROR: filename not given")
			return
		
		runs = []		
		with open(filename) as f:
			lines = f.readlines()
			for line in lines:
				run = eval(line.strip("\n"))
				runs.append(run)
		return runs
		
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
	'''
