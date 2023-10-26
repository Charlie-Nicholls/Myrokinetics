import os
from numpy import full, real, imag, array, loadtxt, transpose, savez
from .ncdf2dict import ncdf2dict as readnc
from .equilibrium import equilibrium
from .templates import systems
from .inputs import scan_inputs
import f90nml
import glob

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
		
		if self['ideal']:
			self.make_ideal_files(directory = run_path)
		if self['gyro']:
			self.make_gyro_files(directory = run_path, group_runs = group_runs)
		self.run_jobs(n_jobs = n_jobs, n_par = n_par, n_sim = n_sim)
	
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
			compile_modules = systems[self['system']]['compile']
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
			if n_jobs is None or n_jobs > len(self._input_files):
				n_jobs = len(self._input_files)
			while n_jobs > 0:
				for input_file in self._input_files:
					input_file = self._input_files[i]
					job_name = input_file[:-3] + ".job"
					directory = "".join([f"{x}/" for x in input_file.split('/')[:-1]])
					os.system(f"cd {directory}")
					jobfile = open(f"{job_name}",'w')
					jobfile.write(f"""{sbatch}

{compile_modules}

which gs2
gs2 --build-config

echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"
echo \"Input: {sub_dir}/{filename}.in\"
gs2 \"{sub_dir}/{filename}.in\"
echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"""")
					jobfile.close()
					
					os.system(f"sbatch \"{job_name}\"")
					os.system(f"cd {cwd}")
					self._input_files.remove(input_file)
					n_jobs -= 1
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
	os.system(f"srun --nodes=1 --ntasks={self.inputs['sbatch']['ntasks-per-node']} gs2 \\\"{{run}}\\\"")
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
			
		for psiN in runs:
			sub_dir = f"{directory}/ideal/psin = {psiN:.4g}"
			os.makedirs(sub_dir,exist_ok=True)
			
			existing_inputs = [] 
			for f in glob.glob(r'itteration_*.in'):
				existing_inputs.append([x for x in f if x.isdigit()])
			itt = max([eval("".join(x)) for x in existing_inputs],default=-1) + 1
			filename = f"itteration_{itt}"
			
			nml = self.eqbm.get_surface_input(psiN = psiN)
			nml.write(f"{sub_dir}/{filename}.in", force=True)
			if self['system'] == 'ypi_server':
				self._input_files.append(f"ideal_ball \"{sub_dir}/{filename}.in\"")
			else:
				compile_modules = systems[self['system']]['compile']
				jobfile = open(f"{sub_dir}/{filename}.job",'w')
				jobfile.write(f"""#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --job-name={self.inputs['run_name']}
#SBATCH --ntasks=1
#SBATCH --output={sub_dir}/{filename}.slurm

{compile_modules}

which gs2
gs2 --build-config

ideal_ball \"{sub_dir}/{filename}.in\"""")
				jobfile.close()
				self._input_files.append(f"sbatch \"{sub_dir}/{filename}.job\"")
	
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
		sub_dir = f"{self.inputs['data_path']}/run_files/" + "/".join([f"{name} = {run[name]:.4g}" for name in self.inputs.dim_order])
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
			for psiN in self.dimensions['psin'].values:
				if os.path.exists(f"{directory}/ideal/psin = {psiN}/itteration_0.ballstab_2d"):
					finished_ideal.append(psiN)
				else:
					unfinished_ideal.append(psiN)
		
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
			directory = self.path
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
			directory = self.inputs['data_path']
		
		if not self['gyro'] and not self['ideal']:
			print("Error: Both Gyro and Ideal are False")
			return
		
		if self['system'] in ['viking','archer2'] and not SlurmSave:
			save_modules = systems[self['system']]['save']
			os.chdir(f"{directory}")
			self._save_nml_diff(directory = directory)
			job = open(f"save_out.job",'w')
			if self['system'] == 'viking':
				job.write(f"""#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --job-name={self.inputs['run_name']}
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --output={directory}/save_out.slurm

{save_modules}

python {directory}/save_out.py""")
			if self['system'] == 'archer2':
				job.write(f"""#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --job-name={self.inputs['run_name']}
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --output={directory}/save_out.slurm
#SBATCH --account={self.inputs['sbatch']['account']}

{save_modules}

python {directory}/save_out.py &
wait""")
			job.close()
			pyth = open(f"save_out.py",'w')
			pyth.write(f"""from Myrokinetics import myro_scan
from numpy import load
with load(\"{directory}/nml_diffs.npz\",allow_pickle = True) as obj:
	nd = obj['name_diffs']
	run = myro_scan(input_file = \"{self.inputs.input_name}\", directory = \"{self.path}\")
	run.namelist_diffs = nd
	run.save_out(filename = \"{filename}\", directory = \"{directory}\",SlurmSave = True,QuickSave = {QuickSave})""")
			pyth.close()
			os.system(f"sbatch \"save_out.job\"")
			os.chdir(f"{self.path}")
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
			only = set({'omega'})
			if not QuickSave:
				only = only | set({'phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob'})
			all_keys = ['omega','phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob']
			
			run_keys = {}
			for dim in self.dimensions.values():
				run_keys[dim.name] = {}
				for val in dim.values:
					run_keys[dim.name][val] = set()
				
			runs = self.get_all_runs()
			for run in runs:
				sub_dir = self.get_run_directory(run)
				try:
					existing_inputs = [] 
					for f in glob.glob(r'itteration_*.in'):
			         			existing_inputs.append([x for x in f if x.isdigit()])
					itt = max([eval("".join(x)) for x in existing_inputs],default=0)
					run_data = readnc(f"{sub_dir}/itteration_{itt}.out.nc",only=only)	
					
					run_key = run_data['attributes']['id']
					for key in run:
						run_keys[key][run[key]].add(run_key)
					
					gyro_data[run_key] = run
					#gyro_data['nml_diffs'] = self.namelist_diffs[?]
					for key in all_keys:
						gyro_data[run_key][key] = None
						
					for key in only:
						try:
							key_data = run_data[key]
							
							if key == 'omega':
								om = key_data[-1,0,0]
								if type(om) != complex:
									om = key_data[-2,0,0]
								gyro_data[run_key]['growth_rate'] = imag(om)
								gyro_data[run_key]['mode_frequency'] = real(om)
								gyro_data[run_key]['omega'] = key_data[:,0,0].tolist()
							elif key in ['phi','apar','bpar']:
								gyro_data[run_key][key] = key_data[0,0,:].tolist()
								if key == 'phi':
									symsum = sum(abs(key_data[0,0,:] + key_data[0,0,::-1]))/sum(abs(key_data[0,0,:]))
									if  symsum > 1.5:
										gyro_data[run_key]['parity'] = 1
									elif symsum < 0.5:
										gyro_data[run_key]['parity'] = -1
									else:
										gyro_data[run_key]['parity'] = 0
							elif key in ['phi2','t','theta', 'gds2', 'jacob']:
								gyro_data[run_key][key] = key_data.tolist()
							elif key in ['epar']:
								epar_path = f"{sub_dir}/itteration_{itt}.epar"
						
								bpar = key_data['bpar'][0,0,:]
								epar_data = loadtxt(epar_path)
								epar = []
								for l in range(len(epar_data[:,3])):
									epar.append(complex(epar_data[l,3],epar_data[l,4]))
								epar = array(epar)
								gyro_data[run_key]['epar'] = epar
						except:
							print(f"Save Error in {sub_dir}/itteration_{itt}: {key}")
							if key == 'omega':
								gyro_data[run_key]['growth_rate'] = nan
								gyro_data[run_key]['mode_frequency'] = nan
								
				except Exception as e:
					print(f"Save Error {sub_dir}/itteration_{itt}: {e}")
		else:
			gyro_data = None

		if self['ideal']:
			ideal_data = {}
			for psiN in psi_itt:
				try:
					ideal_data[psiN] = {}
					sub_dir = f"{directory}/ideal/psin = {psiN:.4g}"
					existing_inputs = [] 
					for f in glob.glob(r'itteration_*.in'):
						existing_inputs.append([x for x in f if x.isdigit()])
					itt = max([eval("".join(x)) for x in existing_inputs],default=0)

					shear = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_shat")
					bp = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_bp")
					stab = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_2d")
					
					ideal_data[psiN]['beta_prime'] = [abs(x) for x in bp]
					ideal_data[psiN]['shear'] = shear.tolist()
					ideal_data[psiN]['stabilities'] = transpose(stab).tolist()
				except:
					ideal_data[psiN]['beta_prime'] = None
					ideal_data[psiN]['shear'] = None
					ideal_data[psiN]['stabilities'] = None
					print(f"Save Error for ideal psiN = {psiN}")
		else:
			ideal_data = None
		
		data = {'gyro': gyro_data,
			'ideal': ideal_data,
			'equilibrium': equilibrium,
			'_run_keys': run_keys,
			}
		
		self.file_lines = {'eq_file': self.eqbm._eq_lines, 'kin_file': self.eqbm._kin_lines, 'template_file': self.eqbm._template_lines}
		savez(f"{self.path}/{filename}", inputs = self.inputs.inputs, data = data, files = self.file_lines)
	'''
	def check_cancelled(self, directory = None, doPrint = True):
		if self.inputs['data_path'] is None:
			self.create_run_info()		
		if directory is None:
			directory = self.inputs['data_path']
		if not self['gyro']:
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
						if len(ids) == 4:
							i, j, k, t = [eval(x) for x in ids]
							cancelled.add((p,i,j,k,t))
						elif len(ids) == 3:
							f = open(line.strip("\n"))
							lins = f.readlines()
							f.close()
							for l in lins:
								if ".in" in l:
									inp = l.split("/")[-1].split(".")[0]
							i, j, k, t = [eval(x) for x in inp.split("_")]
							for ti in range(t, self['n_theta0']):
								cancelled.add((p,i,j,k,ti))
		if doPrint:		
			print(f"{len(cancelled)} Cancelled Runs")
			return
		else:
			return cancelled
	
	def rerun_cancelled(self, directory = None, checkSetup = True, group_runs = None):
		if self.inputs['data_path'] is None:
			self.create_run_info()
		if directory is None:
			directory = self.inputs['data_path']
		cancelled = self.check_cancelled(directory = directory, doPrint = False)
		if len(cancelled) == 0:
			print("No Cancelled Runs")
			return
		for p, i, j, k, t in cancelled:
			if os.path.exists(f"{directory}/{self['psiNs'][p]}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc"):
				os.remove(f"{directory}/{self['psiNs'][p]}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc")
		self.make_gyro_files(directory = directory, specificRuns = cancelled, group_runs = group_runs)
		self.run_jobs()
	
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
	
	def load_info(self, directory = None, filename = "save_info.npz"):
		from numpy import load
		if directory is None:
			directory = os.path.join(self.path, self.inputs['run_name'])
		with load(f"{directory}/{filename}",allow_pickle = True) as obj:
			nd = obj['name_diffs']
			info = obj['info'].item()
			self.inputs = info
			self.namelist_diffs = nd
