import os
from numpy import full, real, imag, array, loadtxt, transpose, savez
from .ncdf2dict import ncdf2dict as readnc
from .equillibrium import equillibrium
from .templates import modules, save_modules
from .inputs import scan_inputs
import f90nml
import glob

'''
GYROKINETIC SCAN PERFORMER
'''

class myro_scan(object):
	def __init__(self, input_file = None, eq_file = None, kin_file = None, template_file = None, kinetics_type = "PEQDSK", directory = "./", run_name = None):
		if directory == "./":
			directory = os.getcwd()
		self.path = directory
		self.template_name = template_file
		self._template_path = self.info = self.dat = self.file_lines = self.verify = self.dimensions = self.namelist_diffs = self.eqbm =  None
		self.jobs = []
		if run_name:
			self.run_name = run_name
		elif input_file:
			self.run_name = input_file.split("/")[-1].split(".")[0]
		
		self.load_inputs(input_file = input_file, directory = directory)
		self.eqbm = self.equillibrium = equillibrium(eq_file = eq_file, kin_file = kin_file, kinetics_type = kinetics_type, template_file = template_file, directory = directory, inputs = self.inputs)
	
	def __getitem__(self, key):
		if key == "inputs":
			self.inputs()
		else:
        		return self.inputs[key]

	def print_inputs(self):
        	self.inputs.print_inputs()
        	
	def keys(self):
		return self.inputs.keys()
	
	def load_geqdsk(self, eq_file, directory = None):
		self.eqbm.load_geqdsk(eq_file = eq_file, directory = directory)
	
	def load_kinetics(self, kin_file, kinetics_type = None, directory = None):
		self.eqbm.load_kinetics(self, kin_file = kin_file, kinetics_type = kinetics_type, directory = directory)
	
	def load_pyro(self, template_file = None, directory = None):
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
	
	def run_scan(self, gyro = None, ideal = None, directory = None, group_runs = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		
		if not self.info:
			self._create_run_info()
		else:
			self.info['itteration'] += 1
		run_path = self.info['data_path']
		
		if not os.path.exists(run_path):
			os.mkdir(run_path)
		
		if not self._check_setup(ideal = ideal, gyro = gyro):
			return
		if self['Ideal']:
			self._make_ideal_files(directory = run_path)
		if self['Gyro']:
			self._make_gyro_files(directory = run_path, group_runs = group_runs)
		self._run_jobs()
	
	def _check_setup(self, ideal = None, gyro = None):
		if not self.inputs.check_scan():
			return False
		
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
			if not self.eqbm.kin_name:
				print("ERROR: No kin_file loaded")
			return False
		
		if gyro and self.namelist_diffs is None:
			self.namelist_diffs = full(self.inputs.shape,{})
		
		if not os.path.exists(self.info['data_path']):
			os.mkdir(self.info['data_path'])

		if not self.eqbm.pyro:
			self.eqbm.load_pyro()

		if self.template_name is not None:
			os.system(f"cp \"{self.eqbm._template_path}/{self.template_name}\" \"{self.info['data_path']}/{self.template_name}\"")
		os.system(f"cp \"{self.eqbm._kin_path}/{self.eqbm.kin_name}\" \"{self.info['data_path']}/{self.eqbm.kin_name}\"")
		os.system(f"cp \"{self.eqbm._eq_path}/{self.eqbm.eq_name}\" \"{self.info['data_path']}/{self.eqbm.eq_name}\"")
		if not self.inputs.input_name:
			self.inputs.input_name = f"{self.run_name}.in"
		self.write_scan_input(filename = self.inputs.input_name, directory = self.info['data_path'],doPrint=False)

		return True
	
	def _run_jobs(self, n = None):
		if n is None or n > len(self.jobs):
			n = len(self.jobs)
		for i in range(n):
			os.system(self.jobs[i])
		self.jobs = self.jobs[n:]
	
	def _make_ideal_files(self, directory = None, specificRuns = None, checkSetup = True):
		if checkSetup:
			if not self._check_setup():
				return
		if directory is None:
			directory = self.info['data_path']
		if specificRuns:
			runs = specificRuns
		else:
			'''
			check = self.check_complete(directory = directory, doPrint = False, gyro = False, ideal = True)
			if check['ideal_complete']:
				print(f"{len(check['ideal_complete'])} Existing Ideal Runs Detected")
			runs = check['ideal_incomplete']
			'''
			runs = self.dimensions['psin'].values
			
		for psiN in runs:
			sub_dir = f"{directory}/ideal/psin = {psiN:.2g}"
			os.makedirs(sub_dir,exist_ok=True)
			
			existing_inputs = [] 
			for f in glob.glob(r'itteration_*.in'):
				existing_inputs.append([x for x in f if x.isdigit()])
			itt = max([eval("".join(x)) for x in existing_inputs],default=-1) + 1
			filename = f"itteration_{itt}"
			
			nml = self.eqbm.get_surface_input(psiN = psiN)
			nml.write(f"{sub_dir}/{filename}.in", force=True)
			if self['System'] == 'plasma':
				self.jobs.append(f"ideal_ball \"{sub_dir}/{filename}.in\"")
			else:
				jobfile = open(f"{sub_dir}/{filename}.job",'w')
				jobfile.write(f"""#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --job-name={self.info['run_name']}
#SBATCH --ntasks=1
#SBATCH --output={sub_dir}/{filename}.slurm

{modules}

which gs2
gs2 --build-config

ideal_ball \"{sub_dir}/{filename}.in\"""")
				jobfile.close()
				self.jobs.append(f"sbatch \"{sub_dir}/{filename}.job\"")
	
	def get_all_runs(self):
		def loop(n,variables={},runs=[]):
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
				
	def _make_gyro_files(self, directory = None, checkSetup = True, specificRuns = None, group_runs = None):
		if checkSetup:
			if not self._check_setup():
				return
		if directory is None:
			directory = self.info['data_path']
		if not specificRuns:
			'''
			check = self.check_complete(directory = directory, doPrint = False, gyro = True, ideal = False)
			if check['gyro_complete']:
				print(f"{len(check['gyro_complete'])} Existing Gyro Runs Detected")
			runs = check['gyro_incomplete']
			'''
			runs = self.get_all_runs()
		else:
			runs = specificRuns
		
		for run in runs:
			sub_dir = f"{directory}/" + "/".join([f"{name} = {run[name]:.2g}" for name in self.inputs.dim_order])
			os.makedirs(sub_dir,exist_ok=True)
			existing_inputs = [] 
			for f in glob.glob(r'itteration_*.in'):
				existing_inputs.append([x for x in f if x.isdigit()])
			itt = max([eval("".join(x)) for x in existing_inputs],default=-1) + 1
			filename = f"itteration_{itt}"
			
			subnml = self.eqbm.get_gyro_input(run = run)
			subnml.write(f"{sub_dir}/{filename}.in", force=True)
						
			if self['System'] == 'plasma':
				self.jobs.append(f"mpirun -np 8 gs2 \"{sub_dir}/{filename}.in\"")
			else:
				jobfile = open(f"{sub_dir}/{filename}.job",'w')
				jobfile.write(f"""#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name={self.info['run_name']}
#SBATCH --ntasks=1
#SBATCH --output={sub_dir}/{filename}.slurm

{modules}

which gs2
gs2 --build-config

echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"
echo \"Input: {sub_dir}/{filename}.in\"
gs2 \"{sub_dir}/{filename}.in\"
echo \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"""")
				jobfile.close()
				self.jobs.append(f"sbatch \"{sub_dir}/{filename}.job\"")
	
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
		self.info = {'run_name': self.run_name, 'run_uuid': str(ID), 'data_path': run_path, 'input_file': self.inputs.input_name, 'eq_file_name': self.eqbm.eq_name, 'template_file_name': self.template_name, 'kin_file_name': self.eqbm.kin_name, 'kinetics_type': self.eqbm.kinetics_type, 'run_data': date, '_eq_file_path': self.eqbm._eq_path, '_kin_file_path': self.eqbm._kin_path, '_template_file_path': self._template_path, 'itteration': 0}
	
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
				for b in range(self['n_beta']):
					for s in range(self['n_shat']):
						for k in range(self['n_aky']):
							for t in range(self['n_theta0']):
								if os.path.exists(f"{directory}/{psiN}/{b}/{s}/{k}/{b}_{s}_{k}_{t}.out.nc"):
									finished_gyro.add((p,b,s,k,t))
								else:
									unfinished_gyro.add((p,b,s,k,t))

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
	
	def quick_save(self, filename = None, directory = None, SlurmSave = False):
		self.save_out(filename = filename, directory = directory, SlurmSave = SlurmSave, QuickSave = True)
	def save_out(self, filename = None, directory = None, SlurmSave = False, QuickSave = False):
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
		
		if self['System'] in ['viking','archer'] and not SlurmSave:
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
	run = myro_scan(eq_file = \"{self.eqbm.eq_name}\", kin_file = \"{self.eqbm.kin_name}\", input_file = \"{self.inputs.input_name}\", kinetics_type = \"{self.eqbm.kinetics_type}\", template_file = \"{self.template_name}\", directory = \"{self.path}\", run_name = \"{self.run_name}\")
	run.info = info
	run.namelist_diffs = nd
	run.save_out(filename = \"{filename}\", directory = \"{directory}\",SlurmSave = True,QuickSave = {QuickSave})""")
			pyth.close()
			os.system(f"sbatch \"save_out.job\"")
			os.chdir(f"{self.path}")
			return
			
		if not self._check_setup():
			return
			
		psiNs = self['psiNs']
		beta_prime_values = full((self['n_psiN']),None)
		shear_values = full((self['n_psiN']),None)
		
		if self['Gyro']:
			#beta_prime_axis = full((self['n_psiN'],self['n_beta']),None).tolist()
			#shear_axis = full((self['n_psiN'],self['n_shat']),None).tolist()
			#akys = full((self['n_psiN'],self['n_beta'],self['n_shat']),None).tolist()
			grs = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
			mfs = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
			syms = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
			
			if self['Epar']:
				eparNs = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
			else:
				eparNs = None
			
			if not QuickSave:
				omega = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				phi = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				apar = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				bpar = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				phi2 = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				time = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				theta = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				jacob = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
				gds2 = full((self['n_psiN'],self['n_beta'],self['n_shat'],self['n_aky'],self['n_theta0']),None).tolist()
			else:
				omega = phi = apar = bpar = phi2 = time = theta = jacob = gds2 = None
                #beta_prime_axis = shear_axis = akys = None
		else:
				grs = mfs = syms = self.inputs['n_shat'] = self.inputs['n_beta'] = self.inputs['aky_values'] = eparNs = omega = phi = apar = bpar = phi2 = time = theta = jacob = gds2 = None
                
		if self['Ideal']:
			beta_prime_axis_ideal = full((self['n_psiN'],self['n_beta_ideal']),None).tolist()
			shear_axis_ideal = full((self['n_psiN'],self['n_shat_ideal']),None).tolist()
			stabilities = full((self['n_psiN'],self['n_beta_ideal'],self['n_shat_ideal']),None).tolist()
		else:
			beta_prime_axis_ideal = shear_axis_ideal = stabilities = None#self.inputs['n_shat_ideal'] = self.inputs['n_beta_ideal'] = 
		
		for p, psiN in enumerate(psiNs):
			run_path = os.path.join(directory,str(psiN))
			nml = f90nml.read(f"{run_path}/{psiN}.in")
			shear = nml['theta_grid_eik_knobs']['s_hat_input']
			beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
			shear_values[p] = shear
			beta_prime_values[p] = beta_prim
			
			if self['Gyro']:
				for i, bp in enumerate(self['betas']):
					for j, sh in enumerate(self['shats']):
						for k, ky in enumerate(self['akys']):
							for t, t0 in enumerate(self['theta0s']):
								try:
									if self['Epar'] and not QuickSave:
										data = readnc(f"{run_path}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc",only=['omega','phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob'])
									elif self['Epar']:
										data = readnc(f"{run_path}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc",only=['omega','phi','bpar'])
									elif not QuickSave:
										data = readnc(f"{run_path}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc",only=['omega','phi','apar','bpar','phi2','t','theta', 'gds2', 'jacob'])
									else:
										data = readnc(f"{run_path}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc",only=['omega','phi'])	
									
									om = data['omega'][-1,0,0]
									if type(om) != complex:
										om = data['omega'][-2,0,0]
									grs[p][i][j][k][t] = imag(om)
									mfs[p][i][j][k][t] = real(om)
									
									symsum = sum(abs(data['phi'][0,0,:] + data['phi'][0,0,::-1]))/sum(abs(data['phi'][0,0,:]))
									if  symsum > 1.9:
										syms[p][i][j][k][t] = 1
									elif symsum < 1:
										syms[p][i][j][k][t] = -1
									else:
										syms[p][i][j][k][t] = 0
									
									if not QuickSave:
										try:
											omega[p][i][j][k] = data['omega'][:,0,0].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: omega")
										try:
											phi[p][i][j][k] = data['phi'][0,0,:].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: phi")
										try:
											apar[p][i][j][k] = data['apar'][0,0,:].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: apar")
										try:
											bpar[p][i][j][k] = data['bpar'][0,0,:].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: bpar")
										try:
											phi2[p][i][j][k] = data['phi2'].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: phi2")
										try:
											time[p][i][j][k] = data['t'].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: time")
										try:
											theta[p][i][j][k] = data['theta'].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: theta")
										try:
											jacob[p][i][j][k] = data['jacob'].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: jacob")
										try:
											gds2[p][i][j][k] = data['gds2'].tolist()
										except: 
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: gds2")
									
									if self['Epar']:
										epar_path = f"{run_path}/{i}/{j}/{k}/{i}_{j}_{k}.epar"
										try:
											bpar = data['bpar'][0,0,:]
											epar_data = loadtxt(epar_path)
											epar = []
											for l in range(len(epar_data[:,3])):
												epar.append(complex(epar_data[l,3],epar_data[l,4]))
											epar = array(epar)
											epar_norm = max(abs(epar))/max(abs(bpar))
											
											eparNs[p][i][j][k][t] = epar_norm
										except:
											print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}: epar")
									
								except Exception as e:
									print(f"Save Error {psiN}/{i}/{j}/{k}/{i}_{j}_{k}_{t}: {e}")
									grs[p][i][j][k][t] = None
									mfs[p][i][j][k][t] = None
									syms[p][i][j][k][t] = None
									if self['Epar']:
										eparNs[p][i][j][k][t] = None

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
		'beta_prime_axis': [self['betas']]*self['n_psi'],
		'shear_axis': [self['shears']]*self['n_psi'],
		'beta_prime_axis_ideal': beta_prime_axis_ideal,
		'shear_axis_ideal': shear_axis_ideal,
		'growth_rates_all': grs,
		'mode_frequencies_all': mfs,
		'parities_all': syms,
		'ideal_stabilities': stabilities,
		'eparN_all': eparNs,
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
		if self.info is None:
			self._create_run_info()
		if directory is None:
			directory = self.info['data_path']
		cancelled = self.check_cancelled(directory = directory, doPrint = False)
		if len(cancelled) == 0:
			print("No Cancelled Runs")
			return
		for p, i, j, k, t in cancelled:
			if os.path.exists(f"{directory}/{self['psiNs'][p]}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc"):
				os.remove(f"{directory}/{self['psiNs'][p]}/{i}/{j}/{k}/{i}_{j}_{k}_{t}.out.nc")
		self._make_gyro_files(directory = directory, specificRuns = cancelled, group_runs = group_runs)
		self._run_jobs()
	
	def rerun(self, runs = None, nml = None, directory = None, group_runs = None):
		if runs is None:
			print("ERROR: runs not given")
			return
		if nml is None:
			print("ERROR: nml not given, if you wish to rerun with no changes please use nml = {}")
			return
			
		self._check_setup()
		
		if type(nml) == str:
			nml = f90nml.read(nml)
		for p,i,j,k,t in runs:
			self.namelist_diffs[p][i][j][k][t] = nml
		if self.info:
			self.info['itteration'] += 1
		self._make_gyro_files(specificRuns = runs, directory = directory, group_runs = group_runs)
		self._run_jobs()
	
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
			directory = os.path.join(self.path, self.run_name)
		with load(f"{directory}/{filename}",allow_pickle = True) as obj:
			nd = obj['name_diffs']
			info = obj['info'].item()
			self.info = info
			self.namelist_diffs = nd
