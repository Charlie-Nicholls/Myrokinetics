import os
from numpy import *
from ncdf2dict import ncdf2dict as readnc
from .equillibrium import equillibrium
import f90nml

'''
DIAGNOSTIC SET SCAN PERFORMER
'''

class myro_set_scan(object):

	def __init__(self, eq_file = None, kin_file = None, input_file = None, template_file = None, kinetics_type = "PEQDSK", directory = "./", run_name = None):
		self._create_empty_inputs()
		self.template_name = template_file
		self.input_name = input_file
		self.run_name = run_name
		self._template_path = self.info = self.pyro = self._template_lines = None

		if directory == "./":
			directory = os.getcwd()
		self.path = directory
		
		self.eqbm = self.equillibrium = equillibrium(eq_file = eq_file, kin_file = kin_file, kinetics_type = kinetics_type, directory = directory)
		
		if input_file is not None:
			self.load_inputs()
	
	def _create_empty_inputs(self):
		self.inputs = {
		'variable': None,
		'values' : None,
		'key': None,
		'shat': None,
		'beta_prime': None,
		'psiNs': None,
		'aky_values': None,
		'Miller': True,
		'Viking': False,
		'Fixed_delt': False}
		
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
			if key == 'psiNs' or key == 'aky_values':
				if self.inputs[key] is None:
					arr_in = input(f"Input values for {key} (Current value: None): ")
				else:
					arr_in = input(f"Input values for {key} (Current value: {self.inputs[key]}): ")
				if arr_in != '':
					if arr_in[0] != "[":
						arr_in = "[" + arr_in
					if arr_in[-1] != "]":
						arr_in = arr_in + "]"
					self.inputs[key] = list(eval(arr_in))
					self.inputs[key].sort()
			else:
				if self.inputs[key] is None:
					inp = input(f"Input value for {key} (Current value: None): ")
					if inp != '':
						self.inputs[key] = eval(inp)
				else:
					inp = input(f"Input value for {key} (Current value: {self.inputs[key]}): ")
					if inp != '':
						self.inputs[key] = eval(inp)
	
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
				if key == 'shear':
					key = 'shat'
				if key in self.inputs.keys():
					if key == 'psiNs' or key == 'aky_values':
						if val[0] != "[":
							val = "[" + val
						if val[-1] != "]":
							val = val + "]"
					if key == 'variable' or 'key':
						if val[0] != "\'":
							val = "\'" + val
						if val[-1] != "\'":
							val = val + "\'"
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
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		if template_file is not None:
				self.template_name = template_file
		if self.template_name:
			with open(os.path.join(directory,self.template_name)) as tfile:
				self._template_lines = tfile.readlines()
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		if directory == "./":
			directory = os.getcwd()
		self._template_path = directory
		
		self.pyro = self.eqbm.load_pyro(template_file = self.template_name, directory = directory)
			
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
		self.info = {'run_name': self.run_name, 'run_uuid': str(ID), 'data_path': run_path, 'input_file': self.input_name, 'eq_file_name': self.eqbm.eq_name, 'template_file_name': self.template_name, 'kin_file_name': self.eqbm.kin_name, 'kinetics_type': self.eqbm.kinetics_type, 'run_data': date, '_eq_file_path': self.eqbm._eq_path, '_kin_file_path': self.eqbm._kin_path, '_template_file_path': self._template_path}

	def check_complete(self, directory = None, doPrint = True):
		if self.info is None:
			self._create_run_info()		
		if directory is None:
			directory = self.info['data_path']
			
		unfinished = []
		finished = []
		for psiN in self.inputs['psiNs']:
			for k, aky in enumerate(self.inputs['aky_values']):
				for v, val in enumerate(self.inputs['values']):
					if os.path.exists(f"{directory}/{psiN}/{v}_{k}.out.nc"):
						finished.append([psiN,v,k])
					else:
						unfinished.append([psiN,v,k])
		
		if doPrint:
			print(f"Runs Complete: {len(finished)} | Incomplete : {len(unfinished)}")
			return
		else:
			return {'complete': finished, 'incomplete': unfinished}
	
	def _check_setup(self):
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
		
		if type(self.inputs['aky_values']) in [int,float]:
			self.inputs['aky_values'] = [self.inputs['aky_values']]
		elif self.inputs['aky_values'] is None:
			self.inputs['aky_values'] = [0.1]
		else:
			self.inputs['aky_values'].sort()

		if type(self.inputs['values']) in [int,float]:
			self.inputs['values'] = [self.inputs['values']]
		elif self.inputs['values'] is None:
			empty_elements.append('values')
		else:
			self.inputs['values']
			
		if self.inputs['variable'] is None:
			empty_elements.append('variable')
		if self.inputs['shat'] is None:
			empty_elements.append('shat')
		if self.inputs['beta_prime'] is None:
			empty_elements.append('beta_prime')

		if empty_elements:
			print(f"ERROR: the following inputs are empty: {empty_elements}")
			return False
		
		return True
	
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
		
		if self.inputs['Miller']:
			nml['theta_grid_eik_knobs']['iflux'] = 0
			nml['theta_grid_eik_knobs']['local_eq'] = True
		else:
			nml['theta_grid_eik_knobs']['eqfile'] = os.path.join(self.path,self.eqbm.eq_name)
			nml['theta_grid_eik_knobs']['efit_eq'] =  True
			nml['theta_grid_eik_knobs']['iflux'] = 1
			nml['theta_grid_eik_knobs']['local_eq'] = False

		nml.write(file_name, force=True)
		return shear, beta_prim, tprim, fprim, beta, nml
	
	def run_scan(self, directory = None):
		if not self._check_setup():
			return
		check = self.check_complete(directory = directory, doPrint = False)
		if check['complete']:
			print(f"{len(check['complete'])} Existing Runs Detected")
		
		self._create_run_info()
		
		if directory is not None:
			self.info['data_path'] = directory
		else:
			directory = self.info['data_path']
		
		if not os.path.exists(directory):
			os.mkdir(directory)
	
		for psiN in self.inputs['psiNs']:
			run_path = os.path.join(directory, str(psiN))
			if not os.path.exists(run_path):
				os.mkdir(run_path)
		
			shear, beta_prim, tprim, fprim, beta, nml = self._make_fs_in(run_path=run_path, psiN=psiN)
			if self.inputs['key'] is None:
				for key1 in nml.keys():
					for key2 in nml[key1].keys():
						if key2 == self.inputs['variable']:
							self.inputs['key'] = key1
				if self.inputs['key'] is None:
					print("ERROR: Could not find namelist key")
				return
			for k, aky in enumerate(self.inputs['aky_values']):
				for v, val in enumerate(self.inputs['values']):
						if [psiN,v,k] in check['incomplete']:
							subnml = nml
							subnml['theta_grid_eik_knobs']['s_hat_input'] = self.inputs['shat']
							subnml['theta_grid_eik_knobs']['beta_prime_input'] = self.inputs['beta_prime']
							subnml['kt_grids_single_parameters']['aky'] = aky
							for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
								mul = self.inputs['beta_prime']/(-2*(nml[spec]['tprim'] + nml[spec]['fprim'])*beta)
								subnml[spec]['tprim'] = nml[spec]['tprim']*mul
								subnml[spec]['fprim'] = nml[spec]['fprim']*mul
							if self.inputs['Fixed_delt'] is False:
								subnml['knobs']['delt'] = 0.04/aky
							subnml[self.inputs['key']][self.inputs['variable']] = val
							subnml.write(f"{run_path}/{v}_{k}.in", force=True)
							
							if not self.inputs['Viking']:
								os.system(f"mpirun -np 8 gs2 \"{run_path}/{v}_{k}.in\"")
							else:
								jobfile = open(f"{run_path}/{v}_{k}.job",'w')
								jobfile.write(f"#!/bin/bash\n#SBATCH --time=05:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n\nmodule purge\nmodule load tools/git\nmodule load compiler/ifort\nmodule load mpi/impi\nmodule load numlib/FFTW\nmodule load data/netCDF/4.6.1-intel-2018b\nmodule load data/netCDF-Fortran/4.4.4-intel-2018b\nmodule load numlib/imkl/2018.3.222-iimpi-2018b\nmodule load lang/Python/3.7.0-intel-2018b\nexport GK_SYSTEM=viking\nexport MAKEFLAGS=-IMakefiles\nexport PATH=$PATH:$HOME/gs2/bin\n\necho \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"\n\"Input: {run_path}/{v}_{k}.in\"\necho \"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\"\n\nwhich gs2\n\ngs2 --build-config\n\ngs2 \"{run_path}/{v}_{k}.in\"")
								
								jobfile.close()
								os.chdir(f"{run_path}")
								os.system(f"sbatch \"{run_path}/{v}_{k}.job\"")
								os.chdir(f"{self.path}")
								
	def save_out(self, filename = None, directory = None, VikingSave = False):
		if filename is None and self.run_name is None:
			filename = input("Output File Name: ")
			filename = filename.split(".")[0]
		elif filename is None:
			filename = self.run_name
			
		if self.info is None:
			self._create_run_info()
		if directory is None:
			directory = self.info['data_path']
			
		if self.inputs['Viking'] and not VikingSave:
			os.chdir(f"{directory}")
			job = open(f"save_out.job",'w')
			job.write(f"#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n#SBATCH --mem=10gb\n\nmodule load lang/Python/3.7.0-intel-2018b\nmodule swap lang/Python lang/Python/3.10.4-GCCcore-11.3.0\n\nsource $HOME/pyroenv2/bin/activate\n\npython {directory}/save_out.py")
			job.close()
			pyth = open(f"save_out.py",'w')
			pyth.write(f"from myrokinetics import myro_set_scan\n\nrun = myro_set_scan(eq_file = \"{self.eqbm.eq_name}\", kin_file = \"{self.eqbm.kin_name}\", input_file = \"{self.input_name}\", kinetics_type = \"{self.eqbm.kinetics_type}\", template_file = \"{self.template_name}\", directory = \"{self.path}\", run_name = \"{self.run_name}\")\nrun.save_out(filename = \"{filename}\", directory = \"{directory}\",VikingSave = True)")
			pyth.close()
			os.system(f"sbatch \"save_out.job\"")
			os.chdir(f"{self.path}")
			return
		
		psiNs = self.inputs['psiNs']
		values = self.inputs['values']
		akys = self.inputs['aky_values']
		in_lines = full((len(psiNs),len(akys),len(values),None)).tolist()
		out_dict = full((len(psiNs),len(akys),len(values),None)).tolist()
		
		for p, psiN in enumerate(psiNs):
			for k, aky in enumerate(values):
				for v, val in enumerate(akys):
					try:
						in_lines[p][k][v] = f90nml.read(f"{self.info['data_path']}/{psiN}/{k}_{v}.in")
					except Exception as e:
						print(f"Input Save Error {psiN}/{k}_{v}: {e}")
						in_lines[p][k][v] = None
					try:
						out_dict[p][k][v] = readnc(f"{self.info['data_path']}/{psiN}/{k}_{v}.out.nc")
					except Exception as e:
						print(f"Output Save Error {psiN}/{k}_{v}: {e}")
						out_dict[p][k][v] = None
		file_lines = {'eq_file': self.eqbm._eq_lines, 'kin_file': self.eqbm._kin_lines, 'template_file': self._template_lines}
		savez(f"{self.path}/{filename}", inputs = self.inputs, run_info = self.info, input_namelists = in_nmls, output_dicts = out_dict, files = file_lines)
