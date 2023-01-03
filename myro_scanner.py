import os
import sys
from numpy import *
from ast import literal_eval
from pyrokinetics import Pyro
from pathlib import Path
from ncdf2dict import ncdf2dict as readnc
from geqdsk_reader import geqdsk

from .peqdsk_reader import peqdsk
from .ammend_peqdsk import AmmendPEQDSK

'''
GYROKINETIC SCAN PERFORMER
'''

class myro_scan(object):
	def __init__(self, eq_file = None, kin_file = None, input_file = None, template_file = None, kinetics_type = "PEQDSK", directory = "./", run_name = None):
		self._create_empty_inputs()
		self.kinetics_type = kinetics_type
		self.eq_name = eq_file
		self.kin_name = kin_file
		self.template_name = template_file
		self.input_name = input_file
		self.run_name = run_name
		self._eq_path = self._kin_path = self._template_path = self.gdat = self.kdat = self.info = self.pyro = self._eq_lines = self._kin_lines = self._template_lines = None

		if directory == "./":
			directory = os.getcwd() 
		self.path = directory
		if eq_file is not None:
			self.load_geqdsk()
		if kin_file is not None:
			self.load_kinetics()
		if input_file is not None:
			self.load_inputs()
	
	def __getitem__(self, key):
		if key == "inputs":
			self.inputs()
		else:
        		return self.inputs[key]

	def inputs(self):
        	for key, val in self.inputs.items():
        		print(f"{key} = {val}")
        	
	def keys(self):
		print(self.inputs.keys())
		
	def _create_empty_inputs(self):
		self.inputs = {
		'shat_min': None,
		'shat_max': None,
		'beta_mul': None,
		'beta_div': None,
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
		'Epar': False}

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
					arr_in = input(f"Input values for {key}: ")
				else:
					arr_in = input(f"Input values for {key} (Current value {self.inputs[key]}): ")
				if arr_in != '':
					if arr_in[0] != "[":
						arr_in = "[" + arr_in
					if arr_in[-1] != "]":
						arr_in = arr_in + "]"
					self.inputs[key] = list(literal_eval(arr_in))
					self.inputs[key].sort()
			else:
				if self.inputs[key] is None:
					inp = input(f"Input value for {key}: ")
					if inp != '':
						self.inputs[key] = literal_eval(inp)
				else:
					inp = input(f"Input value for {key} (Current value: {self.inputs[key]}): ")
					if inp != '':
						self.inputs[key] = literal_eval(inp)
	
	def load_geqdsk(self, eq_file = None, directory = "./"):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		self._eq_path = directory
		if self.eq_name is None and eq_file is None:
			raise ValueError("ERROR: No GEQDSK file given")
		elif self.eq_name is None:
			self.eq_name = eq_file
		with open(os.path.join(directory,self.eq_name)) as gfile:
			self._eq_lines = gfile.readlines()
		self.gdat = geqdsk(filename = self.eq_name, directory = directory)
	
	def load_kinetics(self, kin_file = None, kinetics_type = None, directory = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		self._kin_path = directory
		if self.kin_name is None and kin_file is None:
			raise ValueError("ERROR: No Kinetics file given")
		elif self.kin_name is None:
			self.kin_name = kin_file
		with open(os.path.join(directory,self.kin_name)) as kfile:
			self._kin_lines = kfile.readlines()
		if kinetics_type is not None:
			self.kinetics_type = kinetics_type
		if self.kinetics_type.upper() == "SCENE":
			import xarray as xr
			self.kdat = xr.open_dataset(os.path.join(self._kin_path,self.kin_name))
		elif self.kinetics_type.upper() == "PEQDSK":
			self.kdat = peqdsk(self.kin_name, directory)
			try:
				self.kdat['rhonorm']
			except:
				if self.gdat is None:
					self.load_geqdsk()
				AmmendPEQDSK(peq_file = os.path.join(directory,self.kin_name), geq = self.gdat)
				self.kdat = peqdsk(self.kin_name, directory)
				try:
					self.kdat['rhonorm']
				except:
					raise ValueError("ERROR: Could not load rho data from PEQDSK file")
					
			
		else:
			print(f"ERROR: Kinetics type {self.kinetics_type} not recognised. Currently supported: SCENE, PEQDSK")
	
	def load_pyro(self, template_file = None, directory = None):
		if self.eq_name is None or self.kin_name is None:
			if self.eq_name is None:
				print("ERROR: No equillibrium file given")
			if self.kin_name is None:
				print("ERROR: No kinetics file given")
			return
		if self.kdat is None:
			self.load_kinetics()
		if self.gdat is None:
			self.load_geqdsk()
		eq_file = Path(self._eq_path) / self.eq_name
		kin_file = Path(self._kin_path) / self.kin_name
		
		if self.template_name is None and template_file is None:
			self.pyro = Pyro(
				eq_file=eq_file,
			 	eq_type="GEQDSK",
			 	kinetics_file=kin_file,
			 	kinetics_type=self.kinetics_type)
		else:
			if self.template_name is None:
				self.template_name = template_file
			if directory is None and self.path is None:
				directory = "./"
			if directory is None:
				directory = self.path
			self._template_path = directory
			with open(os.path.join(directory,self.template_name)) as tfile:
				self._template_lines = tfile.readlines()
			self.pyro = Pyro(
				eq_file=eq_file,
			 	eq_type="GEQDSK",
			 	kinetics_file=kin_file,
			 	kinetics_type=self.kinetics_type,
			 	gk_file = Path(self.template_name))
	
		self.local_geometry = "Miller"
		self.pyro.gk_code = "GS2"	
	
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
				if key in self.inputs.keys() or key == 'akys':
					if key == 'akys':
						key = 'aky_values'
					if key == 'psiNs' or key == 'aky_values':
						if val[0] != "[":
							val = "[" + val
						if val[-1] != "]":
							val = val + "]"
					self.inputs[key] = literal_eval(val)
			
	def write_input_file(self, filename = None, directory = "./"):
		if self.input_name is None and filename is None:
			filename = input("Input File Name: ")
			if "." not in filename:
				filename = filename + ".in"
		elif self.input_name is None:
			self.input_name = filename
		
		with open(self.input_name,'w') as in_file:
			for key in self.inputs.keys():
				in_file.write(f"{key} = {self.inputs[key]}\n")

	def _make_fs_in(self, run_path = None, psiN = None):
		if self.pyro == None:
			self.load_pyro()
		if run_path == None:
			run_path = self.path
		if psiN == None:
			print("ERROR: please speicify psiN")
			return
		if self.inputs['Miller'] is False:
			eq_dir = os.path.join(self.path,self.eq_name)
			self.pyro.add_flags({"theta_grid_eik_knobs":{"eqfile": f"{eq_dir}", "efit_eq": True}})
		self.pyro.load_local_geometry(psi_n=psiN)
		self.pyro.load_local_species(psi_n=psiN)
		file_name=f"{run_path}/{psiN}.in"
		self.pyro.write_gk_file(file_name)
		with open(file_name) as f:
			lines = f.readlines()
			f_new = open(file_name,'w')
			for line in lines:
				if line.strip("\t\n ").split(" = ")[0] == "s_hat_input":
					shear = literal_eval(line.strip("\t\n").split(" = ")[1])
				if line.strip("\t\n ").split(" = ")[0] == "beta_prime_input":
					beta_prim = literal_eval(line.strip("\t\n").split(" = ")[1])
				if line.strip("\t\n ").split(" = ")[0] == "tprim":
					tprim = literal_eval(line.strip("\t\n").split(" = ")[1])
				if line.strip("\t\n ").split(" = ")[0] == "fprim":
					fprim = literal_eval(line.strip("\t\n").split(" = ")[1])
				if line.strip("\t\n ").split(" = ")[0] == "beta":
					beta = literal_eval(line.strip("\t\n").split(" = ")[1])
					'''
					pref = self.pyro.local_species.nref * self.pyro.local_species.tref * 1.602176634e-19
					bref = self.pyro.local_geometry.B0
					beta = pref / bref**2 * 8 * pi * 1e-7
					'''
					f_new.write(f"    beta = {beta}\n")
				elif self.inputs['Miller'] is False and line.strip("\t\n ").split(" = ")[0] == "iflux":
					f_new.write("    iflux = 1\n")
				elif self.inputs['Miller'] is False and line.strip("\t\n ").split(" = ")[0] == "local_eq":
					f_new.write("    local_eq = .false.\n")
				elif line.strip("\t\n ").split(" = ")[0] == "write_final_epar":
					if self.inputs['Epar']:
						f_new.write("    write_final_epar = .true.\n")
					else:
						f_new.write("    write_final_epar = .false.\n")
				else:
					f_new.write(line)
		if self.inputs['Ideal']:
			ballstab_knobs = f"\n&ballstab_knobs\n    n_shat = {self.inputs['n_shat_ideal']}\n    n_beta = {self.inputs['n_beta_ideal']}\n    shat_min = {self.inputs['shat_min']}\n    shat_max = {self.inputs['shat_max']}\n    beta_mul = {self.inputs['beta_mul']}\n    beta_div = {self.inputs['beta_div']}\n/\n"
			f_new.write(ballstab_knobs)
		f_new.close()
		return shear, beta_prim, tprim, fprim, beta
	
	def run_scan(self, gyro = None, ideal = None, miller = None, directory = None):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path
		
		if gyro is None:
			gyro = self.inputs['Gyro']
		elif type(gyro) == bool:	
			self.inputs['Gyro'] = gyro
		if ideal is None:
			ideal = self.inputs['Ideal']
		elif type(ideal) == bool:
			self.inputs['Ideal'] = ideal
		if miller is None:
			miller = self.inputs['Miller']
		elif type(miller) == bool:
			self.inputs['Miller'] = miller
		
		self._create_run_info()
		run_path = self.info['data_path']
		
		try:
			os.mkdir(run_path)
		except:
			pass	
		#Split for if dir in file names
		if self.template_name is not None:
			os.system(f"cp \"{directory}/{self.template_name}\" \"{run_path}/{self.template_name}\"")
		os.system(f"cp \"{self._kin_path}/{self.kin_name}\" \"{run_path}/{self.kin_name}\"")
		os.system(f"cp \"{self._eq_path}/{self.eq_name}\" \"{run_path}/{self.eq_name}\"")
		if self.input_name is not None:		
			os.system(f"cp \"{self.path}/{self.input_name}\" \"{run_path}/{self.input_name}\"")
		
		if self.inputs['Gyro']:
			self.run_gyro(directory = run_path)
		if self.inputs['Ideal']:
			self.run_ideal(directory = run_path)
		if not self.inputs['Viking']:
			self.save_out(directory = run_path)
					
	def run_ideal(self, directory = None):
		if self.info is None:
			self._create_run_info()
		if self.inputs['n_shat_ideal'] is None:
			self.inputs['n_shat_ideal'] = self.inputs['n_shat']
			print("n_shat_ideal is empty, setting to n_shat")
		if self.inputs['n_beta_ideal'] is None:
			self.inputs['n_beta_ideal'] = self.inputs['n_beta']
			print("n_beta_ideal is empty, setting to n_beta")
		if directory is None:
			directory = self.info['data_path']
		empty_elements = []
		for key in self.inputs.keys():
			if self.inputs[key] is None and key not in ['akys','n_beta','n_beta_ideal'] :
				empty_elements.append(key)
		if empty_elements:
			raise ValueError(f"ERROR: the following inputs are empty: {empty_elements}")
			return
		self.inputs['Ideal'] = True
		if not self.pyro:
			self.load_pyro()
		if type(self.inputs['psiNs']) == int or type(self.inputs['psiNs']) == float:
			self.inputs['psiNs'] = [self.inputs['psiNs']]

		check = self.check_complete(directory = directory, doPrint = False)
		if check['ideal_complete']:
			print("Existing Ideal Runs Detected")

		for psiN in self.inputs['psiNs']:
			run_path = os.path.join(directory, str(psiN))
			if psiN not in check['ideal_complete']:
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
					jobfile.write(f"#!/bin/bash\n#SBATCH --time=01:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n\nmodule purge\nmodule load tools/git\nmodule load compiler/ifort\nmodule load mpi/impi\nmodule load numlib/FFTW\nmodule load data/netCDF/4.6.1-intel-2018b\nmodule load data/netCDF-Fortran/4.4.4-intel-2018b\nmodule load numlib/imkl/2018.3.222-iimpi-2018b\nmodule load lang/Python/3.7.0-intel-2018b\nexport GK_SYSTEM=viking\nexport MAKEFLAGS=-IMakefiles\nexport PATH=$PATH:$HOME/gs2/bin\n\nwhich gs2\n\ngs2 --build-config\n\nideal_ball \"{run_path}/{psiN}.in\"")
					jobfile.close()
					os.chdir(f"{run_path}")
					os.system(f"sbatch \"{run_path}/{psiN}.job\"")
					os.chdir(f"{self.path}")
			else:	
				pass
				#print(f"Skipping Ideal {psiN}: Existing Run Detected")

	def run_gyro(self, directory = None):
		if self.info is None:
			self._create_run_info()
		if directory is None:
			directory = self.info['data_path']
		empty_elements = []
		for key in self.inputs.keys():
			if self.inputs[key] is None and key not in ['n_beta_ideal','n_shat_ideal']:
				empty_elements.append(key)
		if empty_elements:
			print(f"ERROR: the following inputs are empty: {empty_elements}")
			return
		self.inputs['Gyro'] = True
		if not self.pyro:
			self.load_pyro()
		if type(self.inputs['psiNs']) == int or type(self.inputs['psiNs']) == float:
			self.inputs['psiNs'] = [self.inputs['psiNs']]
		if type(self.inputs['aky_values']) == int or type(self.inputs['aky_values']) == float:
			self.inputs['aky_values'] = [self.inputs['aky_values']]
		
		check = self.check_complete(directory = directory, doPrint = False)
		if check['gyro_complete']:
			print(f"Existing Gyro Runs Detected")

		if len(check['gyro_incomplete']) > 10000:
			group_runs = True
		else:
			group_runs = False
		
		for psiN in self.inputs['psiNs']:
			run_path = os.path.join(directory, str(psiN))
			try:
				os.mkdir(run_path)
			except:
				pass
				
			shear, beta_prim, tprim, fprim, beta = self._make_fs_in(run_path=run_path, psiN=psiN)
			f = open(f"{run_path}/{psiN}.in")
			lines = f.readlines()
			f.close()
			beta_min = beta_prim/self.inputs['beta_div']
			beta_max = beta_prim*self.inputs['beta_mul']
			
			for i in range(self.inputs['n_beta']):
				for j in range(self.inputs['n_shat']):
					group_kys = []
					fol = str(i) + "_" + str(j)
					sub_path = os.path.join(run_path,fol)
					try:
						os.mkdir(sub_path)
					except:
						pass
					bp = (beta_max - beta_min)*i/(self.inputs['n_beta']-1) + beta_min
					sh = (self.inputs['shat_max'] - self.inputs['shat_min'])*j/(self.inputs['n_shat']-1) + self.inputs['shat_min']
					mul = bp/(-2*(tprim + fprim)*beta)

					for k, aky in enumerate(self.inputs['aky_values']):
						if [psiN,i,j,k] not in check['gyro_complete']:
							infile = open(f"{sub_path}/{fol}_{k}.in",'w')
							for line in lines:
								if line.strip("\t\n ").split(" = ")[0] == "s_hat_input":
									infile.write(f"    s_hat_input = {sh}\n")
								elif line.strip("\t\n ").split(" = ")[0] == "beta_prime_input":
									infile.write(f"    beta_prime_input = {bp}\n")
								elif line.strip("\t\n ").split(" = ")[0] == "tprim":
									tprim_new = mul*tprim
									infile.write(f"    tprim = {tprim_new}\n")
								elif line.strip("\t\n ").split(" = ")[0] == "fprim":
									fprim_new = mul*fprim
									infile.write(f"    fprim = {fprim_new}\n")
								elif line.strip("\t\n ").split(" = ")[0] == "beta":
									infile.write(f"    beta = {beta}\n")
								elif line.strip("\t\n ").split(" = ")[0] == "aky":
									infile.write(f"    aky = {aky}\n")
								
								elif self.inputs['Fixed_delt'] is False and line.strip("\t\n ").split(" = ")[0] == "delt":
									delt = 0.025/aky
									infile.write(f"    delt = {delt}\n")
								else:
									infile.write(line)
							infile.close()
							if not self.inputs['Viking']:
								os.system(f"mpirun -np 8 gs2 \"{sub_path}/{fol}_{k}.in\"")
							elif not group_runs:
								
								jobfile = open(f"{sub_path}/{fol}_{k}.job",'w')
								hours = len(self.inputs['aky_values'])
								jobfile.write(f"#!/bin/bash\n#SBATCH --time={hours}:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n\nmodule purge\nmodule load tools/git\nmodule load compiler/ifort\nmodule load mpi/impi\nmodule load numlib/FFTW\nmodule load data/netCDF/4.6.1-intel-2018b\nmodule load data/netCDF-Fortran/4.4.4-intel-2018b\nmodule load numlib/imkl/2018.3.222-iimpi-2018b\nmodule load lang/Python/3.7.0-intel-2018b\nexport GK_SYSTEM=viking\nexport MAKEFLAGS=-IMakefiles\nexport PATH=$PATH:$HOME/gs2/bin\n\nwhich gs2\n\ngs2 --build-config\n\ngs2 \"{sub_path}/{fol}_{k}.in\"")
								
								jobfile.close()
								os.chdir(f"{sub_path}")
								os.system(f"sbatch \"{sub_path}/{fol}_{k}.job\"")
								os.chdir(f"{self.path}")
							
							else:
								group_kys.append(k)
							
						else:
							pass
							#print(f"Skipping Gyro {psiN}/{fol}/{fol}_{k}: Existing Run Detected")
				
					if self.inputs['Viking'] and group_runs and group_kys:
						jobfile = open(f"{sub_path}/{fol}.job",'w')
						jobfile.write(f"#!/bin/bash\n#SBATCH --time=02:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n\nmodule purge\nmodule load tools/git\nmodule load compiler/ifort\nmodule load mpi/impi\nmodule load numlib/FFTW\nmodule load data/netCDF/4.6.1-intel-2018b\nmodule load data/netCDF-Fortran/4.4.4-intel-2018b\nmodule load numlib/imkl/2018.3.222-iimpi-2018b\nmodule load lang/Python/3.7.0-intel-2018b\nexport GK_SYSTEM=viking\nexport MAKEFLAGS=-IMakefiles\nexport PATH=$PATH:$HOME/gs2/bin\n\nwhich gs2\n\ngs2 --build-config\n\n")

						for k in group_kys:
							jobfile.write(f"gs2 \"{sub_path}/{fol}_{k}.in\"\n")
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
		if self.path == "./":
			path = os.getcwd()
		else:
			path = self.path
		if self.run_name is None:
			run_path = os.path.join(path, str(ID))
		else:	
			run_path = os.path.join(path, self.run_name)
		try:
			from datetime import datetime as dt
			date = dt.now().strftime("%d/%m/%Y, %H:%M:%S")
		except:
			print("ERROR: unable to import datetime module, setting run date to None")
			data = None
		self.info = {'run_name': self.run_name, 'run_uuid': str(ID), 'data_path': run_path, 'input_file': self.input_name, 'eq_file_name': self.eq_name, 'template_file_name': self.template_name, 'kin_file_name': self.kin_name, 'kinetics_type': self.kinetics_type, 'run_data': date, '_eq_file_path': self._eq_path, '_kin_file_path': self._kin_path, '_template_file_path': self._template_path}

	def check_complete(self, directory = None, doPrint = True):
		if self.info is None:
			self._create_run_info()		
		if directory is None:
			directory = self.info['data_path']

		unfinished_gyro = []
		finished_gyro = []
		if self.inputs['Gyro']:
			for psiN in self.inputs['psiNs']:
				for i in range(self.inputs['n_beta']):
					for j in range(self.inputs['n_shat']):
						for k, aky in enumerate(self.inputs['aky_values']):
							if os.path.exists(f"{directory}/{psiN}/{i}_{j}/{i}_{j}_{k}.out.nc"):
								finished_gyro.append([psiN,i,j,k])
							else:
								unfinished_gyro.append([psiN,i,j,k])

		unfinished_ideal = []
		finished_ideal = []
		if self.inputs['Ideal']:
			for psiN in self.inputs['psiNs']:
				if os.path.exists(f"{directory}/{psiN}/{psiN}.ballstab_2d"):
					finished_ideal.append(psiN)
				else:
					unfinished_ideal.append(psiN)
		
		if doPrint:
			print(f"Gyro Runs Complete: {len(finished_gyro)} | Incomplete : {len(unfinished_gyro)}")
			print(f"Ideal Runs Complete: {len(finished_ideal)} | Incomplete : {len(unfinished_ideal)}")
			return
		else:
			return {'gyro_complete': finished_gyro, 'gyro_incomplete': unfinished_gyro, 'ideal_complete': finished_ideal, 'ideal_incomplete': unfinished_ideal}
			
	def detailed_save(self, filename = None, directory = None, VikingSave = False):
		self.save_out(filename = filename, directory = directory, VikingSave = VikingSave, DetailedSave = True)
	def save_out(self, filename = None, directory = None, VikingSave = False, DetailedSave = False):
		if filename is None and self.run_name is None:
			filename = input("Output File Name: ")
			filename = filename.split(".")[0]
		elif filename is None:
			filename = self.run_name
			
		if self.info is None:
			self._create_run_info()
		if directory is None:
			directory = self.info['data_path']

		if self.info is None:
			self._create_run_info()
		
		if not self.inputs['Gyro'] and not self.inputs['Ideal']:
			print("Error: Both Gyro and Ideal are False")
			return
		
		if self.inputs['Viking'] and not VikingSave:
			os.chdir(f"{directory}")
			job = open(f"save_out.job",'w')
			job.write(f"#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --job-name={self.info['run_name']}\n#SBATCH --ntasks=1\n#SBATCH --mem=10gb\n\nmodule load lang/Python/3.7.0-intel-2018b\nmodule swap lang/Python lang/Python/3.10.4-GCCcore-11.3.0\n\nsource $HOME/pyroenv2/bin/activate\n\npython {directory}/save_out.py")
			job.close()
			pyth = open(f"save_out.py",'w')
			pyth.write(f"from myrokinetics import myro_scan\n\nrun = myro_scan(eq_file = \"{self.eq_name}\", kin_file = \"{self.kin_name}\", input_file = \"{self.input_name}\", kinetics_type = \"{self.kinetics_type}\", template_file = \"{self.template_name}\", directory = \"{self.path}\", run_name = \"{self.run_name}\")\nrun.save_out(filename = \"{filename}\", directory = \"{directory}\",VikingSave = True,DetailedSave = {DetailedSave})")
			pyth.close()
			os.system(f"sbatch \"save_out.job\"")
			os.chdir(f"{self.path}")
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
			
			if DetailedSave:
				omega = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				phi = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				apar = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				phi2 = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				time = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
				theta = full((len(psiNs),self.inputs['n_beta'],self.inputs['n_shat'],len(self.inputs['aky_values'])),None).tolist()
			else:
				omega = phi = apar = phi2 = time = theta = None
		
		else:
				gr = mf = grs = sym = syms = beta_prime_axis = shear_axis = akys = self.inputs['n_shat'] = self.inputs['n_beta'] = self.inputs['aky_values'] = eparN = eparNs = omega = phi = apar = phi2 = time = theta = None
		
		if self.inputs['Ideal']:
			beta_prime_axis_ideal = full((len(psiNs),self.inputs['n_beta_ideal']),None).tolist()
			shear_axis_ideal = full((len(psiNs),self.inputs['n_shat_ideal']),None).tolist()
			stabilities = full((len(psiNs),self.inputs['n_beta_ideal'],self.inputs['n_shat_ideal']),None).tolist()
		else:
			beta_prime_axis_ideal = shear_axis_ideal = self.inputs['n_shat_ideal'] = self.inputs['n_beta_ideal'] = stabilities = None
		
		for psiN in psiNs:
			idx = psiNs.index(psiN)
			run_path = os.path.join(directory,str(psiN))
			orig = open(f"{run_path}/{psiN}.in")
			lines = orig.readlines()
			orig.close()
			for line in lines:
				if line.strip("\t\n ").split(" = ")[0] == "s_hat_input":
					shear_values[idx] = literal_eval(line.strip("\t\n").split(" = ")[1])
				if line.strip("\t\n ").split(" = ")[0] == "beta_prime_input":
					beta_prim = abs(literal_eval(line.strip("\t\n").split(" = ")[1]))
					beta_prime_values[idx] = beta_prim
			
			if self.inputs['Gyro']:
				for i in range(self.inputs['n_beta']):
					per = 100*(idx*self.inputs['n_beta']+i)/(self.inputs['n_beta']*len(psiNs))
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
							
								if self.inputs['Epar'] and DetailedSave:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi','bpar','apar','phi2','t','theta'])
								elif self.inputs['Epar']:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi','bpar'])
								elif DetailedSave:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi','apar','phi2','t','theta'])
								else:
									data = readnc(f"{run_path}/{fol}/{fol}_{k}.out.nc",only=['omega','phi'])	
								
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
								
								if DetailedSave:
									try:
										omega[idx][i][j][k] = data['omega'][:,0,0].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: omega")
									try:
										phi[idx][i][j][k] = data['phi'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: phi")
									try:
										apar[idx][i][j][k] = data['apar'][0,0,:].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: apar")
									try:
										phi2[idx][i][j][k] = data['phi2'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: phi2")
									try:
										time[idx][i][j][k] = data['t'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: time")
									try:
										theta[idx][i][j][k] = data['theta'].tolist()
									except: 
										print(f"Save Error {psiN}/{fol}/{fol}_{k}: theta")
								
								if self.inputs['Epar']:
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
					
						aky_idx = gr_aky.index(amax(array(gr_aky)[isfinite(gr_aky)]))

						grs[idx][i][j] = gr_list
						mfs[idx][i][j] = mf_list
						syms[idx][i][j] = sym_list
						akys[idx][i][j] = self.inputs['aky_values'][aky_idx]
						gr[idx][i][j] = gr_list[aky_idx]
						mf[idx][i][j] = mf_list[aky_idx]
						sym[idx][i][j] = sym_list[aky_idx]
						if self.inputs['Epar']:
							eparNs[idx][i][j] = epars
							eparN[idx][i][j] = epars[aky_idx]
							
				beta_min = beta_prim/self.inputs['beta_div']
				beta_max = beta_prim*self.inputs['beta_mul']
				for i in range(self.inputs['n_beta']):
					beta_prime_axis[idx][i] = abs((beta_max - beta_min)*i/(self.inputs['n_beta']-1) + beta_min)
				for i in range(self.inputs['n_shat']):
					shear_axis[idx][i] = (self.inputs['shat_max'] - self.inputs['shat_min'])*i/(self.inputs['n_shat']-1) + self.inputs['shat_min']

			if self.inputs['Ideal']:
				shear = loadtxt(f"{run_path}/{psiN}.ballstab_shat")
				bp = loadtxt(f"{run_path}/{psiN}.ballstab_bp")
				stab = loadtxt(f"{run_path}/{psiN}.ballstab_2d")
				
				bp = [abs(x) for x in bp]
				beta_prime_axis_ideal[idx] = bp
				shear_axis_ideal[idx] = shear
				stabilities[idx] = transpose(stab)
			
		dat = {'beta_prime_values':beta_prime_values,
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
		'phi2': phi2,
		'time': time,
		'theta': theta
		}
		file_lines = {'eq_file': self._eq_lines, 'kin_file': self._kin_lines, 'template_file': self._template_lines}
		
		
		savez(f"{self.path}/{filename}", inputs = self.inputs, data = dat, run_info = self.info, files = file_lines)
