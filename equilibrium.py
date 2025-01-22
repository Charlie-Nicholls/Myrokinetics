import os
import f90nml
from .inputs import scan_inputs
from .templates import template_dir
from .codes import codes
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
			self.inputs.inputs['files']['template_name'] = self.inputs.code.template
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
		
		self._template_lines = self.inputs.code.get_template_lines()

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
		if psiN not in self.surface_namelists.keys():
			if self.pyro is None:
				self.load_pyro()
		
			self.pyro.load_local(psi_n=psiN)
			self.pyro.update_gk_code()
			
			nml = deepcopy(self.pyro.gk_input.data)
			for dim in self.inputs.single_parameters.values():
				nml = dim.single_edit_nml(nml)
				
			nml = self.inputs.code._get_surface_input(nml)
			for dim in self.inputs.single_parameters.values():
				nml = dim.single_edit_nml(nml)
			
			self.surface_namelists[psiN] = nml
			
		return deepcopy(self.surface_namelists[psiN])

	def get_gyro_input(self, run = None, indexes = None, namelist_diff = {}):
		nml = self.inputs.code.get_gyro_input(run = run, indexes=indexes, namelist_diff=namelist_diff)
		return nml
	
	def write_nml(self, nml, directory = ".", filename = None):
		 self.inputs.code.write_nml(nml=nml,directory=directory,filename=filename)

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
