import os

class equillibrium(object):
	
	def __init__(self, eq_file = None, kin_file = None, kinetics_type = None, directory = None):
		self.eq_name = eq_file
		self.kin_name = kin_file
		self.kinetics_type = kinetics_type
		self.path = directory
		self.eq_data = self.kin_data = self._eq_path = self._kin_path = self._eq_lines = self._kin_lines = None
		if eq_file:
			self.load_geqdsk()
		if kin_file:
			self.load_kinetics()

	def load_geqdsk(self, eq_file = None, directory = None):
		from .geqdsk_reader import geqdsk
		if directory is None and self._eq_path is None and self.path is None:
			directory = "./"
		elif directory is None and self._eq_path is None:
			directory = self.path
		elif directory is None:
			directory = self._eq_path
		if directory == "./":
			directory = os.getcwd()
		self._eq_path = directory
		
		if self.eq_name is None and eq_file is None:
			print("ERROR: No GEQDSK file given")
			return
		elif eq_file is None:
			eq_file = self.eq_name
		self.eq_name = eq_file
		
		with open(os.path.join(directory,self.eq_name)) as efile:
			self._eq_lines = efile.readlines()
			
		self.eq_data = geqdsk(filename = self.eq_name, directory = directory)
	
	def load_kinetics(self, kin_file = None, kinetics_type = None, directory = None):
		if directory is None and self._kin_path is None and self.path is None:
			directory = "./"
		elif directory is None and self._kin_path is None:
			directory = self.path
		elif directory is None:
			directory = self._kin_path
		if directory == "./":
			directory = os.getcwd()
		self._kin_path = directory
		
		if self.kin_name is None and kin_file is None:
			print("ERROR: No Kinetics file given")
			return
		elif kin_file is None:
			kin_file = self.kin_name
		self.kin_name = kin_file
		
		with open(os.path.join(directory,self.kin_name)) as kfile:
			self._kin_lines = kfile.readlines()
			
		if kinetics_type is not None:
			self.kinetics_type = kinetics_type
		if self.kinetics_type is None:
			print("ERROR: kinetics_type not specified")
		
		if self.kinetics_type.upper() == "SCENE":
			import xarray as xr
			self.kin_data = xr.open_dataset(os.path.join(self._kin_path,self.kin_name))
		elif self.kinetics_type.upper() == "PEQDSK":
			from .peqdsk_reader import peqdsk
			self.kin_data = peqdsk(filename = self.kin_name, directory = directory)
			if 'rhonorm' not in self.kin_data.keys():
				self.AmmendPEQDSK(peq_file = os.path.join(directory,self.kin_name), geq = self.eq_data)
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

		eq_file = Path(self._eq_path) / self.eq_name
		kin_file = Path(self._kin_path) / self.kin_name
		
		if template_file is None:
			pyro = Pyro(
				eq_file=eq_file,
			 	eq_type="GEQDSK",
			 	kinetics_file=kin_file,
			 	kinetics_type=self.kinetics_type)
		else:
			if directory is None and self.path is None:
				directory = "./"
			elif directory is None:
				directory = self.path
			pyro = Pyro(
				eq_file=eq_file,
			 	eq_type="GEQDSK",
			 	kinetics_file=kin_file,
			 	kinetics_type=self.kinetics_type,
			 	gk_file=Path(directory) / template_file)

		pyro.gk_code = "GS2"
		return pyro
	
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
		
		f = open(self.kin_file,'a')
		f.write(f"{len(rho)+1} psinorm rho rhonorm")
		f.write("\n 0.0000000   0.0000000   0.0000000")
		for i in range(len(rho)):
			f.write(f"\n {psi_n[i]:.7f}   {rho[i]:.7f}   {rhonorm[i]:.7f}")  
		f.close()
		return
