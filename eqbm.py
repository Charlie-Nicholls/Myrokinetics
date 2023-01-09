from numpy import *
from matplotlib.pyplot import *
from pyrokinetics import Pyro
from geqdsk_reader import geqdsk
from .peqdsk_reader import peqdsk

class eqbm(object):

	def __init__(self, eq_file = None, kin_file = None, eq_path = None, kin_path = None, eq_type = None, kin_type = None):
		self.eq_file = eq_file
		self.kin_file = kin_file
		self.eq_path = eq_path
		self.kin_path = kin_path
		self.eq_type = eq_type
		self.kin_type = kin_type
		
		self.eq = self.kin = None
		
		if eq_file:
			self.load_eq()
		if kin_file:
			self.load_kin()
		
	def load_eq(self, eq_file = None, eq_path = None, eq_type = None):
		if eq_file:
			self.eq_file = eq_file
		elif self.eq_file:
			eq_file = self.eq_file
		else:
			print("ERROR: No eq_file given")
		if eq_path:
			self.eq_path = eq_path
		elif self.eq_file:
			eq_path = self.eq_path
		else:
			self.eq_path = eq_path = "./"
		if eq_type:
			self.eq_type = eq_type
		elif self.eq_type:
			eq_type = self.eq_type
		else:
			self.eq_type = eq_type = "GEQDSK"
			
		if eq_type.upper() = "GEQDSK":
			self.eq = geqdsk(filename = eq_file, directory = eq_path)
		else:
			print(f"Equillibrium type {eq_type} not supported")
	
	def load_kin(self, kin_file = None, kin_path = None, kin_type = None):
		if kin_file:
			self.kin_file = kin_file
		elif self.kin_file:
			kin_file = self.kin_file
		else:
			print("ERROR: No kin_file given")
		if kin_path:
			self.kin_path = kin_path
		elif self.kin_file:
			kin_path = self.kin_path
		else:
			self.kin_path = kin_path = "./"
		if kin_type:
			self.kin_type = kin_type
		elif self.kin_type:
			kin_type = self.kin_type
		else:
			self.kin_type = kin_type = "PEQDSK"
			
		if kin_type.upper() = "PEQDSK":
			self.kin = peqdsk(filename = kin_file, directory = kin_path)
		elif kin_type.upper() == "SCENE":
			import xarray as xr
			self.kin = xr.open_dataset(os.path.join(kin_path,kin_file))
		else:
			print(f"Equillibrium type {kin_type} not supported")
			
	'''
	def getsb(eq, kin):
		pyro = Pyro(eq_file=eq,eq_type="GEQDSK",kinetics_file=kin,kinetics_type="PEQDSK",gk_file="../templates/template.gs2")
		pyro.local_geometry = "Miller"
		pyro.gk_code = "GS2"
		
		bp = []
		sh = []
		for psiN in psiNs:
			pyro.load_local_geometry(psi_n=psiN)
			bp.append(pyro.local_geometry['beta_prime'])
			sh.append(pyro.local_geometry['shat'])
			
		return bp, sh
		
	psiNs = linspace(0.01,1,100)

	b190, s190 = getsb("g045272.190","p045272.190")
	b210, s210 = getsb("g045272.210","p045272.210")
	b475, s475 = getsb("g045272.475","p045272.475")

	fig, ax = subplots(2,1)

	ax[0].plot(psiNs, [-x for x in b190], 'b', label="190ms")
	ax[0].plot(psiNs, [-x for x in b210], 'r', label="210ms")
	ax[0].plot(psiNs, [-x for x in b475], 'k', label="475ms")
	ax[0].plot(0.01, -b190[0], 'b.')
	ax[0].plot(0.01, -b210[0], 'r.')
	ax[0].plot(0.01, -b475[0], 'k.')
	ax[1].plot(psiNs, s190, 'b', label="190ms")
	ax[1].plot(psiNs, s210, 'r', label="210ms")
	ax[1].plot(psiNs, s475, 'k', label="475ms")
	ax[1].plot(0.01, s190[0], 'b.')
	ax[1].plot(0.01, s210[0], 'r.')
	ax[1].plot(0.01, s475[0], 'k.')
	ax[0].legend(loc=0)
	ax[1].legend(loc=0)
	ax[0].set_xlabel("psiN")
	ax[1].set_xlabel("psiN")
	ax[0].set_ylabel("beta_prime")
	ax[1].set_ylabel("shear")
	show()

	def AmmendPEQDSK(peq_file = None, geq_file = None, peq = None, geq = None):
		if geq_file is None and geq is None:
			print("ERROR: No GEQDSK file or object provided")
			return
		if geq is None:
			from geqdsk_reader import geqdsk
			geq = geqdsk(geq_file)
		if peq is None:
			from .peqdsk_reader import peqdsk
			peq = peqdsk(peq_file)
			
		psi_n = peq['psinorm']
		rho = []
		for psiN in psi_n:
			if psiN < 1.0e-3:
				psi_n = psi_n[psi_n != psiN]	
			else:
				fs = geq.flux_surface(psiN = psiN)
				rho.append((max(fs['R']) - min(fs['R']))/2)
		rhonorm = rho/max(rho)
		
		f = open(peq_file,'a')
		f.write(f"{len(rho)+1} psinorm rho rhonorm")
		f.write("\n 0.0000000   0.0000000   0.0000000")
		for i in range(len(rho)):
			f.write(f"\n {psi_n[i]:.7f}   {rho[i]:.7f}   {rhonorm[i]:.7f}")  
		f.close()
		return
