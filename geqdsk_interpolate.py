from .geqdsk_reader import geqdsk as readg
from matplotlib.pyplot import *
from scipy.interpolate import RectBivariateSpline, InterpolatedUnivariateSpline, RegularGridInterpolator
from matplotlib.widgets import TextBox
from numpy import linspace, zeros, diff, insert, zeros
from pyrokinetics import Pyro

class geqdsk_interpolator(object):
	def __init__(self, file1 = None, file2 = None):
		if file1 is None or file2 is None:
			print("ERROR: file1 and/or file2 is None")
		self.geq1 = readg(file1)
		self.geq2 = readg(file2)
		self.data = {}
		self._variables_sin = ['rmaxis','zmaxis','simag','sibry','bcentr','current']
		self._variables_arr = ['fpol','pres','ffprim','pprime','qpsi','psi']
		self._constants = ['nr','nw','nz','nh','rdim','zdim','rcentr','rleft','zmid','rlim','zlim','R','Z','psiN','rbbbs','zbbbs']
		'''
		psirz needs unique handling, descriptor can be set when making file
		NOT IMPLIMENTED: 'rbbbs','zbbbs' vary between files but are not yet interpolated between,
		for now just uses geq1's values as an approximation, they should be similar anyway
		'''
		self.create_interpolation()
	
	def create_interpolation(self, res = 300):
		dat = {}
		psiNs = self.geq1['psiN']
		self._psiNs = psiNs
		times = linspace(0,1,res)

		for cvar in self._constants:
			self.data[cvar] = self.geq1[cvar]

		for svar in self._variables_sin:
			dat[svar] = zeros((res))
			for t, time in enumerate(times):
				dat[svar][t] = (self.geq1[svar]-self.geq2[svar])*time + self.geq1[svar]
			self.data[svar] = InterpolatedUnivariateSpline(times,dat[svar])

		for avar in self._variables_arr:
			dat[avar] = zeros((res,len(psiNs)))
			for t, time in enumerate(times):
		  		dat[avar][t] = (self.geq1[avar]-self.geq2[avar])*time + self.geq1[avar]
			self.data[avar] = RectBivariateSpline(times,psiNs,dat[avar])
		
		dat['psirz'] = zeros((res,len(psiNs),len(psiNs)))
		for t, time in enumerate(times):
	  		dat['psirz'][t] = (self.geq1['psirz']-self.geq2['psirz'])*time + self.geq1['psirz']
		self.data['psirz'] = RegularGridInterpolator((times,self.data['R'],self.data['Z']),dat['psirz'])
	
	def var(self, var = None, t = None):
		if var is None:
			print("ERROR: var must be specified")
			return
		if t is None or t < 0 or t > 1:
			print("ERROR: t must be specified between 0 and 1")
			return
		if var in self._constants:
			return self.data[var]
		elif var in self._variables_sin:
			return self.data[var](t)
		elif var in self._variables_arr:
			return self.data[var](t,self._psiNs)[0]
		else:
			print(f"ERROR: var must be in {self.data.keys()}")
			return
	
	def time(self, t = None, psiNs = None, descriptor = None):
		if t is None or t < 0 or t > 1:
			print("ERROR: t must be specified between 0 and 1")
			return
		if psiNs is None:
			psiNs = self._psiNs
		geq = {}
		for cvar in self._constants:
			geq[cvar] = self.data[cvar]
		for svar in self._variables_sin:
			geq[svar] = self.data[svar](t)
		for avar in self._variables_arr:
			geq[avar] = zeros((len(psiNs)))
			for p, psiN in enumerate(psiNs):
				geq[avar][p] = self.data[avar](t,psiN)
		geq['psirz'] = zeros((len(geq['R']),len(geq['Z'])))
		for r, R in enumerate(geq['R']):
			for z, Z in enumerate(geq['Z']):
				geq['psirz'][r][z] = self.data['psirz']((t,R,Z))
		if descriptor is None:
			descriptor = f"Interpolated GEQDSK for t = {t}"
		geq['descriptor'] = descriptor
		return geq
	
	def write_geqdsk(self, t = None, filename = None):
		if t is None or t < 0 or t > 1:
			print("ERROR: t must be specified between 0 and 1")
			return
		a = self.time(t = t)
		a.write_geqdsk(filename = filename)
	
	def plot(self, var = 'beta_prime'):
		def draw_fig(val = None):
			val = eval(val)
			k = self.var(var = var, t = val)
			ax.cla()
			ax.plot(self._psiNs,k0,'r', label = "eq1")
			ax.plot(self._psiNs,k,'k--', label = "new eq")
			ax.plot(self._psiNs,k1,'b', label = "eq2")
			ax.legend(loc = 0)
			ax.set_title(var)
			fig.canvas.draw_idle()
		fig, ax = subplots()
		k0 = self.var(var = var, t = 0)
		k1 = self.var(var = var, t = 1)
		taxes = fig.add_axes([0.5, 0.01, 0.05, 0.04])
		tbox = TextBox(taxes, 'time:', initial = '0.5')
		tbox.on_submit(draw_fig)
		draw_fig("0.5")
		show()
	
	class geqdsk_time(object):
		def __init__(self, t = None, psiNs = None):
			self._psiNs = psiNs
			self.time = self.t = t 
		
		def write_geqdsk(self, filename = None):
			if filename is None:
				filename = f"p{self.t}.geqdsk"
			'''
			f = open(filename,'w')
			f.write(f"{len(self.psiN)} psinorm ne(10^20/m^3) dne/dpsiN\n")
			for p, psi in enumerate(self.psiN):
				differ = diff(self.ne)/diff(self.psiN)
				differ = insert(differ, 0,0)
				f.write(f" {psi:.7f}   {self.ne[p]:.7f}   {differ[p]:.7f}\n")
			f.write(f"{len(self.psiN)} psinorm te(KeV) dte/dpsiN\n")
			for p, psi in enumerate(self.psiN):
				differ = diff(self.ne)/diff(self.psiN)
				differ = insert(differ, 0,0)
				f.write(f" {psi:.7f}   {self.te[p]:.7f}   {differ[p]:.7f}\n")
			f.close()
			'''


