from matplotlib.pyplot import show, ion, subplots, axes
from matplotlib.widgets import Slider
from numpy import log, nan, zeros, array, polyfit, exp
from numpy import sum as np_sum
from scipy.stats import pearsonr

class plot_nl_phi2:
	def __init__(self, reader, settings = {}):
		self.reader = reader
		self.settings = {'y_axis_type': 'phi2_by_ky', 'plot_type': 'group', 'fit': True, "fontsizes": {"legend": 10,"axis": 11,"title": 13,"suptitle": 20},"suptitle": None,}
		self.settings['run'] = self.reader.get_all_runs()[0]
		for key, val in settings.items():
			if key in self.settings:
				self.settings[key] = val
		self.open_plot()

	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")

	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(8.8,5.8))
		if self.settings['suptitle']:
			self.fig.suptitle(self.settings['suptitle'],fontsize=self.settings['fontsizes']['suptitle'])
		if 'ky' in self['y_axis_type']:
			var = 'ky'
		elif 'kx' in self['y_axis_type']:
			var = 'kx'
		if self['plot_type'] == 'single':
			self.slider = Slider(axes([0.25, 0.01, 0.5, 0.03]), f"{var} index:", 0, len(self.reader[var])-1, valinit = 0, valstep = 1)
			self.slider.on_changed(self.draw_fig)
		self._tmin = 5
		self._tmax = 7.5
		ion()
		show()
		self.draw_fig()
	
	def set_plot_type(self, plot_type):
		if plot_type not in ['single','group','mesh']:
			print(f"ERROR: {plot_type} not valid, valid = single, group, mesh")
			return
		self.settings['plot_type'] = plot_type
		self.draw_fig()

	def toggle_fit(self):
		self.settings['fit'] = not self.settings['fit']
		self.draw_fig()

	def draw_fig(self, val = None):
		self.ax.cla()
		if 'ky' in self['y_axis_type']:
			var = 'ky'
		elif 'kx' in self['y_axis_type']:
			var = 'kx'
		self.x = self.reader('t',self['run'])
		phi2 = array(self.reader(self['y_axis_type'],self['run'])).T
		if self['plot_type'] == 'single':
			sl_id = self.slider.val
			sl = self.reader[var][sl_id]
			self.y = phi2[sl_id,:]
			self.ax.plot(self.x,self.y)
			self.ax.set_yscale('log')
			self.ax.set_title(f"{self['y_axis_type']} | {var} = {sl}",fontsize=self['fontsizes']['title'])
			self.ax.set_ylabel("$\phi^{2}$",fontsize=self['fontsizes']['axis'])
			if self['fit']:
				lt = int(self._tmin/(self.x[1]-self.x[0]))
				ut = int(self._tmax/(self.x[1]-self.x[0])) + 1
				fit = polyfit(self.x[lt:ut],log(self.y[lt:ut]),1)
				self._gr = fit[0]/2
				pr = pearsonr(self.x[lt:ut],log(self.y[lt:ut]))
				self.ax.plot(self.x[lt:ut],exp(array(self.x[lt:ut])*fit[0] + fit[1]),'r',label=f"GR = {self._gr:+.2e}\nR = {pr[0]:.4f}")
				self.ax.legend(loc=2,fontsize=self['fontsizes']['legend'])
		elif self['plot_type'] == 'group':
			self.y = phi2
			self.ax.plot(self.x,array(self.y).T)
			self.ax.set_yscale('log')
			self.ax.set_title(self['y_axis_type'],fontsize=self['fontsizes']['title'])
			self.ax.set_ylabel("$\phi^{2}$",fontsize=self['fontsizes']['axis'])
		elif self['plot_type'] == 'mesh':
			self.y = self.reader[var]
			self.z = log(array(phi2)).tolist()
			self.ax.pcolormesh(self.x,self.y,self.z)
			self.ax.set_title(f"log({self['y_axis_type']})",fontsize=self['fontsizes']['title'])
			self.ax.set_ylabel(f"{var}",fontsize=self['fontsizes']['axis'])
			
		self.ax.set_xlabel(f"Time ({len(self.x)} steps)",fontsize=self['fontsizes']['axis'])
			
class plot_phi2_by_mode:
	def __init__(self, reader, settings = {}):
		self.reader = reader
		self.settings = {'avg': 1}
		self.settings['run'] = self.reader.get_all_runs()[0]
		for key, val in settings.items():
			if key in self.settings:
				self.settings[key] = val
		self.open_plot()

	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")

	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(8.8,5.8))
		self.t = self.reader('t',self['run'])
		self.slider = Slider(axes([0.25, 0.01, 0.5, 0.03]), f"time index:", self['avg']-1, len(self.t)-1, valinit = len(self.t)-1, valstep = 1)
		self.slider.on_changed(self.draw_fig)
		ion()
		show()
		self.draw_fig()
	
	def draw_fig(self, val = None):
		self.ax.cla()
		t = self.slider.val
		kxs = self.reader['kx']
		kys = self.reader['ky']
		phi2s = zeros((len(kxs),len(kys)))
		for kyid, ky in enumerate(kys):
			for kxid, kx in enumerate(kxs):
				run = {'ky': ky, 'kx': kx}
				phi2s[kxid,kyid] = log(sum(self.reader('phi2',run)[int(t-self['avg']):int(t)])/self['avg'])
		self.x = kxs
		self.y = kys
		self.z = phi2s
		self.ax.pcolormesh(kxs,kys,phi2s.T)
		self.ax.set_title(f"log(<phi2>)")
		self.ax.set_xlabel("kx")
		self.ax.set_ylabel("ky")

class plot_hflux:
	def __init__(self, reader, settings = {}):
		self.reader = reader
		self.settings = {}
		self.settings['run'] = self.reader.get_all_runs()[0]
		for key, val in settings.items():
			if key in self.settings:
				self.settings[key] = val
		self.open_plot()

	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")

	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(8.8,5.8))
		ion()
		show()
		self.draw_fig()
	
	def draw_fig(self):
		self.ax.cla()
		self.x = self.reader('t',self['run'])
		self.y = self.reader('heat_flux_tot',self['run'])
		self.ax.plot(self.x,self.y)
		self.ax.set_xlabel(f"Time ({len(self.x)} steps)")
		self.ax.set_ylabel("Heat Flux")
		self.ax.set_yscale('log')

class plot_zonality:
	def __init__(self, reader, settings = {}):
		self.reader = reader
		self.settings = {}
		self.settings['run'] = self.reader.get_all_runs()[0]
		for key, val in settings.items():
			if key in self.settings:
				self.settings[key] = val
		self.open_plot()

	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")

	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(8.8,5.8))
		ion()
		show()
		self.draw_fig()
	
	def draw_fig(self):
		self.ax.cla()
		self.x = self.reader('t',self['run'])
		self.y = self.reader('phi2_by_ky',self['run'])[0]/np_sum(self.reader('phi2_by_ky',self['run'])[1:],axis=0)
		self.y = [y for yi, y in enumerate(self.y) if self.x[yi] not in [None,nan]]
		self.x = [x for x in self.x if x not in [None,nan]]
		self.ax.plot(self.x,self.y)
		self.ax.set_xlabel(f"Time ({len(self.x)} steps)")
		self.ax.set_ylabel("Zonality")
