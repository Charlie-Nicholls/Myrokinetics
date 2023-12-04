from numpy import linspace, pi, zeros, real, imag
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.pyplot import figure, subplots, ion, show
from copy import deepcopy
from .slider_ax import slider_axes

default_settings = {"suptitle": None,
		"var": 'phi',
		"periods": 3,
		"polar": True,
		"run": {},
		"fontsizes": {"axis": 11,"title": 13,"suptitle": 20},
		"visible": {"suptitle": True, "title": True, 'absolute': True, 'real': True, 'imag': True},
}

slider_settings = {"slider_1": {"axis": [0.25, 0.01, 0.5, 0.03]},
		"slider_2": {"axis": [0.25, 0.05, 0.5, 0.03]},
		"slider_3": {"axis": [0.92, 0.2, 0.03, 0.5], "orientation": 'vertical'}, 
		"slider_4": {"axis": [0.96, 0.2, 0.03, 0.5], "orientation": 'vertical'},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "slider_5": False, "slider_6": False, "slider_7": False, "slider_8": False, "slider_9": False},
}

class plot_theta(object):
	def __init__(self, reader, settings = {}, sliders = None):
		self.reader = reader
		self.settings = {}
		self.sliders = sliders
		defaults = deepcopy(default_settings)
		for key in settings:
			if key not in defaults:
				print(f"ERROR: {key} not found")
			else:
				self.settings[key] = settings[key]
		for key in defaults:
			if key not in self.settings:
				self.settings[key] = defaults[key]
			elif type(self.settings[key]) == dict:
				for skey in defaults[key]:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
				
		if self['var'] == 0:
			self.settings['var'] = "phi"
		if self['var'] == 1:
			self.settings['var'] = "apar"
		if self['var'] == 2:
			self.settings['var'] = "bpar"
		if self['var'] == 3:
			self.settings['var'] = "epar"
		if self['var'] not in ['omega','phi','apar','bpar','epar','phi2']:
			print(f"ERROR: variable name/value {self['var']} not supported. supported: phi/0, apar/1, bpar/2, epar/3")
			return
			
		self.open_plot()
	
	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
	
	def save_plot(self, filename = None):
		if filename is None:
			filename = f"Theta_{self['var']}_{self['slider_1']['id']}_{self['slider_2']}_{self['slider_3']['id']}_{self['slider_4']['id']}"
		self.fig.savefig(filename)
	
	def open_plot(self):
		if self['polar']:
			self.fig = figure(figsize=(14.6,7))
			self.fig.add_subplot(121)
			self.fig.add_subplot(122,polar=True)
			self.ax = self.fig.get_axes()
		else:
			self.fig, self.ax = subplots(1,2,figsize=(14.6,7))
		
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])

		self.dims = self.reader.inputs.dim_order		
		for dim in [x for x in self.dims if x not in self.settings['run']]:
			self.settings['run'][dim] = self.reader.dimensions[dim].values[0]
		
		if self.sliders in [None,'seperate']:
			self.fig.subplots_adjust(bottom=0.15)
			if self.sliders == 'seperate':
				self.sliders = slider_axes(reader=self.reader)
			else:
				slider_defaults = deepcopy(slider_settings)
				slider_defaults['dims'] = self.dims
				self.sliders = slider_axes(reader=self.reader,settings=slider_defaults,ax=self.ax)
		self.sliders.add_plot(self)
		
		ion()
		show()
		self.draw_fig()
	
	def set_slider(self, num = None, key = None, dimension_type = None):
		self.sliders.set_slider(num = num, key = key, dimension_type = dimension_type)
	
	def set_visible(self, key, val = None):
		if key not in self['visible'] and key not in self.sliders['visible']:
			print(f"ERROR: key not found, valid keys: {list(self['visible'].keys())} {list(self.sliders['visible'].keys())}")
			return
		if 'slider_' in key:
			self.sliders.set_visible(key=key,val=val)
		else:
			if val not in [True,False]:
				val = not self['visible'][key]
			self.settings['visible'][key] = val
			if key == 'suptitle':
				self.fig._suptitle.set_visible(self['visible']['suptitle'])
			elif key == 'title':
				self.ax.legend_.set_visible(self['visible']['title'])
	
	def set_options_fontsize(self, fontsize):
		self.settings['fontsizes']['ch_box'] = fontsize
		for text in self.options.labels:
			text.set_fontsize(fontsize)
	
	def set_axis_fontsize(self, fontsize):
		self.settings['fontsizes']['axis'] = fontsize
		self.ax.xaxis.label.set_fontsize(fontsize)
		self.ax.yaxis.label.set_fontsize(fontsize)
		self.fig.canvas.draw_idle()
		
	def set_title_fontsize(self, fontsize):
		self.settings['fontsizes']['title'] = fontsize
		self.draw_fig()
		
	def set_suptitle_fontsize(self, fontsize):
		self.settings['fontsizes']['suptitle'] = fontsize
		self.fig._suptitle.set_fontsize(fontsize)
		
	def set_suptitle(self, title):
		self.settings['suptitle'] = title
		self.fig.suptitle(title,fontsize=self.settings['fontsizes']['suptitle'])
	
	def draw_fig(self, val = None):
		for key, sli in self.sliders.sliders.items():
			dim = self.sliders.settings[key]['dimension_type']
			if dim in self.dims:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				
		self.ax[0].cla()
		self.ax[1].cla()
		run = self['run']
		theta = linspace(0, 2*pi, 100)
		rt_fun = zeros((len(theta)))
		it_fun = zeros((len(theta)))
		re_fun = InterpolatedUnivariateSpline(self.reader('theta',run),real(self.reader(self['var'],run)))
		ie_fun = InterpolatedUnivariateSpline(self.reader('theta',run),imag(self.reader(self['var'],run)))
		for t, th in enumerate(theta):
			sampling = True
			i = 0
			while sampling:
				th_p = self.reader('theta',run)[0] + th + 2*pi*i
				if self.reader('theta',run)[0] <= th_p <= self.reader('theta',run)[-1]:
					rt_fun[t] += re_fun(th_p)
					it_fun[t] += ie_fun(th_p)
				elif th_p > self.reader('theta',run)[-1]:
					sampling = False
				i = i+1
		
		if self['polar']:
			if self['visible']['real']:
				self.ax[1].plot(theta, rt_fun,'r--',label="real")
			if self['visible']['imag']:
				self.ax[1].plot(theta, it_fun,'b--',label="imaginary")
			if self['visible']['absolute']:
				self.ax[1].plot(theta,[(x**2 + y**2)**0.5 for x,y in zip(rt_fun,it_fun)],'k--',label="absolute")
		else:
			for i in range(self['periods']):
				if self['visible']['real']:
					self.ax[1].plot(theta+i*2*pi, rt_fun,'b--',label="real")
				if self['visible']['imag']:
					self.ax[1].plot(theta+i*2*pi, it_fun,'r--',label="imaginary")
				if self['visible']['absolute']:
					self.ax[1].plot(theta+i*2*pi, [(x**2 + y**2)**0.5 for x,y in zip(rt_fun,it_fun)],'k--',label="absolute")
			self.ax[1].set_xlabel("Theta")
			self.ax[1].set_ylabel(self['var'])
		if self['visible']['real']:
			self.ax[0].plot(self.reader('theta',run), real(self.reader(self['var'],run)),'r--',label="real")
		if self['visible']['imag']:
			self.ax[0].plot(self.reader('theta',run), imag(self.reader(self['var'],run)),'b--',label="imaginary")
		if self['visible']['absolute']:
			self.ax[0].plot(self.reader('theta',run), [abs(x) for x in self.reader(self['var'],run)],'k--',label="absolute")
		self.ax[0].legend(loc=0)
		self.ax[0].set_xlabel("Ballooning Angle")
		self.ax[0].set_ylabel(self['var'])
		self.fig.canvas.draw_idle()
