from numpy import linspace, pi, zeros, real, imag
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.pyplot import *
from copy import deepcopy

default_settings = {"suptitle": None,
		"var": 'phi',
		"periods": 3,
		"polar": True,
		"slider_1": {"dimension_type": None, "id": 0},
		"slider_2": {"dimension_type": None, "id": 0},
		"slider_3": {"dimension_type": None, "id": 0},
		"slider_4": {"dimension_type": None, "id": 0},
		"run": {},
		"fontsizes": {"axis": 11,"title": 13,"suptitle": 20},
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "suptitle": True, "title": True},
}

class plot_theta(object):
	def __init__(self, reader, settings = {}):
		self.reader = reader
		
		self.settings = {}
		self.init_settings = {}
		defaults = deepcopy(default_settings)
		for key in settings:
			if key not in defaults:
				print(f"ERROR: {key} not found")
			else:
				self.settings[key] = settings[key]
				self.init_settings[key] = settings[key]
		for key in defaults:
			if key not in self.settings:
				self.settings[key] = defaults[key]
				self.init_settings[key] = defaults[key]
			elif type(self.settings[key]) == dict:
				for skey in defaults[key]:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
						self.init_settings[key][skey] = defaults[key][skey]
				
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
		used_dims = [self.settings[key]['dimension_type'] for key in self.settings.keys() if 'slider_' in key]
		unused_dims = [x for x in self.dims if x not in used_dims]
		slider_keys = [x for x in self.settings if 'slider_' in x]
		empty_sliders = [x for x in slider_keys if self.settings[x]['dimension_type'] is None]
		for sli, key in enumerate(empty_sliders):
			if len(unused_dims) > sli:
				self.settings[key]['dimension_type'] = unused_dims[sli]
				used_dims.append(unused_dims[sli])
			else:
				self.settings[key]['dimension_type'] = None
				self.settings['visible'][key] = False
				
		for dim in [x for x in self.dims if x not in self.settings['run']]:
			self.settings['run'][dim] = self.reader.dimensions[dim].values[0]
		
		self._slider_axes = {'slider_1': axes([0.25, 0.01, 0.5, 0.03],visible=self['visible']['slider_1']),
		'slider_2': axes([0.25, 0.05, 0.5, 0.03],visible=self['visible']['slider_2']),
		'slider_3': axes([0.92, 0.2, 0.03, 0.5],visible=self['visible']['slider_3']),
		'slider_4': axes([0.96, 0.2, 0.03, 0.5],visible=self['visible']['slider_4']),
		}
		self.sliders = {}
		for key in [x for x in self._slider_axes.keys() if self.settings[x]['dimension_type'] is not None]:
			dim = self.reader.dimensions[self.settings[key]['dimension_type']]
			if key in ['slider_1','slider_2']:
				orient = 'horizontal'
			elif key in ['slider_3','slider_4']:
				orient = 'vertical'
			self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index", 0, len(dim)-1, valinit = self[key]['id'], valstep = 1, orientation = orient)
			self.sliders[key].on_changed(self.draw_fig)
			self.set_visible(key,val=self['visible'][key])
			if orient == 'vertical':
				self.sliders[key].label.set_rotation(90)
		
		ion()
		show()
		self.draw_fig()
	
	def set_slider(self, slider_num, dimension_type, visible = True):
		if dimension_type not in self.dims:
			print(f"ERROR: invalid dimension type, valid: {self.dims}")
			return
		key = f"slider_{slider_num}"
		self.settings[key]['dimension_type'] = dimension_type
		self.settings[key]['id'] = 0
		dim = self.reader.dimensions[dimension_type]
		self.sliders[key].ax.remove()
		self._slider_axes[key] = axes(slider_axes[key])
		if key in ['slider_1','slider_2']:
			orient = 'horizontal'
		elif key in ['slider_3','slider_4']:
			orient = 'vertical'
		self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index", 0, len(dim)-1, valinit = 0, valstep = 1, orientation = orient)
		self.sliders[key].on_changed(self.draw_fig)
		if orient == 'vertical':
			self.sliders[key].label.set_rotation(90)
		self.set_visible(key,val=visible)
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys: {self.settings['visible'].keys()}")
			return
		if val not in [True,False]:
			val = not self['visible'][key]
		
		self.settings['visible'][key] = val
		
		if 'slider_' in key:
			self._slider_axes[key].set_visible(self['visible'][key])
			
			if self['visible']['slider_1'] == False and self['visible']['slider_2'] == False:
				self.fig.subplots_adjust(bottom=0.11)
			elif self['visible']['slider_2'] == False: 
				self.fig.subplots_adjust(bottom=0.13)
			else:
				self.fig.subplots_adjust(bottom=0.15)
		elif key == 'suptitle':
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
		self.ax[0].cla()
		self.ax[1].cla()
		for key, sli in self.sliders.items():
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				self.settings[key]['id'] = sli.val
		
		run_id = self.reader.get_run_id(run=self['run'])
		data = self.reader.gyro_data[run_id]
		
		theta = linspace(0, 2*pi, 100)
		rt_fun = zeros((len(theta)))
		it_fun = zeros((len(theta)))
		re_fun = InterpolatedUnivariateSpline(data['theta'],real(data[self['var']]))
		ie_fun = InterpolatedUnivariateSpline(data['theta'],imag(data[self['var']]))
		for t, th in enumerate(theta):
			sampling = True
			i = 0
			while sampling:
				th_p = data['theta'][0] + th + 2*pi*i
				if data['theta'][0] <= th_p <= data['theta'][-1]:
					rt_fun[t] += re_fun(th_p)
					it_fun[t] += ie_fun(th_p)
				elif th_p > data['theta'][-1]:
					sampling = False
				i = i+1
		
		if self['polar']:
			self.ax[1].plot(theta, rt_fun,'b--')
			self.ax[1].plot(theta, it_fun,'r--')
		else:
			for i in range(self['periods']):
				self.ax[1].plot(theta+i*2*pi, rt_fun,'b--')
				self.ax[1].plot(theta+i*2*pi, it_fun,'r--')
			self.ax[1].set_xlabel("Theta")
			self.ax[1].set_ylabel(self['var'])
		self.ax[0].plot(data['theta'], real(data[self['var']]),'b--')
		self.ax[0].plot(data['theta'], imag(data[self['var']]),'r--')
		self.ax[0].set_xlabel("Ballooning Angle")
		self.ax[0].set_ylabel(self['var'])
		self.fig.canvas.draw_idle()
