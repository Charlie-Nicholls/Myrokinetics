from matplotlib.pyplot import subplots, ion, show, Line2D
from numpy import full, nan
from copy import deepcopy
from .slider_ax import slider_axes

default_settings = {"suptitle": None,
		"x_axis_type": "beta_prime",
		"y_axis_type": "growth_rate",
		"xscale": "linear",
		"yscale": "linear",
		"run": {},
		"limit": None,
		"ref_line": {"x_axis": [], "y_axis": []},
		"x_lim": [None, None],
		"y_lim": [None, None],
		"fontsizes": {"title": 13, "axis": 17,"suptitle": 20},
		"visible": {"eqbm": True, "suptitle": True, "title": True, "ref_line": False},
		"colours": {"eqbm": 'k', "points": 'k', "line": 'r', "ref_line": 'b', "ref_points": 'k'},
}

slider_settings = {"slider_1": {"axis": [0.25, 0.01, 0.5, 0.03]},
		"slider_2": {"axis": [0.25, 0.05, 0.5, 0.03]},
		"slider_3": {"axis": [0.92, 0.2, 0.03, 0.5], "orientation": 'vertical'},
		"slider_4": {"axis": [0.96, 0.2, 0.03, 0.5], "orientation": 'vertical'},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "slider_5": False, "slider_6": False, "slider_7": False, "slider_8": False, "slider_9": False},
}

class plot_slice(object):
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
			if key not in settings:
				self.settings[key] = defaults[key]
			elif type(defaults[key]) == dict:
				for skey in defaults[key].keys():
					if skey not in settings[key].keys():
						self.settings[key][skey] = defaults[key][skey]
						
		if self.settings['x_axis_type'] not in self.reader.dimensions:
			self.settings['x_axis_type'] = self.reader.inputs.dim_order[0]
		
		self.open_plot()
		
	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
		
	def save_plot(self, filename = None):
		if filename is None:
			filename = f"Slice_{self['slider_1']['id']}_{self['slider_2']['id']}"
		self.fig.savefig(filename)
		
	def open_plot(self, save = False, filename = None):
		self.fig, self.ax = subplots(figsize=(9,7))
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		self.dims = [x for x in self.reader.inputs.dim_order if x not in [self['x_axis_type']]]
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
		
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])

		ion()
		show()
		self.temp = False
		self.settings['visible']['ref_line'] = True
		self.draw_fig()
		
	def set_slider(self, num = None, key = None, dimension_type = None):
		self.sliders.set_slider(num = num, key = key, dimension_type = dimension_type)
	
	def _load_x_axis(self, axis_type):
		if axis_type not in self.reader.dimensions:
			print(f"ERROR: axis_type not found, valid types {self.reader.dimensions.keys()}")
			return
		if self.settings['x_axis_type'] not in self.dims:
			self.dims.append(self.settings['x_axis_type'])
		self.dims = [x for x in self.dims if x != axis_type]
		self.sliders.settings['dims'] = [x for x in self.sliders['dims'] if x != axis_type]
		self.settings['x_axis_type'] = axis_type
		self.x_axis = self.reader.dimensions[axis_type].values
		self._x_axis_label =  self.reader.dimensions[axis_type].axis_label
		for key in [x for x in self.sliders.settings if 'slider_' in x]:
			if self.sliders.settings[key]['dimension_type'] is not None and self.sliders.settings[key]['dimension_type'] not in self.dims:
				self.sliders.remove_slider(key=key)
			
	def set_x_axis_type(self, axis_type):
		self._load_x_axis(axis_type)
		self.draw_fig()
	
	def _load_y_axis(self, axis_type):
		if axis_type not in ['quasilinear','growth_rate','growth_rate_norm','ql_norm','ql_metric','mode_frequency']:
			print("ERROR: axis_type not found, valid types: 'quasilinear','growth_rate','growth_rate_norm','ql_norm','ql_metric','mode_frequency")
			return
			
		self.settings['y_axis_type'] = axis_type
		
		if axis_type == 'quasilinear':
			self._y_axis_label = "Quasilinear Metric"
			self._y_key = '_quasilinear_keys'
			self.dims = [x for x in self.dims if x not in ['ky']]
			if 'ky' in self['run']:
				self.settings['run'].pop('ky')
		elif axis_type in ['growth_rate','growth_rate_norm','ql_norm','ql_metric','mode_frequency']:
			self._y_key = '_gyro_keys'
			if 'ky' in self.reader.dimensions and 'ky' not in self.dims:
				self.dims.append('ky')
			if 'ky' in self.reader.dimensions and 'ky' not in self['run']:
				self.settings['run']['ky'] = self.reader.dimensions['ky'].values[0]
			if 'theta0' in self.reader.dimensions and 'theta0' not in self.dims:
				self.dims.append('theta0')
			if 'theta0' in self.reader.dimensions and 'theta0' not in self['run']:
				self.settings['run']['theta0'] = self.reader.dimensions['theta0'].values[0]
			if axis_type == 'growth_rate':
				self._y_axis_label = "Growth Rate"
			elif axis_type == 'growth_rate_norm':
				self._y_axis_label = "Growth Rate/$(k_{y}\\rho_{0})^{2}$"
			elif axis_type == 'ql_norm':
				self._y_axis_label = "QL Normalised Growth Rate"
			elif axis_type == 'ql_metric':
				self._y_axis_label = "Quasilinear Metric"
			elif axis_type == 'mode_frequency':
				self._y_axis_label = "Mode Frequency"
		
		for key in [x for x in self.settings if 'slider_' in x]:
			if self.sliders.settings[key]['dimension_type'] is not None and self.sliders.settings[key]['dimension_type'] not in self.dims:
				self.sliders.remove_slider(key=key)
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def set_xscale(self, scale):
		self.ax.set_xscale(scale)
		self.settings['xscale'] = scale
		
	def set_yscale(self, scale):
		self.ax.set_yscale(scale)
		self.settings['yscale'] = scale
	
	def set_ylim(self, limit=[None,None]):
		if type(limit) != list or len(limit) != 2:
			print("ERROR: limit must be list of length 2 [lower,upper] (None for no limit)")
			return
		self.settings['y_lim'] = limit
		self.draw_fig()
		
	def set_xlim(self, limit=[None,None]):
		if type(limit) != list or len(limit) != 2:
			print("ERROR: limit must be list of length 2 [lower,upper] (None for no limit)")
			return
		self.settings['x_lim'] = limit
		self.draw_fig()
	
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
		self.draw_fig()
	
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
		handles = []
		for key, sli in self.sliders.sliders.items():
			dim = self.sliders.settings[key]['dimension_type']
			if dim in self.dims:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label=f"{self.reader.dimensions[dim].axis_label} = {self.settings['run'][dim]}",visible = False))
		
		self.ax.cla()
		self.ax.set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax.set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		
		x_vals = self.reader.dimensions[self['x_axis_type']].values
		if self['x_lim'][0] is not None:
			x_vals = [x for x in x_vals if x >= self['x_lim'][0]]
		if self['x_lim'][1] is not None:
			x_vals = [x for x in x_vals if x <= self['x_lim'][1]]
		y_vals = full((len(x_vals)),nan)
		for x_id, x_value in enumerate(x_vals):
			run = self['run'].copy()
			run[self['x_axis_type']] = x_value
			y_vals[x_id] = self.reader(self['y_axis_type'],run)

		self.x_axis = [x for i, x in enumerate(x_vals) if str(y_vals[i]) not in ['nan','inf','-inf']]
		self.y_axis = [y for y in y_vals if str(y) not in ['nan','inf','-inf']]

		self.ax.plot(self.x_axis,self.y_axis,c=self['colours']['line'])
		self.ax.plot(self.x_axis,self.y_axis,'.',c=self['colours']['points'])
		
		if self['visible']['ref_line']:
			if self.temp:
				ge = list(self.reader.data['_ql_fs'].keys())[-1]
				kyid = self.reader['ky'].index(self['run']['ky'])
				t0max = self.reader.data['_ql_fs'][ge]['theta0_maxs'][kyid]
				self.settings['ref_line']['x_axis'] = [x for x in self.x_axis if x <= t0max]
				self.settings['ref_line']['y_axis'] = [self.y_axis[xi] for xi, x in enumerate(self.settings['ref_line']['x_axis'])]
				self.ax.plot(self['ref_line']['x_axis'],self['ref_line']['y_axis'],c=self['colours']['ref_line'])
				self.ax.vlines(t0max,min(self.y_axis),max(self.y_axis))
			else:
				self.ax.plot(self['ref_line']['x_axis'],self['ref_line']['y_axis'],c=self['colours']['ref_line'])
				self.ax.plot(self['ref_line']['x_axis'],self['ref_line']['y_axis'],'.',c=self['colours']['ref_points'])
		
		if self['visible']['eqbm']:
			limits = self.ax.get_ylim()
			low_lim = 0 if min(limits) > 0 else min(limits)
			high_lim = max(limits) if max(limits) > 0 else 0
			if 'psin' in self['run']:
				psiN = self['run']['psin']
			else:
				psiN = self.reader.single_parameters['psin'].values[0]
			if self['x_axis_type'] in self.reader.data['equilibrium'][psiN]:
				eqbm_val = abs(self.reader.data['equilibrium'][psiN][self['x_axis_type']])
				self.ax.vlines(eqbm_val,low_lim,high_lim,self['colours']['eqbm'])
				handles.append(Line2D([0.5,0.5],[0,1],c=self['colours']['eqbm'],label = "Equillibrium"))
		
		if self['y_lim'] != [None,None]:
			limits = list(self.ax.get_ylim())
			if self['y_lim'][0] is not None:
				limits[0] = self['y_lim'][0]
			if self['y_lim'][1] is not None:
				limits[1] = self['y_lim'][1]
			self.ax.set_ylim(limits[0],limits[1])
		self.ax.set_xscale(self['xscale'])
		self.ax.set_yscale(self['yscale'])
		
		self.ax.legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0.5,0.98),loc = "lower center", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax.legend_.set_visible(self['visible']['title'])
		
		self.fig.canvas.draw_idle()
