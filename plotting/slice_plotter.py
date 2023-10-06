from matplotlib.pyplot import *
from matplotlib.widgets import Slider
from numpy import array, full, nan
from copy import deepcopy

default_settings = {"suptitle": None,
		"x_axis_type": "beta_prime",
		"y_axis_type": "quasilinear",
		"slider_1": {"dimension_type": None, "id": 0},
		"slider_2": {"dimension_type": None, "id": 0},
		"slider_3": {"dimension_type": None, "id": 0},
		"run": {},
		"limit": None,
		"fontsizes": {"title": 13, "axis": 17,"suptitle": 20},
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "eqbm": True, "suptitle": True, "title": True},
		"colours": {"eqbm": 'k', "points": 'k', "line": 'r'}
}

class plot_slice(object):
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
			if key not in settings:
				self.settings[key] = defaults[key]
				self.init_settings[key] = defaults[key]
			elif type(self.settings[key]) == dict:
				for skey in defaults[key]:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
						self.init_settings[key][skey] = defaults[key][skey]
		
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
		
		self._axes_vals = {'slider_1':[0.25, 0.01, 0.5, 0.03],
			'slider_2': [0.25, 0.05, 0.5, 0.03],
			'slider_3': [0.92, 0.2, 0.03, 0.5],
		}
		self._slider_axes = {'slider_1': axes(self._axes_vals['slider_1'],visible=self['visible']['slider_1']),
			'slider_2': axes(self._axes_vals['slider_2'],visible=self['visible']['slider_2']),
			'slider_3': axes(self._axes_vals['slider_3'],visible=self['visible']['slider_3']),
		}
		self.dims = [x for x in self.reader.inputs.dim_order if x not in [self['x_axis_type']]]
		
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		
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
		
		self.sliders = {}
		for key in [x for x in self._slider_axes.keys() if self.settings[x]['dimension_type'] is not None]:
			self.set_slider(key = key, dimension_type = self.settings[key]['dimension_type'], _drawfig = False)

		ion()
		show()
		self.draw_fig()
		
	def set_slider(self, num = None, key = None, dimension_type = None, _drawfig = True):
		if dimension_type not in self.dims:
			print(f"ERROR: invalid dimension type, valid: {self.dims}")
			return
		if num is None and key is None:
			print("ERROR: num or key must be given")
		if key is None:	
			key = f"slider_{slider_num}"
		self.settings[key]['dimension_type'] = dimension_type
		self.settings[key]['id'] = 0
		dim = self.reader.dimensions[dimension_type]
		if self._slider_axes[key]:
			self._slider_axes[key].set_visible(False)
			del(self._slider_axes[key])
		self._slider_axes[key] = axes(self._axes_vals[key],visible=self['visible'][key])
		if key in ['slider_1','slider_2']:
			orient = 'horizontal'
		elif key in ['slider_3','slider_4']:
			orient = 'vertical'
		self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index:", 0, len(dim)-1, valinit = 0, valstep = 1, orientation = orient)
		if orient == 'vertical':
			self.sliders[key].label.set_rotation(90)
		self.sliders[key].on_changed(self.draw_fig)
		self.set_visible(key,val=self['visible'][key])
		if _drawfig:
			self.draw_fig()
	
	def _load_x_axis(self, axis_type):
		if axis_type not in self.reader.dimensions:
			print(f"ERROR: axis_type not found, valid types {self.reader.dimensions}")
			return
		
		if self.settings['x_axis_type'] not in self.dims:
			self.dims.append(self.settings['x_axis_type'])
		self.dims = [x for x in self.dims if x != axis_type]
		self.settings['x_axis_type'] = axis_type
		self.x_axis = self.reader.dimensions[axis_type].values
		self._x_axis_label =  self.reader.dimensions[axis_type].axis_label
		for key in [x for x in self.settings if 'slider_' in x]:
			if self.settings[key]['dimension_type'] is not None and self.settings[key]['dimension_type'] not in self.dims:
				self.settings[key]['dimension_type'] = None
				self.settings[key]['id'] = 0
				self.set_visible(key,False)
			
	def set_x_axis_type(self, axis_type):
		self._load_x_axis(axis_type)
		self.draw_fig()
	
	def _load_y_axis(self, axis_type):
		if axis_type not in ['quasilinear','growth_rate','growth_rate_norm','ql_norm']:
			print(f"ERROR: axis_type not found, valid types ['quasilinear','growth_rate','growth_rate_norm',ql_norm]")
			return
			
		self.settings['y_axis_type'] = axis_type
		
		if axis_type == 'quasilinear':
			self._y_axis_label = "Quasilinear Metric"
			self._y_key = '_quasilinear_keys'
			self.dims = [x for x in self.dims if x not in ['theta0','ky']]
			if 'ky' in self['run']:
				self.settings['run'].pop('ky')
			if 'theta0' in self['run']:
				self.settings['run'].pop('theta0')
		elif axis_type in ['growth_rate','growth_rate_norm','ql_norm']:
			self._y_key = '_run_keys'
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
			
		
		for key in [x for x in self.settings if 'slider_' in x]:
			if self.settings[key]['dimension_type'] is not None and  self.settings[key]['dimension_type'] not in self.dims:
				self.settings[key]['dimension_type'] = None
				self.settings[key]['id'] = 0
				self.set_visible(key,False)
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys {self.settings['visible'].keys()}")
			return
		if val not in [True,False]:
			val = not self['visible'][key]
		
		self.settings['visible'][key] = val
		
		if 'slider_' in key:
			self._slider_axes[key].set_visible(self['visible'][key])
		elif key == 'suptitle':
			self.fig._suptitle.set_visible(self['visible']['suptitle'])
		elif key == 'title':
			self.ax.legend_.set_visible(self['visible']['title'])
			
		if not self['visible']['slider_1'] and not self['visible']['slider_2']:
			self.fig.subplots_adjust(bottom=0.11)
		elif not self['visible']['slider_2']: 
			self.fig.subplots_adjust(bottom=0.13)
		else:
			self.fig.subplots_adjust(bottom=0.15)
	
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
		for key in [x for x in self.sliders.keys() if x != 'ql_slider']:
			sli = self.sliders[key]
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				self.settings[key]['id'] = sli.val
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label=f"{self.reader.dimensions[dim].axis_label} = {self.settings['run'][dim]}",visible = False))
		x_axis = list(self.x_axis)
		
		self.ax.cla()
		self.ax.set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax.set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		
		y_vals = full((len(self.reader.dimensions[self['x_axis_type']])),nan)
		for x_id, x_value in enumerate(self.x_axis):
			run = self['run'].copy()
			run[self['x_axis_type']] = x_value
			run_id = self.reader.get_run_id(run = run, keys = self._y_key)
			if self['y_axis_type'] == 'quasilinear':
				y_vals[x_id] = self.reader.data['quasilinear'][run_id]
			else:
				y_vals[x_id] = self.reader.data['gyro'][run_id][self['y_axis_type']]
		self.y_axis = y_vals
		
		
		self.ax.plot(self.x_axis,self.y_axis,c=self['colours']['line'])
		self.ax.plot(self.x_axis,self.y_axis,'.',c=self['colours']['points'])
		
		if self['visible']['eqbm']:
			limits = self.ax.get_ylim()
			if 'psin' in self['run']:
				psiN = self['run']['psin']
			else:
				psiN = self.reader.single_parameters.values[0]
			if self['x_axis_type'] in self.reader.data['equilibrium'][psiN]:
				eqbm_val = abs(self.reader.data['equilibrium'][psiN][self['x_axis_type']])
				self.ax.vlines(eqbm_val,0,limits[1],self['colours']['eqbm'])
				handles.append(Line2D([0.5,0.5],[0,1],c=self['colours']['eqbm'],label = "Equillibrium"))
		
		if self['limit']:
			self.ax.set_ylim(0,self['limit'])
		
		self.ax.legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0.5,0.98),loc = "lower center", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax.legend_.set_visible(self['visible']['title'])
		
		self.fig.canvas.draw_idle()

