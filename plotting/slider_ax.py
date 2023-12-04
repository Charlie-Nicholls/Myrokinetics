from matplotlib.pyplot import subplots, axes, ion, show
from matplotlib.widgets import Slider
from copy import deepcopy

default_settings = {"slider_1": {"dimension_type": None, "id": 0, "axis": [0.2,0.9,0.7,0.05], "orientation": 'horizontal'},
		"slider_2": {"dimension_type": None, "id": 0, "axis": [0.2,0.8,0.7,0.05], "orientation": 'horizontal'},
		"slider_3": {"dimension_type": None, "id": 0, "axis": [0.2,0.7,0.7,0.05], "orientation": 'horizontal'},
		"slider_4": {"dimension_type": None, "id": 0, "axis": [0.2,0.6,0.7,0.05], "orientation": 'horizontal'},
		"slider_5": {"dimension_type": None, "id": 0, "axis": [0.2,0.5,0.7,0.05], "orientation": 'horizontal'},
		"slider_6": {"dimension_type": None, "id": 0, "axis": [0.2,0.4,0.7,0.05], "orientation": 'horizontal'},
		"slider_7": {"dimension_type": None, "id": 0, "axis": [0.2,0.3,0.7,0.05], "orientation": 'horizontal'},
		"slider_8": {"dimension_type": None, "id": 0, "axis": [0.2,0.2,0.7,0.05], "orientation": 'horizontal'},
		"slider_9": {"dimension_type": None, "id": 0, "axis": [0.2,0.1,0.7,0.05], "orientation": 'horizontal'},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "slider_5": True, "slider_6": True, "slider_7": True, "slider_8": True, "slider_9": True},
}

class slider_axes(object):
	def __init__(self, reader, settings = {}, ax = None):
		self.reader = reader
		self.settings = {}
		self.ax = ax
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
		self.plots = set()
		self.open_plot()
		
	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
		
	def open_plot(self, save = False, filename = None):
		if self.ax is None:
			self.fig, self.ax = subplots()
			self.ax.set_axis_off()
		self.axes = {'slider_1': axes(self.settings['slider_1']['axis'],visible=self['visible']['slider_1']),
			'slider_2': axes(self.settings['slider_2']['axis'],visible=self['visible']['slider_2']),
			'slider_3': axes(self.settings['slider_3']['axis'],visible=self['visible']['slider_3']),
			'slider_4': axes(self.settings['slider_4']['axis'],visible=self['visible']['slider_4']),
			'slider_5': axes(self.settings['slider_5']['axis'],visible=self['visible']['slider_5']),
			'slider_6': axes(self.settings['slider_6']['axis'],visible=self['visible']['slider_6']),
			'slider_7': axes(self.settings['slider_7']['axis'],visible=self['visible']['slider_7']),
			'slider_8': axes(self.settings['slider_8']['axis'],visible=self['visible']['slider_8']),
			'slider_9': axes(self.settings['slider_9']['axis'],visible=self['visible']['slider_9']),
		}
		
		if self['dims'] is None:
			self.settings['dims'] = [x for x in self.reader.inputs.dim_order]
		used_dims = [self.settings[key]['dimension_type'] for key in self.settings.keys() if 'slider_' in key]
		unused_dims = [x for x in self['dims'] if x not in used_dims]
		slider_keys = [x for x in self.settings if 'slider_' in x]
		empty_sliders = [x for x in slider_keys if self.settings[x]['dimension_type'] is None]
		for sli, key in enumerate(empty_sliders):
			if len(unused_dims) > sli:
				self.settings[key]['dimension_type'] = unused_dims[sli]
				used_dims.append(unused_dims[sli])
			else:
				self.settings[key]['dimension_type'] = None
				self.set_visible(key,False)
		
		self.sliders = {}
		for key in [x for x in self.axes.keys() if self.settings[x]['dimension_type'] is not None]:
			self.set_slider(key = key, dimension_type = self.settings[key]['dimension_type'], _drawfig = False)

		ion()
		show()
		self.draw_fig()
		
	def set_slider(self, num = None, key = None, dimension_type = None, _drawfig = True):
		if dimension_type not in self['dims']:
			print(f"ERROR: invalid dimension type, valid: {self['dims']}, given: {dimension_type}")
			return
		if num is None and key is None:
			print("ERROR: num or key must be given")
		if key is None:	
			key = f"slider_{num}"
		self.settings[key]['dimension_type'] = dimension_type
		self.settings[key]['id'] = 0
		dim = self.reader.dimensions[dimension_type]
		if self.axes[key]:
			self.axes[key].set_visible(False)
			del(self.axes[key])
		self.axes[key] = axes(self[key]['axis'],visible=self['visible'][key])
		self.sliders[key] = Slider(self.axes[key], f"{dim.axis_label} index:", 0, len(dim)-1, valinit = 0, valstep = 1, orientation = self[key]['orientation'])
		if self[key]['orientation'] == 'vertical':
			self.sliders[key].label.set_rotation(90)
		self.sliders[key].on_changed(self.draw_fig)
		self.set_visible(key,val=self['visible'][key])
		if _drawfig:
			self.draw_fig()
	
	def remove_slider(self, num = None, key = None, _draw_fig = True):
		if num is None and key is None:
			print("ERROR: num or key must be given")
		if key is None:	
			key = f"slider_{num}"
		self.settings[key]['dimension_type'] = None
		self.settings[key]['id'] = 0
		self.set_visible(key,False)
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys {self.settings['visible'].keys()}")
			return
		if val not in [True,False]:
			val = not self['visible'][key]
		self.settings['visible'][key] = val
		self.axes[key].set_visible(self['visible'][key])
	
	def add_plot(self, plot):
		self.plots.add(plot)
	
	def draw_fig(self, val = None):
		for key in [x for x in self.sliders.keys()]:
			sli = self.sliders[key]
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings[key]['id'] = sli.val
		for plot in self.plots:
			plot.draw_fig()

