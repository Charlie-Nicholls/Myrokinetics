from numpy import real, imag, array, pi
from matplotlib.pyplot import *
from matplotlib.widgets import Slider
from copy import deepcopy

default_settings = {"suptitle": None,
		"var": 'phi',
		"slider_1": {"dimension_type": None, "id": 0},
		"slider_2": {"dimension_type": None, "id": 0},
		"slider_3": {"dimension_type": None, "id": 0},
		"slider_4": {"dimension_type": None, "id": 0},
		"run": {},
		"normalisation": "True",
		"fontsizes": {"legend": 10,"ch_box": 8,"axis": 11,"title": 13,"suptitle": 20, "verify": 8},
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "suptitle": True, "title": True, "legend": True, 'absolute': True, 'real': True, 'imag': True, 'divider': True},
		"colours": {"real": 'r', "imag": 'b', "absolute": 'k', "divider": 'g'},
}

class plot_box_diag(object):
	def __init__(self, reader, settings = {}):
		self.reader = reader
		self.verify = reader.verify
		
		self.settings = {}
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
		
		if self['var'] == 1:
			self.settings['var'] = "phi"
		if self['var'] == 2:
			self.settings['var'] = "apar"
		if self['var'] == 3:
			self.settings['var'] = "bpar"
		if self['var'] == 4:
			self.settings['var'] = "epar"
		if self['var'] not in ['phi','apar','bpar','epar']:
			print(f"ERROR: variable name/value {self['var']} not supported. supported: phi/1, apar/2, bpar/3, epar/4")
			return
			
		self.open_plot()


	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
	
	def save_plot(self, filename = None):
		if filename is None:
			filename = f"{self['var']}_{self['slider_1']['id']}_{self['slider_2']['id']}_{self['slider_3']['id']}_{self['slider_4']['id']}"
		self.fig.savefig(filename)
	
	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(8.8,5.8))	
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
				
		self.set_normalisation(self['normalisation'])
		
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
		self.draw_fig()
	
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
	
	def set_normalisation(self, norm):
		self.settings['normalisation'] = norm
		self.draw_fig()
		
	def draw_fig(self, val = None):
		self.ax.cla()
		for key, sli in self.sliders.items():
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				self.settings[key]['id'] = sli.val		
		
		data = self.reader('data',self['run'])
		ky_id = self.reader.dimensions['ky'].values.index(data['ky'])
		if ky_id != 0:
			if 'jtwist' in self['run']:
				jtwist = self['run']['jtwist']
			else:
				jtwist = self.reader.single_parameters['jtwist'].values[0]
			kx_id = self.reader.dimensions['kx'].values.index(data['kx'])
			kx_list = set()
			offset = jtwist*ky_id
			new_id = kx_id%offset
			while new_id < len(self.reader.dimensions['kx']):
				kx_list.add(new_id)
				new_id += offset
			kx_list = list(kx_list)
			kx_list.sort()
			kx_list.reverse()
			field = []
			theta = []
			for i, kxid in enumerate(kx_list):
				run = deepcopy(self['run'])
				run['kx'] = self.reader.dimensions['kx'].values[kxid]
				field += self.reader(self['var'],run)
				theta += list(2*pi*i+array(self.reader('theta',run)))
			theta = list(array(theta)-(len(kx_list)-1)*pi)
			
			title = "".join([f"{self.reader.dimensions[x].axis_label}: {data[x]:.2g} | " for x in [y for y in self.reader.inputs.dim_order if y != 'kx']])
			kxs = [f"{self.reader.dimensions['kx'].values[i]:.2g}" for i in kx_list]
			title += "kx: "
			title += "".join([f"{kx}, " for kx in kxs])[:-2]
			self.ax.set_title(title,fontsize=self['fontsizes']['title'])
			
			if self['var'] == 'phi':
				ylabel = "$\phi$"
			elif self['var'] == 'apar':
				ylabel = "$A_{\parallel}$"
			elif self['var'] == 'bpar':
				ylabel = "$B_{\parallel}$"
			elif self['var'] == 'epar':
				ylabel = "$E_{\parallel}$"
			
			if self['normalisation'] == True:
				norm = max(field)
				ylabel += f" / max({ylabel})"
			else:
				norm = 1
			field_norm = array(field)/norm
			if self['visible']['real']:
				self.ax.plot(theta,real(field_norm),color=self['colours']['real'],label="real")
			if self['visible']['imag']:
				self.ax.plot(theta,imag(field_norm),color=self['colours']['imag'],label="imaginary")
			if self['visible']['absolute']:
				self.ax.plot(theta,[abs(x) for x in field_norm],color=self['colours']['absolute'],linestyle='--',label="absolute")
			if self['visible']['divider']:
				vls = [min(theta) + 2*pi*(i+1) for i in range(len(kx_list)-1)]
				ylim = self.ax.get_ylim()
				self.ax.vlines(vls,ylim[0],ylim[1],color=self['colours']['divider'],linestyle='--')
				
			self.ax.legend(loc=0)
				
			self.ax.set_xlabel("Ballooning Angle",fontsize=self['fontsizes']['axis'])
			self.ax.set_ylabel(ylabel,fontsize=self['fontsizes']['axis'])
		self.fig.canvas.draw_idle()
		return
