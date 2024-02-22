from numpy import real, imag, array, pi, amax
from matplotlib.pyplot import subplots, ion, show
from copy import deepcopy
from .slider_ax import slider_axes

default_settings = {"suptitle": None,
		"var": 'phi',
		"run": {},
		"normalisation": "highest",
		"fontsizes": {"legend": 10,"ch_box": 8,"axis": 11,"title": 13,"suptitle": 20, "verify": 8},
		"visible": {"suptitle": True, "title": True, "legend": True, 'absolute': True, 'real': True, 'imag': True, 'divider': True},
		"colours": {"real": 'r', "imag": 'b', "absolute": 'k', "divider": 'g'},
}

slider_settings = {"slider_1": {"axis": [0.25, 0.01, 0.5, 0.03]},
		"slider_2": {"axis": [0.25, 0.05, 0.5, 0.03]},
		"slider_3": {"axis": [0.92, 0.2, 0.03, 0.5], "orientation": 'vertical'}, 
		"slider_4": {"axis": [0.96, 0.2, 0.03, 0.5], "orientation": 'vertical'},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "slider_5": False, "slider_6": False, "slider_7": False, "slider_8": False, "slider_9": False},
}

class plot_box_diag(object):
	def __init__(self, reader, settings = {}, sliders = None):
		self.reader = reader
		self.verify = reader.verify
		self.sliders = sliders
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
				
		self.set_normalisation(self['normalisation'])
		
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
		if norm.lower() not in ['none','highest','phi','apar','bpar','epar']:
			print("ERROR: invalid normalisation, allowed: ['none','highest','phi','apar','bpar','epar']")
			return
			
		self.settings['normalisation'] = norm.lower()
		self.draw_fig()
		
	def draw_fig(self, val = None):
		for key, sli in self.sliders.sliders.items():
			dim = self.sliders.settings[key]['dimension_type']
			if dim in self.dims:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]	
		
		self.ax.cla()
		run = self['run']
		ky_id = self.reader.dimensions['ky'].values.index(self.reader('ky',run))
		if ky_id != 0:
			if 'jtwist' in self['run']:
				jtwist = self['run']['jtwist']
			else:
				jtwist = self.reader.single_parameters['jtwist'].values[0]
			kx_id = self.reader.dimensions['kx'].values.index(self.reader('kx',run))
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
				th = array(self.reader('theta',run))
				th_width = max(th)
				th_new = th_width + th
				theta += list(max(th_new)*(i+1) + th_new)
			theta = list(array(theta)-max(theta)/2)
			
			if self['visible']['title']:
				title = "".join([f"{self.reader.dimensions[x].axis_label}: {self.reader(x,run):.2g} | " for x in [y for y in self.reader.inputs.dim_order if y != 'kx']])
				kxs = [f"{self.reader.dimensions['kx'].values[i]:.2g}" for i in kx_list]
				title += f"{self.reader.dimensions['kx'].axis_label}: "
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
			
			if self['normalisation'] == 'highest':
					norms = []
					for kxid in kx_list:
						run['kx'] = self.reader.dimensions['kx'].values[kxid]
						norm = max([amax([abs(i) for i in self.reader(x,run)]) for x in [y for y in ['phi','apar','bpar','epar'] if (self.reader(y,run) is not None)]])
						norms.append(norm)
					norm = max(norms)
					ylabel += " / max($\phi,A_{\parallel},B_{\parallel},E_{\parallel}$)"
			elif self['normalisation'] in ['phi','apar','bpar','epar']:
				norms = []
				for kxid in kx_list:
					run['kx'] = self.reader.dimensions['kx'].values[kxid]
					norm = amax([abs(i) for i in self.reader(self['normalisation'],run)])
					norms.append(norm)
				norm = max(norms)
				ylabel += f" / max({self['normalisation']})"
			else:
				norm = 1
			field_norm = array(field)/norm
			abs_line_style = 'solid'
			if self['visible']['real']:
				abs_line_style = 'dashed'
				self.ax.plot(theta,real(field_norm),color=self['colours']['real'],label="real")
			if self['visible']['imag']:
				abs_line_style = 'dashed'
				self.ax.plot(theta,imag(field_norm),color=self['colours']['imag'],label="imaginary")
			if self['visible']['absolute']:
				self.ax.plot(theta,[abs(x) for x in field_norm],color=self['colours']['absolute'],linestyle=abs_line_style,label="absolute")
			if self['visible']['legend']:
				self.ax.legend(loc=0)
			if self['visible']['divider']:
				vls = [min(theta) + 2*th_width*(i+1) for i in range(len(kx_list)-1)]
				ylim = self.ax.get_ylim()
				self.ax.vlines(vls,ylim[0],ylim[1],color=self['colours']['divider'],linestyle='dotted')
				
			self.ax.set_xlabel("Ballooning Angle",fontsize=self['fontsizes']['axis'])
			self.ax.set_ylabel(ylabel,fontsize=self['fontsizes']['axis'])
		self.fig.canvas.draw_idle()
		return
