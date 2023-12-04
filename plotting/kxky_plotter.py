from numpy import transpose, array, amax, isfinite, full, nan
from matplotlib.pyplot import subplots, ion, show, Normalize, colorbar, axes, Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider
from matplotlib.colors import LinearSegmentedColormap
from copy import deepcopy
from .slider_ax import slider_axes

default_settings = {"suptitle": None,
		"eqbm_style": "title",
		"contour_type": 0,
		"x_axis_type": "kx",
		"y_axis_type": "ky",
		"z_axis_type": "growth_rate",
		"xscale": "linear",
		"yscale": "linear",
		"gr_slider": {"scale": 100, "max": None},
		"mf_slider": {"scale": 100, "max": None},
		"run": {},
		"fontsizes": {"title": 13, "ch_box": 8, "axis": 17, "suptitle": 20, "legend": 13},
		"visible": {"gr_slider": True, "mf_slider": True, "op_box": True, "suptitle": True, "title": True, "legend": True},
		"cdict": {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
			'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
			'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))},
}

slider_settings = {"slider_1": {"axis": [0.15, 0.01, 0.5, 0.03]},
		"slider_2": {"axis": [0.15, 0.05, 0.5, 0.03]},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": False, "slider_4": False, "slider_5": False, "slider_6": False, "slider_7": False, "slider_8": False, "slider_9": False},
}

class plot_kxky(object):
	def __init__(self, reader = None, settings = {}, sliders = None):
		self.reader = reader
		self.verify = reader.verify
		self.data = reader.data['gyro']
		self.sliders = sliders
		if self.reader.data['gyro'] is None:
			print("Error: No Gyrokinetic Data")
			return
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
			elif type(self.settings[key]) == dict and key != 'cdict':
				for skey in defaults[key]:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		
		if self['eqbm_style'] not in self._valid_eqbm_styles:
			print("ERROR: eqbm_style not found, valid styles = {self._valid_eqbm_styles}")
			self.settings['eqbm_style'] = "title"
		self.open_plot()
	
	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
	
	def save_plot(self, filename = None):
		if filename is None:
			filename = f"kxky_{self['slider_1']['id']}_{self['slider_2']['id']}"
		self.fig.savefig(filename)
		
	def open_plot(self):
		self.fig, self.ax = subplots(1,2,figsize=(14.6,7))
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		self._load_z_axis(self['z_axis_type'])
		
		self.dims = [x for x in self.reader.inputs.dim_order if x not in [self['x_axis_type'],self['y_axis_type'],'ky','theta0']]
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
		
		self.cmap = LinearSegmentedColormap('GnRd', self['cdict'])
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar_mf = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[1])
		self.cbar_gr = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[0])

		self._sliders = {}
		self.gr_axes = axes([0.93, 0.15, 0.01, 0.73])
		self._sliders['gr_slider'] = Slider(self.gr_axes, 'GR', 0, 100, valinit = self['gr_slider']['scale'], valstep = 1, orientation = 'vertical')
		self._sliders['gr_slider'].on_changed(self.draw_fig)
		
		self.mf_axes = axes([0.97, 0.15, 0.01, 0.73])
		self._sliders['mf_slider'] = Slider(self.mf_axes, 'MF', 0, 100, valinit = self['mf_slider']['scale'], valstep = 1, orientation = 'vertical')
		self._sliders['mf_slider'].on_changed(self.draw_fig)

		self.draw_fig()
		ion()
		show()
	
	def _load_x_axis(self, axis_type):
		if axis_type not in ['ky','theta0','kx']:
			print("ERROR: axis_type not found, valid types ['ky','theta0','kx']")
			return
			
		self.settings['x_axis_type'] = axis_type
		
		if axis_type in ['ky']:
			self.x_axis = self.reader.dimensions['ky'].values
			self._x_axis_label = self.reader.dimensions['ky'].axis_label
		elif axis_type in ['theta0']:
			self.x_axis = self.reader.dimensions['theta0'].values
			self._x_axis_label = self.reader.dimensions['theta0'].axis_label
		elif axis_type in ['kx']:
			self.x_axis = self.reader.dimensions['kx'].values
			self._x_axis_label = r'$k_{x}\rho_{0}$'
		
			
	def set_x_axis_type(self, axis_type):
		self._load_x_axis(axis_type)
		self.draw_fig()
	
	def _load_y_axis(self, axis_type):
		if axis_type not in ['ky','theta0','kx']:
			print("ERROR: axis_type not found, valid types ['ky','theta0','kx']")
			return
			
		self.settings['y_axis_type'] = axis_type
		
		if axis_type in ['ky']:
			self.y_axis = self.reader.dimensions['ky'].values
			self._y_axis_label = self.reader.dimensions['ky'].axis_label
		elif axis_type in ['theta0']:
			self.y_axis = self.reader.dimensions['theta0'].values
			self._y_axis_label = self.reader.dimensions['theta0'].axis_label
		elif axis_type in ['kx']:
			self.y_axis = self.reader.dimensions['theta0'].values
			self._y_axis_label = r'$\rho_{0}$'
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def _load_z_axis(self, axis_type):
		if axis_type not in ['growth_rate','growth_rate_norm','ql_norm']:
			print("ERROR: axis_type not found, valid types ['growth_rate','growth_rate_norm','ql_norm']")
			return
			
		self.settings['z_axis_type'] = axis_type
		
		if axis_type in ['growth_rate']:
			self._z_axis_label = r'$\gamma$'
		elif axis_type in ['growth_rate_norm']:
			self._z_axis_label = r'$\gamma/k_{y}^{2}$'
		elif axis_type in ['ql_norm']:
			self._z_axis_label = "QL metric"
			
	def set_z_axis_type(self, axis_type):
		self._load_z_axis(axis_type)
		self.draw_fig()
	
	def set_xscale(self, scale):
		self.ax[0].set_xscale(scale)
		self.ax[1].set_xscale(scale)
		self.settings['xscale'] = scale
		
	def set_yscale(self, scale):
		self.ax[0].set_yscale(scale)
		self.ax[1].set_yscale(scale)
		self.settings['yscale'] = scale
	
	def set_gr_max(self, val):
		self.settings['gr_slider']['max'] = val
		self.draw_fig()
		
	def set_mf_max(self, val):
		self.settings['mf_slider']['max'] = val
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys: {list(self['visible'].keys())} {list(self.sliders['visible'].keys())}")
			return
		if 'slider_' in key:
			self.sliders.set_visible(key=key,val=val)
		else:
			if val not in [True,False]:
				val = not self['visible'][key]
			self.settings['visible'][key] = val
			if key == 'op_box':
				self.ch_axes.set_visible(self['visible']['op_box'])
			elif key == 'suptitle':
				self.fig._suptitle.set_visible(self['visible']['suptitle'])
			elif key == 'title':
				self.ax.legend_.set_visible(self['visible']['title'])
	
	def set_axis_fontsize(self, fontsize):
		self.settings['fontsizes']['axis'] = fontsize
		self.ax.xaxis.label.set_fontsize(fontsize)
		self.ax.yaxis.label.set_fontsize(fontsize)
		self.fig.canvas.draw_idle()
		
	def set_title_fontsize(self, fontsize):
		self.settings['fontsizes']['title'] = fontsize
		self.draw_fig()
	
	def set_legend_fontsize(self, fontsize):
		self.settings['fontsizes']['legend'] = fontsize
		self.draw_fig()
		
	def set_suptitle_fontsize(self, fontsize):
		self.settings['fontsizes']['suptitle'] = fontsize
		self.fig._suptitle.set_fontsize(fontsize)
		
	def set_suptitle(self, title):
		self.settings['suptitle'] = title
		self.fig.suptitle(title)
	
	def set_eqbm_style(self, eqbm_style):
		if eqbm_style in self._valid_eqbm_styles:
			self.settings['eqbm_style'] = eqbm_style
			self.draw_fig()
		else:
			print(f"ERROR: eqbm_style not found, valid styles = {self._valid_eqbm_styles}")
	
	def set_contour_type(self, contour_type):
		if contour_type in [0,1]:
			self.settings['contour_type'] = contour_type
			self.draw_fig()
		else:
			print(f"ERROR: eqbm_style not found, valid styles = {0,1}")
	
	def set_cdict(self, cdict):
		self.settings['cdict'] = cdict
		self.cmap = LinearSegmentedColormap('User', self.settings['cdict'])
		self.draw_fig()
	
	def draw_fig(self, val = None):
		handles = []
		for key, sli in self.sliders.sliders.items():
			dim = self.sliders.settings[key]['dimension_type']
			if dim in self.dims:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label=f"{self.reader.dimensions[dim].axis_label} = {self.settings['run'][dim]}",visible = False))
		
		x_axis = list(self.x_axis)
		y_axis = list(self.y_axis)
			
		self.ax[0].cla()
		self.ax[0].set_facecolor('grey')
		self.ax[0].set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[0].set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].cla()
		self.ax[1].set_facecolor('grey')
		self.ax[1].set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].set_title("Mode Frequency",fontsize=self['fontsizes']['title'])
		
		z_gr = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		z_mf = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		for x_id, x_value in enumerate(self.x_axis):
			for y_id, y_value in enumerate(self.y_axis):
				run = self['run'].copy()
				run[self['x_axis_type']] = x_value
				run[self['y_axis_type']] = y_value
				z_mf[x_id][y_id] = self.reader('mode_frequency',run)
				z_gr[x_id][y_id] = self.reader(self['z_axis_type'],run)
		self.z_axis_gr = transpose(z_gr)
		self.z_axis_mf = transpose(z_mf)
		
		self.ax[0].set_title(self._z_axis_label,fontsize=self['fontsizes']['title'])
		
		if self['gr_slider']['max']:
			grmax = self._sliders['gr_slider'].val * self['gr_slider']['max']/100
		else:
			grmax = self._sliders['gr_slider'].val * amax(abs(array(self.z_axis_gr)[isfinite(self.z_axis_gr)]))/100

		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		self.settings['gr_slider']['scale'] = self._sliders['gr_slider'].val
		self.cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = self.cmap))
		if self['contour_type'] == 1:
			self.ax[0].contourf(x_axis, y_axis, self.z_axis_gr, cmap = self.cmap, norm=norm_gr)
		else:
			self.ax[0].pcolormesh(x_axis, y_axis, self.z_axis_gr, cmap = self.cmap, norm=norm_gr)
		
		if self['mf_slider']['max']:
			mfmax = self._sliders['mf_slider'].val * self['mf_slider']['max']/100
		else:
			mfmax = self._sliders['mf_slider'].val * amax(abs(array(self.z_axis_mf)[isfinite(self.z_axis_mf)]))/100

		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		self.settings['mf_slider']['scale'] = self._sliders['mf_slider'].val
		self.cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		self.ax[1].pcolormesh(x_axis, y_axis, self.z_axis_mf, norm = norm_mf)
		
		self.ax[0].legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0,1.02),loc = "lower left", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax[0].legend_.set_visible(self['visible']['legend'])
		self.ax[0].set_xscale(self['xscale'])
		self.ax[0].set_yscale(self['yscale'])
		self.ax[1].set_xscale(self['xscale'])
		self.ax[1].set_yscale(self['yscale'])
		
		self.fig.canvas.draw_idle()
