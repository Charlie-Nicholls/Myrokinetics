from numpy import transpose, array, amax, amin, isfinite, linspace, where, full, nan, diff
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap
from copy import deepcopy

default_settings = {"suptitle": None,
		"eqbm_style": "title",
		"contour_type": 0,
		"x_axis_type": "kx",
		"y_axis_type": "ky",
		"z_axis_type": "growth_rate",
		"xscale": "linear",
		"yscale": "linear",
		"slider_1": {"dimension_type": None, "id": 0},
		"slider_2": {"dimension_type": None, "id": 0},
		"gr_slider": {"scale": 100, "max": None},
		"mf_slider": {"scale": 100, "max": None},
		"run": {},
		"fontsizes": {"title": 13, "ch_box": 8, "vr_box": 8, "axis": 17, "suptitle": 20, "legend": 13},
		"visible": {"slider_1": True, "slider_2": True, "gr_slider": True, "mf_slider": True, "op_box": True, "vr_box": True, "suptitle": True, "title": True, "legend": True},
		"cdict": {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
			'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
			'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))},
}

class plot_kxky(object):
	def __init__(self, reader = None, settings = {}):
		self.reader = reader
		self.verify = reader.verify
		self.data = reader.data['gyro']
		if self.data is None:
			print("Error: No Gyrokinetic Data")
			return
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
			elif type(self.settings[key]) == dict and key != 'cdict':
				for skey in defaults[key]:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
						self.init_settings[key][skey] = defaults[key][skey]
		
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
		
		self.cmap = LinearSegmentedColormap('GnRd', self['cdict'])
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar_mf = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[1])
		self.cbar_gr = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[0])
		
		self._slider_axes = {'slider_1': axes([0.15, 0.01, 0.5, 0.03],visible=self['visible']['slider_1']),
			'slider_2': axes([0.15, 0.05, 0.5, 0.03],visible=self['visible']['slider_2']),
		}
		self.sliders = {}
		for key in [x for x in self._slider_axes.keys() if self.settings[x]['dimension_type'] is not None]:
			dim = self.reader.dimensions[self.settings[key]['dimension_type']]
			self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index:", 0, len(dim)-1, valinit = self[key]['id'], valstep = 1)
			self.sliders[key].on_changed(self.draw_fig)
		
		self.gr_axes = axes([0.93, 0.15, 0.01, 0.73])
		self.sliders['gr_slider'] = Slider(self.gr_axes, 'GR', 0, 100, valinit = self['gr_slider']['scale'], valstep = 1, orientation = 'vertical')
		self.sliders['gr_slider'].on_changed(self.draw_fig)
		
		self.mf_axes = axes([0.97, 0.15, 0.01, 0.73])
		self.sliders['mf_slider'] = Slider(self.mf_axes, 'MF', 0, 100, valinit = self['mf_slider']['scale'], valstep = 1, orientation = 'vertical')
		self.sliders['mf_slider'].on_changed(self.draw_fig)
		
		if not self['visible']['slider_1'] and not self['visible']['slider_2'] and not self['visible']['vr_box']:
			self.fig.subplots_adjust(bottom=0.11)
		elif not self['visible']['slider_2'] and not self['visible']['vr_box']: 
			self.fig.subplots_adjust(bottom=0.13)
		else:
			self.fig.subplots_adjust(bottom=0.15)

		ion()
		self.draw_fig()
		if self['gr_slider']['max']:
			self.set_gr_max(self.init_settings['gr_max'])
		if self['mf_slider']['max']:
			self.set_mf_max(self.init_settings['mf_max'])
		
		show()
	
	def _load_x_axis(self, axis_type):
		if axis_type not in ['ky','theta0','kx']:
			print(f"ERROR: axis_type not found, valid types ['ky','theta0','kx']")
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
			print(f"ERROR: axis_type not found, valid types ['ky','theta0','kx']")
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
			print(f"ERROR: axis_type not found, valid types ['growth_rate','growth_rate_norm','ql_norm']")
			return
			
		self.settings['z_axis_type'] = axis_type
		
		if axis_type in ['growth_rate']:
			self._z_axis_label = r'$\gamma$'
		elif axis_type in ['growth_rate_norm']:
			self._z_axis_label = r'$\gamma/k_{y}^{2}$'
		elif axis_type in ['ql_norm']:
			self._z_axis_label = "QL"
			
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
		self.settings['gr_max'] = val
		self.draw_fig()
		
	def set_mf_max(self, val):
		self.settings['mf_max'] = val
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
		elif key == 'op_box':
			self.ch_axes.set_visible(self['visible']['op_box'])
		elif key == 'vr_box':
			self.vr_axes.set_visible(self['visible']['vr_box'])
		elif key == 'suptitle':
			self.fig._suptitle.set_visible(self['visible']['suptitle'])
		elif key == 'title':
			self.ax.legend_.set_visible(self['visible']['title'])
			
		if not self['visible']['slider_1'] and not self['visible']['slider_2'] and not self['visible']['vr_box']:
			self.fig.subplots_adjust(bottom=0.11)
		elif not self['visible']['slider_2'] and not self['visible']['vr_box']: 
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
		for key in [x for x in self.sliders.keys() if x not in ['gr_slider','mf_slider']]:
			sli = self.sliders[key]
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				self.settings[key]['id'] = sli.val
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label=f"{self.reader.dimensions[dim].axis_label} = {self.settings['run'][dim]}",visible = False))
		
		if 'psin' in self['run']:
			psiN = self['run']['psin']
		else:
			psiN = self.reader.single_parameters['psin'].values[0]
		
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
		run_ids = []
		for x_id, x_value in enumerate(self.x_axis):
			for y_id, y_value in enumerate(self.y_axis):
				run = self['run'].copy()
				run[self['x_axis_type']] = x_value
				run[self['y_axis_type']] = y_value
				run_id = self.reader.get_run_id(run = run)
				if run_id is not None:
					run_ids.append(run_id)
					z_mf[x_id][y_id] = self.data[run_id]['mode_frequency']
					z_gr[x_id][y_id] = self.data[run_id][self['z_axis_type']]
				else:
					z_gr[x_id][y_id] = nan
					z_mf[x_id][y_id] = nan
		z_gr = transpose(z_gr)
		z_mf = transpose(z_mf)
		
		self.ax[0].set_title(self._z_axis_label,fontsize=self['fontsizes']['title'])
		
		if self['gr_slider']['max']:
			grmax = self.sliders['gr_slider'].val * self['gr_slider']['max']/100
		else:
			grmax = self.sliders['gr_slider'].val * amax(abs(array(z_gr)[isfinite(z_gr)]))/100

		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		self.settings['gr_slider']['scale'] = self.sliders['gr_slider'].val
		self.cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = self.cmap))
		if self['contour_type'] == 1:
			self.ax[0].contourf(x_axis, y_axis, z_gr, cmap = self.cmap, norm=norm_gr)
		else:
			self.ax[0].pcolormesh(x_axis, y_axis, z_gr, cmap = self.cmap, norm=norm_gr)
		
		if self['mf_slider']['max']:
			mfmax = self.sliders['mf_slider'].val * self['mf_slider']['max']/100
		else:
			mfmax = self.sliders['mf_slider'].val * amax(abs(array(z_mf)[isfinite(z_mf)]))/100

		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		self.settings['mf_slider']['scale'] = self.sliders['mf_slider'].val
		self.cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		self.ax[1].pcolormesh(x_axis, y_axis, z_mf, norm = norm_mf)
		
		self.ax[0].legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0,1.02),loc = "lower left", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax[0].legend_.set_visible(self['visible']['legend'])
		self.ax[0].set_xscale(self['xscale'])
		self.ax[0].set_yscale(self['yscale'])
		self.ax[1].set_xscale(self['xscale'])
		self.ax[1].set_yscale(self['yscale'])
		
		self.fig.canvas.draw_idle()
