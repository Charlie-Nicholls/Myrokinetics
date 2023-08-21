from numpy import transpose, array, amax, amin, isfinite, linspace, where, full, nan
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons, TextBox
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator

default_settings = {"suptitle": None,
		"eqbm_style": "title",
		"contour_type": 0,
		"x_axis_type": "beta_prime",
		"y_axis_type": "shear",
		"slider_1": {"dimension_type": None, "id": 0},
		"ql_slider": {"scale": 100, "max": None},
		"run": {},
		"options": [False,False,True,True,False,False],
		"fontsizes": {"title": 13, "ch_box": 8,"axis": 17,"suptitle": 20},
		"visible": {"slider_1": True, "ql_slider": True, "op_box": True, "suptitle": True, "title": True},
		"cdict": {"red": ((0.0, 1, 1),(1.0, 0.8, 0.8)),
			"green": ((0.0, 1, 1),(1.0, 0.0, 0.0)),
			"blue": ((0.0, 1, 1),(1.0, 0.0, 0.0))},
}

class plot_ql(object):
	def __init__(self, reader, settings = {}):	
		self.reader = reader
		if self.reader['quasilinear'] is None:
			print("Error: No QuasiLinear Data")
			return
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		self._options = ["Show ID","Global Axis Limits","Global Colourbar","Show Equillibrium","Show Ideal","Draw EQBM Contour"]

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
				for skey in self.defaults:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
						self.init_settings[key][skey] = defaults[key][skey]
		
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
			filename = f"QL_{self['slider_1']['id']}"
		self.fig.savefig(filename)
		
	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(9,7))
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
	
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		
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

		self.cmap = LinearSegmentedColormap('WhRd', self.settings['cdict'])
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax)
		
		self.ch_axes = axes([0.75, 0.01, 0.09, 0.09],frame_on = False,visible=self['visible']['op_box'])
		self.options = CheckButtons(self.ch_axes, self._options, self['options'])
		self.options.on_clicked(self.draw_fig)
		self.set_options_fontsize(self['fontsizes']['ch_box'])
		
		self._slider_axes = {'slider_1': axes([0.15, 0.01, 0.5, 0.03],visible=self['visible']['slider_1']),
		}
		self.sliders = {}
		for key in [x for x in self._slider_axes.keys() if self.settings[x]['dimension_type'] is not None]:
			dim = self.reader.dimensions[self.settings[key]['dimension_type']]
			self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index:", 0, len(dim)-1, valinit = self[key]['id'], valstep = 1)
			self.sliders[key].on_changed(self.draw_fig)
		
		self.ql_axes = axes([0.9, 0.13, 0.01, 0.73],visible=self['visible']['ql_slider'])
		self.sliders['ql_slider'] = Slider(self.ql_axes, 'Scale', 0, 100, valinit = self['ql_slider']['scale'], valstep = 1, orientation = 'vertical')
		self.sliders['ql_slider'].on_changed(self.draw_fig)
			
		if self['visible']['slider_1'] == True:	
			self.fig.subplots_adjust(bottom=0.15)
		elif self['visible']['slider_1'] == False:
			self.fig.subplots_adjust(bottom=0.11)
		
		ion()
		self.draw_fig()
		if self['ql_slider']['max']:
			self.set_ql_max(self.settings['ql_slider']['max'])
		show()
	
	def set_slider(self, slider_num, dimension_type, visible = True):
		if dimension_type not in self.dims:
			print(f"ERROR: invalid dimension type, valid: {self.dims}")
			return
		key = f"slider_{slider_num}"
		self.settings[key]['dimension_type'] = dimension_type
		self.settings[key]['id'] = 0
		dim = self.reader.dimensions[dimension_type]
		self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index:", 0, len(dim)-1, valinit = 0, valstep = 1)
		self.set_visible(key,val=visible)
		self.draw_fig()
		
	def _load_x_axis(self, axis_type):
		if axis_type not in ['beta_prime','alpha']:
			print(f"ERROR: axis_type not found, valid types ['beta_prime','alpha']")
			return
			
		self.settings['x_axis_type'] = axis_type
		
		if axis_type in ['alpha']:
			if not self.data['alpha_axis']:
				print("ERROR: Alpha not calculated, use calculate_alpha()")
			else:
				#self.x_axis = self.data['alpha_axis'] needs updating
				self._x_axis_label = r'$\alpha$'
		else:
			self.x_axis = self.reader.dimensions['beta_prime'].values
			self._x_axis_label = r'$\beta^{\prime}$'
			
	def set_x_axis_type(self, axis_type):
		self._load_x_axis(axis_type)
		self.draw_fig()
	
	def _load_y_axis(self, axis_type):
		if axis_type not in ['shear','current']:
			print(f"ERROR: axis_type not found, valid types ['shear','current']")
			return
			
		self.settings['y_axis_type'] = axis_type
		
		if axis_type in ['current']:
			print("Not yet implimented")
		else:
			self.y_axis = self.reader.dimensions['shear'].values
			self._y_axis_label = r'$\hat{s}$'
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def set_ql_max(self, val = None):
		self.settings['ql_slider']['max'] = val
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
			if self['visible']['slider_1'] == True:	
				self.fig.subplots_adjust(bottom=0.15)
			elif self['visible']['slider_1'] == False:
				self.fig.subplots_adjust(bottom=0.11)
				
		elif key == 'ql_slider':
			self.ql_axes.set_visible(self['visible']['ql_slider'])
		elif key == 'op_box':
			self.ch_axes.set_visible(self['visible']['op_box'])
		elif key == 'suptitle':
			self.fig._suptitle.set_visible(self['visible']['suptitle'])
		elif key == 'title':
			self.ax.legend_.set_visible(self['visible']['title'])
			
	def set_options_fontsize(self, fontsize):
		self.settings['fontsizes']['cbox'] = fontsize
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
	
	def set_option(self, option, value):
		if type(option) == str:
			if option not in self._options:
				print(f"ERROR: option not found, valid = {self._options}")
				return
			else:
				option = self._options.index(option)
		if option >= len(self._options):
			print("ERROR: index out of bounds")
			return
		if self.options.get_status()[option] != value:
			self.settings['options'][option] = value
			self.options.set_active(option)
	
	def set_options(self, values):
		if len(values) >= len(self._options):
			print("ERROR: list too long")
		for i, val in enumerate(values):
			if self.options.get_status()[i] != val:
				self.options.set_active(i)
				self.settings['options'][i] = val
	
	def set_cdict(self, cdict):
		self.settings['cdict'] = cdict
		self.cmap = LinearSegmentedColormap('User', self.settings['cdict'])
		self.draw_fig()

	def draw_fig(self, val = None):
		handles = []
		for key in [x for x in self.sliders.keys() if x != 'ql_slider']:
			sli = self.sliders[key]
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				self.settings[key]['id'] = sli.val
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label=f"{self.reader.dimensions[dim].axis_label} = {self.settings['run'][dim]}",visible = False))
		
		if 'psin' in self['run']:
			psiN = self['run']['psin']
		else:
			psiN = self.reader.single_parameters.values[0]
				
		x_axis = list(self.x_axis)
		y_axis = list(self.y_axis)
		
		x_val = abs(self.reader.data['equilibrium'][psiN][self['x_axis_type']])
		y_val = self.reader.data['equilibrium'][psiN][self['y_axis_type']]
		
		self.ax.cla()
		self.ax.set_facecolor('grey')
		self.ax.set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax.set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		
		status = self.options.get_status()
		self.settings['options'] = status
		
		z = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		for x_id, x_value in enumerate(self.x_axis):
			for y_id, y_value in enumerate(self.y_axis):
				run = self['run'].copy()
				run[self['x_axis_type']] = x_value
				run[self['y_axis_type']] = y_value
				run_id = self.reader.get_run_id(run = run, keys = '_quasilinear_keys')
				z[x_id][y_id] = self.reader.data['quasilinear'][run_id]
		
		if self['ql_slider']['max']:
			qlmax = self.sliders['ql_slider'].val * self['ql_slider']['max']/100
		#elif status[2]:
			#qlmax = self.sliders['ql_slider'].val * abs(amax(array(self.data['quasilinear'])[isfinite(self.data['quasilinear'])]))/100
		else:
			try:
				qlmax = self.sliders['ql_slider'].val * amax(abs(array(z)[isfinite(z)]))/100
				if qlmax < 10e-10:
					qlmax = 10e-10
			except:
				qlmax = 1
		
		norm = Normalize(vmin=0,vmax=qlmax)	
		self.settings['ql_slider']['scale'] = self.sliders['ql_slider'].val
		self.cbar.update_normal(ScalarMappable(norm = norm, cmap = self.cmap))
		if self['contour_type'] == 1:
			self.ax.contourf(x_axis,y_axis,z, cmap = self.cmap, norm = norm)
		else:
			self.ax.pcolormesh(x_axis, y_axis, transpose(z), cmap = self.cmap, norm=norm)
			
		if status[5]:
			grid = RegularGridInterpolator((array(x_axis),array(y_axis)),transpose(z))
			eqbm_ql = grid((x_val,y_val))
			self.ax.contourf(x_axis, y_axis, z, [eqbm_ql-eqbm_ql/100,eqbm_ql+eqbm_ql/100], colors = ('b'))
			
		if status[0]:
			dx = x_axis[1] - x_axis[0]
			dy = y_axis[1] - y_axis[0]
			self.ax.set_xticks(array(x_axis)-dx/2, minor=True)
			self.ax.set_yticks(array(y_axis)-dy/2, minor=True)
			
			self.ax.set_xticks(array(x_axis), minor=False)
			self.ax.set_yticks(array(y_axis), minor=False)
			
			self.ax.set_xticklabels([])
			self.ax.set_yticklabels([])
			xlabels = [str(i) for i in range(len(x_axis))]
			ylabels = [str(i) for i in range(len(y_axis))]
			self.ax.set_xticklabels(xlabels,minor=False)
			self.ax.set_yticklabels(ylabels,minor=False)
			
			self.ax.grid(which="minor",color='k')
		
		if status[1]:
			self.ax.set_xlim(amin(self.x_axis),amax(self.x_axis))
			self.ax.set_ylim(amin(self.y_axis),amax(self.y_axis))

		if status[3]:
			self.ax.plot(x_val,y_val,'kx')
			eqbm_pos = f"{x_val:.2f},{y_val:.2f}"
			if status[5]:
				eqbm_pos += f",{eqbm_ql:.2f}"
			if self['eqbm_style'] in ["point",2,"point numless",4]:
				self.ax.annotate("Eqbm",(x_val,y_val),textcoords = "offset points",xytext = (0,7), ha = "center")
			if self['eqbm_style'] in ["split",1,"point",2]:
				self.ax.annotate(eqbm_pos,(x_val,y_val),textcoords = "offset points",xytext = (0,-13), ha = "center")
			if self['eqbm_style'] in ["title",0]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label=f"Equillibrium ({eqbm_pos})",linewidth=0))
			if self['eqbm_style'] in ["split",1,"title numless",3]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label="Equillibrium",linewidth=0))
		
		if status[4]:
			if self.reader.data['ideal_data'] is not None and self.reader.data['ideal_data'][psiN] is not None:
				idata = self.reader.data['ideal_data'][psiN]
				self.ax.contourf(idata[self['x_axis_type']], idata[self['y_axis_type']], idata['stabilities'], [0.01,0.99], colors = ('k'))
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label="Ideal Boundary"))
			else:
				self.ax.text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=self.ax.transAxes,color='k')
		
		self.ax.legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0.5,0.98),loc = "lower center", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax.legend_.set_visible(self['visible']['title'])

		self.fig.canvas.draw_idle()
