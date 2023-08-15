from numpy import transpose, array, amax, amin, isfinite, linspace, where, full, nan
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons, TextBox
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator

default_settings = {"suptitle": None,
		"psi_id": 0,
		"ql_scale": 100,
		"ql_max": None,
		"eqbm_style": "title",
		"contour_type": 0,
		"x_axis_type": 'beta_prime',
		"y_axis_type": 'shear',
		"options": [False,False,True,True,False,False],
		"fontsizes": {"title": 13, "ch_box": 8,"axis": 17,"suptitle": 20},
		"visible": {"psi_sli": True, "ql_sli": True, "op_box": True, "suptitle": True, "title": True},
		"cdict": {"red": ((0.0, 1, 1),(1.0, 0.8, 0.8)),
			"green": ((0.0, 1, 1),(1.0, 0.0, 0.0)),
			"blue": ((0.0, 1, 1),(1.0, 0.0, 0.0))},
}

class plot_ql(object):
	def __init__(self, reader, settings = {}):
			
		self.reader = reader
		self.inputs = reader.inputs
		if 'psin' in self.reader.dimensions:
			self.psiNs = self.reader.dimensions['psin'].values
		else:
			self.psiNs = [self.reader.single_parameters['psin'].values
		
		if self.reader['quasilinear'] is None:
			print("Error: No QuasiLinear Data")
			return
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		self._options = ["Show ID","Global Axis Limits","Global Colourbar","Show Equillibrium","Show Ideal","Draw EQBM Contour"]

		self.settings = {}
		self.init_settings = {}
		for key in settings:
			if key not in default_settings:
				print(f"ERROR: {key} not found")
			else:
				self.settings[key] = settings[key]
				self.init_settings[key] = settings[key]
		for key in default_settings:
			if key not in self.settings:
				self.settings[key] = default_settings[key]
				self.init_settings[key] = default_settings[key]
			elif type(self.settings[key]) == dict and key != 'cdict':
				for skey in self.default_settings:
					if skey not in self.settings[key]:
						self.settings[key][skey] = default_settings[key][skey]
						self.init_settings[key][skey] = default_settings[key][skey]
		
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
			filename = f"QL_{self['psi_id']}"
		self.fig.savefig(filename)
		
	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(9,7))
		self.fig.subplots_adjust(bottom=0.15)
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
	
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		
		self.dims = [x for x in self.reader.inputs.dim_order if x not in [self['x_axis_type'],self['y_axis_type']]]
		
		self.cmap = LinearSegmentedColormap('WhRd', self.settings['cdict'])
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax)
		
		self.ch_axes = axes([0.75, 0.01, 0.09, 0.09],frame_on = False,visible=self['visible']['op_box'])
		self.options = CheckButtons(self.ch_axes, self._options, self['options'])
		self.options.on_clicked(self.draw_fig)
		self.set_options_fontsize(self['fontsizes']['ch_box'])
		
		self.psi_axes = axes([0.15, 0.01, 0.5, 0.03],visible=self['visible']['psi_sli'])
		self.slider = Slider(self.psi_axes, '\u03A8\u2099 index:', 0, len(self.psiNs)-1, valinit = self['psi_id'], valstep = 1)
		if self['visible']['psi_sli']:
			self.fig.subplots_adjust(bottom=0.15)
		self.slider.on_changed(self.draw_fig)
		
		self.ql_axes = axes([0.9, 0.13, 0.01, 0.73],visible=self['visible']['ql_sli'])
		self.ql_slider = Slider(self.ql_axes, 'Scale', 0, 100, valinit = self['ql_scale'], valstep = 1, orientation = 'vertical')
		self.ql_slider.on_changed(self.draw_fig)
		
		if len(self.psiNs) == 1:
			self.set_visible(key = 'psi_sli', val = False)
		
		ion()
		self.draw_fig()
		if self['ql_max']:
			self.set_ql_max(self.init_settings['ql_max'])
		show()
	
	def _load_x_axis(self, axis_type):
		if axis_type not in ['beta_prime','alpha']:
			print(f"ERROR: axis_type not found, valid types ['beta_prime','alpha']")
			return
			
		self.settings['x_axis_type'] = axis_type
		
		if axis_type in ['alpha']:
			if not self.data['alpha_axis']:
				print("ERROR: Alpha not calculated, use calculate_alpha()")
			else:
				self.x_axis = self.data['alpha_axis']
				self.x_axis_ideal = self.data['alpha_axis_ideal']
				self.x_values = self.data['alpha_values']
				self._x_axis_label = "\u03B1"
		else:
			self.x_axis = self.reader.dimensions['beta_prime'].values
			self._x_axis_label = "-\u03B2'"
			
	def set_x_axis_type(self, axis_type):
		self._load_x_axis(axis_type)
		self.draw_fig()
	
	def _load_y_axis(self, axis_type):
		if axis_type not in ['shear','current']:
			print(f"ERRORL axis_type not found, valid types ['shear','current']")
			return
			
		self.settings['y_axis_type'] = axis_type
		
		if axis_type in ['current']:
			print("Not yet implimented")
		else:
			self.y_axis = self.reader.dimensions['shear'].values
			self._y_axis_label = "Shear"
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def set_ql_max(self, val = None):
		self.settings['ql_max'] = val
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys {self.settings.keys()}")
			return
		if val not in [True,False]:
			val = not self['visible'][key]
		
		self.settings['visible'][key] = val
		
		if key == 'psi_sli' and val is True:
			self.fig.subplots_adjust(bottom=0.15)
			self.psi_axes.set_visible(self['visible']['psi_sli'])
		elif key == 'psi_sli' and val is False:
			self.fig.subplots_adjust(bottom=0.11)
			self.psi_axes.set_visible(self['visible']['psi_sli'])
		elif key == 'ql_sli':
			self.ql_axes.set_visible(self['visible']['ql_sli'])
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
		idx = self.slider.val
		psiN = self.psiNs[idx]
		self.settings['psi_id'] = idx
		
		x = list(self.x_axis)
		y = list(self.y_axis)
		
		x_val = abs(self.reader.data['equilibrium'][psiN]['beta_prime'])
		y_val = self.reader.data['equilibrium'][psiN]['shear']
		
		psi_line = Line2D([0,1],[0.5,0.5],color='k',label=f"\u03A8\u2099 = {psiN}",visible = False)
		ideal_line =  None
		eqbm_line = None
		
		self.ax.cla()
		self.ax.set_facecolor('grey')
		self.ax.set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax.set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		
		status = self.options.get_status()
		self.settings['options'] = status
		
		z = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		for x_id, x_value in enumerate(self.x_axis):
			for y_id, y_value in enumerate(self.y_axis):
				run = {'psin': psiN, self['x_axis_type']: x_value, self['y_axis_type']: y_value}
				run_id = self.reader.get_run_id(run = run, keys = '_quasilinear_keys')
				z[x_id][y_id] = self.reader.data['quasilinear'][run_id]
		
		if self['ql_max']:
			qlmax = self.ql_slider.val * self['ql_max']/100
		#elif status[2]:
			#qlmax = self.ql_slider.val * abs(amax(array(self.data['quasilinear'])[isfinite(self.data['quasilinear'])]))/100
		else:
			try:
				qlmax = self.ql_slider.val * amax(abs(array(z)[isfinite(z)]))/100
				if qlmax < 10e-10:
					qlmax = 10e-10
			except:
				qlmax = 1
		
		norm = Normalize(vmin=0,vmax=qlmax)	
		self.settings['ql_scale'] = self.ql_slider.val
		self.cbar.update_normal(ScalarMappable(norm = norm, cmap = self.cmap))
		if self['contour_type'] == 1:
			self.ax.contourf(x,y,z, cmap = self.cmap, norm = norm)
		else:
			self.ax.pcolormesh(x, y, transpose(z), cmap = self.cmap, norm=norm)
			
		if status[5]:
			grid = RegularGridInterpolator((array(x),array(y)),transpose(z))
			eqbm_ql = grid((x_val,y_val))
			self.ax.contourf(x, y, z, [eqbm_ql-eqbm_ql/100,eqbm_ql+eqbm_ql/100], colors = ('b'))
			
		if status[0]:
			dx = x[1] - x[0]
			dy = y[1] - y[0]
			self.ax.set_xticks(array(x)-dx/2, minor=True)
			self.ax.set_yticks(array(y)-dy/2, minor=True)
			
			self.ax.set_xticks(array(x), minor=False)
			self.ax.set_yticks(array(y), minor=False)
			
			self.ax.set_xticklabels([])
			self.ax.set_yticklabels([])
			xlabels = [str(i) for i in range(len(x))]
			ylabels = [str(i) for i in range(len(y))]
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
				eqbm_line = Line2D([0.5],[0.5],marker='x',color='k',label=f"Equillibrium ({eqbm_pos})",linewidth=0)
			if self['eqbm_style'] in ["split",1,"title numless",3]:
				eqbm_line = Line2D([0.5],[0.5],marker='x',color='k',label="Equillibrium",linewidth=0)
		
		if status[4]:
			if self.reader.data['ideal_data'] is not None and self.reader.data['ideal_data'][psiN] is not None:
				self.ax.contourf(self.reader.data['ideal_data'][psiN]['beta_prime'], self.reader.data['ideal_data'][psiN]['shear'], self.reader.data['ideal_data'][psiN]['stabilities'], [0.01,0.99], colors = ('k'))
				ideal_line = Line2D([0,1],[0.5,0.5],color='k',label="Ideal Boundary")
			else:
				self.ax.text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=self.ax.transAxes,color='k')
		
		handles = [line for line in [psi_line,ideal_line,eqbm_line] if line is not None]
		self.ax.legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0.5,0.98),loc = "lower center", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax.legend_.set_visible(self['visible']['title'])

		self.fig.canvas.draw_idle()
