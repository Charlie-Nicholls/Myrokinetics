from numpy import transpose, array, amax, amin, isfinite, linspace, where, full, nan
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap

default_settings = {"suptitle": None,
		"eqbm_style": "title",
		"aky": False,
		"contour_type": 0,
		"x_axis_type": "beta_prime",
		"y_axis_type": "shear",
		"slider_1": {"dimension_type": None, "id": 0},
		"slider_2": {"dimension_type": None, "id": 0},
		"gr_slider": {"scale": 100, "max": None},
		"mf_slider": {"scale": 100, "max": None},
		"run": {},
		"options": [False,True,True,True,False],
		"vr_options": [False,False,False,False,False],
		"fontsizes": {"title": 13, "ch_box": 8, "vr_box": 8, "axis": 17, "suptitle": 20, "legend": 13},
		"visible": {"slider_1": True, "slider_2": True, "gr_slider": True, "mf_slider": True, "op_box": True, "vr_box": True, "suptitle": True, "title": True, "legend": True},
		"cdict": {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
			'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
			'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))},
}

class plot_scan(object):
	def __init__(self, reader = None, settings = {}):
		self.reader = reader
		self.verify = reader.verify
		if self.reader.gyro_data is None:
			print("Error: No Gyrokinetic Data")
			return
		
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
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		self._options = ["Show Parities","Normalised Growth Rate","Global Colorbar","Show Equillibrium","Show Ideal"]
		self._vr_options = ["Show ID","Show Convergence","Show Bad nstep","Show bad fields","Show Omega -> nan"]
		
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
		if filename is None and self['aky']:
			filename = f"Scan_{self['psi_id']}_{self['ky_id']}"
		elif filename is None:
			filename = f"Scan_{self['psi_id']}"
		self.fig.savefig(filename)
		
	def open_plot(self):
		self.fig, self.ax = subplots(1,2,figsize=(14.6,7))
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		
		if self['aky']:
			self.dims = [x for x in self.reader.inputs.dim_order if x not in [self['x_axis_type'],self['y_axis_type']]]
		else:
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
		
		if 'ky' not in self.reader.dimensions and 'ky' not in self.reader.single_parameters:
			self.settings['options'][1] = False
		self.ch_axes = axes([0.85, 0.01, 0.09, 0.1],frame_on = False,visible=self['visible']['op_box'])
		self.options = CheckButtons(self.ch_axes, self._options, self['options'])
		self.options.on_clicked(self.draw_fig)
		
		if self.verify is None:
			self.settings['visible']['vr_box'] = False
		self.vr_axes = axes([0.72, 0.01, 0.09, 0.1],frame_on = False,visible=self['visible']['vr_box'])
		self.vroptions = CheckButtons(self.vr_axes, self._vr_options, self['vr_options'])
		self.vroptions.on_clicked(self.draw_fig)
		
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
	
	def set_gr_max(self, val):
		self.settings['gr_max'] = val
		self.draw_fig()
		
	def set_mf_max(self, val):
		self.settings['mf_max'] = val
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys {self.settings.keys()}")
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
			
	def set_options_fontsize(self, fontsize):
		self.settings['fontsizes']['ch_box'] = fontsize
		self.settings['fontsizes']['vr_box'] = fontsize
		for text in self.options.labels:
			text.set_fontsize(fontsize)
		for text in self.vr_options.labels:
			text.set_fontsize(fontsize)
	
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
	
	def set_vr_options(self, values):
		if len(values) >= len(self._vr_options):
			print("ERROR: list too long")
		for i, val in enumerate(values):
			if self.vr_options.get_status()[i] != val:
				self.vr_options.set_active(i)
				self.settings['vr_options'][i] = val
	
	def set_vr_option(self, option, value):
		if type(option) == str:
			if option not in self._vr_options:
				print(f"ERROR: option not found, valid = {self._vr_options}")
				return
			else:
				option = self._vr_options.index(option)
		if option >= len(self._options):
			print("ERROR: index out of bounds")
			return
		if self.vr_options.get_status()[option] != value:
			self.settings['vr_options'][option] = value
			self.vr_options.set_active(option)
	
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
			psiN = self.reader.single_parameters.values[0]
		
		x_axis = list(self.x_axis)
		y_axis = list(self.y_axis)
		
		x_val = abs(self.reader.data['equilibrium'][psiN][self['x_axis_type']])
		y_val = self.reader.data['equilibrium'][psiN][self['y_axis_type']]
			
		self.ax[0].cla()
		self.ax[0].set_facecolor('grey')
		self.ax[0].set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[0].set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].cla()
		self.ax[1].set_facecolor('grey')
		self.ax[1].set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].set_title("Mode Frequency",fontsize=self['fontsizes']['title'])
		self.ax[0].set_ylim(amin(self.y_axis),amax(self.y_axis))
		self.ax[0].set_xlim(amin(self.x_axis),amax(self.x_axis))
		self.ax[1].set_ylim(amin(self.y_axis),amax(self.y_axis))
		self.ax[1].set_xlim(amin(self.x_axis),amax(self.x_axis))
		
		status = self.options.get_status()
		self.settings['options'] = status
		vr_status = self.vroptions.get_status()
		self.settings['vr_options'] = vr_status
		
		z_gr = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		z_mf = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		if self['aky']:
			for x_id, x_value in enumerate(self.x_axis):
				for y_id, y_value in enumerate(self.y_axis):
					run = self['run'].copy()
					run[self['x_axis_type']] = x_value
					run[self['y_axis_type']] = y_value
					run_id = self.reader.get_run_id(run = run, keys = '_run_keys')
					gr = self.reader.gyro_data[run_id]['growth_rate']
					z_mf[x_id][y_id] = self.reader.gyro_data[run_id]['mode_frequency']
					if status[1]:
						if 'ky' in self.reader.gyro_data[run_id]:
							ky = self.reader.gyro_data[run_id]['ky']
						else:
							ky = self.reader.single_parameters['ky'].values[0]
						z_gr[x_id][y_id] = gr/ky**2
					else:
						z_gr[x_id][y_id] = gr
		else:
			for x_id, x_value in enumerate(self.x_axis):
				for y_id, y_value in enumerate(self.y_axis):
					run = self['run'].copy()
					run[self['x_axis_type']] = x_value
					run[self['y_axis_type']] = y_value
					
					if status[1]:
						run_id = self.reader.get_run_id(run = run, keys = '_norm_gr_keys')
						gr = self.reader.gyro_data[run_id]['growth_rate']
						z_mf[x_id][y_id] = self.reader.gyro_data[run_id]['mode_frequency']
						if 'ky' in self.reader.gyro_data[run_id]:
							ky = self.reader.gyro_data[run_id]['ky']
						else:
							ky = self.reader.single_parameters['ky'].values[0]
						z_gr[x_id][y_id] = gr/ky**2
					else:
						run_id = self.reader.get_run_id(run = run, keys = '_abs_gr_keys')
						gr = self.reader.gyro_data[run_id]['growth_rate']
						z_gr[x_id][y_id] = gr
		z_gr = transpose(z_gr)
		z_mf = transpose(z_mf)
		
		if status[1]:
			self.ax[0].set_title("Growth Rate/$(k_{y}\\rho_{0})^{2}$",fontsize=self['fontsizes']['title'])
		else:
			self.ax[0].set_title("Growth Rate",fontsize=self['fontsizes']['title'])
		
		if self['gr_slider']['max']:
			grmax = self.sliders['gr_slider'].val * self['gr_slider']['max']/100
		#elif status[2]:
			#grmax = self.sliders['gr_slider'].val * abs(amax(array(self.data['growth_rates'])[isfinite(self.data['growth_rates'])]))/100
		else:
			try:
				grmax = self.sliders['gr_slider'].val * amax(abs(array(z_gr)[isfinite(z_gr)]))/100
				if grmax < 10e-10:
					grmax = 10e-10
			except:
				grmax = 1
		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		self.settings['gr_slider']['scale'] = self.sliders['gr_slider'].val
		self.cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = self.cmap))
		if self['contour_type'] == 1:
			self.ax[0].contourf(x_axis, y_axis, z_gr, cmap = self.cmap, norm=norm_gr)
		else:
			self.ax[0].pcolormesh(x_axis, y_axis, z_gr, cmap = self.cmap, norm=norm_gr)
		
		if self['mf_slider']['max']:
			mfmax = self.sliders['mf_slider'].val * self['mf_slider']['max']/100
		#elif status[2]:
			#mfmax = self.sliders['mf_slider'].val * abs(amax(array(self.data['mode_frequencies'])[isfinite(self.data['mode_frequencies'])]))/100
		else:
			try:
				mfmax = self.sliders['mf_slider'].val * amax(abs(array(z_mf)[isfinite(z_mf)]))/100
				if mfmax == 0:
					mfmax = 10e-10
			except:
				mfmax = 1
		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		self.settings['mf_slider']['scale'] = self.sliders['mf_slider'].val
		self.cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		self.ax[1].pcolormesh(x_axis, y_axis, z_mf, norm = norm_mf)
		
		if status[0]:
			if self['aky']:
				parities = array(self.data['parities_all'])[idx,:,:,ky_idx]
			else:
				parities = array(self.data['parities'][idx])
			sym_ids = where(parities == 1)
			anti_ids = where(parities == -1)
			self.ax[0].plot(array(x_axis)[sym_ids[0]], array(y_axis)[sym_ids[1]], '+', color = 'purple', label = 'par')
			self.ax[1].plot(array(x_axis)[sym_ids[0]], array(y_axis)[sym_ids[1]], '+', color = 'purple', label = 'par')
			self.ax[0].plot(array(x_axis)[anti_ids[0]], array(y_axis)[anti_ids[1]], '_', color = 'cyan', label = 'par')
			self.ax[1].plot(array(x_axis)[anti_ids[0]], array(y_axis)[anti_ids[1]], '_', color = 'cyan', label = 'par')
		
		if status[3]:
			self.ax[0].plot(x_val,y_val,'kx')
			self.ax[1].plot(x_val,y_val,'kx')
			if self['eqbm_style'] in ["point",2,"point numless",4]:
				self.ax[0].annotate("Eqbm",(x_val,y_val),textcoords = "offset points",xytext = (0,7), ha = "center")
				self.ax[1].annotate("Eqbm",(x_val,y_val),textcoords = "offset points",xytext = (0,7), ha = "center")
			if self['eqbm_style'] in ["split",1,"point",2]:
				self.ax[0].annotate(f"{x_val:.2f},{y_val:.2f}",(x_val,y_val),textcoords = "offset points",xytext = (0,-13), ha = "center")
				self.ax[1].annotate(f"{x_val:.2f},{y_val:.2f}",(x_val,y_val),textcoords = "offset points",xytext = (0,-13), ha = "center")
			if self['eqbm_style'] in ["title",0]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label=f"Equillibrium ({x_val:.2f},{y_val:.2f})",linewidth=0))
			if self['eqbm_style'] in ["split",1,"title numless",3]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label="Equillibrium",linewidth=0))
		
		if status[4]:
			if self.reader.data['ideal_data'] is not None and self.reader.data['ideal_data'][psiN] is not None:
				idata = self.reader.data['ideal_data'][psiN]
				self.ax[0].contourf(idata[self['x_axis_type']], idata[self['y_axis_type']], idata['stabilities'], [0.01,0.99], colors = ('k'))
				self.ax[1].contourf(idata[self['x_axis_type']], idata[self['y_axis_type']], idata['stabilities'], [0.01,0.99], colors = ('k'))
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label="Ideal Boundary"))
			else:
				self.ax.text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=self.ax.transAxes,color='k')
		
		if vr_status[0]:
			dx = x_axis[1] - x_axis[0]
			dy = y_axis[1] - y_axis[0]
			self.ax[0].set_xticks(array(x_axis)-dx/2, minor=True)
			self.ax[1].set_xticks(array(x_axis)-dx/2, minor=True)
			self.ax[0].set_yticks(array(y_axis)-dy/2, minor=True)
			self.ax[1].set_yticks(array(y_axis)-dy/2, minor=True)
			
			self.ax[0].set_xticks(array(x_axis), minor=False)
			self.ax[1].set_xticks(array(x_axis), minor=False)
			self.ax[0].set_yticks(array(y_axis), minor=False)
			self.ax[1].set_yticks(array(y_axis), minor=False)
			
			self.ax[0].set_xticklabels([])
			self.ax[1].set_xticklabels([])
			self.ax[0].set_yticklabels([])
			self.ax[1].set_yticklabels([])
			xlabels = [str(i) for i in range(len(x_axis))]
			ylabels = [str(i) for i in range(len(y_axis))]
			self.ax[0].set_xticklabels(xlabels,minor=False)
			self.ax[1].set_xticklabels(xlabels,minor=False)
			self.ax[0].set_yticklabels(ylabels,minor=False)
			self.ax[1].set_yticklabels(ylabels,minor=False)
			
			self.ax[0].grid(which="minor",color='k')
			self.ax[1].grid(which="minor",color='k')
			
		if vr_status[1]:
			if self['aky']:
				un_x = [i for p,i,j,k in self.verify['unconverged'] if p == idx and k == ky_idx]
				un_y = [j for p,i,j,k in self.verify['unconverged'] if p == idx and k == ky_idx]
				uns_x = [i for p,i,j,k in self.verify['unconverged_stable'] if p == idx and k == ky_idx]
				uns_y = [j for p,i,j,k in self.verify['unconverged_stable'] if p == idx and k == ky_idx]
				co_x = [i for p,i,j,k in self.verify['converged'] if p == idx and k == ky_idx]
				co_y = [j for p,i,j,k in self.verify['converged'] if p == idx and k == ky_idx]
				cof_x = [i for p,i,j,k in self.verify['converged_fit'] if p == idx and k == ky_idx]
				cof_y = [j for p,i,j,k in self.verify['converged_fit'] if p == idx and k == ky_idx]
			else:
				un_x = [i for p,i,j,k in self.verify['unconverged'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				un_y = [j for p,i,j,k in self.verify['unconverged'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				uns_x = [i for p,i,j,k in self.verify['unconverged_stable'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				uns_y = [j for p,i,j,k in self.verify['unconverged_stable'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				co_x = [i for p,i,j,k in self.verify['converged'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				co_y = [j for p,i,j,k in self.verify['converged'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				cof_x = [i for p,i,j,k in self.verify['converged_fit'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				cof_y = [j for p,i,j,k in self.verify['converged_fit'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
			
			self.ax[0].plot(array(x_axis)[un_x], array(y_axis)[un_y], 'c^', label = 'U')
			self.ax[1].plot(array(x_axis)[un_x], array(y_axis)[un_y], 'c^', label = 'U')
			self.ax[0].plot(array(x_axis)[uns_x], array(y_axis)[uns_y], 'bv', label = 'Us')
			self.ax[1].plot(array(x_axis)[uns_x], array(y_axis)[uns_y], 'bv', label = 'Us')
			self.ax[0].plot(array(x_axis)[co_x], array(y_axis)[co_y], 'k>', label = 'C')
			self.ax[1].plot(array(x_axis)[co_x], array(y_axis)[co_y], 'k>', label = 'C')
			self.ax[0].plot(array(x_axis)[cof_x], array(y_axis)[cof_y], 'm<', label = 'Cf')
			self.ax[1].plot(array(x_axis)[cof_x], array(y_axis)[cof_y], 'm<', label = 'Cf')
			self.ax[0].legend(loc = 'center', bbox_to_anchor=(0.5,1.1), ncol = 4)
		
		if vr_status[2]:
			if self['aky']:
				ns_x = [i for p,i,j,k in self.verify['nstep'] if p == idx and k == ky_idx]
				ns_y = [j for p,i,j,k in self.verify['nstep'] if p == idx and k == ky_idx]
			else:
				ns_x = [i for p,i,j,k in self.verify['nstep'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				ns_y = [j for p,i,j,k in self.verify['nstep'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
			self.ax[0].plot(array(x_axis)[ns_x], array(y_axis)[ns_y], 'kX')
			self.ax[1].plot(array(x_axis)[ns_x], array(y_axis)[ns_y], 'kX')
		
		if vr_status[3]:
			for bpid, bp in enumerate(x_axis):
				for shid, sh in enumerate(y_axis):
					s = ''
					if not self['aky']:
						ky_idx = self.inputs['aky_values'].index(self.data['akys'][idx][bpid][shid])
					if (idx, bpid, shid, ky_idx) in self.verify['phi']:
						s += 'p,'
					if (idx, bpid, shid, ky_idx) in self.verify['apar']:
						s += 'a,'
					if (idx, bpid, shid, ky_idx) in self.verify['bpar']:
						s += 'b,'
					if s != '':
						self.ax[0].text(bp, sh, s[:-1], color = 'k',ha='center',va='center',size=7, rotation = 45)
						self.ax[1].text(bp, sh, s[:-1], color = 'k',ha='center',va='center',size=7, rotation = 45)
		
		if vr_status[4]:
			if self['aky']:
				na_x = [i for p,i,j,k in self.verify['nan'] if p == idx and k == ky_idx]
				na_y = [j for p,i,j,k in self.verify['nan'] if p == idx and k == ky_idx]
			else:
				na_x = [i for p,i,j,k in self.verify['nan'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				na_y = [j for p,i,j,k in self.verify['nan'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
			self.ax[0].plot(array(x_axis)[na_x], array(y_axis)[na_y], 'kX')
			self.ax[1].plot(array(x_axis)[na_x], array(y_axis)[na_y], 'kX')
		
		self.ax[0].legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0,1.02),loc = "lower left", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax[0].legend_.set_visible(self['visible']['legend'])
		
		self.fig.canvas.draw_idle()
