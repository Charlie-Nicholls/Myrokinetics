from numpy import transpose, array, amax, isfinite, full, nan, diff
from matplotlib.pyplot import subplots, ion, show, Normalize, colorbar, axes, Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap
from copy import deepcopy
from .slider_ax import slider_axes

default_settings = {"suptitle": None,
		"eqbm_style": "title",
		"aky": False,
		"contour_type": 0,
		"x_axis_type": "beta_prime",
		"y_axis_type": "shear",
		"z_axis_type": "growth_rate",
		"gr_slider": {"scale": 100, "max": None},
		"mf_slider": {"scale": 100, "max": None},
		"run": {},
		"options": [False,False,True,False],
		"vr_options": [False,False,False,False,False],
		"fontsizes": {"title": 13, "ch_box": 8, "vr_box": 8, "axis": 17, "suptitle": 20, "legend": 13},
		"visible": {"slider_1": True, "slider_2": True, "gr_slider": True, "mf_slider": True, "op_box": True, "vr_box": True, "suptitle": True, "title": True, "legend": True},
		"cdict": {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
			'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
			'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))},
}

slider_settings = {"slider_1": {"axis": [0.25, 0.01, 0.5, 0.03]},
		"slider_2": {"axis": [0.25, 0.05, 0.5, 0.03]},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": False, "slider_4": False, "slider_5": False, "slider_6": False, "slider_7": False, "slider_8": False, "slider_9": False},
}

class plot_scan(object):
	def __init__(self, reader = None, settings = {}, sliders = None):
		self.reader = reader
		self.verify = reader.verify
		self.data = reader.data['gyro']
		self.sliders = sliders
		if self.data is None:
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
		self._options = ["Show Parities","Global Colorbar","Show Equillibrium","Show Ideal"]
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
		if filename is None:
			filename = f"Scan_{self['slider_1']['id']}_{self['slider_2']['id']}"
		self.fig.savefig(filename)
		
	def open_plot(self):
		self.fig, self.ax = subplots(1,2,figsize=(14.6,7))
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		self._load_z_axis(self['z_axis_type'])
		
		if self['aky']:
			self.dims = [x for x in self.reader.inputs.dim_order if x not in [self['x_axis_type'],self['y_axis_type']]]
		else:
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
		
		
		self._sliders = {}
		self.gr_axes = axes([0.93, 0.15, 0.01, 0.73])
		self._sliders['gr_slider'] = Slider(self.gr_axes, 'GR', 0, 100, valinit = self['gr_slider']['scale'], valstep = 1, orientation = 'vertical')
		self._sliders['gr_slider'].on_changed(self.draw_fig)
		
		self.mf_axes = axes([0.97, 0.15, 0.01, 0.73])
		self._sliders['mf_slider'] = Slider(self.mf_axes, 'MF', 0, 100, valinit = self['mf_slider']['scale'], valstep = 1, orientation = 'vertical')
		self._sliders['mf_slider'].on_changed(self.draw_fig)

		mfs = []
		for run in self.reader.get_all_runs():
			val = self.reader('mode_frequency',run)
			if str(val) not in ['-inf','inf','nan']:
				mfs.append(val)
		self._mf_max = max(mfs,default=1)
		self.draw_fig()
		ion()
		show()
	
	def set_slider(self, num = None, key = None, dimension_type = None):
		self.sliders.set_slider(num = num, key = key, dimension_type = dimension_type)
	
	def _load_x_axis(self, axis_type):
		if axis_type not in ['beta_prime','alpha']:
			print("ERROR: axis_type not found, valid types ['beta_prime','alpha']")
			return
			
		self.settings['x_axis_type'] = axis_type
		
		if axis_type in ['alpha']:
			if not self.data['alpha_axis']:
				print("ERROR: Alpha not calculated, use calculate_alpha()")
			else:
				#self.x_axis = self.data['alpha_axis'] needs updating
				self._x_axis_label = r'$\alpha$'
		elif axis_type in ['beta_prime']:
			self.x_axis = self.reader.dimensions['beta_prime'].values
			self._x_axis_label = self.reader.dimensions['beta_prime'].axis_label	
			
	def set_x_axis_type(self, axis_type):
		self._load_x_axis(axis_type)
		self.draw_fig()
	
	def _load_y_axis(self, axis_type):
		if axis_type not in ['shear','current']:
			print("ERROR: axis_type not found, valid types ['shear','current']")
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
	
	def _load_z_axis(self, axis_type):
		if axis_type not in ['growth_rate','growth_rate_norm','ql_norm']:
			print("ERROR: axis_type not found, valid types ['growth_rate','growth_rate_norm','ql_norm']")
			return
		if axis_type == 'ql_norm' and not self['aky']:
			print("ERROR: ql_norm axis type only allowed on seperate ky plots")
			
		self.settings['z_axis_type'] = axis_type
		
		zs = []
		for run in self.reader.get_all_runs():
			val = self.reader(self['z_axis_type'],run)
			if str(val) not in ['-inf','inf','nan']:
				zs.append(val)
		self._z_max = max(zs,default=100)
		
		if axis_type in ['growth_rate']:
			self._z_axis_label = r'$\gamma$'
		elif axis_type in ['growth_rate_norm']:
			self._z_axis_label = r'$\gamma/k_{y}^{2}$'
		elif axis_type in ['ql_norm']:
			self._z_axis_label = "QL Metric"
			
	def set_z_axis_type(self, axis_type):
		self._load_z_axis(axis_type)
		self.draw_fig()
	
	def set_gr_max(self, val):
		self.settings['gr_slider']['max'] = val
		self.draw_fig()
		
	def set_mf_max(self, val):
		self.settings['mf_slider']['max'] = val
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
			if key == 'op_box':
				self.ch_axes.set_visible(self['visible']['op_box'])
			elif key == 'vr_box':
				self.vr_axes.set_visible(self['visible']['vr_box'])
			elif key == 'suptitle':
				self.fig._suptitle.set_visible(self['visible']['suptitle'])
			elif key == 'title':
				self.ax.legend_.set_visible(self['visible']['title'])
			elif key == 'gr_slider':
				self.gr_axes.set_visible(self['visible']['gr_slider'])
			elif key == 'mf_slider':
				self.mf_axes.set_visible(self['visible']['mf_slider'])
		self.draw_fig()				

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
		for key, sli in self.sliders.sliders.items():
			dim = self.sliders.settings[key]['dimension_type']
			if dim in self.dims:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label=f"{self.reader.dimensions[dim].axis_label} = {self.settings['run'][dim]}",visible = False))
		
		if 'psin' in self['run']:
			psiN = self['run']['psin']
		else:
			psiN = self.reader.single_parameters['psin'].values[0]
		
		x_axis = list(self.x_axis)
		y_axis = list(self.y_axis)
		
		x_val = abs(self.reader.data['equilibrium'][psiN][self['x_axis_type']])
		y_val = self.reader.data['equilibrium'][psiN][self['y_axis_type']]
		
		status = self.options.get_status()
		self.settings['options'] = status
		vr_status = self.vroptions.get_status()
		self.settings['vr_options'] = vr_status
			
		self.ax[0].cla()
		self.ax[0].set_facecolor('grey')
		self.ax[0].set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[0].set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].cla()
		self.ax[1].set_facecolor('grey')
		self.ax[1].set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax[1].set_title("Mode Frequency",fontsize=self['fontsizes']['title'])
		
		#self.ax[0].set_ylim(amin(self.y_axis),amax(self.y_axis))
		#self.ax[0].set_xlim(amin(self.x_axis),amax(self.x_axis))
		#self.ax[1].set_ylim(amin(self.y_axis),amax(self.y_axis))
		#self.ax[1].set_xlim(amin(self.x_axis),amax(self.x_axis))
		
		
		
		z_gr = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		z_mf = full((len(self.reader.dimensions[self['x_axis_type']]),len(self.reader.dimensions[self['y_axis_type']])),nan)
		run_ids = []
		for x_id, x_value in enumerate(self.x_axis):
			for y_id, y_value in enumerate(self.y_axis):
				run = self['run'].copy()
				run[self['x_axis_type']] = x_value
				run[self['y_axis_type']] = y_value
				run_ids.append(self.reader.get_run_id(run))
				z_type = self['z_axis_type']
				if not self['aky'] and z_type == 'growth_rate':
					z_type = 'abs_gr'
				elif not self['aky'] and z_type == 'growth_rate_norm':
					z_type = 'norm_gr'
				z_mf[x_id][y_id] = self.reader('mode_frequency',run) if self.reader('mode_frequency',run) is not None else nan
				z_gr[x_id][y_id] = self.reader(z_type,run) if self.reader(z_type,run) is not None else nan
		
		self.z_axis_gr = z_gr = transpose(z_gr)
		self.z_axis_mf = z_mf = transpose(z_mf)
		self.ax[0].set_title(self._z_axis_label,fontsize=self['fontsizes']['title'])
		
		if self['gr_slider']['max']:
			grmax = self._sliders['gr_slider'].val * self['gr_slider']['max']/100
		elif status[1]:
			grmax = self._sliders['gr_slider'].val * self._z_max/100
		else:
			try:
				grmax = self._sliders['gr_slider'].val * amax(abs(array(z_gr)[isfinite(z_gr)]))/100
				if grmax < 10e-10:
					grmax = 10e-10
			except:
				grmax = 1
		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		self.settings['gr_slider']['scale'] = self._sliders['gr_slider'].val
		self.cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = self.cmap))
		if self['contour_type'] == 1:
			self.ax[0].contourf(x_axis, y_axis, z_gr, cmap = self.cmap, norm=norm_gr)
		else:
			self.ax[0].pcolormesh(x_axis, y_axis, z_gr, cmap = self.cmap, norm=norm_gr)
		
		if self['mf_slider']['max']:
			mfmax = self._sliders['mf_slider'].val * self['mf_slider']['max']/100
		elif status[1]:
			mfmax = self._sliders['mf_slider'].val * self._mf_max/100
		else:
			try:
				mfmax = self._sliders['mf_slider'].val * amax(abs(array(z_mf)[isfinite(z_mf)]))/100
				if mfmax == 0:
					mfmax = 10e-10
			except:
				mfmax = 1
		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		self.settings['mf_slider']['scale'] = self._sliders['mf_slider'].val
		self.cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		self.ax[1].pcolormesh(x_axis, y_axis, z_mf, norm = norm_mf)
		
		if status[0]:
			sym_x = []
			sym_y = []
			ant_x = []
			ant_y = []
			for run_id in run_ids:
				run = self.data[run_id]
				if 'parity' in run:
					if run['parity'] == 1:
						sym_x.append(run[self['x_axis_type']])
						sym_y.append(run[self['y_axis_type']])
					elif run['parity'] == -1:
						ant_x.append(run[self['x_axis_type']])
						ant_y.append(run[self['y_axis_type']])
			self.ax[0].plot(sym_x, sym_y, '+', color = 'purple', label = 'par')
			self.ax[1].plot(sym_x, sym_y, '+', color = 'purple', label = 'par')
			self.ax[0].plot(ant_x, ant_y, '_', color = 'cyan', label = 'par')
			self.ax[1].plot(ant_x, ant_y, '_', color = 'cyan', label = 'par')
		
		if status[2]:
			self.ax[0].plot(x_val,y_val,'kx')
			self.ax[1].plot(x_val,y_val,'kx')
			if self['eqbm_style'] in ["point",2,"point numless",4]:
				self.ax[0].annotate("Eqbm",(x_val,y_val),textcoords = "offset points",xytext = (0,7), ha = "center")
				self.ax[1].annotateself.ax[0].set_title(self._z_axis_label,fontsize=self['fontsizes']['title'])("Eqbm",(x_val,y_val),textcoords = "offset points",xytext = (0,7), ha = "center")
			if self['eqbm_style'] in ["split",1,"point",2]:
				self.ax[0].annotate(f"{x_val:.2f},{y_val:.2f}",(x_val,y_val),textcoords = "offset points",xytext = (0,-13), ha = "center")
				self.ax[1].annotate(f"{x_val:.2f},{y_val:.2f}",(x_val,y_val),textcoords = "offset points",xytext = (0,-13), ha = "center")
			if self['eqbm_style'] in ["title",0]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label=f"Equillibrium ({x_val:.2f},{y_val:.2f})",linewidth=0))
			if self['eqbm_style'] in ["split",1,"title numless",3]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label="Equillibrium",linewidth=0))
		
		if status[3]:
			if 'theta0' in self['run']:
				theta0 = self['run']['theta0']
			else:
				theta0 = list(self.reader.data['_ideal_keys']['theta0'].keys())[0]
			run_id = self.reader.get_run_id(run={'psin': psiN,'theta0': theta0},keys='_ideal_keys')
			if run_id is not None:
				idata = self.reader.data['ideal'][run_id]
				self.ax[0].contourf(idata[self['x_axis_type']], idata[self['y_axis_type']], idata['stabilities'], [0.01,0.99], colors = ('k'))
				self.ax[1].contourf(idata[self['x_axis_type']], idata[self['y_axis_type']], idata['stabilities'], [0.01,0.99], colors = ('k'))
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label="Ideal Boundary"))
			else:
				self.ax.text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=self.ax.transAxes,color='k')
		
		if vr_status[0]:
			xticks = [x+dx/2 for dx,x in zip(diff(x_axis),x_axis)]
			yticks = [y+dy/2 for dy,y in zip(diff(y_axis),y_axis)]
			self.ax[0].set_xticks(xticks, minor=True)
			self.ax[1].set_xticks(xticks, minor=True)
			self.ax[0].set_yticks(yticks, minor=True)
			self.ax[1].set_yticks(yticks, minor=True)
			
			self.ax[0].set_xticks(x_axis, minor=False)
			self.ax[1].set_xticks(x_axis, minor=False)
			self.ax[0].set_yticks(y_axis, minor=False)
			self.ax[1].set_yticks(y_axis, minor=False)
			
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
			
		if vr_status[1] and self.verify is not None:
			un_x = [self.data[run_id][self['x_axis_type']] for run_id in run_ids if run_id in self.verify['unconverged']]
			un_y = [self.data[run_id][self['y_axis_type']] for run_id in run_ids if run_id in self.verify['unconverged']]
			uns_x = [self.data[run_id][self['x_axis_type']] for run_id in run_ids if run_id in self.verify['unconverged_stable']]
			uns_y = [self.data[run_id][self['y_axis_type']] for run_id in run_ids if run_id in self.verify['unconverged_stable']]
			co_x = [self.data[run_id][self['x_axis_type']] for run_id in run_ids if run_id in self.verify['converged']]
			co_y = [self.data[run_id][self['y_axis_type']] for run_id in run_ids if run_id in self.verify['converged']]
			cof_x = [self.data[run_id][self['x_axis_type']] for run_id in run_ids if run_id in self.verify['converged_fit']]
			cof_y = [self.data[run_id][self['y_axis_type']] for run_id in run_ids if run_id in self.verify['converged_fit']]
			
			self.ax[0].plot(un_x, un_y, 'c^', label = 'U')
			self.ax[1].plot(un_x, un_y, 'c^', label = 'U')
			self.ax[0].plot(uns_x, uns_y, 'bv', label = 'Us')
			self.ax[1].plot(uns_x, uns_y, 'bv', label = 'Us')
			self.ax[0].plot(co_x, co_y, 'k>', label = 'C')
			self.ax[1].plot(co_x, co_y, 'k>', label = 'C')
			self.ax[0].plot(cof_x, cof_y, 'm<', label = 'Cf')
			self.ax[1].plot(cof_x, cof_y, 'm<', label = 'Cf')
			self.ax[0].legend(loc = 'center', bbox_to_anchor=(0.5,1.1), ncol = 4)
		
		if vr_status[2]:
			ns_x = [self.data[run_id][self['x_axis_type']] for run_id in run_ids if run_id in self.verify['nstep']]
			ns_y = [self.data[run_id][self['y_axis_type']] for run_id in run_ids if run_id in self.verify['nstep']]

			self.ax[0].plot(ns_x, ns_y, 'kX')
			self.ax[1].plot(ns_x, ns_y, 'kX')
		
		if vr_status[3]:
			for run_id in run_ids:
				s = ''
				if run_id in self.verify['phi']:
					s += 'p,'
				if run_id in self.verify['apar']:
					s += 'a,'
				if run_id in self.verify['bpar']:
					s += 'b,'
				if s != '':
					self.ax[0].text(self.data[run_id][self['x_axis_type']], self.data[run_id][self['y_axis_type']], s[:-1], color = 'k',ha='center',va='center',size=7, rotation = 45)
					self.ax[1].text(self.data[run_id][self['x_axis_type']], self.data[run_id][self['y_axis_type']], s[:-1], color = 'k',ha='center',va='center',size=7, rotation = 45)
		
		if vr_status[4]:
			na_x = [self.data[run_id][self['x_axis_type']] for run_id in run_ids if run_id in self.verify['nan']]
			na_y = [self.data[run_id][self['y_axis_type']] for run_id in run_ids if run_id in self.verify['nan']]
			self.ax[0].plot(na_x, na_y, 'kX')
			self.ax[1].plot(na_x, na_y, 'kX')
		
		self.ax[0].legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0,1.02),loc = "lower left", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax[0].legend_.set_visible(self['visible']['legend'])
		
		self.fig.canvas.draw_idle()
