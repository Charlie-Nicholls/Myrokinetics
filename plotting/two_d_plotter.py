from numpy import transpose, array, amax, amin, isfinite, linspace, where, full, nan, diff
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons, TextBox
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import RegularGridInterpolator
from copy import deepcopy
from .slider_ax import slider_axes

default_settings = {"suptitle": None,
		"eqbm_style": "title",
		"contour_type": 0,
		"x_axis_type": "beta_prime",
		"y_axis_type": "shear",
		"z_slider": {"scale": 100, "max": None},
		"run": {},
		"options": [False,False,True,True,False,False],
		"fontsizes": {"title": 13, "ch_box": 8,"axis": 17,"suptitle": 20},
		"visible": {"z_slider": True, "op_box": True, "suptitle": True, "title": True},
		"cdict": {"red": ((0.0, 1, 1),(1.0, 0.8, 0.8)),
			"green": ((0.0, 1, 1),(1.0, 0.0, 0.0)),
			"blue": ((0.0, 1, 1),(1.0, 0.0, 0.0))},
}

slider_settings = {"slider_1": {"axis": [0.15, 0.01, 0.5, 0.03]},
		"slider_2": {"axis": [0.15, 0.05, 0.5, 0.03]},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": False, "slider_4": False, "slider_5": False, "slider_6": False, "slider_7": False, "slider_8": False, "slider_9": False},
}

class plot_2d(object):
	def __init__(self, reader, settings = {}, sliders = None):	
		self.reader = reader
		self.sliders = slider
		if self.reader['quasilinear'] is None:
			print("Error: No QuasiLinear Data")
			return
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		self._options = ["Show ID","Global Axis Limits","Global Colourbar","Show Equillibrium","Show Ideal","Draw EQBM Contour"]

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
			filename = f"self['z_axis_type']_{self['slider_1']['id']}_{self['slider_2']['id']}"
		self.fig.savefig(filename)
		
	def open_plot(self):
		self.fig, self.ax = subplots(figsize=(9,7))
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
	
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		self._load_z_axis(self['z_axis_type'])
		
		self.dims = [x for x in self.reader.inputs.dim_order if x not in [self['x_axis_type'],self['y_axis_type'],self['z_axis_type']]]
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

		self.cmap = LinearSegmentedColormap('WhRd', self.settings['cdict'])
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax)
		
		self.ch_axes = axes([0.75, 0.01, 0.09, 0.09],frame_on = False,visible=self['visible']['op_box'])
		self.options = CheckButtons(self.ch_axes, self._options, self['options'])
		self.options.on_clicked(self.draw_fig)
		self.set_options_fontsize(self['fontsizes']['ch_box'])
		
		
		self._sliders = {}
		self.z_axes = axes([0.9, 0.13, 0.01, 0.73],visible=self['visible']['z_slider'])
		self._sliders['z_slider'] = Slider(self.z_axes, 'Scale', 0, 100, valinit = self['z_slider']['scale'], valstep = 1, orientation = 'vertical')
		self._sliders['z_slider'].on_changed(self.draw_fig)
			
		zs = [x for x in self.reader.data['quasilinear'].values() if str(x) not in ['-inf','inf','nan']]
		self._z_max = max(zs,default=1)
		
		ion()
		self.draw_fig()
		if self['z_slider']['max']:
			self.set_z_max(self.settings['z_slider']['max'])
		show()
	
	def set_slider(self, num = None, key = None, dimension_type = None):
		self.sliders.set_slider(num = num, key = key, dimension_type = dimension_type)
		
	def _load_x_axis(self, axis_type):
		if axis_type not in ['beta_prime','alpha']:
			print(f"ERROR: axis_type not found, valid types ['beta_prime','alpha']")
			return
			
		self.settings['x_axis_type'] = axis_type
		
		if axis_type in ['alpha']:
			print("Not yet implimented")
			self._x_axis_label = r'$\alpha$'
		else:
			self.x_axis = self.reader.dimensions['beta_prime'].values
			self._x_axis_label = self.reader.dimensions['beta_prime'].axis_label
			
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
			self._y_axis_label = self.reader.dimensions['shear'].axis_label
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def _load_z_axis(self, axis_type):
		if axis_type not in ['growth_rate','growth_rate_norm','mode_frequency','quasilinear','ql_norm','ql_metric']:
			print(f"ERROR: axis_type not found, valid types ['growth_rate','growth_rate_norm','ql_norm','ql_metric','mode_frequency']")
			return
			
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
		elif axis_type in ['mode_frequency']:
			self._z_axis_label = "Mode Frequency"
		elif axis_type in ['ql_norm']:
			self._z_axis_label = "QL norm"
		elif axis_type in ['ql_metric']:
			self._z_axis_label = "QL norm (gs2)"
		elif axis_type in ['quasilinear']:
			self._z_axis_label = "Quasilinear Metric"
			self.dims = [x for x in self.dims if x not in ['ky']]#&theta0
			if 'ky' in self['run']:
				self.settings['run'].pop('ky')
			#remove ky sliders
			#if 'theta0' in self['run']:
				#self.settings['run'].pop('theta0')
			
	def set_z_axis_type(self, axis_type):
		self._load_z_axis(axis_type)
		self.draw_fig()
	
	def set_max(self, val = None):
		self.settings['z_slider']['max'] = val
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
			if key == 'z_slider':
				self.z_axes.set_visible(self['visible']['z_slider'])
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
		for key in self.sliders.sliders:
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
				z[x_id][y_id] = self.reader(self['z_axis_type'],run)
		z = transpose(z)
		self.z_axis = z
		
		if self['z_slider']['max']:
			z_max = self._sliders['z_slider'].val * self['z_slider']['max']/100
		elif status[2]:
			z_max = self._sliders['z_slider'].val * self._z_max/100
		else:
			try:
				z_max = self._sliders['z_slider'].val * amax(abs(array(z)[isfinite(z)]))/100
				if z_max < 10e-10:
					z_max = 10e-10
			except:
				z_max = 1
		
		norm = Normalize(vmin=0,vmax=z_max)	
		self.settings['z_slider']['scale'] = self._sliders['z_slider'].val
		self.cbar.update_normal(ScalarMappable(norm = norm, cmap = self.cmap))
		if self['contour_type'] == 1:
			self.ax.contourf(x_axis,y_axis,z, cmap = self.cmap, norm = norm)
		else:
			self.ax.pcolormesh(x_axis, y_axis, z, cmap = self.cmap, norm=norm)
			
		if status[5]:
			grid = RegularGridInterpolator([x_axis,y_axis],transpose(z))
			if min(x_axis) > x_val or max(x_axis) < x_val or min(y_axis) > y_val or max(y_axis) < y_val:
				print("ERRROR: eqbm point outside scan range")
			else:
				eqbm_val = grid((x_val,y_val))
				self.ax.contourf(x_axis, y_axis, z, [eqbm_val-eqbm_val/100,eqbm_val+eqbm_val/100], colors = ('b'))
			
		if status[0]:
			xticks = [x+dx/2 for dx,x in zip(diff(x_axis),x_axis)]
			yticks = [y+dy/2 for dy,y in zip(diff(y_axis),y_axis)]
			self.ax.set_xticks(xticks, minor=True)
			self.ax.set_yticks(yticks, minor=True)
			
			self.ax.set_xticks(x_axis, minor=False)
			self.ax.set_yticks(y_axis, minor=False)
			
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
				eqbm_pos += f",{eqbm_val:.2f}"
			if self['eqbm_style'] in ["point",2,"point numless",4]:
				self.ax.annotate("Eqbm",(x_val,y_val),textcoords = "offset points",xytext = (0,7), ha = "center")
			if self['eqbm_style'] in ["split",1,"point",2]:
				self.ax.annotate(eqbm_pos,(x_val,y_val),textcoords = "offset points",xytext = (0,-13), ha = "center")
			if self['eqbm_style'] in ["title",0]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label=f"Equillibrium ({eqbm_pos})",linewidth=0))
			if self['eqbm_style'] in ["split",1,"title numless",3]:
				handles.append(Line2D([0.5],[0.5],marker='x',color='k',label="Equillibrium",linewidth=0))
		
		if status[4]:
			if 'theta0' in self['run']:
				theta0 = self['run']['theta0']
			else:
				theta0 = self.reader.data['_ideal_keys']['theta0'].keys()[0]
			run_id = self.reader.get_run_id(run={'psin': psiN,'theta0': theta0},key='_ideal_keys')
			if run_id is not None:
				idata = self.reader.data['ideal'][run_id]
				self.ax.contourf(idata[self['x_axis_type']], idata[self['y_axis_type']], idata['stabilities'], [0.01,0.99], colors = ('k'))
				handles.append(Line2D([0,1],[0.5,0.5],color='k',label="Ideal Boundary"))
			else:
				self.ax.text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=self.ax.transAxes,color='k')
		
		self.ax.legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0.5,0.98),loc = "lower center", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax.legend_.set_visible(self['visible']['title'])

		self.fig.canvas.draw_idle()
