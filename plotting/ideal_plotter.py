from numpy import *
from matplotlib.pyplot import *
import matplotlib.patches as pt
from matplotlib.widgets import Slider, CheckButtons
from time import sleep
from copy import deepcopy

default_settings = {"suptitle": None,
		"eqbm_style": "split",
		"x_axis_type": 'beta_prime',
		"y_axis_type": 'shear',
		"slider_1": {"dimension_type": None, "id": 0},
		"run": {},
		"options": [True,False,True,False],
		"fontsizes": {"title": 13, "ch_box": 8,"axis": 17,"suptitle": 20},
		"visible": {"slider_1": True, "op_box": True, "suptitle": True, "title": True},
		"colours": {"unstable": 'r', "stable": 'g', "boundary": 'k'}
}

class plot_ideal(object):
	def __init__(self, reader, settings = {}):
		self.reader = reader
		if self.reader.data['ideal_data'] is None:
			print("Error: No ideal_ball data")
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		self._options = ["Colour","Global Axis Limits","Show Equillibrium","Show Legend"]
		
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
			if key not in settings:
				self.settings[key] = defaults[key]
				self.init_settings[key] = defaults[key]
			elif type(self.settings[key]) == dict:
				for skey in self.defaults[key]:
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
	
	def open_plot(self, save = False, filename = None):
		self.fig, self.ax = subplots(figsize=(10,7))
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		
		self.dims = [x for x in self.reader.inputs.dim_order if x == 'psin']

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
		
		self._slider_axes = {'slider_1': axes([0.15, 0.01, 0.5, 0.03],visible=self['visible']['slider_1']),
		}
		self.sliders = {}
		for key in [x for x in self._slider_axes.keys() if self.settings[x]['dimension_type'] is not None]:
			dim = self.reader.dimensions[self.settings[key]['dimension_type']]
			self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index:", 0, len(dim)-1, valinit = self[key]['id'], valstep = 1)
			self.sliders[key].on_changed(self.draw_fig)
		
		if self['visible']['slider_1'] == True:	
			self.fig.subplots_adjust(bottom=0.15)
		elif self['visible']['slider_1'] == False:
			self.fig.subplots_adjust(bottom=0.11)

		self.ch_axes = axes([0.72, 0.01, 0.09, 0.1],frame_on=False)
		self.options = CheckButtons(self.ch_axes, self._options, self['options'])
		self.options.on_clicked(self.draw_fig)
		
		

		ion()
		show()
		self.draw_fig()
	
	def _load_x_axis(self, axis_type):
		if axis_type not in ['beta_prime','alpha']:
			print(f"ERROR: axis_type not found, valid types ['beta_prime','alpha']")
			return
			
		self.settings['x_axis_type'] = axis_type
		
		if axis_type in ['alpha']:
			if not self.data['alpha_axis_ideal']:
				print("ERROR: Alpha not calculated, use calculate_alpha()")
			else:
				self._x_axis_label = r'$\alpha$'
		else:
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
			self._y_axis_label = r'$\hat{s}$'
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys {self.settings['visible'].keys()}")
			return
		if val not in [True,False]:
			val = not self['visible'][key]
		
		self.settings['visible'][key] = val
		
		if key == 'psi_sli':
			self.psi_axes.set_visible(self['visible']['psi_sli'])
		elif key == 'op_box':
			self.ch_axes.set_visible(self['visible']['op_box'])
		elif key == 'suptitle':
			self.fig._suptitle.set_visible(self['visible']['suptitle'])
		elif key == 'title':
			self.ax.legend_.set_visible(self['visible']['title'])
		if self['visible']['op_box'] == False and self['visible']['psi_sli'] == False:
			self.fig.subplots_adjust(bottom=0.11)
	
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
		
	def set_suptitle_fontsize(self, fontsize):
		self.settings['fontsizes']['suptitle'] = fontsize
		self.fig._suptitle.set_fontsize(fontsize)
		
	def set_suptitle(self, title):
		self.settings['suptitle'] = title
		self.fig.suptitle(title, fontsize = self.settings['fontsizes']['suptitle'])
	
	def set_eqbm_style(self, eqbm_style):
		if eqbm_style in self._valid_eqbm_styles:
			self.settings['eqbm_style'] = eqbm_style
			self.draw_fig()
		else:
			print(f"ERROR: eqbm_style not found, valid styles = {self._valid_eqbm_styles}")
	
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
	
	def draw_fig(self, val = None):
		for key, sli in self.sliders.items():
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				self.settings[key]['id'] = sli.val
		
		if 'psin' in self['run']:
			psiN = self['run']['psin']
		else:
			psiN = self.reader.single_parameters.values[0]
		
		x_axis = self.reader.data['ideal_data'][psiN][self['x_axis_type']]
		y_axis = self.reader.data['ideal_data'][psiN][self['y_axis_type']]
		stab = self.reader.data['ideal_data'][psiN]['stabilities']
		x_val = abs(self.reader.data['equilibrium'][psiN][self['x_axis_type']])
		y_val = self.reader.data['equilibrium'][psiN][self['y_axis_type']]
		
		self.ax.cla()
		self.ax.set_facecolor('lightgrey')
		self.ax.set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax.set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		
		psi_line = Line2D([0,1],[0.5,0.5],color='k',label=f"{self.reader.dimensions['psin'].axis_label} = {psiN}",visible = False)
		ideal_line = None
		eqbm_line = None
		s_patch = None
		u_patch = None
		
		status = self.options.get_status()
		if status[0]:
			self.ax.contourf(x_axis, y_axis, stab, [0,0.01,0.99,1,], colors = (self['colours']['stable'],self['colours']['boundary'],self['colours']['unstable']), extend = "both")
		else:
			self.ax.contourf(x_axis, y_axis, stab, [0.01,0.99], colors = (self['colours']['boundary']))
			
		if status[1]:
			self.ax.set_xlim(amin(self.x_axis_ideal),amax(self.x_axis_ideal))
			self.ax.set_ylim(amin(self.y_axis_ideal),amax(self.y_axis_ideal))
		
		if status[2]:
			self.ax.plot(x_val,y_val,'kx')
			if self['eqbm_style'] in ["point",2,"point numless",4]:
				self.ax.annotate("Eqbm",(x_val,y_val),textcoords = "offset points",xytext = (0,7), ha = "center")
			if self['eqbm_style'] in ["split",1,"point",2]:
				self.ax.annotate(f"{x_val:.2f},{y_val:.2f}",(x_val,y_val),textcoords = "offset points",xytext = (0,-13), ha = "center")
			if self['eqbm_style'] in ["title",0]:
				eqbm_line = Line2D([0.5],[0.5],marker='x',color='k',label=f"Equillibrium ({x_val:.2f},{y_val:.2f})",linewidth=0)
			if self['eqbm_style'] in ["split",1,"title numless",3]:
				eqbm_line = Line2D([0.5],[0.5],marker='x',color='k',label="Equillibrium",linewidth=0)
		
		if status[3]:
			ideal_line = Line2D([0,1],[0.5,0.5],color=self['colours']['boundary'],label="Ideal Boundary")
			if status[0]:
				s_patch = pt.Patch(color=self['colours']['stable'], label='Stable')
				u_patch = pt.Patch(color=self['colours']['unstable'], label='Unstable')
				
		handles = [line for line in [psi_line,eqbm_line,ideal_line,s_patch,u_patch] if line is not None]
		self.ax.legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0.5,0.98),loc = "lower center", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax.legend_.set_visible(self['visible']['title'])
		
		self.fig.canvas.draw_idle()
