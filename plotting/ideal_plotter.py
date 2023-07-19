from numpy import *
from matplotlib.pyplot import *
import matplotlib.patches as pt
from matplotlib.widgets import Slider, CheckButtons
from time import sleep

default_settings = {"suptitle": None,
		"psi_id": 0,
		"eqbm_style": "split",
		"x_axis_type": 0,
		"y_axis_type": 0,
		"options": [True,False,True,False],
		"fontsizes": {"title": 13, "ch_box": 8,"axis": 17,"suptitle": 20},
		"visible": {"psi_sli": True, "op_box": True, "suptitle": True, "title": True},
		"colours": {"unstable": 'r', "stable": 'g', "boundary": 'k'}
}

class plot_ideal(object):
	def __init__(self, data = None, inputs = None, settings = {}):
		if data is None:
			print("ERROR: data not given")
			return
		if inputs is None:
			print("ERROR: input not given")
			return
			
		self.data = data
		self.inputs = inputs
		self.psiNs = self.inputs['psiNs']
		if self.data['ideal_stabilities'] is None:
			print("Error: No ideal_ball data")
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		self._options = ["Colour","Global Axis Limits","Show Equillibrium","Show Legend"]
		
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
			elif type(self.settings[key]) == dict:
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
	
	def open_plot(self, save = False, filename = None):
		self.fig, self.ax = subplots(figsize=(10,7))
		self.fig.subplots_adjust(bottom=0.15)
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		self._load_x_axis(self['x_axis_type'])
		self._load_y_axis(self['y_axis_type'])
		
		self.psi_axes = axes([0.15, 0.01, 0.5, 0.03])
		self.slider = Slider(self.psi_axes, '\u03A8\u2099 index:', 0, len(self.psiNs)-1, valinit = self['psi_id'], valstep = 1)
		self.slider.on_changed(self.draw_fig)
		
		self.ch_axes = axes([0.72, 0.01, 0.09, 0.1],frame_on=False)
		self.options = CheckButtons(self.ch_axes, self._options, self['options'])
		self.options.on_clicked(self.draw_fig)
		
		if len(self.psiNs) == 1:
			self.set_visible(key = 'psi_sli', val = False)

		ion()
		show()
		self.draw_fig()
	
	def _load_x_axis(self, axis_type):
		if axis_type not in [0,'beta_prime',1,'alpha']:
			print(f"ERROR: axis_type not found, valid types [0,'beta_prime',1,'alpha']")
			return
			
		self.settings['x_axis_type'] = axis_type
		
		if axis_type in [1,'alpha']:
			if not self.data['alpha_axis_ideal']:
				print("ERROR: Alpha not calculated, use calculate_alpha()")
			else:
				self.x_axis_ideal = self.data['alpha_axis_ideal']
				self.x_values = self.data['alpha_values']
				self._x_axis_label = "\u03B1"
		else:
			self.x_axis_ideal = self.data['beta_prime_axis_ideal']
			self.x_values = self.data['beta_prime_values']
			self._x_axis_label = "-\u03B2'"
			
	def set_x_axis_type(self, axis_type):
		self._load_x_axis(axis_type)
		self.draw_fig()
	
	def _load_y_axis(self, axis_type):
		if axis_type not in [0,'shear',1,'current']:
			print(f"ERRORL axis_type not found, valid types [0,'shear',1,'current']")
			return
			
		self.settings['y_axis_type'] = axis_type
		
		if axis_type in [1,'current']:
			print("Not yet implimented")
		else:
			self.y_axis_ideal = self.data['shear_axis_ideal']
			self.y_values = self.data['shear_values']
			self._y_axis_label = "Shear"
			
	def set_y_axis_type(self, axis_type):
		self._load_y_axis(axis_type)
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys {self.settings.keys()}")
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
		idx = self.slider.val
		psiN = self.psiNs[idx]
		self.settings['psi'] = idx
		
		x = self.x_axis_ideal[idx]
		y = self.y_axis_ideal[idx]
		stab = self.data['ideal_stabilities'][idx]
		x_val = abs(self.x_values[idx])
		y_val = self.y_values[idx]
		
		self.ax.cla()
		self.ax.set_facecolor('lightgrey')
		self.ax.set_ylabel(self._y_axis_label,fontsize=self['fontsizes']['axis'])
		self.ax.set_xlabel(self._x_axis_label,fontsize=self['fontsizes']['axis'])
		
		psi_line = Line2D([0,1],[0.5,0.5],color='k',label=f"\u03A8\u2099 = {psiN}",visible = False)
		ideal_line = None
		eqbm_line = None
		s_patch = None
		u_patch = None
		
		status = self.options.get_status()
		if status[0]:
			self.ax.contourf(x, y, stab, [0,0.01,0.99,1,], colors = (self['colours']['stable'],self['colours']['boundary'],self['colours']['unstable']), extend = "both")
		else:
			self.ax.contourf(x, y, stab, [0.01,0.99], colors = (self['colours']['boundary']))
			
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
