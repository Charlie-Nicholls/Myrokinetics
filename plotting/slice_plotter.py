from matplotlib.pyplot import *
from matplotlib.widgets import Slider
from numpy import array

default_settings = {"suptitle": None,
		"psi_id": 0,
		"sh_id": 0,
		"limit": None,
		"fontsizes": {"title": 13, "axis": 17,"suptitle": 20},
		"visible": {"psi_sli": True, "sh_sli": True, "eqbm": True, "suptitle": True, "title": True},
		"colours": {"eqbm": 'k', "points": 'k', "line": 'r'}
}

class plot_slice(object):
	def __init__(self, scan = None, settings = {}):
		if scan is None:
			print("ERROR: no scan dictionary given")
			return
		self.info = scan['info']
		self.data = scan['data']
		self.inputs = scan['inputs']
		self.psiNs = self.inputs['psiNs']
		
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
		
		self.open_plot()
		
	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
		
	def save_plot(self, filename = None):
		if filename is None:
			filename = f"Slice_{self['psi_id']}_{self['sh_id']}"
		self.fig.savefig(filename)
		
	def open_plot(self, save = False, filename = None):
		self.fig, self.ax = subplots(figsize=(9,7))
		self.fig.subplots_adjust(bottom=0.15)
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		self.psi_axes = axes([0.25, 0.01, 0.5, 0.03])
		self.psi_slider = Slider(self.psi_axes, '\u03A8\u2099 index:', 0, len(self.psiNs)-1, valinit = self['psi_id'], valstep = 1)
		self.psi_slider.on_changed(self.draw_fig)
		
		self.sh_axes = axes([0.25, 0.05, 0.5, 0.03])
		self.sh_slider = Slider(self.sh_axes, 'shear index:', 0, len(self.data['shear_axis'][0])-1, valinit = self['sh_id'], valstep = 1)
		self.sh_slider.on_changed(self.draw_fig)
		
		if len(self.psiNs) == 1:
			self.set_visible(key = 'psi_sli', val = False)
		if len(self.data['shear_axis'][0]) == 1:
			self.set_visible(key = 'sh_sli', val = False)
		
		ion()
		show()
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
		elif key == 'sh_sli':
			self.sh_axes.set_visible(self['visible']['sh_sli'])
		elif key == 'suptitle':
			self.fig._suptitle.set_visible(self['visible']['suptitle'])
		elif key == 'title':
			self.ax.legend_.set_visible(self['visible']['title'])
		if self['visible']['sh_sli'] == False and self['visible']['psi_sli'] == False:
			self.fig.subplots_adjust(bottom=0.11)
		elif self['visible']['sh_sli'] == False: 
			self.fig.subplots_adjust(bottom=0.13)
	
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
	
	def draw_fig(self, val = None):
		self.ax.cla()
		psi_id = self.psi_slider.val
		sh_id = self.sh_slider.val
		self.settings['psi_id'] = psi_id
		self.settings['sh_id'] = psi_id
		
		self.ax.set_ylabel("Quasilinear Metric",fontsize=self['fontsizes']['axis'])
		self.ax.set_xlabel(f"-\u03B2'",fontsize=self['fontsizes']['axis'])
		
		psi_line = Line2D([0,1],[0.5,0.5],color='k',label=f"\u03A8\u2099 = {self.psiNs[psi_id]}",visible = False)
		shear_line = Line2D([0,1],[0.5,0.5],color='k',label=f"shear = {self.data['shear_axis'][psi_id][sh_id]}",visible = False)
		eqbm_line = None
		
		qls = array(self.data['quasilinear'])[psi_id,:,sh_id]
		
		self.ax.plot(self.data['beta_prime_axis'][psi_id],qls,c=self['colours']['line'])
		self.ax.plot(self.data['beta_prime_axis'][psi_id],qls,'.',c=self['colours']['points'])
		
		if self['visible']['eqbm']:
			limits = self.ax.get_ylim()
			self.ax.vlines(abs(self.data['beta_prime_values'][psi_id]),0,limits[1],self['colours']['eqbm'])
			eqbm_line = Line2D([0.5,0.5],[0,1],c=self['colours']['eqbm'],label = "Equillibrium")
		
		if self['limit']:
			self.ax.set_ylim(0,self['limit'])
		
		handles = [line for line in [psi_line,shear_line,eqbm_line] if line is not None]
		self.ax.legend(ncol = len(handles), handles = handles, bbox_to_anchor= (0.5,0.98),loc = "lower center", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax.legend_.set_visible(self['visible']['title'])
		
		self.fig.canvas.draw_idle()

