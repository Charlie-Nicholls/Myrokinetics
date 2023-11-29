from numpy import real, imag, log, polyfit, array, exp, amax
from matplotlib.pyplot import *
from matplotlib.widgets import Slider, CheckButtons, TextBox
from scipy.stats import pearsonr
from copy import deepcopy
from .slider_ax import slider_axes

default_settings = {"suptitle": None,
		"var": 'omega',
		"run": {},
		"normalisation": "highest",
		"options": [True],
		"fontsizes": {"legend": 10,"ch_box": 8,"axis": 11,"title": 13,"suptitle": 20, "verify": 8},
		"visible": {"op_box": True, "suptitle": True, "title": True, "legend": True, "verify": True, 'absolute': True, 'real': True, 'imag': True},
		"colours": {"real": 'r', "imag": 'b', "absolute": 'k', "divider": 'g'},
}

slider_settings = {"slider_1": {"axis": [0.25, 0.01, 0.5, 0.03]},
		"slider_2": {"axis": [0.25, 0.05, 0.5, 0.03]},
		"slider_3": {"axis": [0.92, 0.2, 0.03, 0.5], "orientation": 'vertical'}, 
		"slider_4": {"axis": [0.96, 0.2, 0.03, 0.5], "orientation": 'vertical'},
		"dims": None,
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "slider_5": False, "slider_6": False, "slider_7": False, "slider_8": False, "slider_9": False},
}

class plot_diag(object):
	def __init__(self, reader, settings = {}, sliders = None):
		self.reader = reader
		self.verify = reader.verify
		self.settings = {}
		self.sliders = sliders
		defaults = deepcopy(default_settings)
		for key in settings:
			if key not in defaults:
				print(f"ERROR: {key} not found")
			else:
				self.settings[key] = settings[key]
		for key in defaults:
			if key not in self.settings:
				self.settings[key] = defaults[key]
			elif type(self.settings[key]) == dict:
				for skey in defaults[key]:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
		
		if self['var'] == 0:
			self.settings['var'] = "omega"
		if self['var'] == 1:
			self.settings['var'] = "phi"
		if self['var'] == 2:
			self.settings['var'] = "apar"
		if self['var'] == 3:
			self.settings['var'] = "bpar"
		if self['var'] == 4:
			self.settings['var'] = "epar"
		if self['var'] == 5:
			self.settings['var'] = "phi2"
		if self['var'] == 6:
			self.settings['var'] = "phi2_avg"
		if self['var'] == 7:
			self.settings['var'] = "jacob"
		if self['var'] not in ['omega','phi','apar','bpar','epar','phi2','phi2_avg','jacob']:
			print(f"ERROR: variable name/value {self['var']} not supported. supported: omega/0, phi/1, apar/2, bpar/3, epar/4, phi2/5, phi2_avg/6, jacob/7")
			return
		self.open_plot()

	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
	
	def save_plot(self, filename = None):
		if filename is None:
			filename = f"{self['var']}_{self['slider_1']['id']}_{self['slider_2']['id']}_{self['slider_3']['id']}_{self['slider_4']['id']}"
		self.fig.savefig(filename)
	
	def open_plot(self):
		if self['var'] == 'omega':
			self.fig, [self.ax, self.ax2] = subplots(2,1,figsize=(10,10))
		else:
			self.fig, self.ax = subplots(figsize=(8.8,5.8))
		
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		if self['var'] == "phi2":
			chaxes = axes([0.8, 0.01, 0.09, 0.1],frame_on = False)
			self.options = CheckButtons(chaxes, ["Show Fit"],self['options'])
			self.options.on_clicked(self.draw_fig)
		
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])

		self.dims = self.reader.inputs.dim_order		
		for dim in [x for x in self.dims if x not in self.settings['run']]:
			self.settings['run'][dim] = self.reader.dimensions[dim].values[0]
		
		if self.sliders in [None,'seperate']:
			slider_defaults = deepcopy(slider_settings)
			slider_defaults['dims'] = self.dims
			self.fig.subplots_adjust(bottom=0.15)
			if self.sliders == 'seperate':
				sl_ax = None
			else:
				sl_ax = self.ax
			self.sliders = slider_axes(reader=self.reader,settings=slider_defaults,ax=sl_ax)
		self.sliders.add_plot(self)
				
		self.set_normalisation(self['normalisation'])
		
		ion()
		show()
		self.draw_fig()
	
	def set_slider(self, num = None, key = None, dimension_type = None):
		self.sliders.set_slider(num = num, key = key, dimension_type = dimension_type)
	
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
			elif key == 'suptitle':
				self.fig._suptitle.set_visible(self['visible']['suptitle'])
			elif key == 'title':
				self.ax.legend_.set_visible(self['visible']['title'])
		self.draw_fig()
	
	def set_options_fontsize(self, fontsize):
		self.settings['fontsizes']['ch_box'] = fontsize
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
	
	def set_normalisation(self, norm):
		if norm.lower() not in ['none','highest','phi','apar','bpar','epar']:
			print("ERROR: invalid normalisation, allowed: ['none','highest','phi','apar','bpar','epar']")
			return
			
		self.settings['normalisation'] = norm.lower()
		self.draw_fig()
		
	def draw_fig(self, val = None):
		for key, sli in self.sliders.sliders.items():
			dim = self.sliders.settings[key]['dimension_type']
			if dim in self.dims:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
					
		self.ax.cla()
		run=self['run']
		run_id = self.reader.get_run_id(run)
		title = "".join([f"{self.reader.dimensions[x].axis_label}: {self.reader(x,run):.2g} | " for x in self.reader.inputs.dim_order])[:-3]
		self.ax.set_title(title,fontsize=self['fontsizes']['title'])
		
		if self['var'] == 'omega':
			if self.reader('omega',run) is None or self.reader('t',run) is None:
				print("ERROR: data not found")
				return
			else:
				omega = self.reader('omega',run)
				t = self.reader('t',run)
			
			self.ax2.cla()
			
			self.ax2.plot(t,real(omega),color=self['colours']['real'],label="mode frequency")
			self.ax.plot(t,imag(omega),color=self['colours']['imag'],label="growth rate")
			
			self.ax.text(0.01,0.99,f"GR: {self.reader('growth_rate',run):+.2e}\nOmega[-1]: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax.transAxes,fontsize=self['fontsizes']['legend'])
			self.ax2.text(0.01,0.99,f"MF: {real(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax2.transAxes,fontsize=self['fontsizes']['legend'])
			self.ax.set_ylabel("Growth Rate",fontsize=self['fontsizes']['axis'])
			self.ax.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			self.ax2.set_ylabel("Mode_Frequency",fontsize=self['fontsizes']['axis'])
			self.ax2.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			self.ax.legend(loc=1)
			self.ax2.legend(loc=1)
			
		elif self['var'] in ['phi','apar','bpar','epar','jacob']:
			if self.reader(self['var'],run) is None or self.reader('theta',run) is None:
				print("ERROR: data not found")
				return
			else:
				field = self.reader(self['var'],run)
				theta = self.reader('theta',run)
				
			if self['var'] == 'phi':
				ylabel = "$\phi$"
			elif self['var'] == 'apar':
				ylabel = "$A_{\parallel}$"
			elif self['var'] == 'bpar':
				ylabel = "$B_{\parallel}$"
			elif self['var'] == 'epar':
				ylabel = "$E_{\parallel}$"
			elif self['var'] == 'jacob':
				ylabel = "Jacobian"
			
			if self['var'] in ['phi','apar','bpar','epar']:
				if self['normalisation'] == 'highest':
					norm = max([amax([abs(i) for i in self.reader(x,run)]) for x in [y for y in ['phi','apar','bpar','epar'] if (self.reader(y,run) is not None)]])
					ylabel += " / max($\phi,A_{\parallel},B_{\parallel},E_{\parallel}$)"
				elif self['normalisation'] in ['phi','apar','bpar','epar']:
					norm = amax([abs(i) for i in self.reader(self['normalisation'],run)])
					ylabel += f" / max({self['normalisation']})"
				else:
					norm = 1
				field_norm = array(field)/norm
				if self['visible']['real']:
					self.ax.plot(theta,real(field_norm),color=self['colours']['real'],label="real")
				if self['visible']['imag']:
					self.ax.plot(theta,imag(field_norm),color=self['colours']['imag'],label="imaginary")
				if self['visible']['absolute']:
					self.ax.plot(theta,[abs(x) for x in field_norm],color=self['colours']['absolute'],linestyle='--',label="absolute")
				self.ax.legend(loc=0)
			else:
				self.ax.plot(theta,field,'r')
				
			self.ax.set_xlabel("Ballooning Angle",fontsize=self['fontsizes']['axis'])
			self.ax.set_ylabel(ylabel,fontsize=self['fontsizes']['axis'])
			
		elif self['var'] in ['phi2','phi2_avg']:
			if self.reader(self['var'],run) is None or self.reader('t',run) is None:
				print("ERROR: data not found")
				return
			else:
				phi2 = [x for x in self.reader(self['var'],run) if x !=0]
				t = [x for xi, x in enumerate(self.reader('t',run)) if self.reader(self['var'],run)[xi] != 0]
			if len(phi2) > 0:
				self.ax.plot(t,phi2,'k')
				if max(phi2) > 1e267:
					self.ax.set_ylim(min(phi2)/100, 1e267) #Display error occurs on log scale above 1e267
				self.ax.set_ylabel("$\phi^{2}$",fontsize=self['fontsizes']['axis'])
				self.ax.set_yscale('log')
				self.ax.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
				
				if self['var'] == 'phi2' and self.options.get_status()[0] and self.verify is not None:
					if run_id in self.verify.nts:
						nt = self.verify.nts[run_id]
					else:
						nt = None
					if nt is not None:
						fit = polyfit(t[-nt:],log(phi2[-nt:]),1)
						fitgr = fit[0]/2
						pr = pearsonr(t[-nt:],log(phi2[-nt:]))
						gradgr = log(phi2[-1]/phi2[0])/(t[-1]-t[0])/2
						omega = self.reader('omega',run)
						
						self.ax.plot(t[-nt:],exp(array(t[-nt:])*fit[0] + fit[1]),'r',label=f"GR = {fitgr:+.2e}\nR = {pr[0]:.4f}")
						self.ax.plot([t[0],t[-1]],[phi2[0],phi2[-1]],'b',label=f"AVG GR = {gradgr:+.2e}")
						self.ax.text(0.01,0.99,f"GR: {self.reader('growth_rate',run):+.2e}\nOmega[-1]: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax.transAxes,fontsize=self['fontsizes']['legend'])
						self.ax.set_xlabel(f"Time ({len(t)} steps) | nt = {nt}")
						self.ax.legend(loc=1,fontsize=self['fontsizes']['legend'])
		
		if self.verify is not None and self['visible']['verify']:
			bad = []
			if run_id in self.verify['nstep']:
				bad.append('nstep')
			if run_id in self.verify['nan']:
				bad.append('nan')
			if run_id in self.verify['unconv']:
				bad.append('unconverged')
			if run_id in self.verify['negative']:
				bad.append('negative')
			if run_id in self.verify['order']:
				bad.append('order')
			if self['var'] == 'phi' and run_id in self.verify['phi']:
				bad.append('phi')
			if self['var'] == 'apar' and run_id in self.verify['apar']:
				bad.append('apar')
			if self['var'] == 'bpar' and run_id in self.verify['apar']:
				bad.append('bpar')
			if bad:
				self.ax.text(0.01,0.01,f"BAD RUN: {str(bad)[1:-1]}",ha='left',va='bottom',transform=self.ax.transAxes,color='r')
		self.fig.canvas.draw_idle()
		return
