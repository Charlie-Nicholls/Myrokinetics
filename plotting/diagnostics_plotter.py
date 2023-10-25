from numpy import real, imag, log, polyfit, array, exp, amax
from matplotlib.pyplot import *
from matplotlib.widgets import Slider, CheckButtons, TextBox
from scipy.stats import pearsonr
from copy import deepcopy

default_settings = {"suptitle": None,
		"var": 'omega',
		"slider_1": {"dimension_type": None, "id": 0},
		"slider_2": {"dimension_type": None, "id": 0},
		"slider_3": {"dimension_type": None, "id": 0},
		"slider_4": {"dimension_type": None, "id": 0},
		"run": {},
		"options": [True],
		"fontsizes": {"legend": 10,"ch_box": 8,"axis": 11,"title": 13,"suptitle": 20, "verify": 8},
		"visible": {"slider_1": True, "slider_2": True, "slider_3": True, "slider_4": True, "op_box": True, "suptitle": True, "title": True, "legend": True, "verify": True, 'absolute': True, 'real': True, 'imag': True},
}

class plot_diag(object):
	def __init__(self, reader, settings = {}):
		self.reader = reader
		self.verify = reader.verify
		
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
			elif type(self.settings[key]) == dict:
				for skey in defaults[key]:
					if skey not in self.settings[key]:
						self.settings[key][skey] = defaults[key][skey]
						self.init_settings[key][skey] = defaults[key][skey]
		
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
			self.settings['var'] = "jacob"
		if self['var'] not in ['omega','phi','apar','bpar','epar','phi2','jacob']:
			print(f"ERROR: variable name/value {self['var']} not supported. supported: omega/0, phi/1, apar/2, bpar/3, epar/4, phi2/5, jacob/6")
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
		
		self._slider_axes = {'slider_1': axes([0.25, 0.01, 0.5, 0.03],visible=self['visible']['slider_1']),
		'slider_2': axes([0.25, 0.05, 0.5, 0.03],visible=self['visible']['slider_2']),
		'slider_3': axes([0.92, 0.2, 0.03, 0.5],visible=self['visible']['slider_3']),
		'slider_4': axes([0.96, 0.2, 0.03, 0.5],visible=self['visible']['slider_4']),
		}
		self.sliders = {}
		for key in [x for x in self._slider_axes.keys() if self.settings[x]['dimension_type'] is not None]:
			dim = self.reader.dimensions[self.settings[key]['dimension_type']]
			if key in ['slider_1','slider_2']:
				orient = 'horizontal'
			elif key in ['slider_3','slider_4']:
				orient = 'vertical'
			self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index", 0, len(dim)-1, valinit = self[key]['id'], valstep = 1, orientation = orient)
			self.sliders[key].on_changed(self.draw_fig)
			self.set_visible(key,val=self['visible'][key])
			if orient == 'vertical':
				self.sliders[key].label.set_rotation(90)
		
		ion()
		show()
		self.draw_fig()
	
	def set_slider(self, slider_num, dimension_type, visible = True):
		if dimension_type not in self.dims:
			print(f"ERROR: invalid dimension type, valid: {self.dims}")
			return
		key = f"slider_{slider_num}"
		self.settings[key]['dimension_type'] = dimension_type
		self.settings[key]['id'] = 0
		dim = self.reader.dimensions[dimension_type]
		self.sliders[key].ax.remove()
		self._slider_axes[key] = axes(slider_axes[key])
		if key in ['slider_1','slider_2']:
			orient = 'horizontal'
		elif key in ['slider_3','slider_4']:
			orient = 'vertical'
		self.sliders[key] = Slider(self._slider_axes[key], f"{dim.axis_label} index", 0, len(dim)-1, valinit = 0, valstep = 1, orientation = orient)
		self.sliders[key].on_changed(self.draw_fig)
		if orient == 'vertical':
			self.sliders[key].label.set_rotation(90)
		self.set_visible(key,val=visible)
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys: {self.settings['visible'].keys()}")
			return
		if val not in [True,False]:
			val = not self['visible'][key]
		
		self.settings['visible'][key] = val
		
		if 'slider_' in key:
			self._slider_axes[key].set_visible(self['visible'][key])
			
			if self['visible']['slider_1'] == False and self['visible']['slider_2'] == False:
				self.fig.subplots_adjust(bottom=0.11)
			elif self['visible']['slider_2'] == False: 
				self.fig.subplots_adjust(bottom=0.13)
			else:
				self.fig.subplots_adjust(bottom=0.15)
			
		elif key == 'op_box':
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
		
	def draw_fig(self, val = None):
		self.ax.cla()
		for key, sli in self.sliders.items():
			dim = self[key]['dimension_type']
			if dim is not None:
				self.settings['run'][dim] = self.reader.dimensions[dim].values[sli.val]
				self.settings[key]['id'] = sli.val		
		
		run_id = self.reader.get_run_id(run=self['run'])
		data = self.reader.data['gyro'][run_id]
		
		title = "".join([f"{self.reader.dimensions[x].axis_label}: {data[x]:.2g} | " for x in self.reader.inputs.dim_order])[:-3]
		self.ax.set_title(title,fontsize=self['fontsizes']['title'])
		
		if self['var'] == 'omega':
			if data['omega'] is None or data['t'] is None:
				print("ERROR: data not found")
			else:
				omega = data['omega']
				t = data['t']
			
			self.ax2.cla()
			
			self.ax2.plot(t,real(omega),'r',label="mode frequency")
			self.ax.plot(t,imag(omega),'b',label="growth rate")
			
			self.ax.text(0.01,0.99,f"GR: {data['growth_rate']:+.2e}\nOmega[-1]: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax.transAxes,fontsize=self['fontsizes']['legend'])
			self.ax2.text(0.01,0.99,f"MF: {real(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax2.transAxes,fontsize=self['fontsizes']['legend'])
			self.ax.set_ylabel("Growth Rate",fontsize=self['fontsizes']['axis'])
			self.ax.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			self.ax2.set_ylabel("Mode_Frequency",fontsize=self['fontsizes']['axis'])
			self.ax2.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			self.ax.legend(loc=1)
			self.ax2.legend(loc=1)
			
		elif self['var'] in ['phi','apar','bpar','epar','jacob']:
			if data[self['var']] is None or data['theta'] is None:
				print("ERROR: data not found")
			else:
				field = data[self['var']]
				theta = data['theta']
			
			if self['var'] in ['phi','apar','bpar','epar']:
				norm = max([amax([abs(i) for i in data[x]]) for x in ['phi','apar','bpar','epar'] if (x in data and data[x] is not None)])
				field_norm = field/norm
				if self['visible']['real']:
					self.ax.plot(theta,real(field_norm),'r',label="real")
				if self['visible']['imag']:
					self.ax.plot(theta,imag(field_norm),'b',label="imaginary")
				if self['visible']['absolute']:
					self.ax.plot(theta,[abs(x) for x in field_norm],'k--',label="absolute")
				self.ax.legend(loc=0)
			else:
				self.ax.plot(theta,field,'r')
			self.ax.set_xlabel("Ballooning Angle",fontsize=self['fontsizes']['axis'])
			
			if self['var'] == 'phi':
				ylabel = "Electric Potential"
			elif self['var'] == 'apar':
				ylabel = "Parallel Mangetic Potential"
			elif self['var'] == 'bpar':
				ylabel = "Parallel Mangetic Field"
			elif self['var'] == 'epar':
				ylabel = "Parallel Electric Field"
			elif self['var'] == 'jacob':
				ylabel = "Jacobian"
			self.ax.set_ylabel(ylabel,fontsize=self['fontsizes']['axis'])
			
		elif self['var'] == 'phi2':
			if data['phi2'] is None or data['t'] is None:
				print("ERROR: data not found")
			else:
				phi2 = data['phi2']
				t = data['t']
			self.ax.plot(t,phi2,'k')
			if max(phi2) > 1e267:
				self.ax.set_ylim(min(phi2)/100, 1e267) #Display error occurs on log scale above 1e267
			self.ax.set_ylabel("$\phi^{2}$",fontsize=self['fontsizes']['axis'])
			self.ax.set_yscale('log')
			self.ax.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			
			if self.options.get_status()[0] and self.verify is not None:
				if run_id in self.verify.nts:
					nt = self.verify.nts[run_id]
				else:
					nt = None
				if nt is not None:
					fit = polyfit(t[-nt:],log(phi2[-nt:]),1)
					fitgr = fit[0]/2
					pr = pearsonr(t[-nt:],log(phi2[-nt:]))
					gradgr = log(phi2[-1]/phi2[0])/(t[-1]-t[0])/2
					omega = data['omega']
					
					self.ax.plot(t[-nt:],exp(array(t[-nt:])*fit[0] + fit[1]),'r',label=f"GR = {fitgr:+.2e}\nR = {pr[0]:.4f}")
					self.ax.plot([t[0],t[-1]],[phi2[0],phi2[-1]],'b',label=f"AVG GR = {gradgr:+.2e}")
					self.ax.text(0.01,0.99,f"GR: {data['growth_rate']:+.2e}\nOmega[-1]: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax.transAxes,fontsize=self['fontsizes']['legend'])
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
