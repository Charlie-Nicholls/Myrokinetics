from numpy import transpose, array, amax, amin, isfinite, linspace, where
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap

default_settings = {"suptitle": None,
		"psi_id": 0,
		"ky_id": 0,
		"gr_scale": 100,
		"mf_scale": 100,
		"gr_max": None,
		"mf_max": None,
		"eqbm_style": "title",
		"gr_type": "Normalised",
		"aky": False,
		"options": [False,False,True,True,False],
		"vr_options": [False,False,False,False,False],
		"fontsizes": {"title": 13, "ch_box": 8, "vr_box": 8, "axis": 17,"suptitle": 20},
		"visible": {"psi_sli": True, "ky_sli": True, "gr_sli": True, "mf_sli": True, "op_box": True, "vr_box": True, "suptitle": True, "title": True},
		"cdict": {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
			'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
			'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))},
}

class plot_scan(object):
	def __init__(self, scan = None, verify = None, settings = {}):
		if scan is None:
			print("ERROR: no scan dictionary given")
			return
		self.scan = scan
		self.data = scan['data']
		self.inputs = scan['inputs']
		self.psiNs = self.inputs['psiNs']
		self.verify = verify
		if self.data['growth_rates'] is None:
			print("Error: No Gyrokinetic Data")
			return
		
		self.bpmin = amin(self.data['beta_prime_axis'])
		self.bpmax = amax(self.data['beta_prime_axis'])
		self.shmin = amin(self.data['shear_axis'])
		self.shmax = amax(self.data['shear_axis'])
		
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
		self.cmap = LinearSegmentedColormap('GnRd', self['cdict'])
		
		self._valid_eqbm_styles = ["title",0,"split",1,"point",2,"title numless",3,"point numless",4]
		self._options = ["Show Parities","Global Axis Limits","Global Colorbar","Show Equillibrium","Show Ideal"]
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
		self.fig.subplots_adjust(bottom=0.15)   
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
	
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar_mf = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[1])
		self.cbar_gr = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[0])
		
		self.ch_axes = axes([0.72, 0.01, 0.09, 0.1],frame_on = False)
		self.options = CheckButtons(self.ch_axes, self._options, self['options'])
		self.options.on_clicked(self.draw_fig)
		
		self.vr_axes = axes([0.85, 0.01, 0.09, 0.1],frame_on = False)
		self.vroptions = CheckButtons(self.vr_axes, self._vr_options, self['vr_options'])
		self.vroptions.on_clicked(self.draw_fig)
		
		self.psi_axes = axes([0.15, 0.01, 0.5, 0.03])
		self.slider = Slider(self.psi_axes, '\u03A8\u2099 index:', 0, len(self.psiNs)-1, valinit = self['psi_id'], valstep = 1)
		self.slider.on_changed(self.draw_fig)
		
		self.ky_axes = axes([0.15, 0.05, 0.5, 0.03])
		if self['aky']:
			self.ky_slider = Slider(self.ky_axes, 'ky index:', 0, len(self.inputs['aky_values'])-1, valinit = self['ky_id'], valstep = 1)
			self.ky_slider.on_changed(self.draw_fig)
		else:
			self.set_visible('ky_sli',False)
		
		self.gr_axes = axes([0.93, 0.15, 0.01, 0.73])
		self.gr_slider = Slider(self.gr_axes, 'GR', 0, 100, valinit = self['gr_scale'], valstep = 1, orientation = 'vertical')
		self.gr_slider.on_changed(self.draw_fig)
		
		self.mf_axes = axes([0.97, 0.15, 0.01, 0.73])
		self.mf_slider = Slider(self.mf_axes, 'MF', 0, 100, valinit = self['mf_scale'], valstep = 1, orientation = 'vertical')
		self.mf_slider.on_changed(self.draw_fig)
		
		if len(self.psiNs) == 1:
			self.set_visible(key = 'psi_sli', val = False)
		
		ion()
		self.draw_fig()
		if self['gr_max']:
			self.set_gr_max(self.init_settings['gr_max'])
		if self['mf_max']:
			self.set_mf_max(self.init_settings['mf_max'])
		show()
		
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
		
		if key == 'psi_sli':
			self.psi_axes.set_visible(self['visible']['psi_sli'])
		elif key == 'ky_sli':
			self.ky_axes.set_visible(self['visible']['ky_sli'])
		elif key == 'op_box':
			self.ch_axes.set_visible(self['visible']['op_box'])
		elif key == 'vr_box':
			self.vr_axes.set_visible(self['visible']['vr_box'])
		elif key == 'suptitle':
			self.fig._suptitle.set_visible(self['visible']['suptitle'])
		elif key == 'title':
			self.ax.legend_.set_visible(self['visible']['title'])
		if not self['visible']['ky_sli'] and not self['visible']['psi_sli'] and not self['visible']['op_box']:
			self.fig.subplots_adjust(bottom=0.11)
		elif not self['visible']['ky_sli'] and not self['visible']['op_box']: 
			self.fig.subplots_adjust(bottom=0.13)
			
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
		self.fig.suptitle(title)
	
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
		idx = self.slider.val
		psiN = self.psiNs[idx]
		self.settings['psi'] = idx
		
		x = list(self.data['beta_prime_axis'][idx])
		y = list(self.data['shear_axis'][idx])
		
		beta_prime = abs(self.data['beta_prime_values'][idx])
		shear = self.data['shear_values'][idx]
		
		psi_line = Line2D([0,1],[0.5,0.5],color='k',label=f"\u03A8\u2099 = {psiN}",visible = False)
		ky_line = None
		ideal_line =  None
		eqbm_line = None
		
		if self['aky']:
			ky_idx = self.ky_slider.val
			ky = self.inputs['aky_values'][ky_idx]
			self.settings['ky_id'] = ky_idx
			ky_line = Line2D([0,1],[0.5,0.5],color='k',label=f"ky = {ky}", visible = False)
		else:
			self.settings['ky_id'] = None
			
		
		self.ax[0].cla()
		self.ax[0].set_facecolor('grey')
		self.ax[0].set_ylabel("Shear")
		self.ax[0].set_xlabel("-\u03B2'")
		self.ax[1].cla()
		self.ax[1].set_facecolor('grey')
		self.ax[1].set_ylabel("Shear")
		self.ax[1].set_xlabel("-\u03B2'")
		self.ax[1].set_title("Mode Frequency")
		
		status = self.options.get_status()
		self.settings['options'] = status
		vr_status = self.vroptions.get_status()
		self.settings['vr_options'] = vr_status
		
		if self['aky']:
			z_gr = transpose(array(self.data['growth_rates_all'])[idx,:,:,ky_idx]).tolist()
			z_mf = transpose(array(self.data['mode_frequencies_all'])[idx,:,:,ky_idx]).tolist()
			self.ax[0].set_title(f"Growth Rate")
		else:
			z_gr = transpose(self.data['growth_rates'][idx]).tolist()
			z_mf = transpose(self.data['mode_frequencies'][idx]).tolist()
			if self['gr_type'] == "Normalised":
				self.ax[0].set_title(f"Growth Rate/ky\u00b2")
			else:
				self.ax[0].set_title(f"Growth Rate")
		
		if self['gr_max']:
			grmax = self.gr_slider.val * self['gr_max']/100
		elif status[2]:
			grmax = self.gr_slider.val * abs(amax(array(self.data['growth_rates'])[isfinite(self.data['growth_rates'])]))/100
		else:
			try:
				grmax = self.gr_slider.val * amax(abs(array(z_gr)[isfinite(z_gr)]))/100
				if grmax < 10e-10:
					grmax = 10e-10
			except:
				grmax = 1
		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		self.settings['gr_scale'] = self.gr_slider.val
		self.cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = self.cmap))
		self.ax[0].pcolormesh(x, y, z_gr, cmap = self.cmap, norm=norm_gr)
		
		if self['mf_max']:
			mfmax = self.mf_slider.val * self['mf_max']/100
		elif status[2]:
			mfmax = self.mf_slider.val * abs(amax(array(self.data['mode_frequencies'])[isfinite(self.data['mode_frequencies'])]))/100
		else:
			try:
				mfmax = self.mf_slider.val * amax(abs(array(z_mf)[isfinite(z_mf)]))/100
				if mfmax == 0:
					mfmax = 10e-10
			except:
				mfmax = 1
		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		self.settings['mf_scale'] = self.mf_slider.val
		self.cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		self.ax[1].pcolormesh(x, y, z_mf, norm = norm_mf)
		
		if status[0]:
			if self['aky']:
				parities = array(self.data['parities_all'])[idx,:,:,ky_idx]
			else:
				parities = array(self.data['parities'][idx])
			sym_ids = where(parities == 1)
			anti_ids = where(parities == -1)
			self.ax[0].plot(array(x)[sym_ids[0]], array(y)[sym_ids[1]], '+', color = 'purple', label = 'par')
			self.ax[1].plot(array(x)[sym_ids[0]], array(y)[sym_ids[1]], '+', color = 'purple', label = 'par')
			self.ax[0].plot(array(x)[anti_ids[0]], array(y)[anti_ids[1]], '_', color = 'cyan', label = 'par')
			self.ax[1].plot(array(x)[anti_ids[0]], array(y)[anti_ids[1]], '_', color = 'cyan', label = 'par')
		
		if status[1]:
			self.ax[0].set_ylim(self.shmin,self.shmax)
			self.ax[0].set_xlim(self.bpmin,self.bpmax)
			self.ax[1].set_ylim(self.shmin,self.shmax)
			self.ax[1].set_xlim(self.bpmin,self.bpmax)

		if status[3]:
			self.ax[0].plot(beta_prime,shear,'kx')
			self.ax[1].plot(beta_prime,shear,'kx')
			if self['eqbm_style'] in ["point",2,"point numless",4]:
				self.ax[0].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
				self.ax[1].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			if self['eqbm_style'] in ["split",1,"point",2]:
				self.ax[0].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
				self.ax[1].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
			if self['eqbm_style'] in ["title",0]:
				eqbm_line = Line2D([0.5],[0.5],marker='x',color='k',label=f"Equillibrium ({round(beta_prime,2)},{round(shear,2)})",linewidth=0)
			if self['eqbm_style'] in ["split",1,"title numless",3]:
				eqbm_line = Line2D([0.5],[0.5],marker='x',color='k',label="Equillibrium",linewidth=0)
		
		if status[4]:
			if self.data['ideal_stabilities'] is not None and self.data['ideal_stabilities'][idx] is not None:
				self.ax[0].contourf(self.data['beta_prime_axis_ideal'][idx], self.data['shear_axis_ideal'][idx], self.data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
				self.ax[1].contourf(self.data['beta_prime_axis_ideal'][idx], self.data['shear_axis_ideal'][idx], self.data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
				ideal_line = Line2D([0,1],[0.5,0.5],color='k',label="Ideal Boundary")
			else:
				self.ax[0].text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=self.ax[0].transAxes,color='k')
		
		if vr_status[0]:
			dx = x[1] - x[0]
			dy = y[1] - y[0]
			self.ax[0].set_xticks(array(x)-dx/2, minor=True)
			self.ax[1].set_xticks(array(x)-dx/2, minor=True)
			self.ax[0].set_yticks(array(y)-dy/2, minor=True)
			self.ax[1].set_yticks(array(y)-dy/2, minor=True)
			
			self.ax[0].set_xticks(array(x), minor=False)
			self.ax[1].set_xticks(array(x), minor=False)
			self.ax[0].set_yticks(array(y), minor=False)
			self.ax[1].set_yticks(array(y), minor=False)
			
			self.ax[0].set_xticklabels([])
			self.ax[1].set_xticklabels([])
			self.ax[0].set_yticklabels([])
			self.ax[1].set_yticklabels([])
			xlabels = [str(i) for i in range(len(x))]
			ylabels = [str(i) for i in range(len(y))]
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
			
			self.ax[0].plot(array(x)[un_x], array(y)[un_y], 'c^', label = 'U')
			self.ax[1].plot(array(x)[un_x], array(y)[un_y], 'c^', label = 'U')
			self.ax[0].plot(array(x)[uns_x], array(y)[uns_y], 'bv', label = 'Us')
			self.ax[1].plot(array(x)[uns_x], array(y)[uns_y], 'bv', label = 'Us')
			self.ax[0].plot(array(x)[co_x], array(y)[co_y], 'k>', label = 'C')
			self.ax[1].plot(array(x)[co_x], array(y)[co_y], 'k>', label = 'C')
			self.ax[0].plot(array(x)[cof_x], array(y)[cof_y], 'm<', label = 'Cf')
			self.ax[1].plot(array(x)[cof_x], array(y)[cof_y], 'm<', label = 'Cf')
			self.ax[0].legend(loc = 'center', bbox_to_anchor=(0.5,1.1), ncol = 4)
		
		if vr_status[2]:
			if self['aky']:
				ns_x = [i for p,i,j,k in self.verify['nstep'] if p == idx and k == ky_idx]
				ns_y = [j for p,i,j,k in self.verify['nstep'] if p == idx and k == ky_idx]
			else:
				ns_x = [i for p,i,j,k in self.verify['nstep'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				ns_y = [j for p,i,j,k in self.verify['nstep'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
			self.ax[0].plot(array(x)[ns_x], array(y)[ns_y], 'kX')
			self.ax[1].plot(array(x)[ns_x], array(y)[ns_y], 'kX')
		
		if vr_status[3]:
			for bpid, bp in enumerate(x):
				for shid, sh in enumerate(y):
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
			self.ax[0].plot(array(x)[na_x], array(y)[na_y], 'kX')
			self.ax[1].plot(array(x)[na_x], array(y)[na_y], 'kX')
		
		handles = [line for line in [psi_line,ky_line,ideal_line,eqbm_line] if line is not None]
		self.ax[0].legend(ncol = 4, handles = handles, bbox_to_anchor= (0,1.02),loc = "lower left", fontsize = self['fontsizes']['title'], frameon = False)
		self.ax[0].legend_.set_visible(self['visible']['title'])
		
		self.fig.canvas.draw_idle()
