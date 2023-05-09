from numpy import transpose, array, amax, amin, isfinite, linspace, where
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap

class plot_scan(object):
	def __init__(self, scan = None, verify = None, aky = False, init = [0,0], gr_type = "Normalised"):
		if scan is None:
			print("ERROR: no scan dictionary given")
			return
		self.scan = scan
		self.data = scan['data']
		self.inputs = scan['inputs']
		self.psiNs = self.inputs['psiNs']
		self.verify = verify
		self.aky = aky
		self.init = init
		if self.data['growth_rates'] is None:
			print("Error: No Gyrokinetic Data")
			return
		
		self.bpmin = amin(self.data['beta_prime_axis'])
		self.bpmax = amax(self.data['beta_prime_axis'])
		self.shmin = amin(self.data['shear_axis'])
		self.shmax = amax(self.data['shear_axis'])
		cdict = {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
			'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
			'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))}
		self.cmap = LinearSegmentedColormap('GnRd', cdict)
		
		self.open_plot()
		
	def save_plot(self, filename = None):
		self.open_plot(save = True, filename = None)
		
	def open_plot(self, save = False, filename = "Scan"):
		self.fig, self.ax = subplots(1,2,figsize=(14.6,7))
		self.fig.subplots_adjust(bottom=0.15)   
		self.fig.suptitle(self.scan['info']['run_name'])
	
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar_mf = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[1])
		self.cbar_gr = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax[0])
		
		try:
			op_init = self.options.get_status()
		except:
			op_init = [False,False,False,True,False]
		chaxes = axes([0.72, 0.01, 0.09, 0.1],frame_on = False)
		self.options = CheckButtons(chaxes, ["Show Parities","Global Axis Limits","Global Colorbar","Show Equillibrium","Show Ideal"], op_init)
		self.options.on_clicked(self.draw_fig)
		
		try:
			vr_init = self.vroptions.get_status()
		except:
			vr_init = [False,False,False,False,False]
		vraxes = axes([0.85, 0.01, 0.09, 0.1],frame_on = False)
		self.vroptions = CheckButtons(vraxes, ["Show ID","Show Convergence","Show Bad nstep","Show bad fields","Show Omega -> nan"], vr_init)
		self.vroptions.on_clicked(self.draw_fig)
		
		try:
			sl_init = self.slider.val
		except:
			sl_init = self.init[0]
		slaxes = axes([0.15, 0.01, 0.5, 0.03])
		self.slider = Slider(slaxes, 'psiN index:', 0, len(self.psiNs)-1, valinit = sl_init, valstep = 1)
		self.slider.on_changed(self.draw_fig)
		
		if self.aky:
			try:
				ky_init = self.ky_slider.val
			except:
				ky_init = self.init[1]
			kyaxes = axes([0.15, 0.05, 0.5, 0.03])
			self.ky_slider = Slider(kyaxes, 'ky index:', 0, len(self.inputs['aky_values'])-1, valinit = ky_init, valstep = 1)
			self.ky_slider.on_changed(self.draw_fig)
		
		try:
			gr_init = self.gr_slider.val
		except:
			gr_init = 100
		gr_axes = axes([0.93, 0.15, 0.01, 0.73])
		self.gr_slider = Slider(gr_axes, 'GR', 0, 100, valinit = 100, valstep = 1, orientation = 'vertical')
		self.gr_slider.on_changed(self.draw_fig)
		
		try:
			mf_init = self.mf_slider.val
		except:
			mf_init = 100
		mf_axes = axes([0.97, 0.15, 0.01, 0.73])
		self.mf_slider = Slider(mf_axes, 'MF', 0, 100, valinit = 100, valstep = 1, orientation = 'vertical')
		self.mf_slider.on_changed(self.draw_fig)
		
		if save:
			if filename is None and aky:
				filename = f"Scan_{self.slider.val}"
			elif filename is None:
				filename = f"Scan_{self.slider.val}_{self.ky_slider.val}"
			self.draw_fig()
			self.fig.savefig(filename)
		else:
			ion()
			show()
			self.draw_fig()
	
	def draw_fig(self, val = None):
		if len(self.psiNs) > 1:
			idx = self.slider.val
			psiN = self.psiNs[idx]
			
		else:
			idx = 0
			psiN = self.psiNs[0]
		
		x = list(self.data['beta_prime_axis'][idx])
		y = list(self.data['shear_axis'][idx])
		
		beta_prime = abs(self.data['beta_prime_values'][idx])
		shear = self.data['shear_values'][idx]
		
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
		vr_status = self.vroptions.get_status()
		
		if self.aky:
			ky_idx = self.ky_slider.val
			ky = self.inputs['aky_values'][ky_idx]
			z_gr = transpose(array(self.data['growth_rates_all'])[idx,:,:,ky_idx]).tolist()
			z_mf = transpose(array(self.data['mode_frequencies_all'])[idx,:,:,ky_idx]).tolist()
			if gr_type == "Normalised":
				self.ax[0].set_title(f"Growth Rate/ky\u00b2 | PsiN: {psiN} | ky: {ky}")
			else:
				self.ax[0].set_title(f"Growth Rate | PsiN: {psiN} | ky: {ky}")
		else:
			z_gr = transpose(self.data['growth_rates'][idx]).tolist()
			z_mf = transpose(self.data['mode_frequencies'][idx]).tolist()
			if gr_type == "Normalised":
				self.ax[0].set_title(f"Growth Rate/ky\u00b2 | PsiN: {psiN}")
			else:
				self.ax[0].set_title(f"Growth Rate | PsiN: {psiN}")
		
		if status[2]:
			mfmax = self.mf_slider.val * abs(amax(array(self.data['mode_frequencies'])[isfinite(self.data['mode_frequencies'])]))/100
			grmax = self.gr_slider.val * abs(amax(array(self.data['growth_rates'])[isfinite(self.data['growth_rates'])]))/100
		else:
			try:
				grmax = self.gr_slider.val * amax(abs(array(z_gr)[isfinite(z_gr)]))/100
				if grmax < 10e-10:
					grmax = 10e-10
			except:
				grmax = 1

			try:
				mfmax = self.mf_slider.val * amax(abs(array(z_mf)[isfinite(z_mf)]))/100
				if mfmax == 0:
					mfmax = 10e-10
			except:
				mfmax = 1
		
		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		self.cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		
		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		self.cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = self.cmap))
		self.ax[0].pcolormesh(x, y, z_gr, cmap = self.cmap, norm=norm_gr)
		self.ax[1].pcolormesh(x, y, z_mf, norm = norm_mf)
		
		if status[0]:
			if self.aky:
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
			self.ax[0].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			self.ax[0].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
			self.ax[1].plot(beta_prime,shear,'kx')
			self.ax[1].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			self.ax[1].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
		
		if status[4]:
			if self.data['ideal_stabilities'] is not None and self.data['ideal_stabilities'][idx] is not None:
				self.ax[0].contourf(self.data['beta_prime_axis_ideal'][idx], self.data['shear_axis_ideal'][idx], self.data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
				self.ax[1].contourf(self.data['beta_prime_axis_ideal'][idx], self.data['shear_axis_ideal'][idx], self.data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
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
			if self.aky:
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
			if self.aky:
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
			if self.aky:
				na_x = [i for p,i,j,k in self.verify['nan'] if p == idx and k == ky_idx]
				na_y = [j for p,i,j,k in self.verify['nan'] if p == idx and k == ky_idx]
			else:
				na_x = [i for p,i,j,k in self.verify['nan'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
				na_y = [j for p,i,j,k in self.verify['nan'] if p == idx and k == self.inputs['aky_values'].index(self.data['akys'][p][i][j])]
			self.ax[0].plot(array(x)[na_x], array(y)[na_y], 'kX')
			self.ax[1].plot(array(x)[na_x], array(y)[na_y], 'kX')

		self.fig.canvas.draw_idle()
