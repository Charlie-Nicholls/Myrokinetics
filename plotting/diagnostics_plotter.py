from numpy import real, imag, log, polyfit, array, exp
from matplotlib.pyplot import *
from matplotlib.widgets import Slider, CheckButtons, TextBox
from scipy.stats import pearsonr

default_settings = {"suptitle": None,
		"psi_id": 0,
		"bp_id": 0,
		"sh_id": 0,
		"ky_id": 0,
		"var": 0,
		"aky": True,
		"absolute": False,
		"options": [True],
		"fontsizes": {"legend": 10,"ch_box": 8,"axis": 11,"title": 13,"suptitle": 20, "verify": 8},
		"visible": {"psi_sli": True, "bp_sli": True, "sh_sli": True, "ky_sli": True, "op_box": True, "suptitle": True, "title": True, "legend": True, "verify": True},
}

class plot_diag(object):
	def __init__(self, data = None, inputs = None, info = None, verify = None, settings = {}):
		if data is None:
			print("ERROR: data not given")
			return
		if inputs is None:
			print("ERROR: input not given")
			return
			
		self.data = data
		self.inputs = inputs
		self.psiNs = self.inputs['psiNs']
		self.verify = verify
		self.info = info
		
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
		
		if self['var'] == "omega":
			self.settings['var'] = 0
		if self['var'] == "phi":
			self.settings['var'] = 1
		if self['var'] == "apar":
			self.settings['var'] = 2
		if self['var'] == "bpar":
			self.settings['var'] = 3
		if self['var'] == "phi2":
			self.settings['var'] = 4
		if self['var'] not in [0,1,2,3,4]:
			print(f"ERROR: variable name/value {self['var']} not supported. supported: omega/0, phi/1, apar/2, bpar/3, phi2/4")
			return
			
		self.open_plot()


	def __getitem__(self, key):
		if key in self.settings:
			return self.settings[key]
		else:
			print(f"ERROR: {key} not found")
	
	def save_plot(self, filename = None):
		if filename is None:
			filename = f"{self['var']}_{self['psi_id']}_{self['bp_id']}_{self['sh_id']}_{self['ky_id']}"
		self.fig.savefig(filename)
	
	def open_plot(self):
		if self['var'] == 0:
			self.fig, [self.ax, self.ax2] = subplots(2,1,figsize=(10,10))
		else:
			self.fig, self.ax = subplots(figsize=(8.8,5.8))
		
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])
		
		if self['var'] == 4:
			chaxes = axes([0.8, 0.01, 0.09, 0.1],frame_on = False)
			self.options = CheckButtons(chaxes, ["Show Fit"],self['options'])
			self.options.on_clicked(self.draw_fig)
		
		self.fig.subplots_adjust(bottom=0.15)
		if self['suptitle']:
			self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes']['suptitle'],visible=self['visible']['suptitle'])

		self.psi_axes = axes([0.25, 0.01, 0.5, 0.03])
		self.bp_axes = axes([0.25, 0.05, 0.5, 0.03])
		self.sh_axes = axes([0.02, 0.25, 0.021, 0.5])

		self.psi_slider = Slider(self.psi_axes, '\u03A8\u2099 index:', 0, len(self.psiNs)-1, valinit = self['psi_id'], valstep = 1)
		self.psi_slider.on_changed(self.draw_fig)
		
		self.bp_slider = Slider(self.bp_axes, "-\u03B2' index:", 0, len(self.data['beta_prime_axis'][0])-1, valinit = self['bp_id'], valstep = 1)
		self.bp_slider.on_changed(self.draw_fig)
		
		self.sh_slider = Slider(self.sh_axes, 'shear\nindex:', 0, len(self.data['shear_axis'][0])-1, valinit = self['sh_id'], valstep = 1,orientation = 'vertical')
		self.sh_slider.on_changed(self.draw_fig)
		
		
		self.ky_axes = axes([0.93, 0.25, 0.021, 0.5])
		self.ky_slider = Slider(self.ky_axes, 'ky index:', 0, len(self.inputs['aky_values'])-1, valinit = self['ky_id'], valstep = 1,orientation = 'vertical')
		self.ky_slider.on_changed(self.draw_fig)
		
		if not self['aky']:
			self.set_visible('ky_sli',False)
		
		ion()
		show()
		self.draw_fig()
	
	def set_visible(self, key, val = None):
		if key not in self['visible']:
			print(f"ERROR: key not found, valid keys: {self.settings['visible'].keys()}")
			return
		if val not in [True,False]:
			val = not self['visible'][key]
		
		self.settings['visible'][key] = val
		
		if key == 'psi_sli':
			self.psi_axes.set_visible(self['visible']['psi_sli'])
		elif key == 'bp_sli':
			self.bp_axes.set_visible(self['visible']['bp_sli'])
		elif key == 'sh_sli':
			self.sh_axes.set_visible(self['visible']['sh_sli'])
		elif key == 'ky_sli':
			self.ky_axes.set_visible(self['visible']['ky_sli'])
		elif key == 'op_box':
			self.ch_axes.set_visible(self['visible']['op_box'])
		elif key == 'suptitle':
			self.fig._suptitle.set_visible(self['visible']['suptitle'])
		elif key == 'title':
			self.ax.legend_.set_visible(self['visible']['title'])
		elif key == 'verify':
			self.draw_fig()
		if self['visible']['bp_sli'] == False and self['visible']['psi_sli'] == False:
			self.fig.subplots_adjust(bottom=0.11)
		elif self['visible']['bp_sli'] == False: 
			self.fig.subplots_adjust(bottom=0.13)
	
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
		psi_idx = self.psi_slider.val
		bp_idx = self.bp_slider.val
		sh_idx = self.sh_slider.val
		psiN = self.psiNs[psi_idx]
		
		bp = self.data['beta_prime_axis'][psi_idx][bp_idx]
		sh = self.data['shear_axis'][psi_idx][sh_idx]
		
		if not self['aky']:
			ky = self.data['akys'][psi_idx][bp_idx][sh_idx]
			aky_idx = self.inputs['aky_values'].index(ky)
			ky.slider.set_val(aky_idx)
		else:
			aky_idx = self.ky_slider.val
			ky = self.inputs['aky_values'][aky_idx]
		
		self.settings['psi_id'] = psi_idx
		self.settings['bp_id'] = bp_idx
		self.settings['sh_id'] = sh_idx
		self.settings['ky_id'] = aky_idx

		if self.data['time'] is None and self['var'] in [0,3]:
			path = self.info['data_path']
			try:
				from ..ncdf2dict import ncdf2dict as readnc
				run = readnc(f"{path}/{psiN}/{bp_idx}_{sh_idx}/{psi_idx}_{bp_idx}_{sh_idx}_{aky_idx}.out.nc")
				t = run['t']
			except:
				print(f"ERROR: Unable to read data, either use Detailed Save or ensure data is at specified folder {path}")
				return
		else:
			t = self.data['time'][psi_idx][bp_idx][sh_idx][aky_idx]
			
		if self.data['theta'] is None and self['var'] in [1,2]:
			path = self.info['data_path']
			try:
				from ..ncdf2dict import ncdf2dict as readnc
				run = readnc(f"{path}/{psiN}/{bp_idx}_{sh_idx}/{psi_idx}_{bp_idx}_{sh_idx}_{aky_idx}.out.nc")
				t = run['theta']
			except:
				print(f"ERROR: Unable to read data, either use Detailed Save or ensure data is at specified folder {path}")
				return
		else:
			theta = self.data['theta'][psi_idx][bp_idx][sh_idx][aky_idx]
		
		self.ax.set_title(f"psiN: {psiN} | -\u03B2': {round(bp,3)} | shear: {round(sh,3)} | aky: {ky} | file: {psiN}/{bp_idx}_{sh_idx}/{psi_idx}_{bp_idx}_{sh_idx}_{aky_idx}",fontsize=self['fontsizes']['title'])
		
		if self['var'] == 0:
			if self.data['omega'] is None:
				omega = run['omega'][:,0,0]
			else:
				omega = self.data['omega'][psi_idx][bp_idx][sh_idx][aky_idx]
			
			self.ax2.cla()
			
			self.ax2.plot(t,real(omega),'r',label="mode frequency")
			self.ax.plot(t,imag(omega),'b',label="growth rate")
			
			self.ax.text(0.01,0.99,f"GR: {self.data['growth_rates_all'][psi_idx][bp_idx][sh_idx][aky_idx]:+.2e}\nOmega[-1]: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax.transAxes,fontsize=self['fontsizes']['legend'])
			self.ax2.text(0.01,0.99,f"MF: {real(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax2.transAxes,fontsize=self['fontsizes']['legend'])
			self.ax.set_ylabel("Growth Rate",fontsize=self['fontsizes']['axis'])
			self.ax.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			self.ax2.set_ylabel("Mode_Frequency",fontsize=self['fontsizes']['axis'])
			self.ax2.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			self.ax.legend(loc=1)
			self.ax2.legend(loc=1)
			
		elif self['var'] == 1:
			if self.data['phi'] is None:
				phi = run['phi'][0,0,:]
			else:
				phi = self.data['phi'][psi_idx][bp_idx][sh_idx][aky_idx]
			
			self.ax.plot(theta,real(phi),'r',label="real")
			self.ax.plot(theta,imag(phi),'b',label="imaginary")
			if self['absolute']:
				self.ax.plot(theta,[abs(x) for x in phi],'k--',label="absolute")
				
			self.ax.set_ylabel("Electric Potential",fontsize=self['fontsizes']['axis'])
			self.ax.set_xlabel("Ballooning Angle",fontsize=self['fontsizes']['axis'])
			self.ax.legend(loc=0)
			
		elif self['var'] == 2:
			if self.data['apar'] is None:
				apar = run['apar'][0,0,:]
			else:
				apar = self.data['apar'][psi_idx][bp_idx][sh_idx][aky_idx]
			self.ax.plot(theta,real(apar),'r',label="real")
			self.ax.plot(theta,imag(apar),'b',label="imaginary")
			if self['absolute']:
				self.ax.plot(theta,[abs(x) for x in apar],'k--',label="absolute")
			
			self.ax.set_ylabel("Parallel Mangetic Potential",fontsize=self['fontsizes']['axis'])
			self.ax.set_xlabel("Ballooning Angle",fontsize=self['fontsizes']['axis'])
			self.ax.legend(loc=0)
		
		elif self['var'] == 3:
			if self.data['bpar'] is None:
				bpar = run['bpar'][0,0,:]
			else:
				bpar = self.data['bpar'][psi_idx][bp_idx][sh_idx][aky_idx]
			self.ax.plot(theta,real(bpar),'r',label="real")
			self.ax.plot(theta,imag(bpar),'b',label="imaginary")
			if self['absolute']:
				self.ax.plot(theta,[abs(x) for x in bpar],'k--',label="absolute")
			
			self.ax.set_ylabel("Parallel Mangetic Field",fontsize=self['fontsizes']['axis'])
			self.ax.set_xlabel("Ballooning Angle",fontsize=self['fontsizes']['axis'])
			self.ax.legend(loc=0)
			
		elif self['var'] == 4:
			if self.data['phi2'] is None:
				phi2 = run['phi2']
			else:
				phi2 = self.data['phi2'][psi_idx][bp_idx][sh_idx][aky_idx]
			self.ax.plot(t,phi2,'k')
			if max(phi2) > 1e267:
				self.ax.set_ylim(min(phi2)/100, 1e267) #Display error occurs on log scale above 1e267
			self.ax.set_ylabel("Phi2",fontsize=self['fontsizes']['axis'])
			self.ax.set_yscale('log')
			self.ax.set_xlabel(f"Time ({len(t)} steps)",fontsize=self['fontsizes']['axis'])
			
			if self.options.get_status()[0] and self.verify is not None :
				nt = self.verify.nts[psi_idx][bp_idx][sh_idx][aky_idx]
				if nt is not None:
					fit = polyfit(t[-nt:],log(phi2[-nt:]),1)
					fitgr = fit[0]/2
					pr = pearsonr(t[-nt:],log(phi2[-nt:]))
					gradgr = log(phi2[-1]/phi2[0])/(t[-1]-t[0])/2
					omega = self.data['omega'][psi_idx][bp_idx][sh_idx][aky_idx]
					
					self.ax.plot(t[-nt:],exp(array(t[-nt:])*fit[0] + fit[1]),'r',label=f"GR = {fitgr:+.2e}\nR = {pr[0]:.4f}")
					self.ax.plot([t[0],t[-1]],[phi2[0],phi2[-1]],'b',label=f"AVG GR = {gradgr:+.2e}")
					self.ax.text(0.01,0.99,f"GR: {self.data['growth_rates_all'][psi_idx][bp_idx][sh_idx][aky_idx]:+.2e}\nOmega[-1]: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=self.ax.transAxes,fontsize=self['fontsizes']['legend'])
					self.ax.set_xlabel(f"Time ({len(t)} steps) | nt = {nt}")
					self.ax.legend(loc=1,fontsize=self['fontsizes']['legend'])
		
		if self.verify is not None and self['visible']['verify']:
			bad = []
			if (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['nstep']:
				bad.append('nstep')
			if (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['nan']:
				bad.append('nan')
			if (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['unconv']:
				bad.append('unconverged')
			if (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['negative']:
				bad.append('negative')
			if (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['order']:
				bad.append('order')
			if self['var'] == 1 and (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['phi']:
				bad.append('phi')
			if self['var'] == 2 and (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['apar']:
				bad.append('apar')
			if self['var'] == 3 and (psi_idx,bp_idx,sh_idx,aky_idx) in self.verify['apar']:
				bad.append('bpar')
			if bad:
				self.ax.text(0.01,0.01,f"BAD RUN: {str(bad)[1:-1]}",ha='left',va='bottom',transform=self.ax.transAxes,color='r')
		self.fig.canvas.draw_idle()
		return

	
	
def plot_diag_single(data = None, var = 0, fig = None, ax = None, ax2 = None):
	if data is None:
		print("ERROR: no data dictionary given")
		return
	if var == "omega":
		var = 0
	if var == "phi":
		var = 1
	if var == "apar":
		var = 2
	if var == "bpar":
		var = 3
	if var == "phi2":
		var = 4
	if var not in [0,1,2,3,4]:
		print(f"ERROR: variable name/value {var} not supported. supported: omega/0, phi/1, apar/2, phi2/3")
		return
	
	doShow = False
	if fig is None or ax is None:
		if var == 0:
			fig, [ax, ax2] = subplots(2,1,figsize=(10,10))
		else:
			fig, ax = subplots(figsize=(8.8,5.8))
		doShow = True
	
	if var == 0:
		omega = data['omega']
		t = data['t']
		
		ax.plot(t,imag(omega),'b',label="growth rate")
		ax.text(0.01,0.99,f"GR: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=ax.transAxes)
		ax.set_ylabel("Growth Rate")
		ax.set_xlabel(f"Time ({len(t)} steps)")
		ax.legend(loc=0)
		
		ax2.cla()	
		ax2.plot(t,real(omega),'r',label="mode frequency")
		ax2.text(0.01,0.99,f"MF: {real(omega[-1]):+.2e}",ha='left',va='top',transform=ax2.transAxes)
		ax2.set_ylabel("Mode_Frequency")
		ax2.set_xlabel(f"Time ({len(t)} steps)")
		ax2.legend(loc=1)

	elif var == 1:
		phi = data['phi']
		theta = data['theta']

		ax.plot(theta,real(phi),'r',label="real")
		ax.plot(theta,imag(phi),'b',label="imaginary")
		#ax.plot(theta,[abs(x) for x in phi],'k--',label="absolute")
		
		ax.set_ylabel("Electric Potential")
		ax.set_xlabel("Ballooning Angle")
		ax.legend(loc=0)
		
	elif var == 2:
		apar = data['apar']
		theta = data['theta']

		ax.plot(theta,real(apar),'r',label="real")
		ax.plot(theta,imag(apar),'b',label="imaginary")
		#ax.plot(theta,[abs(x) for x in apar],'k--',label="absolute")
		
		ax.set_ylabel("Parallel Mangetic Potential")
		ax.set_xlabel("Ballooning Angle")
		ax.legend(loc=0)
		
	elif var == 3:
		bpar = data['bpar']
		theta = data['theta']

		ax.plot(theta,real(bpar),'r',label="real")
		ax.plot(theta,imag(bpar),'b',label="imaginary")
		#ax.plot(theta,[abs(x) for x in bpar],'k--',label="absolute")
		
		ax.set_ylabel("Parallel Mangetic Field")
		ax.set_xlabel("Ballooning Angle")
		ax.legend(loc=0)
	
	elif var == 4:
		def draw_fig(val):
			ax.cla()
			if max(phi2) > 1e267:
				ax.set_ylim(min(phi2)/100, 1e267) #Display error occurs on log scale above 1e267
				ax.autoscale(False)
			else:
				ax.autoscale(True)
			ax.plot(t,phi2,'k')
			ax.set_ylabel("Phi2")
			ax.set_yscale('log')
			ax.set_xlabel(f"Time ({len(t)} steps)")
			
			nt = int(eval(val))
			fit = polyfit(t[-nt:],log(phi2[-nt:]),1)
			pr = pearsonr(t[-nt:],log(phi2[-nt:]))
			grad = log(phi2[-1]/phi2[0])/(t[-1]-t[0])
			gr = fit[0]/2
			ax.plot(t[-nt:],exp(array(t[-nt:])*fit[0] + fit[1]),'r',label=f"GR = {gr:+.2e}\nR = {pr[0]:.4f}")
			ax.plot([t[0],t[-1]],[phi2[0],phi2[-1]],'b',label=f"AVG GR = {grad/2:+.2e}")
			ax.text(0.01,0.99,f"Omega[-1]: {imag(data['omega'][-1]):+.2e}",ha='left',va='top',transform=ax.transAxes)
			ax.legend(loc=1)
			
		phi2 = data['phi2']
		t = data['t'] 
		
		taxes = fig.add_axes([0.8, 0.03, 0.05, 0.04])
		tbox = TextBox(taxes, 'nt:')
		tbox.on_submit(draw_fig)
		
		if len(phi2) > 100:
			tbox.set_val(f"{len(phi2)//10}")
		else:
			tbox.set_val("10")
			
	if doShow:
		show()
	return

def plot_diag_set(runs = None, var = 0, init = None):
	def draw_fig(val):
		ax.cla()
		for run in runs.values():
			correct_run = True
			for key in sliders.keys():
				if run[key] != sliders[key]['vals'][sliders[key]['slider'].val]:
					correct_run = False
			if correct_run:
				data = run['data']
				break
		plot_diag_single(data = data, var = var, fig = fig, ax = ax, ax2 = ax2)
		ax.set_title(f"{[x for x in sliders.keys()]} = {[run[x] for x in sliders.keys()]}")
		fig.canvas.draw_idle()
	
	if var == 0:
		fig, [ax, ax2] = subplots(2,1,figsize=(10,10))
	else:
		fig, ax = subplots(figsize=(8.8,5.8))
		ax2 = None
	subplots_adjust(bottom=0.15)
	
	sliders = {}
	for i, key in enumerate([x for x in runs[0].keys() if x != 'data']):
		sliders[key] = {}
		vals = set()
		for run in runs.values():
			vals.add(run[key])
		sliders[key]['vals'] = list(vals)
		sliders[key]['vals'].sort()
		ypos = 0.01 + 0.02*i
		slaxes = axes([0.25, ypos, 0.5, 0.03])
		try:
			ini = init[i]
		except:
			ini = 0
		sliders[key]['slider'] = Slider(slaxes, f'{key} index:', 0, len(sliders[key]['vals'])-1, valinit = ini, valstep = 1)
		sliders[key]['slider'].on_changed(draw_fig)

	draw_fig(init)
	show()
