from numpy import real, imag
from matplotlib.pyplot import *
from matplotlib.widgets import Slider

def plot_diag(scan = None, var = 0, aky = True, init = [0,0,0,0], verify = None):
	if scan is None:
		print("ERROR: no scan dictionary given")
		return
	data = scan['data']
	inputs = scan['inputs']
	psiNs = inputs['psiNs']
	def draw_fig(val = None):
		ax.cla()
		psi_idx = psi_slider.val
		bp_idx = bp_slider.val
		sh_idx = sh_slider.val
		psiN = psiNs[psi_idx]
		
		bp = data['beta_prime_axis'][psi_idx][bp_idx]
		sh = data['shear_axis'][psi_idx][sh_idx]
		if aky:
			aky_idx = ky_slider.val
			ky = inputs['aky_values'][aky_idx]
		else:
			ky = data['akys'][psi_idx][bp_idx][sh_idx]
			aky_idx = inputs['aky_values'].index(ky)

		if data['time'] is None and var in [0,3]:
			path = scan['info']['data_path']
			try:
				from ncdf2dict import ncdf2dict as readnc
				run = readnc(f"{path}/{psiN}/{bp_idx}_{sh_idx}/{psi_idx}_{bp_idx}_{sh_idx}_{aky_idx}.out.nc")
				t = run['t']
			except:
				print(f"ERROR: Unable to read data, either use Detailed Save or ensure data is at specified folder {path}")
				return
		else:
			t = data['time'][psi_idx][bp_idx][sh_idx][aky_idx]
			
		if data['theta'] is None and var in [1,2]:
			path = scan['info']['data_path']
			try:
				from ncdf2dict import ncdf2dict as readnc
				run = readnc(f"{path}/{psiN}/{bp_idx}_{sh_idx}/{psi_idx}_{bp_idx}_{sh_idx}_{aky_idx}.out.nc")
				t = run['theta']
			except:
				print(f"ERROR: Unable to read data, either use Detailed Save or ensure data is at specified folder {path}")
				return
		else:
			theta = data['theta'][psi_idx][bp_idx][sh_idx][aky_idx]
		
		ax.set_title(f"psiN: {psiN} | -\u03B2': {round(bp,3)} | shear: {round(sh,3)} | aky: {ky} | file: {psiN}/{psi_idx}_{bp_idx}_{sh_idx}/{bp_idx}_{sh_idx}_{aky_idx}")
		
		if var == 0:
			if data['omega'] is None:
				omega = run['omega'][:,0,0]
			else:
				omega = data['omega'][psi_idx][bp_idx][sh_idx][aky_idx]
			
			ax2.cla()
			
			ax2.plot(t,real(omega),'r',label="mode frequency")
			ax.plot(t,imag(omega),'b',label="growth rate")
			
			ax.text(0.01,0.99,f"GR: {imag(omega[-1]):+.2e}",ha='left',va='top',transform=ax.transAxes)
			ax2.text(0.01,0.99,f"MF: {real(omega[-1]):+.2e}",ha='left',va='top',transform=ax2.transAxes)
			ax.set_ylabel("Growth Rate")
			ax.set_xlabel(f"Time ({len(t)} steps)")
			ax2.set_ylabel("Mode_Frequency")
			ax2.set_xlabel(f"Time ({len(t)} steps)")
			
		elif var == 1:
			if data['phi'] is None:
				phi = run['phi'][0,0,:]
			else:
				phi = data['phi'][psi_idx][bp_idx][sh_idx][aky_idx]

			ax.plot(theta,real(phi),'r--',label="real")
			ax.plot(theta,imag(phi),'b--',label="imaginary")
			ax.plot(theta,[abs(x) for x in phi],'k',label="absolute")
			
			ax.set_ylabel("Electric Potential")
			ax.set_xlabel("Ballooning Angle")
			ax.legend(loc=0)
			
		elif var == 2:
			if data['apar'] is None:
				apar = run['apar'][0,0,:]
			else:
				apar = data['apar'][psi_idx][bp_idx][sh_idx][aky_idx]
			ax.plot(theta,real(apar),'r--',label="real")
			ax.plot(theta,imag(apar),'b--',label="imaginary")
			ax.plot(theta,[abs(x) for x in apar],'k',label="absolute")
			
			ax.set_ylabel("Parallel Mangetic Potential")
			ax.set_xlabel("Ballooning Angle")
			ax.legend(loc=0)
		elif var == 3:
			if data['phi2'] is None:
				phi2 = run['phi2']
			else:
				phi2 = data['phi2'][psi_idx][bp_idx][sh_idx][aky_idx]
			ax.plot(t,phi2,'k')
			
			ax.set_ylabel("Phi2")
			ax.set_yscale('log')
			ax.set_xlabel(f"Time ({len(t)} steps)")
		
		if verify is not None:
			bad = []
			if verify['nstep'] is not None and [psi_idx,bp_idx,sh_idx,aky_idx] in verify['nstep']:
				bad.append('nstep')
			if verify['other'] is not None and [psi_idx,bp_idx,sh_idx,aky_idx] in verify['other']:
				bad.append('other')
			if verify['unconv'] is not None and [psi_idx,bp_idx,sh_idx,aky_idx] in verify['unconv']:
				bad.append('unconverged')
			if var == 1 and verify['phi'] is not None and [psi_idx,bp_idx,sh_idx,aky_idx] in verify['phi']:
				bad.append('phi')
			if var == 2 and verify['apar'] is not None and [psi_idx,bp_idx,sh_idx,aky_idx] in verify['apar']:
				bad.append('apar')
			if bad:
				ax.text(0.01,0.01,f"BAD RUN: {str(bad)[1:-1]}",ha='left',va='bottom',transform=ax.transAxes,color='r')
		fig.canvas.draw_idle()
		return

	if var == "omega":
		var = 0
	if var == "phi":
		var = 1
	if var == "apar":
		var = 2
	if var == "phi2":
		var = 3
	if var not in [0,1,2,3]:
		print(f"ERROR: variable name/value {var} not supported. supported: omega/0, phi/1, apar/2, phi2/3")
		return
	
	if var == 0:
		fig, [ax, ax2] = subplots(2,1,figsize=(10,10))
	else:
		fig, ax = subplots(figsize=(8.8,5.8))
	subplots_adjust(bottom=0.15)
	fig.suptitle(scan['info']['run_name'])

	psiaxes = axes([0.25, 0.01, 0.5, 0.03])
	bpaxes = axes([0.25, 0.05, 0.5, 0.03])
	shaxes = axes([0.02, 0.25, 0.021, 0.5])

	psi_slider = Slider(psiaxes, 'psiN index:', 0, len(psiNs)-1, valinit = init[0], valstep = 1)
	psi_slider.on_changed(draw_fig)
	
	bp_slider = Slider(bpaxes, "-\u03B2' index:", 0, len(data['beta_prime_axis'][0])-1, valinit = init[1], valstep = 1)
	bp_slider.on_changed(draw_fig)
	
	sh_slider = Slider(shaxes, 'shear\nindex:', 0, len(data['shear_axis'][0])-1, valinit = init[2], valstep = 1,orientation = 'vertical')
	sh_slider.on_changed(draw_fig)
	
	if aky:
		kyaxes = axes([0.93, 0.25, 0.021, 0.5])
		ky_slider = Slider(kyaxes, 'ky index:', 0, len(inputs['aky_values'])-1, valinit = init[3], valstep = 1,orientation = 'vertical')
		ky_slider.on_changed(draw_fig)
		
	draw_fig()
	show()
	
def plot_diag_single(data = None, var = 0, fig = None, ax = None):
	if data is None:
		print("ERROR: no data dictionary given")
		return
	if var == "omega":
		var = 0
	if var == "phi":
		var = 1
	if var == "apar":
		var = 2
	if var == "phi2":
		var = 3
	if var not in [0,1,2,3]:
		print(f"ERROR: variable name/value {var} not supported. supported: omega/0, phi/1, apar/2, phi2/3")
		return
	
	doShow = False
	if fig is None or ax is None:
		fig, ax = subplots()
		doShow = True
	
	if var == 0:
		omega = data['omega'][:,0,0]
		t = data['t']
		
		ax.plot(t,real(omega),'r',label="mode frequency")
		ax.plot(t,imag(omega),'b',label="growth rate")
		
		ax.text(0.01,0.99,f"GR: {imag(omega[-1]):+.2e}\nMF: {real(omega[-1]):+.2e}",ha='left',va='top',transform=ax.transAxes)
		ax.set_ylabel("Omega")
		ax.set_xlabel(f"Time ({len(t)} steps)")
		ax.legend(loc=0)
		
	elif var == 1:
		phi = data['phi'][0,0,:]
		theta = data['theta']

		ax.plot(theta,real(phi),'r--',label="real")
		ax.plot(theta,imag(phi),'b--',label="imaginary")
		ax.plot(theta,[abs(x) for x in phi],'k',label="absolute")
		
		ax.set_ylabel("Electric Potential")
		ax.set_xlabel("Ballooning Angle")
		ax.legend(loc=0)
		
	elif var == 2:
		apar = data['apar'][0,0,:]
		theta = data['theta']

		ax.plot(theta,real(apar),'r--',label="real")
		ax.plot(theta,imag(apar),'b--',label="imaginary")
		ax.plot(theta,[abs(x) for x in apar],'k',label="absolute")
		
		ax.set_ylabel("Parallel Mangetic Potential")
		ax.set_xlabel("Ballooning Angle")
		ax.legend(loc=0)
	elif var == 3:
		phi2 = data['phi2']
		t = data['t']

		ax.plot(t,phi2,'k')
		
		ax.set_ylabel("Phi2")
		ax.set_yscale('log')
		ax.set_xlabel("Time")
	
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
		plot_diag_single(data = data, var = var, fig = fig, ax = ax)
		ax.set_title(f"{[x for x in sliders.keys()]} = {[run[x] for x in sliders.keys()]}")
		fig.canvas.draw_idle()
	
	fig, ax = subplots(figsize=(8.8,6.8))
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
