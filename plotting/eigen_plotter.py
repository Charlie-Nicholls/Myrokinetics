from numpy import *
from matplotlib.pyplot import *
from matplotlib.widgets import Slider

def plot_eigen(scan = None, var = 0, aky = False, init = [0,0,0,0]):
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
				run = readnc(f"{path}/{psiN}/{bp_idx}_{sh_idx}/{bp_idx}_{sh_idx}_{aky_idx}.out.nc")
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
				run = readnc(f"{path}/{psiN}/{bp_idx}_{sh_idx}/{bp_idx}_{sh_idx}_{aky_idx}.out.nc")
				t = run['theta']
			except:
				print(f"ERROR: Unable to read data, either use Detailed Save or ensure data is at specified folder {path}")
				return
		else:
			theta = data['theta'][psi_idx][bp_idx][sh_idx][aky_idx]
		
		ax.set_title(f"psiN: {psiN} | -\u03B2': {round(bp,3)} | shear: {round(sh,3)} | aky: {ky} | file: {psiN}/{bp_idx}_{sh_idx}/{bp_idx}_{sh_idx}_{aky_idx}")
		
		if var == 0:
			if data['omega'] is None:
				omega = run['omega'][:,0,0]
			else:
				omega = data['omega'][psi_idx][bp_idx][sh_idx][aky_idx]
			
			ax.plot(t,real(omega),'r',label="mode frequency")
			ax.plot(t,imag(omega),'b',label="growth rate")
			
			ax.set_ylabel("Omega")
			ax.set_xlabel(f"Time ({len(t)} steps)")
			ax.legend(loc=0)
			
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
			ax.set_xlabel("Time")
		
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
		kys = len(inputs['aky_values'])
		kyaxes = axes([0.93, 0.25, 0.021, 0.5])
		ky_slider = Slider(kyaxes, 'ky index:', 0, len(inputs['aky_values'])-1, valinit = init[3], valstep = 1,orientation = 'vertical')
		ky_slider.on_changed(draw_fig)
		
	draw_fig()
	
	show()
