from numpy import linspace, pi, zeros, real, imag, shape
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.pyplot import *

def plot_theta(run = None, var = 0, fig = None, ax = None):
	if run is None:
		print("Run Not Given")
		return

	if var == 0:
		var = "phi"
	if var == 1:
		var = "apar"
	if var not in ["phi","apar"]:
		print(f"ERROR: variable name/value {var} not supported. supported: phi/0, apar/1")
		return

	b_angle = run['theta']
	e_fun = run[var]
	if len(shape(e_fun)) == 3:
		e_fun = e_fun[0,0,:]
	
	doShow = False
	if fig is None or ax is None:
		fig, ax = subplots(1,2,figsize=(14.6,7))
		doShow = True

	theta = linspace(0, 2*pi, 100)
	rt_fun = zeros((len(theta)))
	it_fun = zeros((len(theta)))
	re_fun = InterpolatedUnivariateSpline(b_angle,real(e_fun))
	ie_fun = InterpolatedUnivariateSpline(b_angle,imag(e_fun))
	for t, th in enumerate(theta):
		sampling = True
		i = 0
		while sampling:
			th_p = b_angle[0] + th + 2*pi*i
			if b_angle[0] <= th_p <= b_angle[-1]:
				rt_fun[t] += re_fun(th_p)
				it_fun[t] += ie_fun(th_p)
			elif th_p > b_angle[-1]:
				sampling = False
			i = i+1
	n=4
	for i in range(n):
		ax[1].plot(theta+i*2*pi, rt_fun,'b--')
		ax[1].plot(theta+i*2*pi, it_fun,'r--')
	ax[0].plot(b_angle, real(e_fun),'b--')
	ax[0].plot(b_angle, imag(e_fun),'r--')
	ax[0].set_xlabel("Ballooning Angle")
	ax[1].set_xlabel("Theta")
	ax[0].set_ylabel(var)
	ax[1].set_ylabel(var)
	if doShow:
		show()
	
def plot_theta_scan(scan = None, var = 0, init = [0,0,0,0], aky = True):
	if scan is None:
		print("Scan Not Given")
		return
	
	if var == 0:
		var = "phi"
	if var == 1:
		var = "apar"
	if var not in ["phi","apar"]:
		print(f"ERROR: variable name/value {var} not supported. supported: phi/0, apar/1")
		return	
		
	data = scan['data']
	inputs = scan['inputs']
	psiNs = inputs['psiNs']
	
	def draw_fig(val):
		ax[0].cla()
		ax[1].cla()
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
		
		run = {'theta': data['theta'][psi_idx][bp_idx][sh_idx][aky_idx], var: data[var][psi_idx][bp_idx][sh_idx][aky_idx]}
		plot_theta(run = run, var = var, fig = fig, ax = ax)
		fig.canvas.draw_idle()
	
	fig, ax = subplots(1,2,figsize=(14.6,7))
	subplots_adjust(bottom=0.15)

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
	draw_fig(0)
	show()

def plot_theta_set(runs = None, var = 0, init = None):
	if var == 0:
		var = "phi"
	if var == 1:
		var = "apar"
	if var not in ["phi","apar"]:
		print(f"ERROR: variable name/value {var} not supported. supported: phi/0, apar/1")
		return
	
	if runs is None:
		print("Runs Not Given")
		return

	def draw_fig(val):
		ax[0].cla()
		ax[1].cla()
		for run in runs.values():
			correct_run = True
			for key in sliders.keys():
				if run[key] != sliders[key]['vals'][sliders[key]['slider'].val]:
					correct_run = False
			if correct_run:
				data = run['data'].run
				break
		plot_theta(run = data, var = var, fig = fig, ax = ax)
		ax[0].set_title(f"{[x for x in sliders.keys()]} = {[run[x] for x in sliders.keys()]}")
		fig.canvas.draw_idle()
	
	fig, ax = subplots(1,2,figsize=(14.6,7))
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
