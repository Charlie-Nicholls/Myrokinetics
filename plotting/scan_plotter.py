from numpy import transpose, array, amax, amin, isfinite, linspace
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap

def plot_scan(scan = None, verify = None, aky = False, init = [0,0]):
	if scan is None:
		print("ERROR: no scan dictionary given")
		return
	data = scan['data']
	inputs = scan['inputs']
	psiNs = inputs['psiNs']
	if data['growth_rates'] is None:
		print("Error: No Gyrokinetic Data")
		return
	def draw_fig(val = None):
		nonlocal cbar_mf, cbar_gr
		if len(psiNs) > 1:
			idx = slider.val
			psiN = psiNs[idx]
			
		else:
			idx = 0
			psiN = psiNs[0]
		
		x = list(data['beta_prime_axis'][idx])
		y = list(data['shear_axis'][idx])
		
		beta_prime = abs(data['beta_prime_values'][idx])
		shear = data['shear_values'][idx]
		
		ax[0].cla()
		ax[0].set_facecolor('grey')
		ax[0].set_ylabel("Shear")
		ax[0].set_xlabel("-\u03B2'")
		ax[1].cla()
		ax[1].set_facecolor('grey')
		ax[1].set_ylabel("Shear")
		ax[1].set_xlabel("-\u03B2'")
		ax[1].set_title("Mode Frequency")
		
		if aky:
			ky_idx = ky_slider.val
			ky = inputs['aky_values'][ky_idx]
			z_gr = transpose(array(data['growth_rates_all'])[idx,:,:,ky_idx]).tolist()
			z_mf = transpose(array(data['mode_frequencies_all'])[idx,:,:,ky_idx]).tolist()
			ax[0].set_title(f"Growth Rate | PsiN: {psiN} | ky: {ky}")
		else:
			ky = data['akys'][idx][bpid][shid]
			ky_idx = inputs['aky_values'].index(ky_idx)
			z_gr = transpose(data['growth_rates'][idx]).tolist()
			z_mf = transpose(data['mode_frequencies'][idx]).tolist()
			ax[0].set_title(f"Growth Rate | PsiN: {psiN}")
		
		if options.get_status()[2]:
			mfmax = mf_slider.val * abs(amax(array(data['mode_frequencies'])[isfinite(data['mode_frequencies'])]))/100
			grmax = gr_slider.val * abs(amax(array(data['growth_rates'])[isfinite(data['growth_rates'])]))/100
		else:
			if [i for x in z_gr for i in x if str(i) != 'nan']:
				grmax = gr_slider.val * amax(abs(array(z_gr)[isfinite(z_gr)]))/100
				if grmax < 10e-10:
					grmax = 10e-10
			else:
				grmax = 1

			if [i for x in z_mf for i in x if str(i) != 'nan']:
				mfmax = mf_slider.val * amax(abs(array(z_mf)[isfinite(z_mf)]))/100
				if mfmax == 0:
					mfmax = 10e-10
			else:
				mfmax = 1
		
		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		
		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = cmap))

		ax[0].pcolormesh(x, y, z_gr, cmap = cmap, norm=norm_gr)
		ax[1].pcolormesh(x, y, z_mf, norm = norm_mf)
		
		if options.get_status()[0]:
			if aky:
				parities = array(data['parities_all'])[idx,:,:,ky_idx].tolist()
			else:
				parities = data['parities'][idx]
			sym_ids = where(parities == 1)
			anti_ids = where(parities == -1)
			ax[0].plot(sym_ids[0], sym_ids[1], '+', color = 'purple', label = 'par')
			ax[1].plot(sym_ids[0], sym_ids[1], '+', color = 'purple', label = 'par')
			ax[0].plot(anti_ids[0], anti_ids[1], '_', color = 'cyan', label = 'par')
			ax[1].plot(anti_ids[0], anti_ids[1], '_', color = 'cyan', label = 'par')
		
		if options.get_status()[1]:
			ax[0].set_ylim(shmin,shmax)
			ax[0].set_xlim(bpmin,bpmax)
			ax[1].set_ylim(shmin,shmax)
			ax[1].set_xlim(bpmin,bpmax)

		if options.get_status()[3]:
			ax[0].plot(beta_prime,shear,'kx')
			ax[0].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			ax[0].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
			ax[1].plot(beta_prime,shear,'kx')
			ax[1].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			ax[1].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
		
		if options.get_status()[4]:
			if data['ideal_stabilities'] is not None and data['ideal_stabilities'][idx] is not None:
				ax[0].contourf(data['beta_prime_axis_ideal'][idx], data['shear_axis_ideal'][idx], data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
				ax[1].contourf(data['beta_prime_axis_ideal'][idx], data['shear_axis_ideal'][idx], data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
			else:
				ax[0].text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=ax[0].transAxes,color='k')
		
		if vroptions.get_status()[0]:
			dx = x[1] - x[0]
			dy = y[1] - y[0]
			ax[0].set_xticks(array(x)-dx/2, minor=True)
			ax[1].set_xticks(array(x)-dx/2, minor=True)
			ax[0].set_yticks(array(y)-dy/2, minor=True)
			ax[1].set_yticks(array(y)-dy/2, minor=True)
			
			ax[0].set_xticks(array(x), minor=False)
			ax[1].set_xticks(array(x), minor=False)
			ax[0].set_yticks(array(y), minor=False)
			ax[1].set_yticks(array(y), minor=False)
			
			ax[0].set_xticklabels([])
			ax[1].set_xticklabels([])
			ax[0].set_yticklabels([])
			ax[1].set_yticklabels([])
			xlabels = [str(i) for i in range(len(x))]
			ylabels = [str(i) for i in range(len(y))]
			ax[0].set_xticklabels(xlabels,minor=False)
			ax[1].set_xticklabels(xlabels,minor=False)
			ax[0].set_yticklabels(ylabels,minor=False)
			ax[1].set_yticklabels(ylabels,minor=False)
			
			ax[0].grid(which="minor",color='k')
			ax[1].grid(which="minor",color='k')
			
		if vroptions.get_status()[1]:
			for bpid, bp in enumerate(x):
				for shid, sh in enumerate(y):
					if (idx, bpid, shid, ky_idx) in verify['unconverged']:
						ax[0].text(bp, sh, 'U', color = 'k',ha='center',va='center',size=7)
						ax[1].text(bp, sh, 'U', color = 'k',ha='center',va='center',size=7)
					elif (idx, bpid, shid, ky_idx) in verify['unconverged_stable']:
						ax[0].text(bp, sh, 'Us', color = 'k',ha='center',va='center',size=7)
						ax[1].text(bp, sh, 'Us', color = 'k',ha='center',va='center',size=7)
					elif (idx, bpid, shid, ky_idx) in verify['converged']:
						ax[0].text(bp, sh, 'C', color = 'k',ha='center',va='center',size=7)
						ax[1].text(bp, sh, 'C', color = 'k',ha='center',va='center',size=7)
					elif (idx, bpid, shid, ky_idx) in verify['converged_fit']:
						ax[0].text(bp, sh, 'Cf', color = 'k',ha='center',va='center',size=7)
						ax[1].text(bp, sh, 'Cf', color = 'k',ha='center',va='center',size=7)
		
		if vroptions.get_status()[2]:
			for bpid, bp in enumerate(x):
				for shid, sh in enumerate(y):
					if (idx, bpid, shid, ky_idx) in verify['nstep']:
						ax[0].text(bp, sh, 'n', color = 'k',ha='center',va='center',size=7)
						ax[1].text(bp, sh, 'n', color = 'k',ha='center',va='center',size=7)
		
		if vroptions.get_status()[3]:
			for bpid, bp in enumerate(x):
				for shid, sh in enumerate(y):
					s = ''
					if (idx, bpid, shid, ky_idx) in verify['phi']:
						s += 'p,'
					if (idx, bpid, shid, ky_idx) in verify['apar']:
						s += 'a,'
					if (idx, bpid, shid, ky_idx) in verify['bpar']:
						s += 'b,'
					if s != '':
						ax[0].text(bp, sh, s[:-1], color = 'k',ha='center',va='center',size=7, rotation = 45)
						ax[1].text(bp, sh, s[:-1], color = 'k',ha='center',va='center',size=7, rotation = 45)
		
		if vroptions.get_status()[4]:
			for bpid, bp in enumerate(x):
				for shid, sh in enumerate(y):
					if (idx, bpid, shid, ky_idx) in verify['nan']:
						ax[0].text(bp, sh, 'n', color = 'k',ha='center',va='center',size=7)
						ax[1].text(bp, sh, 'n', color = 'k',ha='center',va='center',size=7)
		
		fig.canvas.draw_idle()
		return

	bpmin = amin(data['beta_prime_axis'])
	bpmax = amax(data['beta_prime_axis'])
	shmin = amin(data['shear_axis'])
	shmax = amax(data['shear_axis'])
	
	fig, ax = subplots(1,2,figsize=(14.6,7))
	subplots_adjust(bottom=0.15)   
	fig.suptitle(scan['info']['run_name'])
	
	blank_norm = Normalize(vmin=-1,vmax=1)
	cbar_mf = colorbar(ScalarMappable(norm = blank_norm), ax = ax[1])
	cbar_gr = colorbar(ScalarMappable(norm = blank_norm), ax = ax[0])
	
	cdict1 = {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
		'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
		'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))}
	cmap = LinearSegmentedColormap('GnRd', cdict1)

	chaxes = axes([0.72, 0.01, 0.09, 0.1],frame_on = False)
	options = CheckButtons(chaxes, ["Show Parities","Global Axis Limits","Global Colorbar","Show Equillibrium","Show Ideal"],[False,False,False,True,False])
	options.on_clicked(draw_fig)
	
	vraxes = axes([0.85, 0.01, 0.09, 0.1],frame_on = False)
	vroptions = CheckButtons(vraxes, ["Show ID","Show Convergence","Show Bad nstep","Show bad fields","Show Omega -> nan"],[False,False,False,False,False])
	vroptions.on_clicked(draw_fig)
		
	slaxes = axes([0.15, 0.01, 0.5, 0.03])
	slider = Slider(slaxes, 'psiN index:', 0, len(psiNs)-1, valinit = init[0], valstep = 1)
	slider.on_changed(draw_fig)
	
	if aky:
		kyaxes = axes([0.15, 0.05, 0.5, 0.03])
		ky_slider = Slider(kyaxes, 'ky index:', 0, len(inputs['aky_values'])-1, valinit = init[1], valstep = 1)
		ky_slider.on_changed(draw_fig)
	
	gr_axes = axes([0.93, 0.15, 0.01, 0.73])
	gr_slider = Slider(gr_axes, 'GR', 0, 100, valinit = 100, valstep = 1, orientation = 'vertical')
	gr_slider.on_changed(draw_fig)
	
	mf_axes = axes([0.97, 0.15, 0.01, 0.73])
	mf_slider = Slider(mf_axes, 'MF', 0, 100, valinit = 100, valstep = 1, orientation = 'vertical')
	mf_slider.on_changed(draw_fig)
	
	draw_fig()
	show()
