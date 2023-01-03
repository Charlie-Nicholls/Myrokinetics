import os
from numpy import *
from matplotlib.pyplot import *
import matplotlib.patches as pt
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap

'''
GYROKINETIC SCAN ANALYSIS
'''

def plot_scan(scan = None, aky = False, init = [0,0]):
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
		
		beta_prime = data['beta_prime_values'][idx]
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
			z_gr = transpose(data['growth_rates'][idx])
			z_mf = transpose(data['mode_frequencies'][idx])
			ax[0].set_title(f"Growth Rate | PsiN: {psiN}")
		
		if options.get_status()[2]:
			mfmax = mf_slider.val * abs(amax(array(data['mode_frequencies'])[isfinite(data['mode_frequencies'])]))/100
			grmax = gr_slider.val * abs(amax(array(data['growth_rates'])[isfinite(data['growth_rates'])]))/100

		else:
			if aky:
				mfmax = mf_slider.val * abs(amax(array(data['mode_frequencies_all'])[idx,:,:,ky_idx][isfinite(array(data['mode_frequencies_all'])[idx,:,:,ky_idx].tolist())]))/100
				grmax = gr_slider.val * abs(amax(array(data['growth_rates_all'])[idx,:,:,ky_idx][isfinite(array(data['growth_rates_all'])[idx,:,:,ky_idx].tolist())]))/100

			else:
				mfmax = mf_slider.val * abs(amax(array(data['mode_frequencies'][idx])[isfinite(data['mode_frequencies'][idx])]))/100	
				grmax = gr_slider.val * abs(amax(array(data['growth_rates'][idx])[isfinite(data['growth_rates'][idx])]))/100
			if mfmax == 0:
				mfmax = 10e-10
			if grmax == 0:
				grmax = 10e-10
		
		norm_mf = Normalize(vmin=-mfmax,vmax=mfmax)
		cbar_mf.update_normal(ScalarMappable(norm = norm_mf))
		
		norm_gr = Normalize(vmin=-grmax,vmax=grmax)
		cbar_gr.update_normal(ScalarMappable(norm = norm_gr, cmap = cmap))
		
		ax[0].pcolormesh(x, y, z_gr, cmap = cmap, norm=norm_gr)
		ax[1].pcolormesh(x, y, z_mf, norm = norm_mf)
		
		if options.get_status()[0]:
			for bpid, bp in enumerate(x):
				for shid, sh in enumerate(y):
					if data['parities'][idx][bpid][shid] == 1:
						ax[0].plot(bp, sh, '+', color = 'purple')
						ax[1].plot(bp, sh, '+', color = 'purple')
					if data['parities'][idx][bpid][shid] == -1:
						ax[0].plot(bp, sh, '_', color = 'cyan')
						ax[1].plot(bp, sh, '_', color = 'cyan')
		
		if options.get_status()[3]:
			ax[0].plot(beta_prime,shear,'kx')
			ax[0].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			ax[0].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
			ax[1].plot(beta_prime,shear,'kx')
			ax[1].annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			ax[1].annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
		if options.get_status()[1]:
			ax[0].set_ylim(inputs['shat_min'],inputs['shat_max'])
			ax[0].set_xlim(bpmin,bpmax)
			ax[1].set_ylim(inputs['shat_min'],inputs['shat_max'])
			ax[1].set_xlim(bpmin,bpmax)
		
		if options.get_status()[4] and data['ideal_stabilities'] is not None and data['ideal_stabilities'][idx] is not None:
			ax[0].contourf(data['beta_prime_axis_ideal'][idx], data['shear_axis_ideal'][idx], data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
			ax[1].contourf(data['beta_prime_axis_ideal'][idx], data['shear_axis_ideal'][idx], data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
			
		fig.canvas.draw_idle()
		return

	bpmin = amin(data['beta_prime_axis'])
	bpmax = amax(data['beta_prime_axis'])
	grlim = amax(abs(array(data['growth_rates'])))
	mflim = amax(abs(array(data['mode_frequencies'])))
	
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
		
	slaxes = axes([0.15, 0.01, 0.5, 0.03])
	slider = Slider(slaxes, 'psiN index:', 0, len(psiNs)-1, valinit = init[0], valstep = 1)
	slider.on_changed(draw_fig)
	
	if aky:
		kys = len(inputs['aky_values'])
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
