from numpy import amin, amax
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons

def plot_epar(data = None, inputs = None):
	if data is None:
		print("ERROR: no data dictionary given")
		return
	if inputs is None:
		print("ERROR: inputs not given")
		return
	data = data
	psiNs = inputs['psiNs']
	if data['eparN'] is None:
		print("Error: No eparN data")
		return
	def draw_fig(val = None):
		nonlocal cbar
		if len(psiNs) > 1:
			idx = slider.val
			psiN = psiNs[slider.val]
		else:
			idx = 0
			psiN = psiNs[0]
		bp = data['beta_prime_axis'][idx]
		sh = data['shear_axis'][idx]
		eparN = data['eparN'][idx]
		bp_or = data['beta_prime_values'][idx]
		sh_or = data['shear_values'][idx]
		
		ax.cla()
		ax.set_facecolor('lightgrey')
		
		if options.get_status()[0]:
			ax.contourf(bp, sh, eparN, [-1,0,0.1,0.11,1], colors = ('lightgrey','purple','k','y'), extend = "both")
		else:
			ax.pcolormesh(bp,sh,eparN)
		
		if options.get_status()[2]:
			ax.plot(bp_or,sh_or,'kx')
			ax.annotate(f"({bp_or:.2f},{sh_or:.2f})",(bp_or,sh_or),textcoords = "offset points",xytext = (0,10), ha = "center")
		ax.set_ylabel("Shear")
		ax.set_xlabel("-\u03B2'")
		ax.set_title(f"PsiN: {psiN}")
		if options.get_status()[1]:
			ax.set_ylim(amin(data['shear_axis']),amax(data['shear_axis']))
			ax.set_xlim(amin(data['beta_prime_axis']),amax(data['beta_prime_axis']))
		
		fig.canvas.draw_idle()
		return
	
	fig, ax = subplots()
	subplots_adjust(bottom=0.15)
	fig.suptitle(scan['info']['run_name'])

	slaxes = axes([0.15, 0.03, 0.5, 0.03])
	slider = slider = Slider(slaxes, 'psiN index:', 0, len(psiNs)-1, valinit = 0, valstep = 1)
	slider.on_changed(draw_fig)
	
	chaxes = axes([0.72, 0.01, 0.09, 0.1],frame_on=False)
	options = CheckButtons(chaxes, ["Contour","Global Axis Limits","Show Inputs"],[True,False,True])
	options.on_clicked(draw_fig)
	
	bnorm = Normalize(vmin=0,vmax=1)
	cbar = colorbar(ScalarMappable(norm = bnorm), ax=ax)
	draw_fig(psiNs[0])
	show()

