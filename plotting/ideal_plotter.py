from numpy import *
from matplotlib.pyplot import *
import matplotlib.patches as pt
from matplotlib.widgets import Slider, CheckButtons

def plot_ideal(scan = None, init = 0):
	if scan is None:
		print("ERROR: no scan dictionary given")
		return
	init = int(init)
	data = scan['data']
	psiNs = scan['inputs']['psiNs']
	if data['ideal_stabilities'] is None:
		print("Error: No ideal_ball data")
		return
	def draw_fig(val = None):
		if len(psiNs) > 1:
			idx = slider.val
			psiN = psiNs[idx]
			
		else:
			idx = 0
			psiN = psiNs[0]
		bp = data['beta_prime_axis_ideal'][index]
		shear = data['shear_axis_ideal'][index]
		stab = data['ideal_stabilities'][index]
		bp_or = data['beta_prime_values'][index]
		sh_or = data['shear_values'][index]
		
		ax.cla()
		ax.set_facecolor('lightgrey')
		if options.get_status()[0]:
			ax.contourf(bp, shear, stab, [0,0.01,0.99,1,], colors = ('g','y','r'), extend = "both")
		else:
			ax.contourf(bp, shear, stab, [0.01,0.99], colors = ('k'))
		if options.get_status()[2]:
			ax.plot(bp_or,sh_or,'kx')
			ax.annotate("Eqbm",(bp_or,sh_or),textcoords = "offset points",xytext = (0,7), ha = "center")
			ax.annotate(f"{bp_or:.2f},{sh_or:.2f}",(bp_or,sh_or),textcoords = "offset points",xytext = (0,-13), ha = "center")
		ax.set_ylabel("Shear")
		ax.set_xlabel("-\u03B2'")
		if options.get_status()[1]:
			ax.set_ylim(amin(data['shear_axis']),amax(data['shear_axis']))
			ax.set_xlim(amin(data['beta_prime_axis']),amax(data['beta_prime_axis']))
		ax.legend(handles = [green,yellow,red],loc = 'upper center', bbox_to_anchor=(0.5,1.1), ncol = 3)
		ax.set_title(f"PsiN: {psiN}")
		fig.canvas.draw_idle()
		return
		
	def toggles(label):
		draw_fig(slider.val)
	
	fig, ax = subplots()
	subplots_adjust(bottom=0.15)
	fig.suptitle(scan['info']['run_name'])
	
	red = pt.Patch(color='red', label='Unstable')
	green = pt.Patch(color='green', label='Stable')
	yellow = pt.Patch(color='yellow', label='Unknown')
	
	slaxes = axes([0.15, 0.01, 0.5, 0.03])
	slider = Slider(slaxes, 'psiN index:', 0, len(psiNs)-1, valinit = init, valstep = 1)
	slider.on_changed(draw_fig)
	
	chaxes = axes([0.72, 0.01, 0.09, 0.1],frame_on=False)
	options = CheckButtons(chaxes, ["Colour","Global Axis Limits","Show Equillibrium"],[True,False,True])
	options.on_clicked(toggles)
	
	draw_fig(psiNs[0])
	show()
