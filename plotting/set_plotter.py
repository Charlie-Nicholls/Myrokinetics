from matplotlib.pyplot import *
from matplotlib.widgets import Slider, CheckButtons

def plot_set(runs = None, var = None, variables = None, init = 0, aky_axis = False):
	if runs is None:
		print("ERROR: No runs given")
		return
	if var is None:
		var = [x for x in variables.keys() if x not in ['aky','psiN']][0]
	if var not in [x for x in variables.keys() if x not in ['aky','psiN']]:
		print("ERROR: Invalid Key")
		return
	init = int(init)
	
	def draw_fig(val):
		ax[1].cla()
		ax[0].cla()
		if psiNs:
			psiN = psiNs[slider.val]
		for run in runs.values():
			if (not psiNs or run['psiN'] == psiN) and run['data'] is not None:
				if options.get_status()[0] and akys:
					col = colours[variables[var].index(run[var])%len(colours)]
					mar = markers[variables[var].index(run[var])//len(markers)]
					ax[1].plot(run['aky'],run['data'].gr,marker=mar,c=col,label=str(run[var]))
					ax[0].plot(run['aky'],run['data'].mf,marker=mar,c=col,label=str(run[var]))
				else:
					aky = None
					col = colours[0]
					mar = markers[0]
					if akys:
						aky = run['aky']
						col = colours[akys.index(aky)%len(colours)]
						mar = markers[akys.index(aky)//len(markers)]
					ax[1].plot(run[var],run['data'].gr,marker=mar,c=col,label=str(aky))
					ax[0].plot(run[var],run['data'].mf,marker=mar,c=col,label=str(aky))
		if options.get_status()[0]:
			ax[1].set_xlabel("aky")
			if not akys:
				ax[1].text(0,0,"NO AKY DATA")
		else:
			ax[1].set_xlabel(var)
		if options.get_status()[1]:
			ax[0].set_xscale("log")
			ax[1].set_xscale("log")
		ax[1].set_ylabel("Growth Rate")
		ax[0].set_ylabel("Mode Frequency")
		if psiNs:
			fig.suptitle(f"{var} | psiN = {psiN}")
		else:
			fig.suptitle(f"{var}")
		lin, han = ax[0].get_legend_handles_labels()
		leg = dict(zip(han, lin))
		if options.get_status()[0]:
			ax[0].legend(leg.values(),leg.keys(), title = f"{var}", loc='upper left', bbox_to_anchor=(1, 1))
		elif akys:
			ax[0].legend(leg.values(),leg.keys(), title = "aky", loc='upper left', bbox_to_anchor=(1, 1))
		fig.canvas.draw_idle()
	
	fig, ax = subplots(2,1,figsize=(10,7))
	colours = ['r','b','g','k','y','m','c']
	markers = ['o','v','s','D','^','<','>']
	akys = None
	if 'aky' in variables.keys():
		akys = variables['aky']
	psiNs = None
	if 'psiN' in runs[0].keys():
		psiNs = variables['psiN']
	if psiNs:
		subplots_adjust(bottom=0.15)
		slaxes = axes([0.25, 0.01, 0.5, 0.03])
		slider = Slider(slaxes, 'psiN index:', 0, len(psiNs)-1, valinit = init, valstep = 1)
		slider.on_changed(draw_fig)

	chaxes = axes([0.8, 0.01, 0.1, 0.1],frame_on = False)
	options = CheckButtons(chaxes, ["Plot vs aky","Log x axis"],[aky_axis,False])
	options.on_clicked(draw_fig)

	draw_fig(init)
	show()
