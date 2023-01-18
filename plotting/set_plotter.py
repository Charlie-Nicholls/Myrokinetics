from matplotlib.pyplot import *
from matplotlib.widgets import Slider

def plot_set(runs = None, var = None, init = 0):
	if runs is None:
		print("ERROR: No runs given")
		return
	if var is None:
		var = [x for x in runs[0].keys() if x not in ['data','aky','psiN']][0]
	if var not in [x for x in runs[0].keys() if x not in ['data','aky','psiN']]:
		print("ERROR: Invalid Key")
		return
	init = int(init)
	
	def draw_fig(val):
		ax[1].cla()
		ax[0].cla()
		psiN = psiNs[val]
		for run in runs.values():
			if not psiNs or run['psiN'] == psiN:
				aky = None
				col = colours[0]
				if akys:
					aky = run['aky']	
					col = colours[akys.index(aky)%len(colours)]
				if run['data']:
					ax[1].plot(run[var],run['data'].run['gr'],'.',c=col,label=str(aky))
					ax[0].plot(run[var],run['data'].run['mf'],'.',c=col,label=str(aky))
		ax[1].set_xlabel(var)
		ax[1].set_ylabel("Growth Rate")
		ax[0].set_ylabel("Mode Frequency")
		fig.suptitle(f"{var} | psiN = {psiN}")
		if akys:
			lin, han = ax[0].get_legend_handles_labels()
			leg = dict(zip(han, lin))
			ax[0].legend(leg.values(),leg.keys(), title = "aky", loc='upper left', bbox_to_anchor=(1, 1))
		fig.canvas.draw_idle()
	
	fig, ax = subplots(2,1,figsize=(10,7))
	colours = ['r','b','g','k','y','m','c']
	akys = None
	if 'aky' in runs[0].keys():
		akys = set()
		for run in runs.values():
			akys.add(run['aky'])
		akys = list(akys)
		akys.sort()
	psiNs = None
	if 'psiN' in runs[0].keys():
		psiNs = set()
		for run in runs.values():
			psiNs.add(run['psiN'])
		psiNs = list(psiNs)
		psiNs.sort()

	if psiNs:
		slaxes = axes([0.25, 0.01, 0.5, 0.03])
		slider = Slider(slaxes, 'psiN index:', 0, len(psiNs)-1, valinit = init, valstep = 1)
		slider.on_changed(draw_fig)
	draw_fig(init)
	show()
