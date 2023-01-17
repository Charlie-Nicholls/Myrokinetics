from matplotlib.pyplot import *
from matplotlib.widgets import Slider

def plot_set(self, runs = None, var = None, init = 0):
	if runs is None:
		print("ERROR: No runs given")
		return
	if var is None:
		var = [x for x in self.runs[0].keys() if x not in ['data','aky','psiN']][0]
	if var not in [x for x in self.runs[0].keys() if x not in ['data','aky','psiN']]:
		print("ERROR: Invalid Key")
		return
	init = int(init)
	
	def draw_fig(val):
		ax[1].cla()
		ax[0].cla()
		psiN = val
		for run in runs:
			if psiNs and run['psiN'] != psiN:
				break
			aky = None
			if akys:
				aky = str(run['aky'])
			
		col = colours[akys.index(aky)%len(colours)]
		ax[1].plot(run[var],run['data'].run['gr'],'.',c=col,label=aky)
		ax[0].plot(run[var],run['data'].run['mf'],'.',c=col,label=aky)
		ax[1].set_xlabel(var)
		ax[1].set_ylabel("Growth Rate")
		ax[0].set_ylabel("Mode Frequency")
		if akys:
			ax[0].legend(loc=0)
			ax[1].legend(loc=1)
		fig.canvas.draw_idle()
	
	fig, ax = subplots(2,1)
	colours = ['r','b','g','k','y','m','c']
	akys = None
	if 'aky' in runs[0].keys():
		akys = set()
		for run in runs:
			akys.add(run['aky'])
	akys = list(akys).sort()
	psiNs = None
	if 'psiN' in runs[0].keys():
		psiNs = set()
		for run in runs:
			psiNs.add(run['psiN'])
	psiNs = list(psiNs).sort()
	
	if psiNs:
		slaxes = axes([0.25, 0.01, 0.5, 0.03])
		slider = Slider(slaxes, 'psiN index:', 0, len(psiNs)-1, valinit = init, valstep = 1)
		slider.on_changed(draw_fig)
	draw_fig(init)
	show()
