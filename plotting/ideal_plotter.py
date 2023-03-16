from numpy import *
from matplotlib.pyplot import *
import matplotlib.patches as pt
from matplotlib.widgets import Slider, CheckButtons
from time import sleep

class plot_ideal(object):
	def __init__(self, scan = None, init = 0):
		if scan is None:
			print("ERROR: no scan dictionary given")
			return
			
		self.scan = scan
		self.init = int(init)
		self.data = scan['data']
		self.psiNs = scan['inputs']['psiNs']
		if self.data['ideal_stabilities'] is None:
			print("Error: No ideal_ball data")
			
		self.handles = [pt.Patch(color='red', label='Unstable'), pt.Patch(color='green', label='Stable'),pt.Patch(color='yellow', label='Unknown')]
		self.bpmin = amin(self.data['beta_prime_axis_ideal'])
		self.bpmax = amax(self.data['beta_prime_axis_ideal'])
		self.shmin = amin(self.data['shear_axis_ideal'])
		self.shmax = amax(self.data['shear_axis_ideal'])
		
		self.open_plot()
	
	def save_plot(self, filename = None):
		self.open_plot(save = True, filename = filename)
	
	def open_plot(self, save = False, filename = None):
		self.fig, self.ax = subplots(figsize=(10,7))
		self.fig.subplots_adjust(bottom=0.15)
		self.fig.suptitle(self.scan['info']['run_name'])
		
		if len(self.psiNs) > 1:
			slaxes = axes([0.15, 0.01, 0.5, 0.03])
			try:
				valinit = self.slider.val
			except:
				valinit = self.init
			self.slider = Slider(slaxes, 'psiN index:', 0, len(self.psiNs)-1, valinit = valinit, valstep = 1)
			self.slider.on_changed(self.draw_fig)
		
		chaxes = axes([0.72, 0.01, 0.09, 0.1],frame_on=False)
		self.options = CheckButtons(chaxes, ["Colour","Global Axis Limits","Show Equillibrium","Show Legend"],[True,False,True,False])
		self.options.on_clicked(self.toggles)
		
		if save:
			if filename is None:
				filename = f"Ideal_{self.slider.val}"
			self.draw_fig()
			self.fig.savefig(filename)
		else:
			ion()
			show()
			self.draw_fig()
			
	def draw_fig(self, val = None):
		if len(self.psiNs) > 1:
			idx = self.slider.val
			psiN = self.psiNs[idx]
			
		else:
			idx = 0
			psiN = self.psiNs[0]
		bp = self.data['beta_prime_axis_ideal'][idx]
		shear = self.data['shear_axis_ideal'][idx]
		stab = self.data['ideal_stabilities'][idx]
		bp_or = abs(self.data['beta_prime_values'][idx])
		sh_or = self.data['shear_values'][idx]
		
		self.ax.cla()
		self.ax.set_facecolor('lightgrey')
		self.ax.set_ylabel("Shear")
		self.ax.set_xlabel("-\u03B2'")
		status = self.options.get_status()
		
		if status[0]:
			self.ax.contourf(bp, shear, stab, [0,0.01,0.99,1,], colors = ('g','y','r'), extend = "both")
		else:
			self.ax.contourf(bp, shear, stab, [0.01,0.99], colors = ('k'))
			
		if status[1]:
			self.ax.set_ylim(self.shmin,self.shmax)
			self.ax.set_xlim(self.bpmin,self.bpmax)
		
		if status[2]:
			self.ax.plot(bp_or,sh_or,'kx')
			self.ax.annotate("Eqbm",(bp_or,sh_or),textcoords = "offset points",xytext = (0,7), ha = "center")
			self.ax.annotate(f"{bp_or:.2f},{sh_or:.2f}",(bp_or,sh_or),textcoords = "offset points",xytext = (0,-13), ha = "center")
		
		if status[3]:
			self.ax.legend(handles = self.handles,loc = 'upper center', bbox_to_anchor=(0.5,1.1), ncol = 3)
		else:
			self.ax.set_title(f"PsiN: {psiN}")
		self.fig.canvas.draw_idle()
		
	def toggles(self, label):
		self.draw_fig()
		
	def __call__(self, val):
		self.slider.set_val(val)
