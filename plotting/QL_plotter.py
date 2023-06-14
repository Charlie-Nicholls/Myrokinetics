from numpy import transpose, array, amax, amin, isfinite, linspace, where
from matplotlib.pyplot import *
from matplotlib.cm import ScalarMappable
from matplotlib.widgets import Slider, CheckButtons, TextBox
from matplotlib.colors import LinearSegmentedColormap

class plot_QL(object):
	def __init__(self, scan = None, init = [0]):
		if scan is None:
			print("ERROR: no scan dictionary given")
			return
		self.scan = scan
		self.data = scan['data']
		self.inputs = scan['inputs']
		self.psiNs = self.inputs['psiNs']
		self.init = init
		if self.data['quasilinear'] is None:
			print("Error: No QuasiLinear Data")
			return
		
		self.bpmin = amin(self.data['beta_prime_axis'])
		self.bpmax = amax(self.data['beta_prime_axis'])
		self.shmin = amin(self.data['shear_axis'])
		self.shmax = amax(self.data['shear_axis'])
		cdict = {'red':  ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.8, 0.8)),
			'green':  ((0.0, 0.8, 0.8),(0.5, 1, 1),(1.0, 0.0, 0.0)),
			'blue': ((0.0, 0.0, 0.0),(0.5, 1, 1),(1.0, 0.0, 0.0))}
		self.cmap = LinearSegmentedColormap('GnRd', cdict)
		
		self.open_plot()
		
	def save_plot(self, filename = None):
		self.open_plot(save = True, filename = None)
		
	def open_plot(self, save = False, filename = "Scan"):
		self.fig, self.ax = subplots(figsize=(9,7))
		self.fig.subplots_adjust(bottom=0.15)   
		self.fig.suptitle(self.scan['info']['run_name'])
	
		blank_norm = Normalize(vmin=-1,vmax=1)
		self.cbar = colorbar(ScalarMappable(norm = blank_norm), ax = self.ax)
		
		try:
			op_init = self.options.get_status()
		except:
			op_init = [False,False,True,True,False]
		chaxes = axes([0.72, 0.01, 0.09, 0.1],frame_on = False)
		self.options = CheckButtons(chaxes, ["Show ID","Global Axis Limits","Global Colorbar","Show Equillibrium","Show Ideal"], op_init)
		self.options.on_clicked(self.draw_fig)
		
		try:
			sl_init = self.slider.val
		except:
			sl_init = self.init[0]
		if len(self.psiNs) > 1:
			slaxes = axes([0.15, 0.01, 0.5, 0.03])
			self.slider = Slider(slaxes, 'psiN index:', 0, len(self.psiNs)-1, valinit = sl_init, valstep = 1)
			self.slider.on_changed(self.draw_fig)
		
		try:
			ql_init = self.ql_init.val
		except:
			ql_init = 100
		ql_axes = axes([0.97, 0.15, 0.01, 0.73])
		self.ql_slider = Slider(ql_axes, 'Scale', 0, 100, valinit = 100, valstep = 1, orientation = 'vertical')
		self.ql_slider.on_changed(self.draw_fig)
		
		if save:
			if filename is None:
				filename = f"QL_{self.slider.val}"
			self.draw_fig()
			self.fig.savefig(filename)
		else:
			ion()
			show()
			self.draw_fig()
	
	def set_val(self, val):
		if self.options.get_status()[2]:
			sliderval = val/(abs(amax(array(self.data['quasilinear'])[isfinite(self.data['quasilinear'])]))/100)
		else:
			sliderval = val/(amax(abs(array(z)[isfinite(z)]))/100)
		self.ql_slider.set_val(sliderval)
	
	def draw_fig(self, val = None):
		if len(self.psiNs) > 1:
			idx = self.slider.val
			psiN = self.psiNs[idx]
			
		else:
			idx = 0
			psiN = self.psiNs[0]
		
		x = list(self.data['beta_prime_axis'][idx])
		y = list(self.data['shear_axis'][idx])
		
		beta_prime = abs(self.data['beta_prime_values'][idx])
		shear = self.data['shear_values'][idx]
		
		self.ax.cla()
		self.ax.set_facecolor('grey')
		self.ax.set_ylabel("Shear")
		self.ax.set_xlabel("-\u03B2'")
		self.ax.set_title(f"psiN: {psiN}")
		
		status = self.options.get_status()
		
		z = transpose(self.data['quasilinear'][idx]).tolist()
		
		if status[2]:
			qlmax = self.ql_slider.val * abs(amax(array(self.data['quasilinear'])[isfinite(self.data['quasilinear'])]))/100
		else:
			try:
				qlmax = self.ql_slider.val * amax(abs(array(z)[isfinite(z)]))/100
				if qlmax < 10e-10:
					qlmax = 10e-10
			except:
				qlmax = 1

		norm = Normalize(vmin=-qlmax,vmax=qlmax)
		self.cbar.update_normal(ScalarMappable(norm = norm, cmap = self.cmap))
		self.ax.pcolormesh(x, y, z, cmap = self.cmap, norm=norm)
		
		if status[0]:
			dx = x[1] - x[0]
			dy = y[1] - y[0]
			self.ax.set_xticks(array(x)-dx/2, minor=True)
			self.ax.set_yticks(array(y)-dy/2, minor=True)
			
			self.ax.set_xticks(array(x), minor=False)
			self.ax.set_yticks(array(y), minor=False)
			
			self.ax.set_xticklabels([])
			self.ax.set_yticklabels([])
			xlabels = [str(i) for i in range(len(x))]
			ylabels = [str(i) for i in range(len(y))]
			self.ax.set_xticklabels(xlabels,minor=False)
			self.ax.set_yticklabels(ylabels,minor=False)
			
			self.ax.grid(which="minor",color='k')
		
		if status[1]:
			self.ax.set_ylim(self.shmin,self.shmax)
			self.ax.set_xlim(self.bpmin,self.bpmax)

		if status[3]:
			self.ax.plot(beta_prime,shear,'kx')
			self.ax.annotate("Eqbm",(beta_prime,shear),textcoords = "offset points",xytext = (0,7), ha = "center")
			self.ax.annotate(f"{round(beta_prime,2)},{round(shear,2)}",(beta_prime,shear),textcoords = "offset points",xytext = (0,-13), ha = "center")
		
		if status[4]:
			if self.data['ideal_stabilities'] is not None and self.data['ideal_stabilities'][idx] is not None:
				self.ax.contourf(self.data['beta_prime_axis_ideal'][idx], self.data['shear_axis_ideal'][idx], self.data['ideal_stabilities'][idx], [0.01,0.99], colors = ('k'))
			else:
				self.ax.text(0.5,0.5,"No Ideal Data",ha='center',va='center',transform=self.ax.transAxes,color='k')

		self.fig.canvas.draw_idle()
