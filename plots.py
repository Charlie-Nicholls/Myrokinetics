from .plotting import Plotters

class plot(object):
	def __init__(self, settings = {}, reader = None):
		if settings['plot_type'] not in Plotters.keys():
			print("ERROR: plot_type not found")
			return
		self.settings = settings
		self.reader = reader
		if reader is None:
				from .reader import myro_read 
				self.reader = myro_read(self.settings['reader'])
		self.launch_plot()


	def launch_plot(self):
		self.plot = Plotters[self.settings['plot_type']](reader = self.reader, settings = self.settings)
		self.settings = self.plot.settings

	def save_settings(self, filename = None):
		if filename is None:
			filename = f"{reader['run_name']}_{self.settings['plot_type']}"
		if filename[-4:] != ".dat":
			filename += ".dat"

		f = open(filename,'w')
		f.write(f"settings = {settings}")

def open_plot(filename):
	with open(filename,'r') as f:
		settings = eval(f.readlines()[0])
	return Plot(settings=settings)
