viking_modules = ''''''

archer2_modules = ''''''

class code(object):
	def __init__(self):
		self.code_name = 'CODE'
		self.viking_modules = viking_modules
		self.archer2_modules = archer2_modules
		self.template = "template.CODE"
		self.fin_file = "run.fin"
		#from ..dimensions.dimensions_CODE import dimensions_list
		self.dim_list = None #dimensions_list
		if self.dim_list != None:
			self.make_dim_lookup(self.dim_list)
	
	def pass_dependencies(self, scanner):
		self.eqbm = scanner.equilibrium
		self.inputs = scanner.inputs
		self.scanner = self
	
	def make_dim_lookup(self, dim_list):
		dim_lookup = {'_list': [], '_full_list': []}
		for dim in dim_list:
			dim_lookup['_list'].append(dim.name_keys[0].lower())
			for dim_name in dim.name_keys:
				dim_name = dim_name.lower()
				dim_lookup[dim_name] = dim
				dim_lookup['_full_list'].append(dim_name)
		self.dim_lookup = dim_lookup
	
	def check_scan(self, valid = True):
		return valid
	
	def get_template_lines(self):
		template_lines = None
		return template_lines
	
	def get_surface_input(self, psiN):
		self.eqbm.pyro.gk_code = self.code_name
		self.eqbm.pyro.update_gk_code()
		self.eqbm.pyro.load_local(psi_n=psiN)
		self.eqbm.pyro.update_gk_code()
		from copy import deepcopy
		nml = deepcopy(self.eqbm.pyro.gk_input.data)
		for dim in self.inputs.single_parameters.values():
			nml = dim.single_edit_nml(nml)
		nml = self._get_surface_input(nml)
		for dim in self.inputs.single_parameters.values():
			nml = dim.single_edit_nml(nml)
		return nml
		
	def _get_surface_input(self, nml):
		return nml
		
	def get_gyro_input(self, run = None, indexes = None, namelist_diff = {}):
		if run is None and indexes is None:
			print("ERROR: Either indexes or run must be given")
			return None
		
		if run is None:
			if len(indexes) != len(self.inputs.dimensions):
				print(f"ERROR: indexes must be of length {len(self.inputs.dimensions)}, {[self.inputs.dim_order]}")
			run = {}
			for i, dim in zip(indexes,self.inputs.dimensions.values()):
				run[dim.name] = dim.values[i]
		
		if 'psin' in run:	
			psiN = run['psin']
		elif 'psin' in self.inputs.single_parameters:
			psiN = self.inputs.single_parameters['psin'].values[0]
		else:
			print("ERROR: psiN not defined")
			return None
		
		nml = self.get_surface_input(psiN)
		
		for dim_name, dim in self.inputs.dimensions.items():
			nml = dim.edit_nml(nml=nml,val=run[dim_name])
			
		for dim_name, dim in self.inputs.single_parameters.items():
			nml = dim.single_edit_nml(nml)
		
		nml = self._get_gyro_input(run, nml)
		
		for dim_name, dim in self.inputs.dimensions.items():
			nml = dim.edit_nml(nml=nml,val=run[dim_name])
			
		for dim_name, dim in self.inputs.single_parameters.items():
			nml = dim.single_edit_nml(nml)
		
		for key in namelist_diff.keys():
			for skey in namelist_diff[key].keys():
				nml[key][skey] = namelist_diff[key][skey]
		
		return nml
	
	def _get_gyro_input(self, run, nml):
		return nml

	def make_job_files_ypi(self):
		print(f"ERROR: {self.code_name} DOES NOT SUPPORT YPI SERVERS")
		return
	
	jobfile_viking = None
	
	def write_pyth_archer2(self, input_lists):
		print(f"ERROR: {self.code_name} DOES NOT SUPPORT ARCHER2")
		return
	
	def get_non_linear_archer2(self):
		print(f"ERROR: {self.code_name} DOES NOT SUPPORT NON LINEAR ARCHER2")
		run_code = None
		return run_code
		
	def make_gyro_file(self, run, directory):
		input_list_entry = None
		return input_list_entry
	
	def save_out(self, filename = None, directory = None, specificRuns = None, SlurmSave = False, QuickSave = False):
		return
	
	def write_nml(nml, directory = ".", filename = None):
		return
		
		
		
		
