import os
import f90nml

possible_inputs = {
	'knobs': ['gyro','ideal','miller','system','fixed_delt','epar','num_shear_ideal','num_beta_ideal'],
	'dimension_n': ['type','values','min','max','num'],
	}

possible_keys = {'knobs': {
	'miller': ['miller','mill'],
	'ideal': ['ideal'],
	'gyro': ['gyro'],
	'system': ['system','sys','server'],
	'fixed_delt': ['fixed_delt','delt','delt','fix_delt'],
	'epar': ['epar','write_epar'],
	'num_shear_ideal': ['num_shat_ideal','n_shat_ideal','shat_ideal_num','shat_ideal_n','ideal_shat_n','idea_shat_num','num_shear_ideal','n_shear_ideal','shear_ideal_num','shear_ideal_n','ideal_shear_n','idea_shear_num'],
	'num_beta_ideal': ['num_beta_ideal','n_beta_ideal','beta_ideal_num','beta_ideal_n','ideal_beta_n','idea_beta_num','num_beta_prime_ideal','n_beta_prime_ideal','beta_prime_ideal_num','beta_prime_ideal_n','ideal_beta_prime_n','ideal_beta_prime_num'],
	}}

valid_systems = ['ypi_server','viking','archer2']

defaults = {'knobs': {
	'gyro': True,
	'ideal': True,
	'miller': False,
	'system': 'viking',
	'fixed_delt': False,
	'epar': False,
	'num_shear_ideal': None,
	'num_beta_ideal': None,
	}}

class scan_inputs(object):
	def __init__(self, input_file = None, directory = "./", input_dict = None):
		self.input_name = input_file
		if directory == "./":
			directory = os.getcwd()
		self.path = directory
		
		self.inputs = self.dimensions = self.dim_order = None
		
		if input_file is not None:
			self.load_inputs()
		if input_dict is not None:
			self.load_input_dict(input_dict)

	def __getitem__(self, key):
		key = key.lower()
		
		if key in self.inputs.keys():
			return self.inputs[key]

		for pkey in possible_keys:
			for skey in possible_keys[pkey]:
				if key in possible_keys[pkey][skey]:
					for sskey in possible_keys[pkey][skey]:
						if sskey in self.inputs[pkey]:	
							return self.inputs[pkey][sskey]
					if pkey in defaults:
						return defaults[pkey][skey]
		
		if self.dimensions:
			if key in ['order','dim_oder','dimension_order','dimensions_order']:
				return self.dim_order
			for dim in self.dimensions:
				for skey in self.dimensions[dim].possible_keys:
					if key in self.dimensions[dim].possible_keys[skey]:
						return self.dimensions[dim][skey]

		print(f"ERROR: {key} not found")
		return None
		
	def print_inputs(self):
		if inputs is None:
			return
		for key in self.inputs:
			print(f"{key}:")
			for sub_key in self.inputs[key]:
				print(f"\t{sub_key} = {self.inputs[key][sub_key]}")

	def keys(self):
		return self.inputs.keys()
	
	@property
	def shape(self):
		if not self.dimensions:
			return []
		sh = [None] * len(self.dimensions)
		for dim_id, dim_type in enumerate(self.dim_order):
			sh[dim_id] = len(self.dimensions[dim_type])
		return sh
	
	def load_input_dict(self, dic):
		if not self.inputs:
			self.inputs = f90nml.Namelist()
		for key in dic:
			if key not in self.inputs.keys():
				self.inputs[key] = {}
			for sub_key in dic[key]:
				self.inputs[key][sub_key] = dic[key][sub_key]
		self.check_inputs()
		self.load_dimensions()

	def check_inputs(self):
		for key in ['knobs']:
			for skey in self.inputs[key]:
				if skey not in possible_inputs[key]:
					print(f"ERROR: {skey} is not a valid {key} input")
					del(self.inputs[key][skey])
		for key in self.inputs['single_parameters']:
			from .templates import dim_lookup
			if key not in dim_lookup:
				print(f"ERROR: {key} is not a valid parameter, valid parameters: {dim_lookup['_list']}")
				del(self.inputs['single_parameters'][key])
		for key in [x for x in self.inputs if 'dimension_' in x]:
			for skey in self.inputs[key]:
				if skey not in possible_inputs['dimension_n']:
					print(f"ERROR: {skey} is not a valid {key} input")
					del(self.inputs[key][skey])
			for skey in ['type','values','min','max','num']:
				if skey not in self.inputs[key]:
					self.inputs[key][skey] = None
			if type(self.inputs[key]['values']) == list:
				if self.inputs[key]['values'][0] == '[':
					self.inputs[key]['values'] = self.inputs[key]['values'][1:]
				if self.inputs[key]['values'][-1] == ']':
					self.inputs[key]['values'] = self.inputs[key]['values'][:-1]
			if not self.inputs[key]['type']:
				print(f"ERROR: {key} has no type, valid types: {dim_lookup['_list']}")
	
	def check_scan(self):
		valid = True
		if not self.inputs:
			print("ERROR: No inputs loaded")
			return False
		if not self.dimensions:
			self.load_dimensions()
			
		for dim in self.dimensions:
			if self.dimensions[dim].values is None:
				print(f"ERROR: {dim} dimensions invalid")
				valid = False
			
		if 'psin' not in self.dimensions:
			print("ERROR: psiN must be specified as single parameter or scan dimension")
			valid = False

		if self['ideal']:
			if self['num_shear_ideal'] is None and self['n_shat']:
				self.inputs['knobs']['num_shear_ideal'] = self['n_shat']
			elif self['num_shear_ideal'] is None:
				print("ERROR: num_shear_ideal is None")
				valid = False
			if self['num_beta_ideal'] is None and self['n_beta']:
				self.inputs['knobs']['num_beta_ideal'] = self['n_beta']
			elif self['num_beta_ideal'] is None:
				print("ERROR: num_beta_ideal is None")
				valid = False
			if 'beta_prime' not in self.dimensions:
				print("ERROR: beta_prime dimension must be given for ideal scan")
				valid = False
			if 'shear' not in self.dimensions:
				print("ERROR: shear dimension must be given for ideal scan")
				valid = False

		return valid

	def load_inputs(self, filename = None, directory = "./"):
		if self.input_name is None and filename is None:
			filename = input("Input File Name: ")
			if "." not in filename:
				filename = filename + ".in"
		elif self.input_name is None:
			if "." not in filename:
				filename = filename + ".in"
			self.input_name = filename

		self.inputs = f90nml.read(f"{directory}/{self.input_name}")
		self.check_inputs()
		self.load_dimensions()
	
	def load_dimensions(self):
		from .templates import dim_lookup
		dimensions = {}
		single_parameters = {}
		dim_order = []
		for key in [x for x in self.inputs.keys() if 'dimension_' in x]:
			dim_type = self.inputs[key]['type']
			dim_type = dim_type.lower()
			if dim_type and dim_type not in dim_lookup['_full_list']:
				print(f"ERROR: {dim_type} not a valid dimension. Valid = {dim_lookup['_list']}")
			elif dim_type:
				dim = dim_lookup[dim_type](values=self.inputs[key]['values'],mini=self.inputs[key]['min'],maxi=self.inputs[key]['max'],num=self.inputs[key]['num'])
				if dim.name in dimensions:
					print(f"ERROR: {dim_type} defined multiple times")
				else:
					if len(dim) == 1:
						single_parameters[dim.name] = dim
					else:
						dimensions[dim.name] = dim
					
		if 'single_parameters' in self.inputs.keys():
			for dim_type in self.inputs['single_parameters']:
				if dim_type not in dim_lookup['_full_list']:
					print(f"ERROR: {dim_type} not a valid parameter. Must be a dimension, valid = {dim_lookup['_list']}")
				elif dim_type in dimensions:
					print(f"ERROR: {dim_type} defined multiple times. As single parameter and dimension")
				elif dim_type in single_parameters:
					print(f"ERROR: {dim_type} defined multiple times.")
				else:
					dim = dim_lookup[dim_type](values=self.inputs['single_parameters'][dim_type])
					single_parameters[dim.name] = dim
		
		for key in [x for x in self.inputs if 'dimension_' in x]:
			del self.inputs[key]
		self.inputs['single_parameters'] = {}
		
		for dim_type in single_parameters:
			self.inputs['single_parameters'][dim_type] = single_parameters[dim_type].values[0]
				
		for dim in dimensions.values():
			dim_order.append(dim.name)
			dim_id = dim_order.index(dim.name)
			self.inputs[f"dimension_{dim_id+1}"] = {}
			self.inputs[f"dimension_{dim_id+1}"]['type'] = dim.name
			self.inputs[f"dimension_{dim_id+1}"]['values'] = dim.values
			self.inputs[f"dimension_{dim_id+1}"]['min'] = dim.min
			self.inputs[f"dimension_{dim_id+1}"]['max'] = dim.max
			self.inputs[f"dimension_{dim_id+1}"]['num'] = len(dim)
		
		self.dimensions = dimensions
		self.single_parameters = single_parameters
		self.dim_order = dim_order

	def write_scan_input(self, filename = None, directory = "./", doPrint = False):
		if directory is None and self.directory is None:
			directory = "./"
		elif directory is None:
			directory = self.path

		if self.input_name is None and filename is None:
			filename = input("Input File Name: ")
		elif self.input_name is None:
			self.input_name = filename
		if "." not in self.input_name:
			self.input_name = self.input_name + ".in"

		self.inputs.write(f"{directory}/{filename}",force=True)
		
		if doPrint:
			print(f"Created {filename} at {directory}")
			
	def write_default_input(self, filename = None, directory = "./", doPrint = False):
		from .templates import template_dir, inputs_template
		if directory is None and self.directory is None:
			directory = "./"
		elif directory is None:
			directory = self.path

		if self.input_name is None and filename is None:
			filename = "default_input.in"
		elif self.input_name is None:
			self.input_name = filename
		if "." not in self.input_name:
			self.input_name = self.input_name + ".in"
		
		default_in = f90nml.load(f"{template_dir}/{inputs_template}")
		default_in.write(f"{directory}/{filename}",force=True)
		
		if doPrint:
			print(f"Created {filename} at {directory}")

class dimension(object):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		self.values = values
		self.edit_dimension(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = []

	def __getitem__(self, key):
		if type(key) == int:
			return self.values[key]
		key = key.lower()
		if key in ['min']:
			return self.min
		if key in ['max']:
			return self.max
		if key in ['num','n','len']:
			return len(self)
		print(f"{key} not found")
		return None
		

	def sub_validate(self, values):
		return values

	def edit_nml(self, nml, val):
		return nml
		
	def single_edit_nml(self, nml):
		if len(self) != 1:
			print("ERROR: single_edit_nml can only be used if dimension length is 1")
			return None
		return self.edit_nml(nml = nml, val = self.values[0])
	
	@property
	def possible_keys(self):
		keys = {'min': [],'max': [],'num': [],'values': []}
		for name in self.name_keys:
			name = name.lower()
			keys['min'].extend([f'{name}_min',f'min_{name}'])
			keys['max'].extend([f'{name}_max',f'max_{name}'])
			keys['num'].extend([f'{name}_num',f'num_{name}',f'n_{name}',f'{name}_n',f'{name}_number',f'number_{name}',f'len_{name}',f'psi_{name}'])
			keys['values'].extend([f'{name}s',f'{name}_values',f'values_{name}',f'{name}_v',f'{name}_vals'])
		return keys

	@property
	def min(self):
		return min(self.values)

	@property
	def max(self):
		return max(self.values)
		
	@property
	def num(self):
		return len(self.values)
		
	def __len__(self):
		return len(self.values)
		
	@property
	def name(self):
		return self.name_keys[0].lower()


	def edit_dimension(self, values = None, mini = None, maxi = None, num = None):
		if all([x is None for x in [values, mini, maxi, num]]):
			print(f"ERROR: all {self.name} values are None")
			return
		elif values is None and (num == 1 or (mini == maxi != None)):
			if mini:
				values = [mini]
				maxi = mini
			elif maxi:
				values = [maxi]
				maxi = mini
			else:
				print(f"ERROR: too many {self.name} values are None")
				return
		elif values is None and any([mini is None, maxi is None, num is None]):
			print(f"ERROR: too many {self.name} values are None")
			return
		else:
			if type(values) in [int,float]:
				self.values = [values]

			if values is None:
				vals = []
				for i in range(num):
					vals.append((maxi - mini)*i/(num-1) + mini)
				self.values = vals
				
		if self.values is not None:
			self.values.sort()
			self.values = self.sub_validate(self.values)
