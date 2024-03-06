import os
import f90nml
from copy import deepcopy, copy
from .templates import dim_lookup, template_dir, gs2_template, inputs_template, systems

possible_keys = {
	'files': {
	'eq_name': ['eq_name','eq','equilibrium','eq_file'],
	'eq_path': ['eq_path','eq_dir','equilibrium_path','equilibrium_dir','eq_directory'],
	'kin_name': ['kin_name','kin','kinetic','kin_file'],
	'kin_path': ['kin_path','kin_dir','kinetic_path','kinetic_dir','kin_directory'],
	'kin_type': ['kin_type','kinetics_type'],
	'template_name': ['template_name','template','gk_file','template_file'],
	'template_path': ['template_path','template_dir','template_directory'],
	'input_name': ['input_name','input','input_file'],
	'input_path': ['input_path','input_dir','input_directory'],
	},
	'knobs': {
	'miller': ['miller','mill'],
	'ideal': ['ideal'],
	'gyro': ['gyro'],
	'system': ['system','sys','server'],
	'grid_option': ['grid_option', 'grid'],
	'nonlinear': ['nonlinear','non_linear','nl','nonlin','non_lin'],
	'split_nonlinear' : ['split_nonlinear','split_non_linear','split_nl','split_nonlin','split_non_lin'],
	'fixed_delt': ['fixed_delt','delt','delt','fix_delt'],
	'epar': ['epar','write_epar'],
	'num_shear_ideal': ['num_shat_ideal','n_shat_ideal','shat_ideal_num','shat_ideal_n','ideal_shat_n','idea_shat_num','num_shear_ideal','n_shear_ideal','shear_ideal_num','shear_ideal_n','ideal_shear_n','idea_shear_num'],
	'num_beta_ideal': ['num_beta_ideal','n_beta_ideal','beta_ideal_num','beta_ideal_n','ideal_beta_n','idea_beta_num','num_beta_prime_ideal','n_beta_prime_ideal','beta_prime_ideal_num','beta_prime_ideal_n','ideal_beta_prime_n','ideal_beta_prime_num'],
	'force_zero_fs': ['force_zero_fs','force_zero_flow_shear','force_0_fs','force_0_flow_shear']
	},
	'info': {
	'run_name': ['run_name','name'],
	'run_uuid': ['run_uuid','uuid','run_id','id'],
	'run_date': ['run_date','date','time','datetime','timedate'],
	'data_path': ['data_path'],
	'itteration': ['itteration','itt'],
	'description': ['description','desc','notes'],
	},
	'dimension_n': {
	'type': ['type'],
	'values': ['values','vals'],
	'min': ['min','minimum'],
	'max': ['max','maximum'],
	'num': ['num','number','len'],
	'option': ['option','opt','options'],
	},
	}

valid_systems = ['ypi_server','viking','archer2']

default_inputs = {'files': {
	'eq_name': None,
	'eq_path': None,
	'kin_name': None,
	'kin_path': None,
	'kin_type': 'PEQDSK',
	'template_name': None,
	'template_path': None,
	'input_name': None,
	'input_path': None,
	},
	'knobs': {
	'gyro': True,
	'ideal': True,
	'miller': False,
	'system': 'viking',
	'grid_option': 'single',
	'nonlinear': False,
	'split_nonlinear': True,
	'fixed_delt': False,
	'epar': False,
	'num_shear_ideal': None,
	'num_beta_ideal': None,
	'force_zero_fs': True,
	},
	'info': {
	'run_name': None,
	'run_uuid': None,
	'run_date': None,
	'data_path': None,
	'itteration': 0,
	'description': '',
	},
	'single_parameters': {},
	}

class scan_inputs(object):
	def __init__(self, input_file = None, directory = "./", input_dict = None):
		if directory == "./":
			directory = os.getcwd()
		if input_file is None and input_dict is None:
			input_file = inputs_template
			directory = template_dir
		self.path = directory
		self.input_name = input_file
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
		
		if self.dimensions:
			if key in ['order','dim_order','dimension_order','dimensions_order']:
				return self.dim_order
			for dim in self.dimensions:
				for skey in self.dimensions[dim].possible_keys:
					if key in self.dimensions[dim].possible_keys[skey]:
						return self.dimensions[dim][skey]

		print(f"ERROR: {key} not found")
		return None
		
	def print_inputs(self):
		if self.inputs is None:
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
	
	def load_inputs(self, filename = None, directory = None):
		if self.input_name is None and filename is None:
			filename = input("Input File Name: ")
			if "." not in filename:
				filename = filename + ".in"
		elif self.input_name is None:
			if "." not in filename:
				filename = filename + ".in"
			self.input_name = filename
		
		if directory is None and self.path is None:
			directory = os.getcwd()
		elif directory is None:
			directory = self.path

		self.inputs = f90nml.read(f"{directory}/{self.input_name}")
		self.check_inputs()
		self.load_dimensions()
	
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
		defaults = deepcopy(default_inputs)
		for key in defaults:
			if key not in self.inputs:
				self.inputs[key] = {}
			for skey in defaults[key]:
				if skey not in self.inputs[key]:
					self.inputs[key][skey] = defaults[key][skey]
					
		
		for key in ['knobs','files','info']:
			toDelete = []
			toReplace = []
			for skey in self.inputs[key].keys():
				if skey not in default_inputs[key].keys():
					toDelete.append(skey)
					found = False
					for pkey in possible_keys[key].keys():
						if skey in possible_keys[key][pkey]:
							toReplace.append(pkey)
							found = True
					if not found:
						toReplace.append(None)
						
			for old_key, new_key in zip(toDelete,toReplace):
				if new_key is not None:
					self.inputs[key][new_key] = deepcopy(self.inputs[key][old_key])
				else:
					print(f"ERROR: {skey} is not a valid {key} input")
				del(self.inputs[key][old_key])
			
		if self.inputs['files']['input_name'] is None and self.input_name:
			self.inputs['files']['input_name'] = self.input_name
		if self.inputs['files']['input_path'] is None and self.inputs['files']['input_name']:
			self.inputs['files']['input_path'] = self.path
		
		if not self.inputs['files']['template_name']:
			self.inputs['files']['template_name'] = gs2_template
			self.inputs['files']['template_path'] = template_dir
		for key in ['eq','kin']:
			if not self.inputs['files'][f'{key}_name']:
				print(f"ERROR: No {key} file given")
			elif not self.inputs['files'][f'{key}_path']:
				if '/' in self.inputs['files'][f'{key}_name']:
					full = self.inputs['files'][f'{key}_name']
					self.inputs['files'][f'{key}_name'] = full.split('/')[-1]
					self.inputs['files'][f'{key}_path'] = full[:-len(self.inputs['files'][f'{key}_name'])]
				else:
					self.inputs['files'][f'{key}_path'] = self.path
		
		if self['knobs']['system'] in ['archer2','viking']:
			sbatch = copy(systems[self['knobs']['system']]['sbatch'])
			if 'sbatch' not in self.inputs:
				self.inputs['sbatch'] = {}
			for skey in sbatch:
				if skey not in self.inputs['sbatch']:
					if skey == 'job-name':
						self.inputs['sbatch'][skey] = self.inputs['files']['input_name'].split('/')[-1].split('.')[0]
					if skey == 'output':
						self.inputs['sbatch'][skey] = self.inputs['files']['input_name'].split('/')[-1].split('.')[0] + ".slurm"
					if skey == 'error':
						self.inputs['sbatch'][skey] = self.inputs['files']['input_name'].split('/')[-1].split('.')[0] + ".err"
					else:
						self.inputs['sbatch'][skey] = sbatch[skey]
			
			sbatch_save = copy(systems[self['knobs']['system']]['save_sbatch'])
			if 'sbatch_save' not in self.inputs:
				self.inputs['sbatch_save'] = {}
			for skey in sbatch_save:
				if skey not in self.inputs['sbatch_save']:
					if skey == 'job-name':
						self.inputs['sbatch_save'][skey] = self.inputs['files']['input_name'].split('/')[-1].split('.')[0]
					if skey == 'output':
						self.inputs['sbatch_save'][skey] = "save_out.slurm"
					else:
						self.inputs['sbatch_save'][skey] = sbatch_save[skey]
					
		for key in self.inputs['single_parameters']:
			if key not in dim_lookup:
				print(f"ERROR: {key} is not a valid parameter, valid parameters: {dim_lookup['_list']}")
				del(self.inputs['single_parameters'][key])
		for key in [x for x in self.inputs if 'dimension_' in x]:
			for skey in self.inputs[key]:
				if skey not in possible_keys['dimension_n']:
					print(f"ERROR: {skey} is not a valid {key} input")
					del(self.inputs[key][skey])
			for skey in ['type','values','min','max','num','option']:
				if skey not in self.inputs[key]:
					self.inputs[key][skey] = None
			if type(self.inputs[key]['values']) == list:
				if self.inputs[key]['values'][0] == '[':
					self.inputs[key]['values'] = self.inputs[key]['values'][1:]
				if self.inputs[key]['values'][-1] == ']':
					self.inputs[key]['values'] = self.inputs[key]['values'][:-1]
			if not self.inputs[key]['type']:
				print(f"ERROR: {key} has no type, valid types: {dim_lookup['_list']}")
			
			if self['nonlinear'] == True:
				self.inputs['knobs']['grid_option'] = 'box'
	
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
			
		if 'psin' not in self.dimensions and 'psin' not in self.single_parameters:
			print("ERROR: psiN must be specified as single parameter or scan dimension")
			valid = False
		
		if self['grid_option'] == 'box':
			if 'ky' in self.dimensions or 'ky' in self.single_parameters:
				print("ERROR: ky dimension is not compatible with grid_option == box: use nx, ny, y0 and jtwist")
				valid = False
			if 'kx' in self.dimensions or 'kx' in self.single_parameters:
				print("ERROR: kx dimension is not compatible with grid_option == box: use nx, ny, y0 and jtwist")
				valid = False
			if 'theta0' in self.dimensions or 'theta0' in self.single_parameters:
				print("ERROR: theta0 dimension is not compatible with grid_option == box: use nx, ny, y0 and jtwist")
				valid = False
			if not self['fixed_delt']:
				print("grid_option == box runs can only use fixed_delt == True")
				self.inputs['knobs']['fixed_delt'] = True
		
		if self['grid_option'] == 'single':
			if 'ny' in self.dimensions or 'ny' in self.single_parameters:
				print("ERROR: ny dimension is not compatible with grid_option == single: use ky")
				valid = False
			if 'nx' in self.dimensions or 'nx' in self.single_parameters:
				print("ERROR: nx dimension is not compatible with grid_option == single: use ky")
				valid = False
			if 'y0' in self.dimensions or 'y0' in self.single_parameters:
				print("ERROR: nx dimension is not compatible with grid_option == single")
				valid = False
			if 'jtwist' in self.dimensions or 'jtwist' in self.single_parameters:
				print("ERROR: nx dimension is not compatible with grid_option == single")
				valid = False
		
		if ('kx' in self.dimensions or 'kx' in self.single_parameters) and ('theta0' in self.dimensions or 'theta0' in self.single_parameters):
			print("ERROR: cannot define both kx and theta0 dimensions")
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
			if self['nonlinear'] == True:
					print("ERROR: cannot run nonlinear and ideal simultaneously")
					valid = False
			
		if self['nonlinear'] == True and len(self.dimensions) > 0:
			print("ERROR: dimensional scans not currently allowed for Nonlinear runs")
			valid = False

		return valid
		
	def create_run_info(self):
		if self.inputs['info']['run_uuid'] is None:
			try:
				from uuid import uuid4
				self.inputs['info']['run_uuid'] = str(uuid4())
			except:
				print("ERROR: unable to import uuid module, setting ID to None")
				self.inputs['info']['run_uuid'] = None
		if self.inputs['info']['run_name'] is None:
			if self.inputs['files']['input_name']:
				self.inputs['info']['run_name'] = self.inputs['files']['input_name'].split("/")[-1].split(".")[0]
		if self.inputs['info']['data_path'] == './':
			self.inputs['info']['data_path'] = os.getcwd()
		if self.inputs['info']['data_path'] is None:
			path = self.path
			if path[:2] == "./":
				path = os.getcwd() + path[1:]
			if self.inputs['info']['run_name']:
				self.inputs['info']['data_path'] = os.path.join(path, self.inputs['info']['run_name'])
			elif self.inputs['info']['run_uuid']:
				self.inputs['info']['data_path'] = os.path.join(path, self.inputs['info']['run_uuid'])
			else:
				self.inputs['info']['data_path'] = os.path.join(path, "myro_run")
			
		if self.inputs['info']['run_date'] is None:
			try:
				from datetime import datetime as dt
				self.inputs['info']['run_date'] = dt.now().strftime("%d/%m/%Y, %H:%M:%S")
			except:
				print("ERROR: unable to import datetime module, setting run date to None")
				self.inputs['info']['run_date'] = None
		if self.inputs['info']['itteration'] is None:
			self.inputs['info']['itteration'] = 0
	
	def load_dimensions(self):
		dimensions = {}
		single_parameters = {}
		dim_order = []
		for key in [x for x in self.inputs.keys() if 'dimension_' in x]:
			dim_type = self.inputs[key]['type']
			dim_type = dim_type.lower()
			if dim_type and dim_type not in dim_lookup['_full_list']:
				print(f"ERROR: {dim_type} not a valid dimension. Valid = {dim_lookup['_list']}")
			elif dim_type:
				dim = dim_lookup[dim_type](values=self.inputs[key]['values'],mini=self.inputs[key]['min'],maxi=self.inputs[key]['max'],num=self.inputs[key]['num'],option=self.inputs[key]['option'])
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
					print(f"ERROR: {dim_type} not a valid parameter. Must be a dimension, valid = {dim_lookup['_list'].keys()}")
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
			self.inputs[f"dimension_{dim_id+1}"]['option'] = dim.option
		
		self.dimensions = dimensions
		self.single_parameters = single_parameters
		self.dim_order = dim_order

	def write_scan_input(self, filename = None, directory = "./", doPrint = True):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path

		if self.inputs['files']['input_name'] is None and filename is None:
			filename = input("Input File Name: ")
		if self.inputs['files']['input_name'] is None:
			self.inputs['files']['input_name'] = filename
		if filename is None:
			filename = self.inputs['files']['input_name']
		if "." not in filename:
			filename = filename + ".in"
		
		self.inputs.write(f"{directory}/{filename}",force=True)
		
		if doPrint:
			print(f"Created {filename} at {directory}")
			
	def write_blank_input(self, filename = None, directory = "./", doPrint = True):
		if directory is None and self.path is None:
			directory = "./"
		elif directory is None:
			directory = self.path

		if filename is None:
			filename = "default_input.in"
		if "." not in filename:
			filename = filename + ".in"
		
		default_in = f90nml.read(f"{template_dir}/{inputs_template}")
		if default_in['knobs']['system'] in systems and 'sbatch' in systems[default_in['knobs']['system']]:
			default_in['sbatch'] = systems[default_in['knobs']['system']]['sbatch']
		default_in.write(f"{directory}/{filename}",force=True)
		
		if doPrint:
			print(f"Created {filename} at {directory}")
