import os
import f90nml

possible_keys = {'miller': ['miller','mill'],
		'ideal': ['ideal'],
		'gyro': ['gyro'],
		'system': ['system','sys','server'],
		'fixed_delt': ['fixed_delt','delt','delt','fix_delt'],
		'epar': ['epar','write_epar'],
		'num_shear_ideal': ['num_shat_ideal','n_shat_ideal','shat_ideal_num','shat_ideal_n','ideal_shat_n','idea_shat_num','num_shear_ideal','n_shear_ideal','shear_ideal_num','shear_ideal_n','ideal_shear_n','idea_shear_num']
		'num_beta_ideal': ['num_beta_ideal','n_beta_ideal','beta_ideal_num','beta_ideal_n','ideal_beta_n','idea_beta_num','num_beta_prime_ideal','n_beta_prime_ideal','beta_prime_ideal_num','beta_prime_ideal_n','ideal_beta_prime_n','idea_beta_prime_num']
		}

valid_systems = ['plasma','viking','archer']

class scan_inputs(object):
	def __init__(self, input_file = None, directory = "./", input_dict = None):
		from .templates import Dimensions
		self._create_empty_inputs()
		self.input_name = input_file
		if directory == "./":
			directory = os.getcwd()
		self.path = directory

		if input_file is not None:
			self.load_inputs()
		if input_dict is not None:
			self.load_input_dict(input_dict)

	def __getitem__(self, key):
		pkey, skey = self.find_key(key)
		if skey:
			return self.inputs[pkey][skey]
		elif pkey:
			return self.inputs[pkey]

		return None

	def find_key(self, key):
		key = key.lower()
		if key in self.inputs.keys():
			return key, None
		if key in ['psi','psin']:
			return 'psiN', None
		if key == 'ky':
			return 'aky', None

		for pkey in possible_keys:
			if type(possible_keys[pkey]) == dict:
				for skey in possible_keys[pkey]:
					if key in possible_keys[pkey][skey]:
						return pkey, skey
			else:
				if key in possible_keys[pkey]:
					return pkey, None

		print(f"ERROR: {key} not found")
		return None, None

	def print_inputs(self):
		for key in self.inputs:
			if type(self.inputs[key]) == dict:
				print(f"{key}:")
				for sub_key in self.inputs[key]:
					print(f"\t{sub_key} = {self.inputs[key][sub_key]}")
			else:
				print(f"{key} = {self.inputs[key]}")

	def keys(self):
		return self.inputs.keys()

	def _create_empty_inputs(self):
		from .templates import template_dir, inputs_template
		self.inputs = f90nml.Namelist({'knobs':
           {'gyro': True,
			'ideal': True,
            'miller': False,
            'system': 'viking',
            'fixed_delt': False,
            'epar': False,
            'num_shear_ideal': None,
            'num_beta_ideal': None}})

	def edit_inputs(self, key = None, val = 'None'):
		if key and val == 'None':
			print(f"ERROR: please enter val for {key}")
			return
		if key in ['beta_div','beta_mul','shat_div','shat_mul'] and val is None:
			return
		elif key:
			if key.lower() == 'viking':
				key = 'system'
				if val == True:
					val = 'viking'
				elif val == False:
					val = 'plasma'
			pkey, skey = self.find_key(key)
			if skey:
				self.inputs[pkey][skey] = val
			elif pkey:
				self.inputs[pkey] = val
			return

		for key in self.inputs:
			if type(self[key]) == dict:
				for sub_key in self[key]:
					if self[key][sub_key] is None:
						val = input(f"Input value for {key}_{sub_key} (Current value: None): ")
					else:
						val = input(f"Input value for {key}_{sub_key} (Current value: {self[key][sub_key]}): ")
					if sub_key == 'value' and val != "":
						if val[0] != "[":
							val = "[" + val
						if val[-1] != "]":
							val = val+ "]"
					if val != "":
						self.inputs[key][sub_key] = eval(val)
			else:
				if self[key] is None:
					val = input(f"Input value for {key} (Current value: None): ")
				else:
					val = input(f"Input value for {key} (Current value: {self[key]}): ")
				if key == 'System':
					val = val.lower()
					if val != "" and val not in valid_systems:
						print(f"ERROR: {val} is not a valid system, supported: {valid_systems}")
						val = ""
					if val != "":
						self.inputs[key] = val
				if val != "":
						self.inputs[key] = eval(val)

		self.validate_inputs()

	def load_input_dict(self, dic):
		for key in dic:
			if type(dic[key]) == dict:
				for sub_key in dic[key]:
					self.edit_inputs(key = f"{key}_{sub_key}", val = dic[key][sub_key])
			else:
				self.edit_inputs(key = key, val = dic[key])

	def validate_inputs(self):
		valid = True
		if self['gyro']:
			for key in ['psiN','shat','beta','aky','theta0']:
				if all([x is None for x in self[key].values()]):
					print(f"ERROR: all {key} values are None")
					valid = False
				elif self[key]['values'] is None and (self[key]['num'] == 1 or (self[key]['min'] == self[key]['max'] != None)):
					if self[key]['min']:
						self[key]['values'] = [self[key]['min']]
					elif self[key]['max']:
						self[key]['values'] = [self[key]['max']]
					else:
						print(f"ERROR: too many {key} values are None")
						valid = False
				elif self[key]['values'] is None and any([self[key]['min'] is None,self[key]['max'] is None,self[key]['num'] is None]):
					print(f"ERROR: too many {key} values are None")
					valid = False
				else:
					if type(self[key]['values']) in [int,float]:
						self.inputs[key]['values'] = [self[key]['values']]

					if self[key]['values'] is None:
						vals = []
						for i in range(self[key]['num']):
							vals.append((self[key]['max'] - self[key]['min'])*i/(self[key]['num']-1) + self[key]['min'])
						self.inputs[key]['values'] = vals
					if key == 'beta':
						self.inputs[key]['values'] = [abs(x) for x in self.inputs[key]['values']]

					self.inputs[key]['values'].sort()
					self.inputs[key]['min'] = min(self[key]['values'])
					self.inputs[key]['max'] = max(self[key]['values'])
					self.inputs[key]['num'] = len(self[key]['values'])

			if self['shat']['min'] < 1e-4:
				self.inputs['shat']['min'] = 1e-4
				self.inputs['shat']['values'][0] = 1e-4

		return valid

		if self['ideal']:
			if self['n_shat_ideal'] is None and self['n_shat']:
				self.inputs['shat']['n_shat_ideal'] = self['n_shat']
			elif self['n_shat_ideal'] is None:
				print("ERROR: n_shat_ideal is None")
				valid = False
			if self['n_beta_ideal'] is None and self['n_beta']:
				self.inputs['beta']['n_beta_ideal'] = self['n_beta']
			elif self['n_beta_ideal'] is None:
				print("ERROR: n_beta_ideal is None")
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

		inputs = f90nml.read(f"{directory}/{filename}")



		self.validate_inputs()

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

		with open(os.path.join(directory,self.input_name),'w') as in_file:
			for key in self.inputs.keys():
				if type(self.inputs[key]) == dict:
					for skey in self.inputs[key]:
						in_file.write(f"{key}_{skey} = {self.inputs[key][skey]}\n")
				else:
					in_file.write(f"{key} = {self.inputs[key]}\n")
		if doPrint:
			print(f"Created {filename} at {directory}")

class dimension(object):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		self.values = values
		self.edit_dimension(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = []

	def __get__(self, val):
		if val in ['min']:
			return self.min
		if val in ['max']:
			return self.max
		if val in ['num','n','len']:
			return self.len
		return self.values[val]

	def sub_validate(self, values):
		return values

	def edit_nml(self, nml, val):
		return nml
	
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
	def len(self):
		return len(self.values)


	def edit_dimension(self, values, mini, maxi, num):
		if all([x is None for x in values]):
			print(f"ERROR: all {self.name_keys[0]} values are None")
			self.values = []
			return
		elif values is None and (num == 1 or (mini == maxi != None)):
			if mini:
				values = [mini]
				maxi = mini
			elif maxi:
				values = [maxi]
				maxi = mini
			else:
				print(f"ERROR: too many {self.name_keys[0]} values are None")
				self.values = []
				return
		elif values is None and any([mini is None, maxi is None, num is None]):
			print(f"ERROR: too many {self.name_keys[0]} values are None")
			return
		else:
			if type(values) in [int,float]:
				self.values = [values]

			if values is None:
				vals = []
				for i in range(num):
					vals.append((maxi - mini)*i/(num-1) + mini)
				self.values = vals

		self.values.sort()
		self.values = self.sub_validate(values)