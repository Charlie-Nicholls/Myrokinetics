import os

possible_keys = {'psiN': {'min': ['psi_min','min_psi','psin_min','min_psin'],
			'max': ['psi_max','max_psi','psin_max','max_psin'],
			'num': ['psi_num','num_psi','n_psi','psi_n','psi_number','number_psi','len_psi','psi_len','psin_num','num_psin','n_psin','psin_n','psin_number','number_psin','len_psin','psin_len'],
			'values': ['psis','psi_values','values_psi','psi_v','psi_vals','psins','psin_values','values_psin','psin_v','psin_vals']},
		'shat': {'min': ['shat_min','min_shat'],
			'max': ['shat_max','max_shat'],
			'num': ['shat_num','num_shat','n_shat','shat_n','shat_number','number_shat','len_shat','shat_len'],
			'values': ['shats','shat_values','values_shat','shat_v','shat_vals'],
			'ideal_num': ['shat_num_ideal','num_shat_ideal','n_shat_ideal','shat_n_ideal','shat_number_ideal','number_shat_ideal','len_shat_ideal','shat_len_ideal']},
		'beta': {'min': ['beta_min','min_beta'],
			'max': ['beta_max','max_beta'],
			'num': ['beta_num','num_beta','n_beta','beta_n','beta_number','number_beta','len_beta','beta_len'],
			'values': ['betas','beta_values','values_beta','beta_v','beta_vals'],
			'ideal_num': ['beta_num_ideal','num_beta_ideal','n_beta_ideal','beta_n_ideal','beta_number_ideal','number_beta_ideal','len_beta_ideal','beta_len_ideal']},
		'aky': {'min': ['aky_min','min_aky'],
			'max': ['aky_max','max_aky'],
			'num': ['aky_num','num_aky','n_aky','aky_n','aky_number','number_aky','len_aky','aky_len','ky_num','num_ky','n_ky','ky_n','ky_number','number_ky','len_ky','ky_len'],
			'values': ['akys','aky_values','values_aky','aky_v','aky_vals']},
		'theta0': {'min': ['theta_min','min_theta','theta0_min','min_theta0'],
			'max': ['theta_max','max_theta','theta0_max','max_theta0'],
			'num': ['theta_num','num_theta','n_theta','theta_n','theta_number','number_theta','len_theta','theta_len','theta0_num','num_theta0','n_theta0','theta0_n','theta0_number','number_theta0','len_theta0','theta0_len'],
			'values': ['thetas','theta_values','values_theta','theta_v','theta_vals','theta0s','theta0_values','values_theta0','theta0_v','theta0_vals']},
		'Miller': ['miller','mill'],
		'Ideal': ['ideal'],
		'Gyro': ['gyro'],
		'System': ['system','sys','server'],
		'Fixed_delt': ['fixed_delt','delt','delt','fix_delt'],
		'Epar': ['epar','write_epar'],
		}

class inputs(object):
	def __init__(self, input_file = None, directory = "./", input_dict = None):
		self._create_empty_inputs()
		self.input_file = input_file
		self._valid_systems = ['plasma','viking','archer']
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
			return self.inputs['psiN']
		if key == 'ky':
			return self.inputs['aky']
			
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
        	for key, val in self.inputs.items():
        		print(f"{key} = {val}")
        	
	def keys(self):
		return self.inputs.keys()
		
	def _create_empty_inputs(self):
		self.inputs = {
		'psiN': {'min': None, 'max': None, 'num': None, 'values': None},
		'shat': {'min': None, 'max': None, 'num': None, 'values': None, 'ideal_num': None},
		'beta': {'min': None, 'max': None, 'num': None, 'values': None, 'ideal_num': None},
		'aky': {'min': None, 'max': None, 'num': None, 'values': None},
		'theta0': {'min': None, 'max': None, 'num': None, 'values': None},
		'Miller': False,
		'Ideal': True,
		'Gyro': True,
		'System': 'viking',
		'Fixed_delt': False,
		'Epar': False,
		}
	
	def edit_inputs(self, key = None, val = None):
		if key and not val:
			print(f"ERROR: please enter val for {key}")
			return
		elif key:
			pkey, skey = self.find_key(key)
			if skey:
				self.inputs[pkey][skey] = val
			elif pkey:
				self.inputs[pkey] = val
			else:
				print(f"ERROR: {key} not found")
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
					if val != "" and val not in self._valid_systems:
						print(f"ERROR: {val} is not a valid system, supported: {self._valid_systems}")
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
						self.inputs[key]['values'] = [abs(x) for x in vals]
						
					self.inputs[key]['values'].sort()
					self.inputs[key]['min'] = min(self[key]['values'])
					self.inputs[key]['max'] = max(self[key]['values'])
					self.inputs[key]['num'] = len(self[key]['values'])
			
			if self['shat']['min'] < 1e-4:
				self.inputs['shat']['min'] = 1e-4
				self.inputs['shat']['values'][0] = 1e-4
			
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
		if self.input_file is None and filename is None:
			filename = input("Input File Name: ")
			if "." not in filename:
				filename = filename + ".in"
		elif self.input_file is None:
			if "." not in filename:
				filename = filename + ".in"
			self.input_file = filename

		with open(os.path.join(self.path,self.input_file)) as in_file:
			lines = in_file.readlines()
			for line in lines:
				key = line.split(" = ")[0].strip("\t\n ")
				val = line.split(" = ")[1].strip("\t\n ")
				
				if key == 'Viking':
					if val == "True":
						key = "System"
						val = "viking"
					if val == "False":
						key = "System"
						val = "plasma"
				
				pkey, skey = self.find_key(key)
				
				if skey == 'value':
					if val[0] != "[":
						val = "[" + val
					if val[-1] != "]":
						val = val + "]"
						
				if pkey == 'System':
					self.inputs[pkey] = val
				elif skey:
					 self.inputs[pkey][skey] = eval(val)
				elif pkey:
					self.inputs[pkey] = eval(val)
		
		self.validate_inputs()
			
	def write_scan_input(self, filename = None, directory = "./", doPrint = False):
		if directory is None and self.directory is None:
			directory = "./"
		elif directory is None:
			directory = self.path
			
		if self.input_file is None and filename is None:
			filename = input("Input File Name: ")
		elif self.input_file is None:
			self.input_file = filename
		if "." not in self.input_file:
			self.input_file = self.input_file + ".in"
		
		with open(os.path.join(directory,self.input_file),'w') as in_file:
			for key in self.inputs.keys():
				if type(self.inputs[key]) == dict:
					for skey in self.inputs[key]:
						in_file.write(f"{key}_{skey} = {self.inputs[key][skey]}\n")
				else:
					in_file.write(f"{key} = {self.inputs[key]}\n")
		if doPrint:
			print(f"Created {filename} at {directory}")
	
	def print_inputs(self):
		for key in self.inputs:
			if type(dic[key]) == dict:
				print(f"{key}:")
				for sub_key in dic[key]:
					print(f"\t{sub_key} = {self.inputs[key][sub_key]}")
			else:
				print(f"\t{key} = {self.inputs[key]")
