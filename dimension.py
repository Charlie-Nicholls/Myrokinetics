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
