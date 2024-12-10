from ..dimension import dimension

'''
TEMPLATE DIMENSION SUBCLASS | Rememeber to add any new dimensions to dimensions_list
class DimensionType(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['']
	axis_label = ''
	valid_options = []

	def sub_validate(self, values):
		#Any constraints or corrections for your specific dimension
		return values

	def edit_nml(self, nml, val):
		#nml is a gs2 input file, edit and return it with changed dimension value = val
		return nml
'''

class psiN(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['psin','psi','psi_norm']
	axis_label = r'$\psi_{N}$'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 or x > 1 for x in values]):
			print("Error: psiN values outside allowed range (0<x<=1)")
			values = [x for x in values if (0<x<=1)]
		return values

	def edit_nml(self, nml, val):
		return nml

class p_prime(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)
		if self.option is None:
			self.option = 'both'

	name_keys = ['p_prime','p_prime_loc','pprime','pressure_prime','pressure_gradient']
	axis_label = 'p_prime'
	valid_options = ['both','beta']
	
	def sub_validate(self, values):
		values = [abs(x) for x in values]
		values.sort()
		return values

	def edit_nml(self, nml, val):
		from numpy import pi
		nml['p_prime_loc'] = -1*abs(val)
		beta = nml['betae']
		bp_cal = sum((nml[f'rlts_{n}'] + nml[f'rlns_{n}'])*nml[f'as_{n}']*nml[f'taus_{n}'] for n in [1,2,3])*beta
		
		mul = abs(val)/(bp_cal*nml['q_loc']/(nml['rmin_loc']*8*pi))
		if self.option == 'both':
			for n in [1,2,3]:
				nml[f'rlts_{n}'] = nml[f'rlts_{n}']*mul
				nml[f'rlns_{n}'] = nml[f'rlns_{n}']*mul
		elif self.option == 'beta':
			nml['betae'] = beta*mul
		return nml

class q_prime(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['q_prime','q_gradient','q_prime_loc']
	axis_label = 'q_prime'
	valid_options = []
	
	def sub_validate(self, values):
		if any([x<1e-4 for x in values]):
				values = list(set([x if x>1e-4 else 1e-4 for x in values]))
				values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['q_prime_loc'] = val
		return nml

class shear(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['shear','shat','s_hat','sh']
	axis_label = r'$\hat{s}$'
	valid_options = []
	
	def sub_validate(self, values):
		if any([x<1e-4 for x in values]):
				values = list(set([x if x>1e-4 else 1e-4 for x in values]))
				values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['q_prime_loc'] = val * (nml['q_loc']/nml['rmin_loc'])**2
		return nml

class ky_max(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ky_max']
	axis_label = r'$max(k_{y})$'
	valid_options = []

	def sub_validate(self, values):
		if any([x < 0 for x in values]):
			print("Error: ky_max values outside allowed range (x>=0)")
			values = [x for x in values if (x>=0)]
		return values

	def edit_nml(self, nml, val):
		nml['ky'] = val
		return nml

class ky(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ky']
	axis_label = r'$k_{y}$'
	valid_options = []

	def sub_validate(self, values):
		if any([x < 0 for x in values]):
			print("Error: ky values outside allowed range (x>=0)")
			values = [x for x in values if (x>=0)]
		return values

	def edit_nml(self, nml, val):
		return nml

class ky_model(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ky_model','kygrid_model']
	axis_label = 'ky_model'
	valid_options = []

	def sub_validate(self, values):
		if any([x not in [0,1,4] for x in values]):
			print("Error: ky_model value invalid, valid: [0,1,4]")
			values = [x for x in values if (x in [0,1,4])]
		return values

	def edit_nml(self, nml, val):
		nml['kygrid_model'] = val
		return nml

class kx(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['kx','akx','kx_rho0']
	axis_label = r'$k_{x}\rho_{0}$'
	valid_options = []

	def sub_validate(self, values):
		return values

	def edit_nml(self, nml, val):
		print("kx dimension only used for reading data in tglf")
		return nml

class ny(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ny','n_toroidal','ntoroidal','num_toroidal','nky']
	axis_label = 'ny'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: ny values outside the allowed range (x>0)")
			values = [x for x in values if (x>0)]
		if any([x != int(x) for x in values]):
			print("Error: ny values must be integers")
			values = [x for x in values if (x==int(x))]
		return values

	def edit_nml(self, nml, val):
		nml['nky'] = val
		return nml

class sat_rule(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['sat_rule']
	axis_label = 'sat_rule'
	valid_options = []

	def sub_validate(self, values):
		if any([x not in [0,1,2,3] for x in values]):
			print("Error: ky_model value invalid, valid: [0,1,2,3]")
			values = [x for x in values if (x in [0,1,2,3])]
		return values

	def edit_nml(self, nml, val):
		nml['sat_rule'] = val
		return nml

class q_loc(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['q','q_loc']
	axis_label = 'q_loc'
	valid_options = []
	
	def sub_validate(self, values):
		if any([x<0 for x in values]):
				values = [x if x>=0 for x in values]
				values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['q_loc'] = val
		return nml

dimensions_list = [psiN,p_prime,q_prime,shear,ky_max,ky,ky_model,kx,ny,sat_rule,q_loc]
