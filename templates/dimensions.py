from ..inputs import dimension

'''
TEMPLATE DIMENSION SUBCLASS | Rememeber to add any new dimensions to dimensions_list
class DimensionType(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['']
	axis_label = ''

	def sub_validate(self, values):
		#Any constraints or corrections for your specific dimension
		return values

	def edit_nml(self, nml, val):
		#nml is a gs2 input file, edit and return it with changed dimension value = val
		return nml
'''

class psiN(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['psin','psi','psi_norm']
	axis_label = '$\u03C8_{N}$'

	def sub_validate(self, values):
		if any([x <= 0 or x > 1 for x in values]):
			print("Error: psiN values outside allowed range (0<x<=1)")
			values = [x for x in values if (0<x<=1)]
		return values

	def edit_nml(self, nml, val):
		return nml

class beta_prime(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['beta_prime','betap','bp','beta_p','b_p']
	axis_label = '\u03B2\''

	def sub_validate(self, values):
		values = [abs(x) for x in values]
		values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['theta_grid_eik_knobs']['beta_prime_input'] = -1*abs(val)

		bp_cal = sum((nml[spec]['tprim'] + nml[spec]['fprim'])*nml[spec]['dens']*nml[spec]['temp'] for spec in [x for x in nml.keys() if 'species_parameters_' in x])*nml['parameters']['beta']*-1

		mul = -1*abs(val)/bp_cal
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			nml[spec]['tprim'] = nml[spec]['tprim']*mul
			nml[spec]['fprim'] = nml[spec]['fprim']*mul

		return nml

class shear(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['shear','shat','s_hat','sh']
	axis_label = '$\hat{s}$'

	def sub_validate(self, values):
		if any([x<1e-4 for x in values]):
				values = list(set([x if x>1e-4 else 1e-4 for x in values]))
				values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['theta_grid_eik_knobs']['s_hat_input'] = val
		return nml

class ky(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['ky','aky','ky_rho0']
	axis_label = '$k_{y}\rho_{0}$'

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: ky values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		return values

	def edit_nml(self, nml, val):
		nml['kt_grids_single_parameters']['aky'] = val
		return nml

class theta0(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['theta0','theta','t0']
	axis_label = '$\theta_{0}$'

	def sub_validate(self, values):
		from numpy import pi
		if any([x < -pi or x > pi for x in values]):
			print("Error: theta0 values outside allowed range (-pi<=x<=pi)")
			values = [x for x in values if (-pi<=x<=pi)]
		return values

	def edit_nml(self, nml, val):
		nml['kt_grids_single_parameters']['theta0'] = val
		return nml

class nperiod(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['nperiod']
	axis_label = 'nperiod'

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: nperiod values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		return values

	def edit_nml(self, nml, val):
		nml['theta_grid_parameters']['nperiod'] = val
		return nml

class ntheta(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num)

	name_keys = ['ntheta']
	axis_label = 'ntheta'

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: ntheta values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		return values

	def edit_nml(self, nml, val):
		nml['theta_grid_parameters']['ntheta'] = val
		return nml


dimensions_list = [psiN,beta_prime,shear,ky,theta0,nperiod,ntheta]
