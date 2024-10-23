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

class alpha(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['alpha']
	axis_label = r'$\alpha$'
	valid_options = []
	
	def sub_validate(self, values):
		values = [abs(x) for x in values]
		values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['ALPHA_SA'] = val
		return nml

class shear(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['shear','shat','s_hat','sh','s']
	axis_label = r'$\hat{s}$'
	valid_options = []
	
	def sub_validate(self, values):
		if any([x<1e-4 for x in values]):
				values = list(set([x if x>1e-4 else 1e-4 for x in values]))
				values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['SHAT_SA'] = val
		return nml

class ky(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ky','aky','ky_rho0']
	axis_label = r'$k_{y}\rho_{0}$'
	valid_options = []

	def sub_validate(self, values):
		if any([x < 0 for x in values]):
			print("Error: ky values outside allowed range (x>=0)")
			values = [x for x in values if (x>=0)]
		return values

	def edit_nml(self, nml, val):
		nml['KY'] = val
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
		print("kx dimension only used for reading data in cgyro")
		return nml

class theta0(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['theta0','theta','t0']
	axis_label = r'$\theta_{0}$'
	valid_options = []

	def sub_validate(self, values):
		from numpy import pi
		if any([x < -pi or x > pi for x in values]):
			print("Error: theta0 values outside allowed range (-pi<=x<=pi)")
			values = [x for x in values if (-pi<=x<=pi)]
		return values

	def edit_nml(self, nml, val):
		nml['THETA0_SA'] = val
		return nml

class y0(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['y0','ky_min']
	axis_label = 'y0'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: y0 values outside allowed range (x>0)")
			values = [x for x in values if (x>=0)]
		return values

	def edit_nml(self, nml, val):
		nml['KY'] = val
		return nml

class ntheta(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ntheta','n_theta','num_theta']
	axis_label = 'ntheta'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: ntheta values outside the allowed range (x>0)")
			values = [x for x in values if (x>0)]
		if any([x != int(x) for x in values]):
			print("Error: ntheta values must be integers")
			values = [x for x in values if (x==int(x))]
		return values

	def edit_nml(self, nml, val):
		nml['N_THETA'] = val
		nml['THETA_PLOT'] = val #Temporary?
		return nml

class nx(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['nx','n_radial','nradial','num_radial']
	axis_label = 'nx'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: nx values outside the allowed range (x>0)")
			values = [x for x in values if (x>0)]
		if any([x != int(x) for x in values]):
			print("Error: nx values must be integers")
			values = [x for x in values if (x==int(x))]
		return values

	def edit_nml(self, nml, val):
		nml['N_RADIAL'] = val
		return nml

class ny(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ny','n_toroidal','ntoroidal','num_toroidal']
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
		nml['N_TOROIDAL'] = val
		return nml

class delt(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = None)

	name_keys = ['delt','DELTA_T','DELTA']
	axis_label = 'delt'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: delt values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		return values

	def edit_nml(self, nml, val):
		nml['DELTA_T'] = val
		return nml

class jtwist(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = None)

	name_keys = ['jtwist']
	axis_label = 'jtwist'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: jtwist values outside allowed range (x>0)")
		if any([int(x) != x for x in values]):
			print("Error: jtwist must be integers")
			values = [x for x in values if (x>0 and int(x)==x)]
		return values

	def edit_nml(self, nml, val):
		nml['BOX_SIZE'] = val
		return nml

dimensions_list = [psiN,alpha,shear,ky,kx,theta0,y0,ntheta,nx,ny,delt,jtwist]
