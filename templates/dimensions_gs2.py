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

class beta_prime(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['beta_prime','betap','bp','beta_p','b_p']
	axis_label = r'$\beta^{\prime}$'
	valid_options = []
	
	def sub_validate(self, values):
		values = [abs(x) for x in values]
		values.sort()
		return values

	def edit_nml(self, nml, val):
		nml['theta_grid_eik_knobs']['beta_prime_input'] = -1*abs(val)
		
		beta = nml['knobs']['beta'] if 'beta' in nml['knobs'].keys() else nml['parameters']['beta']
		bp_cal = sum((nml[spec]['tprim'] + nml[spec]['fprim'])*nml[spec]['dens']*nml[spec]['temp'] for spec in [x for x in nml.keys() if 'species_parameters_' in x])*beta*-1

		mul = -1*abs(val)/bp_cal
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			nml[spec]['tprim'] = nml[spec]['tprim']*mul
			nml[spec]['fprim'] = nml[spec]['fprim']*mul

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
		nml['theta_grid_eik_knobs']['s_hat_input'] = val
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
		nml['kt_grids_single_parameters']['aky'] = val
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
		nml['kt_grids_single_parameters']['theta0'] = val
		if 'akx' in nml['kt_grids_single_parameters']:
			del(nml['kt_grids_single_parameters']['akx'])
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
		nml['kt_grids_single_parameters']['akx'] = val
		if 'theta0' in nml['kt_grids_single_parameters']:
			del(nml['kt_grids_single_parameters']['theta0'])
		return nml

class nperiod(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['nperiod']
	axis_label = 'nperiod'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: nperiod values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		if any([x != int(x) for x in values]):
			print("Error: nperiod values must be integers")
			values = [x for x in values if (x==int(x))]
		return values

	def edit_nml(self, nml, val):
		nml['theta_grid_parameters']['nperiod'] = val
		return nml

class ntheta(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ntheta']
	axis_label = 'ntheta'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: ntheta values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		if any([x != int(x) for x in values]):
			print("Error: ntheta values must be integers")
			values = [x for x in values if (x==int(x))]
		return values

	def edit_nml(self, nml, val):
		nml['theta_grid_parameters']['ntheta'] = val
		return nml

class bakdif(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['bakdif']
	axis_label = 'bakdif'
	valid_options = []

	def sub_validate(self, values):
		if any([x < 0 for x in values]):
			print("Error: bakdif values outside allowed range (x>=0)")
			values = [x for x in values if (x>0)]
		return values

	def edit_nml(self, nml, val):
		for key in [x for x in nml.keys() if 'dist_fn_species_knobs_' in x]:
			nml[key]['bakdif'] = val
		return nml
		
class fexpr(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['fexpr']
	axis_label = 'fexpr'
	valid_options = []

	def sub_validate(self, values):
		if any([x < 0 or x > 1 for x in values]):
			print("Error: fexpr values outside allowed range (0<=x<=1)")
			values = [x for x in values if (0<=x<=1)]
		return values

	def edit_nml(self, nml, val):
		for key in [x for x in nml.keys() if 'dist_fn_species_knobs_' in x]:
			nml[key]['fexpr'] = val
		return nml

class delt(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = None)

	name_keys = ['delt']
	axis_label = 'delt'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: delt values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		return values

	def edit_nml(self, nml, val):
		nml['knobs']['delt'] = val
		return nml

class vnewk(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)
		if self.option is None:
			self.option = 'all'

	name_keys = ['vnewk']
	axis_label = 'vnewk'
	valid_options = ['all','electron','ion','deuterium','tritium','impurity']

	def sub_validate(self, values):
		if any([x < 0 for x in values]):
			print("Error: delt values outside allowed range (x>=0)")
			values = [x for x in values if (x>=0)]
		return values

	def edit_nml(self, nml, val):
		for key in [x for x in nml.keys() if 'species_parameters_' in x]:
			if nml[key]['type'] == 'electron' and self.option in ['all','electron']:
				nml[key]['vnewk'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['mass'] == 1 and self.option in ['all','ion','deuterium']:
				nml[key]['vnewk'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['z'] == 1 and self.option in ['all','ion','tritium']:
				nml[key]['vnewk'] = val
			elif nml[key]['type'] == 'ion' and self.option in ['all','impurity']:
				nml[key]['vnewk'] = val
		return nml

class tprim(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)
		if self.option is None:
			self.option = 'all'
		self._fix_fprim = True
		if self._fix_fprim:
			print("WARNING: USING EXPERIMENTAL FIX FPRIM SETTINGS")
		else:
			print("WARNING: TPRIM AND FPRIM HAVE NOT BEEN TESTED FOR CONSISTENCY IN LATEST VERSION, PLEASE CHECK BEFORE USE")

	name_keys = ['tprim']
	axis_label = 'tprim'
	valid_options = ['all','electron','ion','deuterium','tritium','impurity']

	def sub_validate(self, values):
		return values

	def edit_nml(self, nml, val):
		for key in [x for x in nml.keys() if 'species_parameters_' in x]:
			if nml[key]['type'] == 'electron' and self.option in ['all','electron']:
				nml[key]['tprim'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['mass'] == 1 and self.option in ['all','ion','deuterium']:
				nml[key]['tprim'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['z'] == 1 and self.option in ['all','ion','tritium']:
				nml[key]['tprim'] = val
			elif nml[key]['type'] == 'ion' and self.option in ['all','impurity']:
				nml[key]['tprim'] = val
			else:
				break
			if self._fix_fprim and self.option != 'all':
				nml[key]['fprim'] = self.fprimcal(nml,val)
		return nml
	
	def fprimcal(nml,tprim):
		bp = abs(nml['theta_grid_eik_knobs']['beta_prime_input'])
		beta = nml['knobs']['beta'] if 'beta' in nml['knobs'].keys() else nml['parameters']['beta']
		sp12 = [(nml[spec]['tprim'] + nml[spec]['fprim'])*nml[spec]['dens']*nml[spec]['temp'] for spec in ['species_parameters_1','species_parameters_2']]
		dp = ((bp/beta) - sp12[0] - sp12[1])/(nml['species_parameters_3']['dens']*nml['species_parameters_3']['temp']) - tprim
		return dp


class fprim(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)
		print("WARNING: TPRIM AND FPRIM HAVE NOT BEEN TESTED FOR CONSISTENCY IN LATEST VERSION, PLEASE CHECK BEFORE USE")
		if self.option is None:
			self.option = 'all'

	name_keys = ['fprim']
	axis_label = 'fprim'
	valid_options = ['all','electron','ion','deuterium','tritium','impurity']

	def sub_validate(self, values):
		return values

	def edit_nml(self, nml, val):
		for key in [x for x in nml.keys() if 'species_parameters_' in x]:
			if nml[key]['type'] == 'electron' and self.option in ['all','electron']:
				nml[key]['fprim'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['mass'] == 1 and self.option in ['all','ion','deuterium']:
				nml[key]['fprim'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['z'] == 1 and self.option in ['all','ion','tritium']:
				nml[key]['fprim'] = val
			elif nml[key]['type'] == 'ion' and self.option in ['all','impurity']:
				nml[key]['fprim'] = val
		return nml

class mass(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)
		if self.option is None:
			self.option = 'electron'

	name_keys = ['mass']
	axis_label = 'mass'
	valid_options = ['all','electron','ion','deuterium','tritium','impurity']

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: mass values outside allowed range (x>0)")
			values = [x for x in values if (x>0)]
		return values

	def edit_nml(self, nml, val):
		for key in [x for x in nml.keys() if 'species_parameters_' in x]:
			if nml[key]['type'] == 'electron' and self.option in ['all','electron']:
				nml[key]['mass'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['mass'] == 1 and self.option in ['all','ion','deuterium']:
				nml[key]['mass'] = val
			elif nml[key]['type'] == 'ion' and nml[key]['z'] == 1 and self.option in ['all','ion','tritium']:
				nml[key]['mass'] = val
			elif nml[key]['type'] == 'ion' and self.option in ['all','impurity']:
				nml[key]['mass'] = val
		return nml

class nx(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['nx']
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
		if 'kt_grids_box_parameters' not in nml.keys():
			nml['kt_grids_box_parameters'] = {}
		nml['kt_grids_box_parameters']['nx'] = val
		return nml

class ny(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['ny']
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
		if 'kt_grids_box_parameters' not in nml.keys():
			nml['kt_grids_box_parameters'] = {}
		nml['kt_grids_box_parameters']['ny'] = val
		return nml

class y0(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['y0']
	axis_label = 'y0'
	valid_options = []

	def sub_validate(self, values):
		if any([x == 0 for x in values]):
			print("Error: y0 values outside the allowed range (x!=0)")
			values = [x for x in values if (x!=0)]
		return values

	def edit_nml(self, nml, val):
		if 'kt_grids_box_parameters' not in nml.keys():
			nml['kt_grids_box_parameters'] = {}
		nml['kt_grids_box_parameters']['y0'] = val
		return nml

class jtwist(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['jtwist']
	axis_label = 'jtwist'
	valid_options = []

	def sub_validate(self, values):
		if any([x <= 0 for x in values]):
			print("Error: jtwist values outside the allowed range (x>0)")
			values = [x for x in values if (x>0)]
		if any([x != int(x) for x in values]):
			print("Error: jtwist values must be integers")
			values = [x for x in values if (x==int(x))]
		return values

	def edit_nml(self, nml, val):
		if 'kt_grids_box_parameters' not in nml.keys():
			nml['kt_grids_box_parameters'] = {}
		nml['kt_grids_box_parameters']['jtwist'] = val
		return nml

class cfl(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['cfl']
	axis_label = 'cfl'
	valid_options = []

	def sub_validate(self, values):
		if any([(x < 0 or x > 1) for x in values]):
			print("Error: cfl values outside the allowed range (0<=x<=1)")
			values = [x for x in values if (0<=x<=1)]
		return values

	def edit_nml(self, nml, val):
		if 'nonlinear_terms_knobs' not in nml.keys():
			nml['nonlinear_terms_knobs'] = {}
		nml['nonlinear_terms_knobs']['cfl'] = val
		return nml

class g_exb(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):
		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['g_exb','flow_shear','flowshear','fs','gexb']
	axis_label = 'g_exb'
	valid_options = []

	def sub_validate(self, values):
		if any([(x < 0) for x in values]):
			print("Error: g_exb values outside the allowed range (x>=0)")
			values = [x for x in values if (x>=0)]
		return values

	def edit_nml(self, nml, val):
		if 'dist_fn_knobs' not in nml.keys():
			nml['dist_fn_knobs'] = {}
		nml['dist_fn_knobs']['g_exb'] = val
		return nml
	
class qinp(dimension):
	def __init__(self, values = None, mini = None, maxi = None, num = None, option = None):

		super().__init__(values = values, mini = mini, maxi = maxi, num = num, option = option)

	name_keys = ['qinp','q_inp','safety_factor']
	axis_label = 'qinp'
	valid_options = []

	def sub_validate(self, values):
		return values

	def edit_nml(self, nml, val):
		if 'theta_grid_parameters' not in nml.keys():
			nml['theta_grid_parameters'] = {}
		nml['theta_grid_parameters']['qinp'] = val
		return nml

dimensions_list = [psiN,beta_prime,shear,ky,theta0,kx,nperiod,ntheta,bakdif,fexpr,delt,vnewk,tprim,fprim,mass,nx,ny,y0,jtwist,cfl,g_exb,qinp]
