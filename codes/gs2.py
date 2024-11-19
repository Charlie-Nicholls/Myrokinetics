from ..code import code
from copy import deepcopy
from numpy import real, imag, array, loadtxt, transpose, nan, trapz
import os
import f90nml

viking_modules = """module purge
module load gompi/2022b
module load OpenMPI/4.1.4-GCC-12.2.0
module load netCDF-Fortran/4.6.0-gompi-2022b
module load FFTW/3.3.10-GCC-12.2.0
module load OpenBLAS/0.3.21-GCC-12.2.0
module load Python/3.10.8-GCCcore-12.2.0
export GK_SYSTEM=viking
export MAKEFLAGS=-IMakefiles
ulimit -s unlimited
export PATH=${PATH}:${HOME}/gs2/bin
which gs2
gs2 --build-config"""

archer2_modules = """module load PrgEnv-gnu
module load cray-hdf5 cray-netcdf cray-fftw cray-python
export GK_SYSTEM=archer2
export MAKEFLAGS=-IMakefiles
ulimit -s unlimited
export PATH=${PATH}:/work/e281/e281/cnicholls/gs2/bin
source /work/e281/e281/cnicholls/pythenv/bin/activate
which gs2
gs2 --build-config"""

class gs2(code):
	def __init__(self):
		self.code_name = 'GS2'
		self.viking_modules = viking_modules
		self.archer2_modules = archer2_modules
		self.template = "template.gs2"
		from ..dimensions.dimensions_gs2 import dimensions_list
		self.dim_list = dimensions_list
		if self.dim_list != None:
			self.make_dim_lookup(self.dim_list)
	
	def check_scan(self, inputs, valid = True):
		if self['grid_option'] == 'box':
			if 'ky' in inputs.dimensions or 'ky' in inputs.single_parameters:
				print("ERROR: ky dimension is not compatible with grid_option == box: use nx, ny, y0 and jtwist")
				valid = False
			if 'kx' in inputs.dimensions or 'kx' in inputs.single_parameters:
				print("ERROR: kx dimension is not compatible with grid_option == box: use nx, ny, y0 and jtwist")
				valid = False
			if 'theta0' in inputs.dimensions or 'theta0' in inputs.single_parameters:
				print("ERROR: theta0 dimension is not compatible with grid_option == box: use nx, ny, y0 and jtwist")
				valid = False
			if not inputs['fixed_delt']:
				print("grid_option == box runs can only use fixed_delt == True")
				self.inputs['knobs']['fixed_delt'] = True
		return valid
	
	def get_template_lines(self, inputs):
		template_lines = f90nml.read(os.path.join(inputs['template_path'],inputs['template_name']))
		return template_lines
		
	def _get_surface_input(self, inputs, nml):
		nml['theta_grid_parameters']['qinp'] = abs(nml['theta_grid_parameters']['qinp'])
		
		if 'parameters' in nml.keys() and 'beta' in nml['parameters'].keys():
			nml['knobs']['beta'] = nml['parameters']['beta']
			del(nml['parameters'])
		
		beta_prim = nml['theta_grid_eik_knobs']['beta_prime_input']
		
		beta = nml['knobs']['beta']
		bp_cal = sum([(nml[spec]['tprim'] + nml[spec]['fprim'])*nml[spec]['dens']*nml[spec]['temp'] for spec in [x for x in nml.keys() if 'species_parameters_' in x]])*beta*-1
		mul = beta_prim/bp_cal
		
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			nml[spec]['tprim'] = nml[spec]['tprim']*mul
			nml[spec]['fprim'] = nml[spec]['fprim']*mul
			
		for spec in [x for x in nml.keys() if 'species_parameters_' in x]:
			#Set bakdif to 0 for Electormagnetic Runs as a default
			nml[spec]['bakdif'] = 0
			#Force Uprim to 0, reccommendation by David
			nml[spec]['uprim'] = 0
			#TEMPORARY UNTIL PYRO FIXES BUG
			if 'bakdif' in nml[spec].keys():
				del(nml[spec]['bakdif'])
		#TEMPORARY UNTIL PYRO FIXES BUG
		if 'normalisations_knobs' in nml.keys():
			if 'qref' in nml['normalisations_knobs'].keys():
				del(nml['normalisations_knobs']['qref'])

		if 'ngauss' in nml['le_grids_knobs'] and 'npassing' in nml['le_grids_knobs']:
			del(nml['le_grids_knobs']['ngauss'])
		
		if inputs['force_zero_fs']:
			nml['dist_fn_knobs']['g_exb'] = 0
			
		nml['theta_grid_eik_knobs']['equal_arc'] = False
		nml['init_g_knobs']['ginit_option'] = 'random_sine'
		nml['dist_fn_knobs']['mach'] = 0
		
		if 'avail_cpu_time' not in nml['knobs'].keys():
				h, m, s = inputs['sbatch']['time'].split(':')
				nml['knobs']['avail_cpu_time'] = (int(h) * 3600) + (int(m) * 60) + int(s)
		
		if inputs['grid_option'] == 'single':
			nml['kt_grids_knobs']['grid_option'] = 'single'
			if 'kt_grids_single_parameters' not in nml:
				nml['kt_grids_single_parameters'] = {'aky': 0.1, 'theta0': 0}
			if 'kt_grids_box_parameters' in nml:
				del(nml['kt_grids_box_parameters'])
			if 'margin_cpu_time' not in nml['knobs'].keys():
				nml['knobs']['margin_cpu_time'] = 300
		elif inputs['grid_option'] == 'box':
			nml['kt_grids_knobs']['grid_option'] = 'box'
			nml['dist_fn_knobs']['boundary_option'] = 'linked'
			nml['dist_fn_knobs']['esv'] = True
			nml['fields_knobs']['field_option'] = 'local'
			if 'kt_grids_box_parameters' not in nml:
				nml['kt_grids_box_parameters'] = {'nx': 50, 'ny': 50, 'y0': -0.05, 'jtwist': 1}
			if 'kt_grids_single_parameters' in nml:
				del(nml['kt_grids_single_parameters'])
			nml['fields_knobs']['dump_response'] = True
			nml['fields_knobs']['response_dir'] = "response"
			nml['init_g_knobs']['restart_dir'] = "restart"
			nml['gs2_diagnostics_knobs']['save_for_restart'] = True
			nml['gs2_diagnostics_knobs']['nc_sync_freq'] = 1
			if nml['gs2_diagnostics_knobs']['nsave'] > 1000:
				nml['gs2_diagnostics_knobs']['nsave'] = 1
			if 'margin_cpu_time' not in nml['knobs'].keys():
				nml['knobs']['margin_cpu_time'] = 2400
			if 'nperiod' not in inputs.dimensions and 'nperiod' not in inputs.single_parameters:
					nml['theta_grid_parameters']['nperiod'] = 1
		else:
			print("ERROR: grid_option is invalid, valid: ['single','box']")
					
		if inputs['non_linear'] == True:
			if 'nonlinear_terms_knobs' not in nml.keys():
				nml['nonlinear_terms_knobs'] = {}
			nml['nonlinear_terms_knobs']['nonlinear_mode'] = 'on'
			if 'cfl' not in nml['nonlinear_terms_knobs'].keys():
				nml['nonlinear_terms_knobs']['cfl']	= 0.5
			if inputs['split_nonlinear'] == True:
				nml['nonlinear_terms_knobs']['split_nonlinear'] = True
				if 'split_nonlinear_terms_knobs' not in nml.keys():
					nml['split_nonlinear_terms_knobs'] = {'show_statistics': True}
			
		if inputs['Miller']:
			nml['theta_grid_eik_knobs']['local_eq'] = True
		else:
			nml['theta_grid_eik_knobs']['eqfile'] = os.path.join(inputs['data_path'],inputs['eq_name'])
			nml['theta_grid_eik_knobs']['efit_eq'] =  True
			nml['theta_grid_eik_knobs']['local_eq'] = False
		
		if inputs['Epar']:
			nml['gs2_diagnostics_knobs']['write_ascii'] = True
			nml['gs2_diagnostics_knobs']['write_final_epar'] = True
		else:
			nml['gs2_diagnostics_knobs']['write_ascii'] = False
			nml['gs2_diagnostics_knobs']['write_final_epar'] = False
		
		if 'ntheta_geometry' not in nml['theta_grid_eik_knobs'].keys():
			nml['theta_grid_eik_knobs']['ntheta_geometry'] = 4096

		if inputs['Ideal']:
			beta_mul = abs(inputs['beta_prime_max']/beta_prim)
			beta_div = abs(beta_prim/inputs['beta_prime_min'])

			nml['ballstab_knobs'] = {'n_shat': inputs['n_shat_ideal'], 'n_beta': inputs['n_beta_ideal'], 'shat_min': inputs['shat_min'], 'shat_max': inputs['shat_max'], 'beta_div': beta_div, 'beta_mul': beta_mul}

		nml['knobs']['wstar_units'] = False
		return nml
	
	def _get_gyro_input(self, inputs, run, nml):
		if inputs['knobs']['fixed_delt'] == False:
			ky = nml['kt_grids_single_parameters']['aky']
			delt = 0.04/ky
			if delt > 0.01:
				delt = 0.01
			nml['knobs']['delt'] = delt
		return nml
	
	def make_job_files_ypi(self, scanner):
		if n_jobs is None or n_jobs > len(scanner._input_files):
			n_jobs = len(scanner._input_files)
		while n_jobs > 0:
			for input_file in scanner._input_files:
				os.system(f"mpirun -np 8 gs2 \"{input_file}.in\"")
				scanner._input_files.remove(input_file)
				n_jobs -= 1
		return
		
	def make_ideal_job_files_ypi(self, scanner):
		if n_jobs is None or n_jobs > len(self._input_files):
			n_jobs = len(self._input_files)
		while n_jobs > 0:
			for input_file in self._ideal_input_files:
				os.system(f"ideal_ball \"{input_file}.in\"")
				self._ideal_input_files.remove(input_file)
				n_jobs -= 1
		return
	
	
	jobfile_viking = '''echo "${{INFILE}}.in"
gs2 "${INFILE}.in"
if test -f "${INFILE}.out.nc"; then
touch "${INFILE}.fin"
fi'''
		
	def write_pyth_archer2(self, scanner, input_list, filename):

		pyth = open(f"{scanner.inputs['data_path']}/submit_files/{filename}/{filename}.py",'w')
		pyth.write(f"""import os
from joblib import Parallel, delayed
from time import sleep

input_files = {input_list}

def start_run(run, run_attempt = 1):
if run_attempt <= 3:
	os.system(f"echo \\\"Input: {{run}}.in\\\"")
	os.system(f"gs2 \\\"{{run}}.in\\\"")
	if os.path.exists(f"\\\"{{run}}.out.nc\\\"'):
		os.system(f"touch \\\"{{run}}.fin\\\"")
	else:
		sleep(60)
		start_run(run, run_attempt = run_attempt+1)
else:
	print(f"ERROR: {{run}} took too many attempts to start, skipping")

Parallel(n_jobs={scanner.inputs['sbatch']['nodes']})(delayed(start_run)(run) for run in input_files)""")
		pyth.close()
		
		return
	
	def get_non_linear_archer2(self, scanner):
		ntasks = scanner.inputs['sbatch']['nodes']*128//scanner.inputs['sbatch']['cpus-per-task']
		run_code = f'''srun --nodes={scanner.inputs['sbatch']['nodes']} --ntasks={ntasks} --cpus-per-task={scanner.inputs['sbatch']['cpus-per-task']} gs2 \"{list(scanner._input_files)[0]}.in\"
if test -f \"{list(scanner._input_files)[0]}.out.nc\"; then
touch \"{list(scanner._input_files)[0]}.fin\"
fi'''
		return run_code
		
	def make_gyro_file(self, eqbm, run, sub_dir, namelist_diff = {}):
		if scanner.inputs['grid_option'] == 'box':
			os.makedirs(sub_dir+'/response',exist_ok=True)
			os.makedirs(sub_dir+'/restart',exist_ok=True)
		
		existing_inputs = [] 
		for f in glob.glob(r'itteration_*.in'):
			existing_inputs.append([x for x in f if x.isdigit()])
		itt = max([eval("".join(x)) for x in existing_inputs],default=-1)
		filename = f"itteration_{scanner.inputs['itteration']}"
		if not os.path.exists(f"{sub_dir}/{filename}"):
			subnml = self.get_gyro_input(eqbm=eqbm,run=run,namelist_diff=namelist_diff)
			self.write_nml(nml=subnml,directory=sub_dir,filename=f"{filename}.in")
		return f"{sub_dir}/{filename}"
	
	def save_out(self, scanner, filename = None, directory = None, specificRuns = None, QuickSave = False):
		
		psi_itt = scanner.single_parameters['psin'].values if 'psin' in scanner.single_parameters else scanner.dimensions['psin'].values
		equilibrium = {}
		for psiN in psi_itt:
			equilibrium[psiN] = {}
			nml = scanner.eqbm.get_surface_input(psiN)
			equilibrium[psiN]['shear'] = nml['theta_grid_eik_knobs']['s_hat_input']
			equilibrium[psiN]['beta_prime'] = nml['theta_grid_eik_knobs']['beta_prime_input']

		if scanner['gyro']:
			gyro_data = {}
			group_data = {}
			only = set({'omega','kx','ky'})
			if not QuickSave:
				only = only | set({'phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob','ql_metric_by_mode', 'phi2_by_mode'})
			#if scanner.inputs['epar']:
				#only = only | set({'epar'}) NOT CURRENTLY WORKING
			if scanner.inputs['grid_option'] == 'box':
				only = only | set({'phi2_by_kx', 'phi2_by_ky'})
			if scanner.inputs['non_linear'] == True:
				only = only | set({'heat_flux_tot'})
			data_keys = ['growth_rate','mode_frequency','omega','phi','bpar','apar','epar','phi2','parity','ql_metric']
			group_keys = ['phi2_avg','t','theta', 'gds2', 'jacob','heat_flux_tot','phi2_by_kx', 'phi2_by_ky']
			gyro_keys = {}
			for dim in scanner.dimensions.values():
				gyro_keys[dim.name] = {}
				for val in dim.values:
					gyro_keys[dim.name][val] = set()
			
			if scanner.inputs['grid_option'] == 'box':
				kxs = set()
				kys = set()
				gyro_keys['ky'] = {}
				gyro_keys['kx'] = {}
			
			if specificRuns:
				runs = list(specificRuns)
			else:
				runs = scanner.get_all_runs() if scanner.inputs['grid_option'] != 'box' else scanner.get_all_runs(excludeDimensions=['kx','ky'])
			
			for run in runs:
				sub_dir = scanner.get_run_directory(run)
				try:
					existing_inputs = [] 
					for f in glob.glob(r'itteration_*.in'):
						existing_inputs.append([x for x in f if x.isdigit()])
					itt = max([eval("".join(x)) for x in existing_inputs],default=0)
					run_data = readnc(f"{sub_dir}/itteration_{itt}.out.nc",only=only)	
					group_key = run_data['attributes']['id']
					group_data[group_key] = {}
					for key in group_keys:
						group_data[group_key][key] = None
					for xi, kx in enumerate(run_data['kx']):
						for yi, ky in enumerate(run_data['ky']):
							from uuid import uuid4
							run_key = str(uuid4())
							gyro_data[run_key] = deepcopy(run)
							for key in run:
								gyro_keys[key][run[key]].add(run_key)
							gyro_data[run_key]['group_key'] = group_key
							if scanner.inputs['grid_option'] == 'box':
								kxs.add(kx)
								kys.add(ky)
								if ky not in gyro_keys['ky']:
									gyro_keys['ky'][ky] = set()
								if kx not in gyro_keys['kx']:
									gyro_keys['kx'][kx] = set()
								gyro_keys['ky'][ky].add(run_key)
								gyro_keys['kx'][kx].add(run_key)
							if 'kx' not in gyro_data[run_key]:
								gyro_data[run_key]['kx'] = kx
							if 'ky' not in gyro_data[run_key]:
								gyro_data[run_key]['ky'] = ky
							#gyro_data['nml_diffs'] = scanner.namelist_diffs[?]
							for key in data_keys:
								gyro_data[run_key][key] = None
							for key in only:
								try:
									key_data = run_data[key]
									
									if key == 'omega':
										om = key_data[-1,yi,xi]
										if type(om) != complex:
											om = key_data[-2,yi,xi]
										gyro_data[run_key]['growth_rate'] = imag(om)
										gyro_data[run_key]['mode_frequency'] = real(om)
										gyro_data[run_key]['omega'] = key_data[:,yi,xi].tolist()
									elif key in ['phi','apar','bpar']:
										gyro_data[run_key][key] = key_data[yi,xi,:].tolist()
										if key == 'phi':
											try:
												absint = abs(trapz(key_data[yi,xi,:],run_data['theta']))
												intabs = trapz(abs(key_data[yi,xi,:]),run_data['theta'])
												par = 1 - absint/intabs #Equation taken from arXiv:2401.14260v1
												gyro_data[run_key]['parity'] = par
											except:
												gyro_data[run_key]['parity'] = None
									elif key in ['t','theta', 'gds2', 'jacob','heat_flux_tot','phi2_by_kx','phi2_by_ky']:
										group_data[group_key][key] = key_data.tolist()
									elif key in ['phi2']:
										group_data[group_key]['phi2_avg'] = key_data.tolist()
									elif key in ['ql_metric_by_mode']:
										gyro_data[run_key]['ql_metric'] = key_data[-1,yi,xi]
									elif key in ['phi2_by_mode']:
										gyro_data[run_key]['phi2'] = key_data[:,yi,xi]
									elif key in ['epar']:
										epar_path = f"{sub_dir}/itteration_{itt}.epar"
										epar_data = loadtxt(epar_path)
										epar = []
										for l in range(len(epar_data[:,3])):
											epar.append(complex(epar_data[l,3],epar_data[l,4]))
										epar = array(epar)
										gyro_data[run_key]['epar'] = epar
								except Exception as e:
									print(f"Save Error in {sub_dir}/itteration_{itt}: {e}")
									if key == 'omega':
										gyro_data[run_key]['growth_rate'] = nan
										gyro_data[run_key]['mode_frequency'] = nan
										
				except Exception as e:
					print(f"Save Error {sub_dir}/itteration_{itt}: {e}")
			if scanner.inputs['grid_option'] == 'box':
				existing_dim_keys = []
				for key in [x for x in scanner.inputs.inputs.keys() if 'dimension_' in x]:
					existing_dim_keys.append([x for x in key if x.isdigit()])
				dim_n = max([eval("".join(x)) for x in existing_dim_keys],default=1) + 1
				kxs = list(kxs)
				kxs.sort()
				scanner.inputs.inputs[f'dimension_{dim_n}'] = {'type': 'kx', 'values': kxs, 'min': min(kxs), 'max': max(kxs), 'num': len(kxs), 'option': None}
				kys = list(kys)
				kys.sort()
				scanner.inputs.inputs[f'dimension_{dim_n+1}'] = {'type': 'ky', 'values': kys, 'min': min(kys), 'max': max(kys), 'num': len(kys), 'option': None}
				scanner.inputs.load_dimensions()
		else:
			gyro_data = None
			gyro_keys = None

		if scanner['ideal']:
			ideal_keys = {}
			if 'theta0' in scanner.single_parameters:
				theta0_itt = scanner.single_parameters['theta0'].values  
			if 'theta0' in scanner.dimensions:
				theta0_itt = scanner.dimensions['theta0'].values
			else:
				theta0_itt = [0]
			
			ideal_keys['psin'] = {}
			ideal_keys['theta0'] = {}
			for val in psi_itt:
				ideal_keys['psin'][val] = set()
			for val in theta0_itt:
				ideal_keys['theta0'][val] = set()

			ideal_data = {}
			for run in scanner.get_all_ideal_runs():
				run_id = str(uuid4())
				for key in run:
					ideal_keys[key][run[key]].add(run_id)
				ideal_data[run_id] = {}
				try:
					sub_dir = scanner.get_ideal_run_directory(run)
					existing_inputs = [] 
					for f in glob.glob(r'itteration_*.in'):
						existing_inputs.append([x for x in f if x.isdigit()])
					itt = max([eval("".join(x)) for x in existing_inputs],default=0)

					shear = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_shat")
					bp = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_bp")
					stab = loadtxt(f"{sub_dir}/itteration_{itt}.ballstab_2d")
					
					ideal_data[run_id]['beta_prime'] = [abs(x) for x in bp]
					ideal_data[run_id]['shear'] = shear.tolist()
					ideal_data[run_id]['stabilities'] = transpose(stab).tolist()
				except:
					ideal_data[run_id]['beta_prime'] = None
					ideal_data[run_id]['shear'] = None
					ideal_data[run_id]['stabilities'] = None
					print(f"Save Error for ideal run: {run}")
		else:
			ideal_data = None
			ideal_keys = None
		
		data = {'gyro': gyro_data,
			'ideal': ideal_data,
			'group': group_data,
			'equilibrium': equilibrium,
			'_gyro_keys': gyro_keys,
			'_ideal_keys': ideal_keys,
			}
		
		return data
	
	def write_nml(self, nml, directory = ".", filename = None):
		if filename is None:
			filename = "itteration_0.in"
		nml.write(f"{directory}/{filename}", force=True)
		return


