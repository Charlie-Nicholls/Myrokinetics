from ..code import code
from copy import deepcopy
from numpy import array, nan, trapz
import os

viking_modules = """"""

archer2_modules = """module load PrgEnv-gnu
module load cray-hdf5 cray-netcdf cray-fftw cray-python
export GACODE_PLATFORM=ARCHER2
export GACODE_ROOT=/work/e281/e281/cnicholls/gacode
. $GACODE_ROOT/shared/bin/gacode_setup
ulimit -s unlimited
source /work/e281/e281/cnicholls/pythenv/bin/activate
which cgyro
cgyro -h"""

class cgyro(code):
	def __init__(self):
		self.code_name = 'CGYRO'
		self.viking_modules = viking_modules
		self.archer2_modules = archer2_modules
		self.template = "template.cgyro"
		from ..dimensions.dimensions_cgyro import dimensions_list
		self.dim_list = dimensions_list
		if self.dim_list != None:
			self.make_dim_lookup(self.dim_list)
			
	verifications = ['omega','phi2','t','phi','apar','bpar','epar']
	input_name = "input.cgyro"
	output_name = "out.cgyro.run"
	
	def check_scan(self, valid = True):
		if self.inputs['ideal'] == True:
			print("ERROR: Ideal runs cannot be performed on CGYRO")
			valid = False
		return valid
	
	def get_template_lines(self):
		with open(os.path.join(self.inputs['template_path'],self.inputs['template_name']),'r') as f:
			template_lines = f.readlines()
		return template_lines
		
	def _get_surface_input(self, nml):
		if self.inputs['non_linear'] == True:
			nml['NONLINEAR_FLAG'] = 1
		else:
			nml['NONLINEAR_FLAG'] = 0
		
		nml['DELTA_T_METHOD'] = 1
		nml['EQUILIBRIUM_MODEL'] = 2
		
		nml['THETA_PLOT'] = nml['N_THETA'] #TEMPORARY UNTIL I FIND A BETTER SOLUTION TO ENSURING ALWAYS A FACTOR OF NTHETA AND/OR THE TEMPLATE LOADING ISSUE WITH PYRO
		return nml
	
	def _get_gyro_input(self, run, nml):
		if self.inputs['knobs']['fixed_delt'] == False:
			ky = nml['KY']
			delt = 0.04/ky
			if delt > 0.01:
				delt = 0.01
			#nml['DELT_T'] = delt
		return nml
	
	def make_job_files_ypi(self):
		print(f"ERROR: {self.code_name} DOES NOT SUPPORT YPI SERVERS")
		return
		
	def write_pyth_archer2(self, dir_list, filename):
		pyth = open(f"{self.inputs['data_path']}/submit_files/{filename}/{filename}.py",'w')
		pyth.write(f"""import os
from joblib import Parallel, delayed
from time import sleep
from numpy import array

input_dirs = {dir_list}
max_cores = {self.inputs["sbatch"]["nodes"]*self.inputs["sbatch"]["ntasks-per-node"]}

def start_run(run, run_attempt = 1):
	if run_attempt <= 3:
		os.system(f"echo \\\"Input: {{run}}\\\"")
		cwd = os.getcwd()
		os.chdir(f"{{run}}")
		if not os.path.exists("input.report"):
			os.system(f"$GACODE_ROOT/cgyro/bin/cgyro -i . &> input.report")
		f = open("input.report")
		lines = f.readlines()
		n_poss = set()
		for line in lines[3:]:
			n_poss.add(int([x for x in line.split(" ") if x.isdigit()][0]))
		poss_cores = [x for x in n_poss if x <= max_cores]
		cores = max(poss_cores)
		os.system(f"$GACODE_ROOT/cgyro/bin/cgyro -e . -n {{cores}} -nomp 1 -numa 8 -mpinuma 16 -p .")
		os.chdir(f"{{cwd}}")
		if os.path.exists(f"{{run}}/out.cgyro.freq"):
			os.system(f"touch {{run}}/run.fin")
		else:
			sleep(60)
			start_run(run, run_attempt = run_attempt+1)
	else:
		print(f"ERROR: {{run}} took too many attempts to start, skipping")

Parallel(n_jobs={self.inputs['sbatch']['nodes']})(delayed(start_run)(run) for run in input_dirs)""")
		pyth.close()
		return
	
	def get_non_linear_archer2(self):
		indir = list(self.scanner._input_dirs)[0]
		ntasks = self.inputs['sbatch']['nodes']*128//self.inputs['sbatch']['cpus-per-task']
		run_code = f'''cd {indir}
srun --nodes={self.inputs['sbatch']['nodes']} --ntasks={ntasks} --cpus-per-task={self.inputs['sbatch']['cpus-per-task']} $GACODE_ROOT/cgyro/bin/cgyro -e . -n {ntasks} -nomp 1 -numa 8 -mpinuma 16 -p .
if test -f \"./{self.output_name}\"; then
touch \"./run.fin\"
fi'''
		return run_code
		
	def make_gyro_file(self, run, sub_dir, namelist_diff = {}):
		if not os.path.exists(f"{sub_dir}/{self.input_name}"):
			subnml = self.get_gyro_input(run=run,namelist_diff=namelist_diff)
			self.write_nml(nml=subnml,directory=sub_dir,filename=self.input_name)
		return sub_dir
	
	def save_out(self, filename = None, directory = None, specificRuns = None, QuickSave = False):
		psi_itt = self.scanner.single_parameters['psin'].values if 'psin' in self.scanner.single_parameters else self.scanner.dimensions['psin'].values
		equilibrium = {}
		for psiN in psi_itt:
			equilibrium[psiN] = {}
			nml = self.eqbm.get_surface_input(psiN)
			equilibrium[psiN]['shear'] = nml['S']
			equilibrium[psiN]['beta_prime'] = nml['BETA_STAR_SCALE']
		
		if self.scanner['gyro']:
			gyro_data = {}
			group_data = {}
			only = set({'growth_rate','mode_frequency','ky','kx'})
			if not QuickSave:
				only = only | set({'phi','bpar','apar','time','theta','heat'})
			data_keys = ['growth_rate','mode_frequency','omega','phi','bpar','apar','epar','phi2','parity','ql_metric']
			group_keys = ['phi2_avg','t','theta', 'gds2', 'jacob','heat_flux_tot','phi2_by_kx','phi2_by_ky']
			gyro_keys = {}
			for dim in self.scanner.dimensions.values():
				gyro_keys[dim.name] = {}
				for val in dim.values:
					gyro_keys[dim.name][val] = set()
			
			if 'ky' not in gyro_keys.keys():
				gyro_keys['ky'] = {}
			kxs = set()
			kys = set()
			gyro_keys['kx'] = {}
			
			runs = self.scanner.get_all_runs() if specificRuns is None else list(specificRuns)
			for run in runs:
				sub_dir = self.scanner.get_run_directory(run)
				try:
					self.eqbm.pyro.load_gk_output(sub_dir)
					run_data = self.eqbm.pyro.gk_output
					group_key = run_data.attrs['object_uuid']
					group_data[group_key] = {}
					for key in group_keys:
						group_data[group_key][key] = None
					for xi, kx in enumerate(run_data['kx'].data):
						for yi, ky in enumerate(run_data['ky'].data):
							from uuid import uuid4
							run_key = str(uuid4())
							gyro_data[run_key] = deepcopy(run)
							for key in run:
								gyro_keys[key][run[key]].add(run_key)
							gyro_data[run_key]['group_key'] = group_key
	
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
							#gyro_data['nml_diffs'] = self.scanner.namelist_diffs[?]
							for key in data_keys:
								gyro_data[run_key][key] = None
							for key in only:
								try:
									key_data = run_data[key]
									if key == 'growth_rate':
										gyro_data[run_key]['growth_rate'] = float(key_data[xi,yi,-1])
									if key == 'mode_frequency':
										gyro_data[run_key]['mode_frequency'] = float(key_data[xi,yi,-1])
									elif key in ['phi','apar','bpar']:
										phi = array(key_data[:,xi,yi,-1])
										gyro_data[run_key][key] = phi.tolist()
										try:
											absint = abs(trapz(phi,array(run_data['theta'])))
											intabs = trapz(abs(phi),array(run_data['theta']))
											par = 1 - absint/intabs
											gyro_data[run_key]['parity'] = par
										except:
											gyro_data[run_key]['parity'] = None
									elif key in ['time']:
										group_data[group_key]['t'] = array(key_data).tolist()
									elif key in ['theta']:
										group_data[group_key][key] = array(key_data).tolist()
									elif key in ['heat']:
										group_data[group_key][key] = array(key_data[:,:,yi,:]).tolist()
								except Exception as e:
									print(f"Save Error in {sub_dir}: {e}")
									if key == 'growth_rate':
										gyro_data[run_key]['growth_rate'] = nan
									elif key == 'mode_frequency':
										gyro_data[run_key]['mode_frequency'] = nan
										
				except Exception as e:
					print(f"Save Error {sub_dir}: {e}")
			
			existing_dim_keys = []
			for key in [x for x in self.inputs.inputs.keys() if 'dimension_' in x]:
				existing_dim_keys.append([x for x in key if x.isdigit()])
			dim_n = max([eval("".join(x)) for x in existing_dim_keys],default=1) + 1
			kxs = list(kxs)
			kxs.sort()
			self.inputs.inputs[f'dimension_{dim_n}'] = {'type': 'kx', 'values': kxs, 'min': min(kxs), 'max': max(kxs), 'num': len(kxs), 'option': None}
			if 'ky' not in self.scanner.dimensions:
				kys = list(kys)
				kys.sort()
				self.inputs.inputs[f'dimension_{dim_n+1}'] = {'type': 'ky', 'values': kys, 'min': min(kys), 'max': max(kys), 'num': len(kys), 'option': None}
			self.inputs.load_dimensions()

		else:
			gyro_data = None
			gyro_keys = None
		
		data = {'gyro': gyro_data,
			'ideal': None,
			'group': group_data,
			'equilibrium': equilibrium,
			'_gyro_keys': gyro_keys,
			'_ideal_keys': None,
			}
		return data
	
	def write_nml(self, nml, directory = ".", filename = None):
		with open(f"{directory}/input.cgyro", "w") as f:
			for key, value in nml.items():
				f.write( f"{key} = {value}\n".upper())


