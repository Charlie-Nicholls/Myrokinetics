from ..code import code
import os

viking_modules = """"""

archer2_modules = """module load PrgEnv-gnu
module load cray-hdf5 cray-netcdf cray-fftw cray-python
export GACODE_PLATFORM=ARCHER2
export GACODE_ROOT=/work/e281/e281/cnicholls/gacode
. $GACODE_ROOT/shared/bin/gacode_setup
ulimit -s unlimited
source /work/e281/e281/cnicholls/pythenv/bin/activate
which tglf
tglf -h"""

class tglf(code):
	def __init__(self):
		self.code_name = 'TGLF'
		self.viking_modules = viking_modules
		self.archer2_modules = archer2_modules
		self.template = "template.gs2"
		from ..dimensions.dimensions_tglf import dimensions_list
		self.dim_list = dimensions_list
		if self.dim_list != None:
			self.make_dim_lookup(self.dim_list)
	
	def check_scan(self, inputs, valid = True):
		if inputs['ideal'] == True:
			print("ERROR: Ideal runs cannot be performed on TGLF")
			valid = False
		return valid
	
	def get_template_lines(self, inputs):
		with open(os.path.join(inputs['template_path'],inputs['template_name']),'r') as f:
				template_lines = f.readlines()
		return template_lines
		
	def _get_surface_input(self, eqbm, nml):
		return nml
	
	def _get_gyro_input(self, eqbm, run, nml):
		return nml
	
	def make_job_files_ypi(self, scanner):
		print(f"ERROR: {self.code_name} DOES NOT SUPPORT YPI SERVERS")
		return
		
	def write_pyth_archer2(self, scanner, input_list, filename):
		pyth = open(f"{scanner.inputs['data_path']}/submit_files/{filename}/{filename}.py",'w')
		pyth.write(f"""import os
from joblib import Parallel, delayed
from time import sleep
from numpy import array

input_files = {input_list}

def start_run(run, run_attempt = 1):
	if run_attempt <= 3:
		os.system(f"echo \\\"Input: {{run}}\\\"")
		cwd = os.getcwd()
		os.chdir(f"{{run}}")
		os.system(f"$GACODE_ROOT/tglf/bin/tglf -n 1 -e .")
		os.chdir(f"{{cwd}}")
		if os.path.exists(f"{{run}}/out.tglf.run"):
			os.system(f"touch {{run}}/out.tglf.fin")
		else:
			sleep(60)
			start_run(run, run_attempt = run_attempt+1)
	else:
		print(f"ERROR: {{run}} took too many attempts to start, skipping")

Parallel(n_jobs={scanner.inputs['sbatch']['ntasks-per-node']})(delayed(start_run)(run) for run in input_files)""")
		pyth.close()
		return
		
	def make_gyro_file(self, eqbm, run, sub_dir, namelist_diff = {}):
		filename = "input.tglf"
		if not os.path.exists(f"{sub_dir}/{filename}"):
			subnml = self.get_gyro_input(eqbm=eqbm,run=run,namelist_diff=namelist_diff)
			self.write_nml(nml=subnml,directory=sub_dir,filename=filename)
		return sub_dir
	
	def save_out(self, scanner, filename = None, directory = None, specificRuns = None, QuickSave = False):
		psi_itt = scanner.single_parameters['psin'].values if 'psin' in scanner.single_parameters else scanner.dimensions['psin'].values
		equilibrium = {}
		for psiN in psi_itt:
			equilibrium[psiN] = {}
			nml = scanner.eqbm.get_surface_input(psiN)
			equilibrium[psiN]['shear'] = nml['q_prime']
			equilibrium[psiN]['beta_prime'] = nml['p_prime']
		
		if scanner['gyro']:
			gyro_data = {}
			group_data = {}
			only = set({'growth_rate','mode_frequency','ky','kx'})
			if not QuickSave:
				only = only | set({'phi','bpar','apar','time','theta','heat'})
			data_keys = ['growth_rate','mode_frequency','omega','phi','bpar','apar','epar','phi2','parity','ql_metric']
			group_keys = ['phi2_avg','t','theta', 'gds2', 'jacob','heat_flux_tot','phi2_by_kx','phi2_by_ky']
			gyro_keys = {}
			for dim in scanner.dimensions.values():
				gyro_keys[dim.name] = {}
				for val in dim.values:
					gyro_keys[dim.name][val] = set()
			
			if 'ky' not in gyro_keys.keys():
				gyro_keys['ky'] = {}
			kxs = set()
			kys = set()
			gyro_keys['kx'] = {}
			
			runs = scanner.get_all_runs() if specificRuns is None else list(specificRuns)
			for run in runs:
				sub_dir = scanner.get_run_directory(run)
				try:
					scanner.eqbm.pyro.load_gk_output(sub_dir)
					run_data = scanner.eqbm.pyro.gk_output
					group_key = run_data.attrs['object_uuid']
					group_data[group_key] = {}
					for key in group_keys:
						group_data[group_key][key] = None
					for xi, kx in enumerate(run_data['kx'].data):
						for yi, ky in enumerate(run_data['ky'].data):
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
							#gyro_data['nml_diffs'] = scanner.namelist_diffs[?]
							for key in data_keys:
								gyro_data[run_key][key] = None
							for key in only:
								try:
									key_data = run_data[key]
									if key == 'growth_rate':
										gyro_data[run_key]['growth_rate'] = float(key_data[xi,yi,-1])
									if key == 'mode_frequency':
										gyro_data[run_key]['mode_frequency'] = float(key_data[xi,yi,-1])
									elif key in ['time']:
										group_data[group_key]['t'] = array(key_data).tolist()
									elif key in ['theta']:
										group_data[group_key][key] = array(key_data).tolist()
									#elif key in ['heat']:
										#group_data[group_key][key] = array(key_data[:,:,yi,:]).tolist()
								except Exception as e:
									print(f"Save Error in {sub_dir}: {e}")
									if key == 'growth_rate':
										gyro_data[run_key]['growth_rate'] = nan
									elif key == 'mode_frequency':
										gyro_data[run_key]['mode_frequency'] = nan
										
				except Exception as e:
					print(f"Save Error {sub_dir}: {e}")
			
			existing_dim_keys = []
			for key in [x for x in scanner.inputs.inputs.keys() if 'dimension_' in x]:
				existing_dim_keys.append([x for x in key if x.isdigit()])
			dim_n = max([eval("".join(x)) for x in existing_dim_keys],default=1) + 1
			kxs = list(kxs)
			kxs.sort()
			scanner.inputs.inputs[f'dimension_{dim_n}'] = {'type': 'kx', 'values': kxs, 'min': min(kxs), 'max': max(kxs), 'num': len(kxs), 'option': None}
			if 'ky' not in scanner.dimensions:
				kys = list(kys)
				kys.sort()
				scanner.inputs.inputs[f'dimension_{dim_n+1}'] = {'type': 'ky', 'values': kys, 'min': min(kys), 'max': max(kys), 'num': len(kys), 'option': None}
			scanner.inputs.load_dimensions()

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
		with open(f"{directory}/input.tglf", "w") as f:
			for key, value in nml.items():
				f.write(f"{key} = {value}\n".upper())	


