import os
from ..codes import codes
	
class archer2():

	save_modules = """module load PrgEnv-gnu
module load cray-python
source /work/e281/e281/cnicholls/pythenv/bin/activate"""

	sbatch = {
		'time': '01:00:00',
		'job-name': 'myro',
		'nodes': 1,
		'output': 'myro.slurm',
		'ntasks-per-node': 128,
		'cpus-per-task': 1,
		'account': 'e281',
		'partition': 'standard',
		'qos': 'standard',
		'distribution': 'block:block',
		'hint': 'nomultithread',
		}

	save_sbatch = {
		'time': '00:10:00',
		'job-name': 'myro',
		'ntasks': 1,
		'mem': '1gb',
		'output': 'save_out.slurm',
		'account': 'e281',
		'partition': 'serial',
		'qos': 'serial',
		}
	
	requires_slurm_save = True
	
	def make_job_files(self, scanner, n_jobs = None, n_par = 1, n_sim = None):
		compile_modules = codes[scanner.inputs['gk_code']].archer2_modules
		sbatch = "#!/bin/bash"
		for key, val in scanner.inputs['sbatch'].items():
			if key == 'output' and '/' not in val:
				val = f"{scanner.inputs['data_path']}/submit_files/{val}"
			if key == 'ntasks-per-node' and 'ntasks' in scanner.inputs['sbatch'].keys():
                                break
			sbatch = sbatch + f"\n#SBATCH --{key}={val}"
			
		if scanner.inputs['grid_option'] == 'single':
		
			if n_sim is None:
				n_sim = n_par if n_par < 8 else 8
			if n_sim > 8:
				print("Archer supports a maximum of n_sim = 8")
				n_sim = 8
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
			input_dirs = {}
			for n in range(n_par):
				input_dirs[n] = []
			if n_jobs is None or n_jobs*n_par > len(scanner._input_dirs):
				total_jobs = len(scanner._input_dirs)
			else:
				total_jobs = n_jobs*n_par
			dir_list = list(scanner._input_dirs)
			for i in range(total_jobs):
				input_dirs[i%n_par].append(dir_list[i])
				scanner._input_dirs.remove(dir_list[i])
			for n in range(n_par):
				filename = f"gyro_{n}"
				os.makedirs(f"{scanner.inputs['data_path']}/submit_files/{filename}",exist_ok=True)
				sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"{filename}/{scanner.inputs['sbatch']['output']}_{n}")
			
				codes[scanner.inputs['gk_code']].write_pyth_archer2(input_dirs[n], filename)
			
				jobfile = open(f"{scanner.inputs['data_path']}/submit_files/{filename}/{filename}.job",'w')
				jobfile.write(f"""{sbatch_n}

{compile_modules}

python {scanner.inputs['data_path']}/submit_files/{filename}/{filename}.py""")
			
				if n_par > n_sim and n + n_sim < n_par:
					jobfile.write(f"\nsbatch {scanner.inputs['data_path']}/submit_files/gyro_{n+n_sim}/gyro_{n+n_sim}.job")
				jobfile.close()
				
			for n in range(n_sim):
				scanner._jobs.add(f"{scanner.inputs['data_path']}/submit_files/gyro_{n}/gyro_{n}.job")
				
		if scanner.inputs['grid_option'] == 'box':
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
			if scanner.inputs['sbatch']['cpus-per-task'] > 1:
				compile_modules += f"\nexport OMP_NUM_THREADS={scanner.inputs['sbatch']['cpus-per-task']}"
			
			run_code = codes[scanner.inputs['gk_code']].get_box_archer2()
			
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/submit.job",'w')
			jobfile.write(f"""{sbatch}

{compile_modules}

{run_code}""")
			jobfile.close()
			scanner._jobs.add(f"{scanner.inputs['data_path']}/submit_files/submit.job")
		
		return
	
	def make_ideal_job_files(self, scanner, n_jobs = None, n_par = 1, n_sim = None):
		compile_modules = codes[scanner.inputs['gk_code']].viking_modules
		sbatch = "#!/bin/bash"
		for key, val in scanner.inputs['sbatch'].items():
			if key == 'output' and '/' not in val:
				val = f"{scanner.inputs['data_path']}/submit_files/{val}"
			sbatch = sbatch + f"\n#SBATCH --{key}={val}"
		if n_sim is None:
			n_sim = n_par if n_par < 8 else 8
		if n_sim > 8:
			print("Archer supports a maximum of n_sim = 8")
			n_sim = 8
		os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
		input_dirs = {}
		for n in range(n_par):
			input_dirs[n] = []
		if n_jobs is None or n_jobs*n_par > len(scanner._ideal_dirs):
			total_jobs = len(scanner._ideal_dirs)
		else:
			total_jobs = n_jobs*n_par
		dir_list = list(scanner._ideal_input_dirs)
		for i in range(total_jobs):
			input_dirs[i%n_par].append(dir_list[i])
			scanner._ideal_input_dirs.remove(dir_list[i])
		for n in range(n_par):
			sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"{scanner.inputs['sbatch']['output']}_ideal_{n}")
			sbatch_n = sbatch_n.replace(f"#SBATCH --nodes = {scanner.inputs['sbatch']['nodes']}","#SBATCH --nodes = 1")
			filename = f"ideal_{n}"
			pyth = open(f"{scanner.inputs['data_path']}/submit_files/{filename}.py",'w')
			pyth.write(f"""import os
from joblib import Parallel, delayed
from time import sleep

input_dirs = {input_dirs[n]}

def start_run(run, run_attempt = 1):
if run_attempt <= 3:
	os.system(f"echo \\\"Ideal Input: {{run}}/ideal.gs2\\\"")
	os.system(f"srun --nodes=1 --ntasks=1 ideal_ball \\\"{{run}}/ideal.gs2\\\"")
	if os.path.exists(f\"{{run}}/ideal.ballstab_2d\"):
		os.system(f"touch \\\"{{run}}/ideal.fin\\\"")
	else:
		sleep(60)
		start_run(run, run_attempt = run_attempt+1)
else:
	print(f"ERROR: {{run}} took too many attempts to start, skipping")

Parallel(n_jobs={scanner.inputs['sbatch']['ntasks-per-node']})(delayed(start_run)(run) for run in input_dirs)""")
			pyth.close()
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/{filename}.job",'w')
			jobfile.write(f"""{sbatch_n}

{compile_modules}

which gs2
gs2 --build-config

python {scanner.inputs['data_path']}/submit_files/{filename}.py""")
			if n_par > n_sim and n + n_sim < n_par:
				jobfile.write(f"\nsbatch {scanner.inputs['data_path']}/submit_files/ideal_{n+n_sim}.job")
			jobfile.close()
		for n in range(n_sim):
			scanner._ideal_jobs.add(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}.job")	
		return			
