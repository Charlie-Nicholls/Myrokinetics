import os
from ..codes import codes

class viking():

	save_modules = """module load Python/3.10.8-GCCcore-12.2.0"""

	sbatch = {
		'job-name': 'myro',
		'partition': 'nodes',
		'time': '01:00:00',
		'cpus-per-task': 1,
		'account': 'pet-gspt-2019',
		'mail-type': None,
		'mail-user': 'cn762@york.ac.uk',
		'output': 'myro.slurm',
		'error': 'myro.err',
		}

	save_sbatch = {
		'time': '00:10:00',
		'job-name': 'myro',
		'partition': 'nodes',
		'ntasks': 1,
		'mem': '1gb',
		'output': 'save_out.slurm',
		'account': 'pet-gspt-2019',
		'mail-type': None,
		'mail-user': 'cn762@york.ac.uk',
		}
	
	requires_slurm_save = True
	
	def make_job_files(self, scanner, n_jobs = None, n_par = None, n_sim = None):
		compile_modules = codes[scanner.inputs['gk_code']].viking_modules
		sbatch = "#!/bin/bash"
		for key, val in scanner.inputs['sbatch'].items():
			if key == 'output' and '/' not in val:
				val = f"{scanner.inputs['data_path']}/submit_files/{val}"
			sbatch = sbatch + f"\n#SBATCH --{key}={val}"

		n_sim = n_par if n_sim is None else n_sim
		os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
		input_lists = {}
		for n in range(n_par):
			input_lists[n] = []	
		if n_jobs is None or n_jobs*n_par > len(scanner._input_files):
			total_jobs = len(scanner._input_files)
		else:
			total_jobs = n_jobs*n_par
		from numpy import ceil
		if ceil(total_jobs/n_par) > 10000:
			print(f"Viking supports a max of 10,000 jobs per array submission (Currently requesting {ceil(total_jobs/n_par)})")
			return
		input_list = list(scanner._input_files)
		for i in range(total_jobs):
			input_lists[i%n_par].append(input_list[i])
			scanner._input_files.remove(input_list[i])
		for n in range(n_par):
			filename = f"gyro_{n}"
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/{filename}",exist_ok=True)
			sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"{filename}/{scanner.inputs['sbatch']['output']}_0")
			sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}",f"{filename}/{scanner.inputs['sbatch']['error']}_0")
			inlist = open(f"{scanner.inputs['data_path']}/submit_files/gyro_{n}/{filename}.txt",'w')
			for infile in input_lists[n]:
				inlist.write(f"{infile[:-3]}\n")
			inlist.close()
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/{filename}/{filename}.job",'w')
			if len(input_lists[n]) > 1:
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['output']}_0",f"{scanner.inputs['sbatch']['output']}_%a")
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}_0",f"{scanner.inputs['sbatch']['error']}_%a")
				jobfile.write(f"{sbatch_n}")
				jobfile.write(f"\n#SBATCH --array=1-{len(input_lists[n])}\n")
			else:
				jobfile.write(f"{sbatch_n}")
			code_jobfile = codes[scanner.inputs['gk_code']].jobfile_viking
			jobfile.write(f"""
{compile_modules}

INFILE=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {scanner.inputs['data_path']}/submit_files/gyro_{n}/gyro_{n}.txt)
echo "${{INFILE}}.in"

{code_jobfile}

""")
			jobfile.close()
		submit = open(f"{scanner.inputs['data_path']}/submit_files/submit.job",'w')
		submit.write(f"""#!/bin/bash
#SBATCH --job-name=submit
#SBATCH --partition=nodes
#SBATCH --time=00:00:30
#SBATCH --ntasks=1
#SBATCH --mem=10MB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --account={scanner.inputs['sbatch']['account']}
#SBATCH --output={scanner.inputs['data_path']}/submit_files/submit.slurm


""")
		for n in range(n_sim):
			submit.write(f"ID_{n}=$(sbatch --parsable \"{scanner.inputs['data_path']}/submit_files/gyro_{n}/gyro_{n}.job\")\n")
		if n_par > n_sim:
			for n in range(n_sim, n_par):
				submit.write(f"ID_{n}=$(sbatch --parsable --dependency=afterany:$ID_{n-n_sim} \"{scanner.inputs['data_path']}/submit_files/gyro_{n}/gyro_{n}.job\")\n")
		submit.close()
		scanner._jobs.add(f"{scanner.inputs['data_path']}/submit_files/submit.job")
	
	def make_ideal_job_files(self, scanner, n_jobs = None, n_par = 1, n_sim = None):
		compile_modules = codes[scanner.inputs['gk_code']].viking_modules
		sbatch = "#!/bin/bash"
		for key, val in scanner.inputs['sbatch'].items():
			if key == 'output' and '/' not in val:
				val = f"{scanner.inputs['data_path']}/submit_files/{val}"
			sbatch = sbatch + f"\n#SBATCH --{key}={val}"
		n_sim = n_par if n_sim is None else n_sim
		os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
		input_lists = {}
		for n in range(n_par):
			input_lists[n] = []
		if n_jobs is None or n_jobs*n_par > len(scanner._ideal_input_files):
			total_jobs = len(scanner._ideal_input_files)
		else:
			total_jobs = n_jobs*n_par
		input_list = list(scanner._ideal_input_files)
		for i in range(total_jobs):
			input_lists[i%n_par].append(input_list[i])
			scanner._ideal_input_files.remove(input_list[i])
		for n in range(n_par):
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}",exist_ok=True)
			sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"ideal_{n}/{scanner.inputs['sbatch']['output']}_0")
			sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}",f"ideal_{n}/{scanner.inputs['sbatch']['error']}_0")
			sbatch_n = sbatch_n.replace(f"--cpus-per-task={scanner.inputs['sbatch']['cpus-per-task']}","--cpus-per-task=1")
			filename = f"ideal_{n}"
			inlist = open(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}/{filename}.txt",'w')
			for infile in input_lists[n]:
				inlist.write(f"{infile[:-3]}\n")
			inlist.close()
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}/{filename}.job",'w')
			if len(input_lists[n]) > 1:
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['output']}_0", f"{scanner.inputs['sbatch']['output']}_%a")
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}_0", f"{scanner.inputs['sbatch']['error']}_%a")
				jobfile.write(f"{sbatch_n}")
				jobfile.write(f"\n#SBATCH --array=1-{len(input_lists[n])}\n")
			else:
				jobfile.write(f"{sbatch_n}")
			jobfile.write(f"""
{compile_modules}

which gs2
gs2 --build-config

INFILE=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {scanner.inputs['data_path']}/submit_files/ideal_{n}/ideal_{n}.txt)
echo "${{INFILE}}.in"
ideal_ball "${{INFILE}}.in"
if test -f "${{INFILE}}.ballstab_2d"; then
touch "${{INFILE}}.fin"
fi""")
			jobfile.close()
		for n in range(n_sim):
			scanner._ideal_jobs.add(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}/ideal_{n}.job")			
			
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
		'account': 'e281-ypimcf',
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
		'account': 'e281-ypimcf',
		'partition': 'serial',
		'qos': 'serial',
		}
	
	requires_slurm_save = True
	
	def make_job_files(self, scanner, codes, n_jobs = None, n_par = 1, n_sim = None):
		compile_modules = codes[scanner.inputs['gk_code']].archer2_modules
		sbatch = "#!/bin/bash"
		for key, val in scanner.inputs['sbatch'].items():
			if key == 'output' and '/' not in val:
				val = f"{scanner.inputs['data_path']}/submit_files/{val}"
			sbatch = sbatch + f"\n#SBATCH --{key}={val}"
			
		if scanner.inputs['non_linear'] == False:
		
			if n_sim is None:
				n_sim = n_par if n_par < 8 else 8
			if n_sim > 8:
				print("Archer supports a maximum of n_sim = 8")
				n_sim = 8
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
			input_lists = {}
			for n in range(n_par):
				input_lists[n] = []
			if n_jobs is None or n_jobs*n_par > len(scanner._input_files):
				total_jobs = len(scanner._input_files)
			else:
				total_jobs = n_jobs*n_par
			input_list = list(scanner._input_files)
			for i in range(total_jobs):
				input_lists[i%n_par].append(input_list[i])
				scanner._input_files.remove(input_list[i])
			for n in range(n_par):
				filename = f"gyro_{n}"
				os.makedirs(f"{scanner.inputs['data_path']}/submit_files/{filename}",exist_ok=True)
				sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"{filename}/{scanner.inputs['sbatch']['output']}_{n}")
			
			codes[scanner.inputs['gk_code']].write_pyth_archer2(scanner, input_list)
			
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/{filename}/{filename}.job",'w')
			jobfile.write(f"""{sbatch_n}

{compile_modules}

python {scanner.inputs['data_path']}/submit_files/{filename}/{filename}.py &

wait""")
			
			if n_par > n_sim and n + n_sim < n_par:
				jobfile.write(f"\nsbatch {scanner.inputs['data_path']}/submit_files/gyro_{n+n_sim}/gyro_{n+n_sim}.job")
			jobfile.close()
				
			for n in range(n_sim):
				scanner._jobs.add(f"{scanner.inputs['data_path']}/submit_files/gyro_{n}/gyro_{n}.job")
				
		if scanner.inputs['non_linear'] == True:
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
			if scanner.inputs['sbatch']['cpus-per-task'] > 1:
				compile_modules += f"\nexport OMP_NUM_THREADS={scanner.inputs['sbatch']['cpus-per-task']}"
			ntasks = scanner.inputs['sbatch']['ntasks'] if 'ntasks' in scanner.inputs['sbatch'] else scanner.inputs['sbatch']['nodes']*scanner.inputs['sbatch']['ntasks-per-node']
			
			run_code = codes[scanner.inputs['gk_code']].get_non_linear_archer2()
			
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/submit.job",'w')
			jobfile.write(f"""{sbatch}

{compile_modules}

{run_code}""")
			jobfile.close()
			scanner._jobs.add(f"{scanner.inputs['data_path']}/submit_files/submit.job")
		
		return
	
	def make_ideal_job_files(self, scanner, codes, n_jobs = None, n_par = 1, n_sim = None):
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
		input_lists = {}
		for n in range(n_par):
			input_lists[n] = []
		if n_jobs is None or n_jobs*n_par > len(scanner._ideal_input_files):
			total_jobs = len(scanner._ideal_input_files)
		else:
			total_jobs = n_jobs*n_par
		input_list = list(scanner._ideal_input_files)
		for i in range(total_jobs):
			input_lists[i%n_par].append(input_list[i])
			scanner._ideal_input_files.remove(input_list[i])
		for n in range(n_par):
			sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"{scanner.inputs['sbatch']['output']}_ideal_{n}")
			sbatch_n = sbatch_n.replace(f"#SBATCH --nodes = {scanner.inputs['sbatch']['nodes']}","#SBATCH --nodes = 1")
			filename = f"ideal_{n}"
			pyth = open(f"{scanner.inputs['data_path']}/submit_files/{filename}.py",'w')
			pyth.write(f"""import os
from joblib import Parallel, delayed
from time import sleep

input_files = {input_lists[n]}

def start_run(run, run_attempt = 1):
if run_attempt <= 3:
	os.system(f"echo \\\"Ideal Input: {{run}}\\\"")
	os.system(f"srun --nodes=1 --ntasks=1 ideal_ball \\\"{{run}}\\\"")
	if os.path.exists(f\"{{run[:-3]}}.ballstab_2d\"):
		os.system(f"touch \\\"{{run[:-3]}}.fin\\\"")
	else:
		sleep(60)
		start_run(run, run_attempt = run_attempt+1)
else:
	print(f"ERROR: {{run}} took too many attempts to start, skipping")

Parallel(n_jobs={scanner.inputs['sbatch']['ntasks-per-node']})(delayed(start_run)(run) for run in input_files)""")
			pyth.close()
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/{filename}.job",'w')
			jobfile.write(f"""{sbatch_n}

{compile_modules}

which gs2
gs2 --build-config

python {scanner.inputs['data_path']}/submit_files/{filename}.py &

wait""")
			if n_par > n_sim and n + n_sim < n_par:
				jobfile.write(f"\nsbatch {scanner.inputs['data_path']}/submit_files/ideal_{n+n_sim}.job")
			jobfile.close()
		for n in range(n_sim):
			scanner._ideal_jobs.add(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}.job")	
		return
			
				
				
				
				
				
				
				
				
				
