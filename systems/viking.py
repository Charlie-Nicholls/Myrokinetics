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
		input_dirs = {}
		for n in range(n_par):
			input_dirs[n] = []	
		if n_jobs is None or n_jobs*n_par > len(scanner._input_dirs):
			total_jobs = len(scanner._input_dirs)
		else:
			total_jobs = n_jobs*n_par
		from numpy import ceil
		if ceil(total_jobs/n_par) > 10000:
			print(f"Viking supports a max of 10,000 jobs per array submission (Currently requesting {ceil(total_jobs/n_par)})")
			return
		dir_list = list(scanner._input_dirs)
		for i in range(total_jobs):
			input_dirs[i%n_par].append(dir_list[i])
			scanner._input_dirs.remove(dir_list[i])
		for n in range(n_par):
			filename = f"gyro_{n}"
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/{filename}",exist_ok=True)
			sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"{filename}/{scanner.inputs['sbatch']['output']}_0")
			sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}",f"{filename}/{scanner.inputs['sbatch']['error']}_0")
			dirlist = open(f"{scanner.inputs['data_path']}/submit_files/gyro_{n}/{filename}.txt",'w')
			for indir in input_dirs[n]:
				dirlist.write(f"{indir}\n")
			dirlist.close()
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/{filename}/{filename}.job",'w')
			if len(input_dirs[n]) > 1:
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['output']}_0",f"{scanner.inputs['sbatch']['output']}_%a")
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}_0",f"{scanner.inputs['sbatch']['error']}_%a")
				jobfile.write(f"{sbatch_n}")
				jobfile.write(f"\n#SBATCH --array=1-{len(input_dirs[n])}\n")
			else:
				jobfile.write(f"{sbatch_n}")
			code_jobfile = codes[scanner.inputs['gk_code']].jobfile_viking
			jobfile.write(f"""
{compile_modules}

INDIR=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {scanner.inputs['data_path']}/submit_files/gyro_{n}/gyro_{n}.txt)

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
		return
	
	def make_ideal_job_files(self, scanner, n_jobs = None, n_par = 1, n_sim = None):
		compile_modules = codes[scanner.inputs['gk_code']].viking_modules
		sbatch = "#!/bin/bash"
		for key, val in scanner.inputs['sbatch'].items():
			if key == 'output' and '/' not in val:
				val = f"{scanner.inputs['data_path']}/submit_files/{val}"
			sbatch = sbatch + f"\n#SBATCH --{key}={val}"
		n_sim = n_par if n_sim is None else n_sim
		os.makedirs(f"{scanner.inputs['data_path']}/submit_files/",exist_ok=True)
		input_dirs = {}
		for n in range(n_par):
			input_dirs[n] = []
		if n_jobs is None or n_jobs*n_par > len(scanner._ideal_input_dirs):
			total_jobs = len(scanner._ideal_input_dirs)
		else:
			total_jobs = n_jobs*n_par
		dir_list = list(scanner._ideal_input_dirs)
		for i in range(total_jobs):
			input_dirs[i%n_par].append(dir_list[i])
			scanner._ideal_input_dirs.remove(dir_list[i])
		for n in range(n_par):
			os.makedirs(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}",exist_ok=True)
			sbatch_n = sbatch.replace(f"{scanner.inputs['sbatch']['output']}",f"ideal_{n}/{scanner.inputs['sbatch']['output']}_0")
			sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}",f"ideal_{n}/{scanner.inputs['sbatch']['error']}_0")
			sbatch_n = sbatch_n.replace(f"--cpus-per-task={scanner.inputs['sbatch']['cpus-per-task']}","--cpus-per-task=1")
			filename = f"ideal_{n}"
			dirlist = open(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}/{filename}.txt",'w')
			for indir in input_dirs[n]:
				dirlist.write(f"{indir}\n")
			dirlist.close()
			jobfile = open(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}/{filename}.job",'w')
			if len(input_dirs[n]) > 1:
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['output']}_0", f"{scanner.inputs['sbatch']['output']}_%a")
				sbatch_n = sbatch_n.replace(f"{scanner.inputs['sbatch']['error']}_0", f"{scanner.inputs['sbatch']['error']}_%a")
				jobfile.write(f"{sbatch_n}")
				jobfile.write(f"\n#SBATCH --array=1-{len(input_dirs[n])}\n")
			else:
				jobfile.write(f"{sbatch_n}")
			jobfile.write(f"""
{compile_modules}

INDIR=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {scanner.inputs['data_path']}/submit_files/ideal_{n}/ideal_{n}.txt)
echo "${{INDIR}}/ideal_ball.in"
ideal_ball "${{INDIR}}/ideal_ball.in"
if test -f "${{INDIR}}/ideal_ball.ballstab_2d"; then
touch "${{INDIR}}/run.fin"
fi""")
			jobfile.close()
		for n in range(n_sim):
			scanner._ideal_jobs.add(f"{scanner.inputs['data_path']}/submit_files/ideal_{n}/ideal_{n}.job")			
		return	
				
