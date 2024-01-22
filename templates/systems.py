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
export PATH=${PATH}:${HOME}/gs2/bin"""

viking_save_modules = """module load Python/3.10.8-GCCcore-12.2.0"""

viking_sbatch = {
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

viking_save_sbatch = {
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

archer2_modules = """module load PrgEnv-gnu
module load cray-hdf5 cray-netcdf cray-fftw cray-python
export GK_SYSTEM=archer2
export MAKEFLAGS=-IMakefiles
ulimit -s unlimited
export PATH=${PATH}:/work/e281/e281/cnicholls/gs2/bin
source /work/e281/e281/cnicholls/pythenv/bin/activate"""

archer2_save_modules = """module load PrgEnv-gnu
module load cray-python
source /work/e281/e281/cnicholls/pythenv/bin/activate"""

archer2_sbatch = {
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

archer2_save_sbatch = {
	'time': '00:10:00',
	'job-name': 'myro',
	'ntasks': 1,
	'mem': '1gb',
	'output': 'save_out.slurm',
	'account': 'e281',
	'partition': 'serial',
	'qos': 'serial',
	}
	
systems = {'viking': {'modules': viking_modules, 'save_modules': viking_save_modules, 'sbatch': viking_sbatch, 'save_sbatch': viking_save_sbatch},
	'archer2': {'modules': archer2_modules, 'save_modules': archer2_save_modules, 'sbatch': archer2_sbatch, 'save_sbatch': archer2_save_sbatch},
	}
