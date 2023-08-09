complile_viking = """module purge
module load tools/git
module load compiler/ifort
module load mpi/impi
module load numlib/FFTW
module load data/netCDF/4.6.1-intel-2018b
module load data/netCDF-Fortran/4.4.4-intel-2018b
module load numlib/imkl/2018.3.222-iimpi-2018b
module load lang/Python/3.7.0-intel-2018b
export GK_SYSTEM=viking
export MAKEFLAGS=-IMakefiles
ulimit -s unlimited"""

complile_archer2 = """module load PrgEnv-gnu
module load cray-hdf5 cray-netcdf cray-fftw
export GK_SYSTEM=archer2
export MAKEFLAGS=-IMakefiles
ulimit -s unlimited"""

save_viking = """module load lang/Python/3.7.0-intel-2018b
module swap lang/Python lang/Python/3.10.4-GCCcore-11.3.0
source $HOME/pyroenv2/bin/activate"""

save_archer2 = ""

system_modules = {'viking': {'compile': complile_viking, 'save': save_viking},
	'archer2': {'compile': complile_archer2, 'save': save_archer2},
	}
