from pathlib import Path

template_dir = Path(__file__).parent / "templates"
template_dir.resolve()
template_dir = str(template_dir)

gs2_template = "template.gs2"

modules = """module purge
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
export PATH=$PATH:$HOME/gs2/bin
ulimit -s unlimited"""

save_modules = """module load lang/Python/3.7.0-intel-2018b
module swap lang/Python lang/Python/3.10.4-GCCcore-11.3.0
source $HOME/pyroenv2/bin/activate"""
