import os
from ..codes import codes

class ypi_server():

	save_modules = """"""

	sbatch = {
		}

	save_sbatch = {
		}
	
	requires_slurm_save = False
	
	def make_job_files(self, scanner, n_jobs = None, n_par = None, n_sim = None):
		codes[scanner.inputs['gk_code']].make_job_files_ypi(scanner)
		return
	
	def make_ideal_job_files(self, scanner, n_jobs = None, n_par = 1, n_sim = None):
		codes[scanner.inputs['gk_code']].make_ideal_job_files_ypi(scanner)
		return	
				
