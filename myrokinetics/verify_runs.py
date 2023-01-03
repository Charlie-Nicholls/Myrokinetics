from numpy import *

class verify_scan(object):
	
	def __init__(self, scan = None):
		self.scan = scan
		self.new_data = {'gra': scan['data']['growth_rates_all'], 'mfa': scan['data']['mode_frequencies_all']}
		self.bad_runs = {}
		self.check_all()
	
	def __getitem__(self, key):
		key = key.lower()
		if key in ["gra","gr_a","gr_all","growth_rates_all"]:
			return self.new_data['gra']
		elif key in ["mfa", "mf_a", "mf_all", "mode_frequencies_all"]:
			return self.new_data['mfa']
		elif key in ["bad_phi2", "badphi2", "phi2"]:
			return self.bad_runs['phi2']
		elif key in ["bad_nstep", "badnstep", "nstep"]:
			return self.bad_runs['nstep']
		elif key in ['old_gra', 'gra_old', 'old_growth_rates_all', 'growth_rates_all_old', 'old_gr_a', 'gr_a_old']:
			return self.scan['growth_rates_all']
		elif key in ['old_mfa', 'mfa_old', 'old_mode_frequencies_all', 'mode_frequencies_all_old', 'old_mf_a', 'mf_a_old']:
			return self.scan['mode_frequencies_all']
		
	def check_all(self):
		self.check_phi2()
		self.check_nstep()
		
	def check_phi2(self):
		if self.scan['data']['phi2'] is None:
			print("No phi2 data: Use detailed_save for this functionality")
			return
		data = self.scan['data']
		sha = shape(array(data['phi2'],dtype=object))
		bad_phi2 = []
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if data['phi2'][i][j][k][l][-1] == 0.0:
							bad_phi2.append([i,j,k,l])
		if bad_phi2:
			print(f"Found {len(bad_phi2)} Bad phi2 Runs")
			for p,i,j,k in bad_phi2:
				self.new_data['mfa'][p][i][j][k] = nan
				grnew = amin(array(self.scan['data']['growth_rates_all'])[p,:,:,k])
				if grnew < 0:
					self.new_data['gra'][p][i][j][k] = grnew
				else:
					self.new_data['gra'][p][i][j][k] = -amax(array(self.scan['data']['growth_rates_all'])[p,:,:,k])
				
		self.bad_runs['phi2'] = bad_phi2
			
	def check_nstep(self):
		if self.scan['data']['omega'] is None:
			print("No omega data: Use detailed_save for this functionality")
			return
		nstep = nwrite = None
		lines = self.scan['files']['template_file']
		if lines:
			for line in lines:
				if line[0] == '&':
					title = line[1:].strip(" \t\n")
				if title == 'knobs' and 'nstep' in line:
					nstep = eval(line.split("=")[1])
				if title == 'gs2_diagnostic_knobs' and 'nwrite' in line:
					nwrite = eval(line.split("=")[1])
		if not nstep:
			nstep = 100000 #Pyro default
		if not nwrite:
			nwrite = 50 #Pyro Default
		lim = int(nstep/nwrite)
		sha = shape(array(self.scan['data']['omega'],dtype=object))
		bad_nstep = []
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if len(self.scan['data']['omega'][i][j][k][l]) == lim:
							bad_nstep.append([i,j,k,l])
		if bad_nstep:
			print(f"Found {len(bad_nstep)} Bad nstep Runs")
			for p,i,j,k in bad_nstep:
				self.new_data['gra'][p][i][j][k] = nan
				self.new_data['mfa'][p][i][j][k] = nan
			self.bad_runs['nstep'] = bad_nstep
		return
