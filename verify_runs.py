from numpy import shape, array, nonzero, imag, real, nan, amax, polyfit, log

class verify_scan(object):
	
	def __init__(self, scan = None):
		self.scan = scan
		self.new_data = {'gra': scan['data']['growth_rates_all'], 'mfa': scan['data']['mode_frequencies_all']}
		self.bad_runs = {'omega': None, 'other': None, 'unconv': None, 'unconv_low': None, 'phi': None, 'apar': None}
		self.save_errors = {'omega': set(), 'phi2': set(), 'time': set(), 'phi': set(), 'apar': set()}
		self.check_all()
	
	def __getitem__(self, key):
		key = key.lower()
		if key in ["gra","gr_a","gr_all","growth_rates_all"]:
			return self.new_data['gra']
		elif key in ["mfa", "mf_a", "mf_all", "mode_frequencies_all"]:
			return self.new_data['mfa']
		elif key in ["bad_nstep", "badnstep", "nstep"]:
			return self.bad_runs['nstep']
		elif key in ["bad_other", "badother", "other"]:
			return self.bad_runs['other']
		elif key in ["bad_phi", "badphi", "phi"]:
			return self.bad_runs['phi']
		elif key in ["bad_apar", "badapar", "apar"]:
			return self.bad_runs['apar']
		elif key in ["unconv", "unconverged"]:
			return self.bad_runs["unconv"]
		elif key in ['unconv_low', 'unconverged_low','unconv_stable','unconverged_stable','low_unconv','low_unconverged']:
			return self.bad_runs['unconv_low']
		elif key in ['old_gra', 'gra_old', 'old_growth_rates_all', 'growth_rates_all_old', 'old_gr_a', 'gr_a_old']:
			return self.scan['growth_rates_all']
		elif key in ['old_mfa', 'mfa_old', 'old_mode_frequencies_all', 'mode_frequencies_all_old', 'old_mf_a', 'mf_a_old']:
			return self.scan['mode_frequencies_all']
		elif key in ['saveerrors','save_errors']:
			return self.save_errors
		elif key in ['saveerrors_all','save_errors_all']:
			return self.save_errors['omega'] | self.save_errors['phi2'] | self.save_errors['time']
		
	def check_all(self):
		if self.scan['data']['phi2'] is None or self.scan['data']['omega'] is None:
			print("Cannot Verify Runs For Quick Save")
			return
		self.check_convergence()
		self.check_nstep()
		self.check_phi()
		self.check_apar()
		self.check_other()
		self.print_verify()
		
	def print_verify(self):

		print(f"Found {len(self.bad_runs['unconv'])} Unconverged Possibly Unstable Runs and {len(self.bad_runs['unconv_low'])} Unconverged Stable Runs")
		print(f"Found {len(self.bad_runs['nstep'])} Runs Hitting The nstep Limit")
		print(f"Found {len(self.bad_runs['other'])} Runs With omega -> nan")
		print(f"Found {len(self.bad_runs['phi'])} Runs Where phi does not -> 0 as n -> +/-inf")
		print(f"Found {len(self.bad_runs['apar'])} Runs Where apar does not -> 0 as n -> +/-inf")
		runs_with_save_errors = self.save_errors['omega'] | self.save_errors['phi2'] | self.save_errors['time']
		print(f"Found {len(runs_with_save_errors)} Runs With Save Errors")
	
	def check_convergence(self):
		if self.scan['data']['omega'] is None:
			print("ERROR: No omega data")
			return
		if self.scan['data']['phi2'] is None:
			print("ERROR: No phi2 data")
			return
		if self.scan['data']['time'] is None:
			print("ERROR: No time data")
			return
		
		
		sha = shape(array(self.scan['data']['omega'],dtype=object))
		save_errors  = set()
		conv = []
		unconv = []
		unconv_low = []
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						phi2 = self.scan['data']['phi2'][i][j][k][l]
						time = self.scan['data']['time'][i][j][k][l]
						omega = self.scan['data']['omega'][i][j][k][l]
						if omega is None or phi2 is None or time is None:
							if omega is None:
								self.save_errors['omega'] | [i,j,k,l]
							if phi2 is None:
								self.save_errors['phi2'] | [i,j,k,l]
							if time is None:
								self.save_errors['time'] | [i,j,k,l]
							break
						if phi2[-1] == 0:
							phi2 = array(phi2)[nonzero(phi2)].tolist()
							omega = omega[:len(phi2)]
							time = time[:len(phi2)]
							self.new_data['gra'][i][j][k][l] = imag(omega[-1])
							self.new_data['mfa'][i][j][k][l] = real(omega[-1])
						if len(phi2) > 20:
							gr = imag(omega[-1])
							fit = polyfit(time[-10:],log(phi2[-10:]),1)
							if abs((fit[0]-gr)/gr) < 0.05:
								conv.append([i,j,k,l])
							else:
								self.new_data['mfa'][i][j][k][l] = nan
								fit2 = polyfit(time,log(phi2),1)
								if fit2[0] < -0.1:
									self.new_data['gra'][i][j][k][l] = fit2[0]
									unconv_low.append([i,j,k,l])
								else:
									self.new_data['gra'][i][j][k][l] = nan
									unconv.append([i,j,k,l])
						else:
							conv.append([i,j,k,l])
		
		self.bad_runs['unconv'] = unconv
		self.bad_runs['unconv_low'] = unconv_low
			
	def check_nstep(self):
		if self.scan['data']['omega'] is None:
			print("ERROR: No omega data")
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
		save_errors  = set()
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['omega'][i][j][k][l] is None:
							save_errors.add((i,j,k,l))
						elif len(self.scan['data']['omega'][i][j][k][l]) == lim:
							bad_nstep.append([i,j,k,l])
							self.new_data['gra'][i][j][k][l] = nan
							self.new_data['mfa'][i][j][k][l] = nan
		
		self.save_errors['omega'] = self.save_errors['omega'] | save_errors
		self.bad_runs['nstep'] = bad_nstep
		
	def check_phi(self):
		if self.scan['data']['phi'] is None:
			print("ERROR: No phi data")
			return
		sha = shape(array(self.scan['data']['phi'],dtype=object))
		bad_phi = []
		save_errors  = set()
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['phi'][i][j][k][l] is None:
							save_errors.add((i,j,k,l))
						phi_max = max([abs(x) for x in self.scan['data']['phi'][i][j][k][l]])
						if abs(self.scan['data']['phi'][i][j][k][l][0])/phi_max > 0.05 or abs(self.scan['data']['phi'][i][j][k][l][-1])/phi_max > 0.05:
							bad_phi.append([i,j,k,l])
							
		self.save_errors['phi'] = self.save_errors['phi'] | save_errors
		self.bad_runs['phi'] = bad_phi
	
	def check_apar(self):
		if self.scan['data']['apar'] is None:
			print("ERROR: No apar data")
			return
		sha = shape(array(self.scan['data']['apar'],dtype=object))
		bad_apar = []
		save_errors  = set()
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['apar'][i][j][k][l] is None:
							save_errors.add((i,j,k,l))
						apar_max = max([abs(x) for x in self.scan['data']['apar'][i][j][k][l]])
						if abs(self.scan['data']['apar'][i][j][k][l][0])/apar_max > 0.05 or abs(self.scan['data']['apar'][i][j][k][l][-1])/apar_max > 0.05:
							bad_apar.append([i,j,k,l])
							
		self.save_errors['apar'] = self.save_errors['apar'] | save_errors
		self.bad_runs['apar'] = bad_apar
	
	def check_other(self):
		if self.scan['data']['omega'] is None:
			print("ERROR: No omega data")
			return
		sha = shape(array(self.scan['data']['omega'],dtype=object))
		bad_other = []
		save_errors  = set()
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['omega'][i][j][k][l] is None:
							save_errors.add((i,j,k,l))
						elif str(self.scan['data']['omega'][i][j][k][l][-1]) == '(nan+nanj)':
							bad_other.append([i,j,k,l])
			
		self.save_errors['omega'] = self.save_errors['omega'] | save_errors
		self.bad_runs['other'] = bad_other
