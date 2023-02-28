from numpy import shape, array, nonzero, imag, real, nan, amax, polyfit, log, inf, full, sign

class verify_scan(object):
	
	def __init__(self, scan = None):
		self.scan = scan
		self.new_data = {'gra': scan['data']['growth_rates_all'], 'mfa': scan['data']['mode_frequencies_all']}
		self.bad_runs = {'nstep': set(), 'nan': set(), 'phi': set(), 'apar': set(), 'bpar': set()}
		self.convergence = {'converged': set(), 'converged_fit': set(), 'unconverged_stable': set(), 'unconverged': set(), 'uncalculated': set()}
		self.nts = full((len(scan['inputs']['psiNs']),scan['inputs']['n_beta'],scan['inputs']['n_shat'],len(scan['inputs']['aky_values'])),None).tolist()
		self.save_errors = {'omega': set(), 'phi2': set(), 'time': set(), 'phi': set(), 'apar': set(), 'bpar': set()}
		self.check_all()
	
	def __getitem__(self, key):
		key = key.lower()
		if key in ["gra","gr_a","gr_all","growth_rates_all"]:
			return self.new_data['gra']
		elif key in ["mfa", "mf_a", "mf_all", "mode_frequencies_all"]:
			return self.new_data['mfa']
		elif key in ["bad_nstep", "badnstep", "nstep"]:
			return self.bad_runs['nstep']
		elif key in ["bad_nan", "badnan", "nan"]:
			return self.bad_runs['nan']
		elif key in ["bad_phi", "badphi", "phi"]:
			return self.bad_runs['phi']
		elif key in ["bad_apar", "badapar", "apar"]:
			return self.bad_runs['apar']
		elif key in ["bad_bpar", "badbpar", "bpar"]:
			return self.bad_runs['bpar']
		elif key in ["unconv", "unconverged"]:
			return self.convergence["unconverged"]
		elif key in ['unconv_low', 'unconverged_low','unconv_stable','unconverged_stable','low_unconv','low_unconverged']:
			return self.convergence['unconverged_stable']
		elif key in ['conv_fit', 'converged_fit','conv_f','convf','convfit']:
			return self.convergence['converged_fit']
		elif key in ['conv', 'converged']:
			return self.convergence['converged']
		elif key in ['uncalculated', 'uncal']:
			return self.convergence['uncalculated']
		elif key in ['old_gra', 'gra_old', 'old_growth_rates_all', 'growth_rates_all_old', 'old_gr_a', 'gr_a_old']:
			return self.scan['growth_rates_all']
		elif key in ['old_mfa', 'mfa_old', 'old_mode_frequencies_all', 'mode_frequencies_all_old', 'old_mf_a', 'mf_a_old']:
			return self.scan['mode_frequencies_all']
		elif key in ['saveerrors','save_errors']:
			return self.save_errors
		elif key in ['saveerrors_all','save_errors_all']:
			return self.runs_with_save_errors
		elif key in ['errors_all','all_errors']:
			return self.runs_with_errors
		else:
			print(f"{key} Not Found")
	
	def _all_keys(self):
		return ["gra","gr_a","gr_all","growth_rates_all",
		"mfa", "mf_a", "mf_all", "mode_frequencies_all",
		"bad_nstep", "badnstep", "nstep",
		"bad_nan", "badnan", "nan",
		"bad_phi", "badphi", "phi",
		"bad_apar", "badapar", "apar",
		"bad_bpar", "badbpar", "bpar",
		"unconv", "unconverged",
		'unconv_low', 'unconverged_low','unconv_stable','unconverged_stable','low_unconv','low_unconverged',
		'conv_fit', 'converged_fit','conv_f','convf','convfit','converged_fitted','conv_fitted',
		'conv', 'converged',
		'uncalculated', 'uncal',
		'old_gra', 'gra_old', 'old_growth_rates_all', 'growth_rates_all_old', 'old_gr_a', 'gr_a_old',
		'old_mfa', 'mfa_old', 'old_mode_frequencies_all', 'mode_frequencies_all_old', 'old_mf_a', 'mf_a_old',
		'saveerrors','save_errors',
		'errors_all','all_errors']
		
	def check_all(self):
		if self.scan['data']['phi2'] is None or self.scan['data']['omega'] is None:
			print("Cannot Verify Runs For Quick Save")
			return
		self.check_nan()
		self.check_convergence()
		self.check_nstep()
		self.check_phi()
		self.check_apar()
		self.check_bpar()
		self.print_verify
	
	@property
	def print_convergence(self):
		print(f"{len(self.convergence['unconverged'])} Unconverged Possibly Unstable | {len(self.convergence['unconverged_stable'])} Unconverged Stable | {len(self.convergence['converged'])} Converged | {len(self.convergence['converged_fit'])} Converged Fitted | {len(self.convergence['uncalculated'])} Not Calculable")
	
	@property
	def runs_with_errors(self):
		return self.convergence['unconverged'] | self.bad_runs['nstep'] | self.bad_runs['nan'] | self.bad_runs['phi'] | self.bad_runs['apar'] | self.bad_runs['apar'] | self.bad_runs['bpar']
	
	@property
	def runs_with_save_errors(self):
		return self.save_errors['omega'] | self.save_errors['phi2'] | self.save_errors['time'] | self.save_errors['phi'] | self.save_errors['apar'] | self.save_errors['bpar']
	
	@property
	def print_verify(self):
		if len(self.convergence['unconverged']) > 0 or len(self.convergence['unconverged_stable']) > 0:
			print(f"Found {len(self.convergence['unconverged'])} Unconverged Possibly Unstable Runs and {len(self.convergence['unconverged_stable'])} Unconverged Stable Runs")
		if len(self.bad_runs['nstep']) > 0:
			print(f"Found {len(self.bad_runs['nstep'])} Runs Hitting The nstep Limit")
		if len(self.bad_runs['nan']) > 0:
			print(f"Found {len(self.bad_runs['nan'])} Runs With phi2/omega -> nan")
		if len(self.bad_runs['phi']) > 0:
			print(f"Found {len(self.bad_runs['phi'])} Runs Where phi does not -> 0 as n -> +/-inf")
		if len(self.bad_runs['apar']) > 0:
			print(f"Found {len(self.bad_runs['apar'])} Runs Where apar does not -> 0 as n -> +/-inf")
		if len(self.bad_runs['bpar']) > 0:
			print(f"Found {len(self.bad_runs['bpar'])} Runs Where bpar does not -> 0 as n -> +/-inf")
		if len(self.runs_with_errors) > 0:
			print(f"Total: {len(self.runs_with_errors)} Unique Runs With Errors")
		if len(self.runs_with_save_errors) > 0:
			print(f"Found {len(self.runs_with_save_errors)} Runs With Save Errors")
	
	def check_nan(self):
		if self.scan['data']['omega'] is None:
			print("ERROR: No omega data")
			return
		sha = shape(array(self.scan['data']['omega'],dtype=object))
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['omega'][i][j][k][l] is None or self.scan['data']['phi2'][i][j][k][l] is None or self.scan['data']['time'][i][j][k][l] is None:
							if self.scan['data']['omega'][i][j][k][l] is None:
								self.save_errors['omega'].add((i,j,k,l))
							if self.scan['data']['phi2'][i][j][k][l] is None:
								self.save_errors['phi2'].add((i,j,k,l))
							if self.scan['data']['time'][i][j][k][l] is None:
								self.save_errors['time'].add((i,j,k,l))
						elif str(self.scan['data']['phi2'][i][j][k][l][-1]) in ['nan','inf','0.0','0']:
							if str(self.scan['data']['phi2'][i][j][k][l][-1]) == 'nan':
								self.bad_runs['nan'].add((i,j,k,l))
							self.scan['data']['phi2'][i][j][k][l] = [x for x in self.scan['data']['phi2'][i][j][k][l] if str(x) not in ['nan','inf','0.0','0']]
							self.scan['data']['omega'][i][j][k][l] = self.scan['data']['omega'][i][j][k][l][:len(self.scan['data']['phi2'][i][j][k][l])]
							self.scan['data']['time'][i][j][k][l] = self.scan['data']['time'][i][j][k][l][:len(self.scan['data']['phi2'][i][j][k][l])]
							self.new_data['gra'][i][j][k][l] = imag(self.scan['data']['omega'][i][j][k][l][-1])
							self.new_data['mfa'][i][j][k][l] = real(self.scan['data']['omega'][i][j][k][l][-1])
						
	def check_convergence(self):
		from scipy.stats import pearsonr
		if self.scan['data']['omega'] is None:
			print("ERROR: No omega data")
			return
		if self.scan['data']['phi2'] is None:
			print("ERROR: No phi2 data")
			return
		if self.scan['data']['time'] is None:
			print("ERROR: No time data")
			return
		
		def find_gr(nt, time, phi2, gr, mf, gradgr):
			fit = polyfit(time[-nt:],log(phi2[-nt:]),1)
			fgr = fit[0]/2
			pr = abs(pearsonr(time[-nt:],log(phi2[-nt:]))[0])
			if abs((fgr-gr)/gr) < 0.05 and pr > 0.99:
				return gr, mf, 'converged'

			elif pr > 0.99:
				#MARKER FOR MF SAYING POTENTIALLY WRONG
				return fgr, mf, 'converged_fit'
			elif gradgr < -0.05 and nt > 10:
				if fgr < 0:
					return fgr, nan, 'unconverged_stable'
				else:
					return gradgr, nan, 'unconverged_stable'
			elif phi2[-1]/phi2[0] < 0.001 and nt > 10:
				return gradgr, nan, 'unconverged_stable'
			else:
				return nan, nan, 'unconverged'
		
		sha = shape(array(self.scan['data']['omega'],dtype=object))
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						phi2 = self.scan['data']['phi2'][i][j][k][l]
						time = self.scan['data']['time'][i][j][k][l]
						omega = self.scan['data']['omega'][i][j][k][l]
						if omega is None or phi2 is None or time is None:
							if omega is None:
								self.save_errors['omega'].add((i,j,k,l))
							if phi2 is None:
								self.save_errors['phi2'].add((i,j,k,l))
							if time is None:
								self.save_errors['time'].add((i,j,k,l))
						else:		
							if phi2[-1] == 0:
								phi2 = array(phi2)[nonzero(phi2)].tolist()
								omega = omega[:len(phi2)]
								time = time[:len(phi2)]
								self.new_data['gra'][i][j][k][l] = imag(omega[-1])
								self.new_data['mfa'][i][j][k][l] = real(omega[-1])
							
							if 'nan' in [str(x) for x in phi2] or inf in phi2:
								self.convergence['uncalculated'].add((i,j,k,l))

							else:
									
								gr = imag(omega[-1])
								mf = real(omega[-1])
								gradgr = log(phi2[-1]/phi2[0])/(time[-1]-time[0])/2
								findingGR = True
								if len(phi2) > 100:
									nt = len(phi2)//10 + 10
									while findingGR:
										nt = nt - 10
										new_gr, new_mf, cat = find_gr(nt = nt, time = time, phi2 = phi2, gr = gr, mf = mf, gradgr = gradgr)
										
										if cat not in ['unconverged','unconverged_stable'] or nt < 110:
											findingGR = False
								else:
									
									nt = 11
									while findingGR:
										nt = nt - 1
										if nt > len(phi2):
											nt = len(phi2)
										new_gr, new_mf, cat =  find_gr(nt = nt, time = time, phi2 = phi2, gr = gr, mf = mf, gradgr = gradgr)
										if cat != 'unconverged' or nt == 2:
											findingGR = False
								
								self.new_data['gra'][i][j][k][l] = new_gr
								self.new_data['mfa'][i][j][k][l] = new_mf
								self.convergence[cat].add((i,j,k,l))
								self.nts[i][j][k][l] = nt

	def check_nstep(self):
		if self.scan['data']['omega'] is None:
			print("ERROR: No omega data")
			return
		nstep = nwrite = None
		if type(self.scan['files']['template_file']) == list:
			lines = self.scan['files']['template_file']
			title = ''
			for line in lines:
				if line[0] == '&':
					title = line[1:].strip(" \t\n")
				if title == 'knobs' and 'nstep' in line:
					nstep = eval(line.split("=")[1])
				if title == 'gs2_diagnostic_knobs' and 'nwrite' in line:
					nwrite = eval(line.split("=")[1])
		else:
			nstep = self.scan['files']['template_file']['knobs']['nstep']
			nwrite = self.scan['files']['template_file']['gs2_diagnostics_knobs']['nwrite']
		
		if not nstep:
			nstep = 100000 #Pyro default
		if not nwrite:
			nwrite = 50 #Pyro Default
		lim = int(nstep/nwrite)
		sha = shape(array(self.scan['data']['omega'],dtype=object))
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['omega'][i][j][k][l] is None:
							self.save_errors['omega'].add((i,j,k,l))
						elif len(self.scan['data']['omega'][i][j][k][l]) == lim:
							self.bad_runs['nstep'].add((i,j,k,l))
							#self.new_data['gra'][i][j][k][l] = nan
							#self.new_data['mfa'][i][j][k][l] = nan
		
	def check_phi(self):
		if self.scan['data']['phi'] is None:
			print("ERROR: No phi data")
			return
		sha = shape(array(self.scan['data']['phi'],dtype=object))
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['phi'][i][j][k][l] is None:
							self.save_errors['phi'].add((i,j,k,l))
						else:
							phi_max = max([abs(x) for x in self.scan['data']['phi'][i][j][k][l]])
							if phi_max != 0 and (abs(self.scan['data']['phi'][i][j][k][l][0])/phi_max > 0.05 or abs(self.scan['data']['phi'][i][j][k][l][-1])/phi_max > 0.05):
								self.bad_runs['phi'].add((i,j,k,l))

	def check_apar(self):
		if self.scan['data']['apar'] is None:
			print("ERROR: No apar data")
			return
		sha = shape(array(self.scan['data']['apar'],dtype=object))
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['apar'][i][j][k][l] is None:
							self.save_errors['apar'].add((i,j,k,l))
						else:
							apar_max = max([abs(x) for x in self.scan['data']['apar'][i][j][k][l]])
							if apar_max != 0 and (abs(self.scan['data']['apar'][i][j][k][l][0])/apar_max > 0.05 or abs(self.scan['data']['apar'][i][j][k][l][-1])/apar_max > 0.05):
								self.bad_runs['apar'].add((i,j,k,l))
		
	def check_bpar(self):
		if self.scan['data']['bpar'] is None:
			print("ERROR: No bpar data")
			return
		sha = shape(array(self.scan['data']['bpar'],dtype=object))
		for i in range(sha[0]):
			for j in range(sha[1]):
				for k in range(sha[2]):
					for l in range(sha[3]):
						if self.scan['data']['bpar'][i][j][k][l] is None:
							self.save_errors['bpar'].add((i,j,k,l))
						else:
							bpar_max = max([abs(x) for x in self.scan['data']['bpar'][i][j][k][l]])
							if bpar_max != 0 and (abs(self.scan['data']['bpar'][i][j][k][l][0])/bpar_max > 0.05 or abs(self.scan['data']['bpar'][i][j][k][l][-1])/bpar_max > 0.05):
								self.bad_runs['bpar'].add((i,j,k,l))
	
	
