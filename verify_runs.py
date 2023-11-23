from numpy import shape, array, nonzero, imag, real, nan, amax, polyfit, log, inf, full, sign, log10

class verify_scan(object):
	
	def __init__(self, reader):
		self.reader = reader
		self.scan = reader.data['gyro']
		try:
			self._nstep = reader.files['template_file']['knobs']['nstep']
		except:
			self._nstep = 100000 #Pyro default
		try:
			self._nwrite = reader.files['template_file']['gs2_diagnostics_knobs']['nwrite']
		except:
			self._nwrite = 50 #Pyro Default
		self.bad_runs = {'nstep': set(), 'nan': set(), 'phi': set(), 'apar': set(), 'bpar': set(), 'epar': set(), 'order': set(), 'negative': set()}
		self.convergence = {'converged': set(), 'converged_fit': set(), 'unconverged_stable': set(), 'unconverged': set(), 'uncalculated': set()}
		self.save_errors = {'omega': set(), 'phi2': set(), 't': set(), 'phi': set(), 'apar': set(), 'bpar': set(), 'epar': set()}
		self.nts = {}
		self.check_all()
	
	def __getitem__(self, key):
		key = key.lower()
		if key in ["bad_nstep", "badnstep", "nstep"]:
			return self.bad_runs['nstep']
		elif key in ["bad_nan", "badnan", "nan"]:
			return self.bad_runs['nan']
		elif key in ["bad_phi", "badphi", "phi"]:
			return self.bad_runs['phi']
		elif key in ["bad_apar", "badapar", "apar"]:
			return self.bad_runs['apar']
		elif key in ["bad_bpar", "badbpar", "bpar"]:
			return self.bad_runs['bpar']
		elif key in ['order']:
			return self.bad_runs['order']
		elif key in ['negative', 'neg']:
			return self.bad_runs['negative']
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
		elif key in ['saveerrors','save_errors']:
			return self.save_errors
		elif key in ['saveerrors_all','save_errors_all']:
			return self.runs_with_save_errors
		elif key in ['errors_all','all_errors']:
			return self.runs_with_errors
		elif key in ['nt','nts']:
			return self.nts
		else:
			print(f"{key} Not Found")
	
	def _all_keys(self):
		return ["bad_nstep", "badnstep", "nstep",
		"bad_nan", "badnan", "nan",
		"bad_phi", "badphi", "phi",
		"bad_apar", "badapar", "apar",
		"bad_bpar", "badbpar", "bpar",
		'order',
		'negative', 'neg',
		"unconv", "unconverged",
		'unconv_low', 'unconverged_low','unconv_stable','unconverged_stable','low_unconv','low_unconverged',
		'conv_fit', 'converged_fit','conv_f','convf','convfit','converged_fitted','conv_fitted',
		'conv', 'converged',
		'uncalculated', 'uncal',
		'saveerrors','save_errors',
		'errors_all','all_errors',
		'nt','nts',
		]
		
	def print_convergence(self):
		print(f"{len(self.convergence['unconverged'])} Unconverged Possibly Unstable | {len(self.convergence['unconverged_stable'])} Unconverged Stable | {len(self.convergence['converged'])} Converged | {len(self.convergence['converged_fit'])} Converged Fitted | {len(self.convergence['uncalculated'])} Not Calculable")
	
	@property
	def runs_with_errors(self):
		return self.convergence['unconverged'] | self.convergence['uncalculated'] | self.bad_runs['phi'] | self.bad_runs['apar'] | self.bad_runs['bpar'] | self.bad_runs['epar']
	
	@property
	def runs_with_save_errors(self):
		return self.save_errors['omega'] | self.save_errors['phi2'] | self.save_errors['t'] | self.save_errors['phi'] | self.save_errors['apar'] | self.save_errors['bpar']
	
	def print_verify(self):
		if len(self.convergence['unconverged']) > 0 or len(self.convergence['unconverged_stable']) > 0:
			print(f"Found {len(self.convergence['unconverged'])} unconverged possibly unstable runs | {len(self.convergence['unconverged_stable'])} unconverged stable runs")
		if len(self.convergence['uncalculated']) > 0:
			print(f"Found {len(self.convergence['uncalculated'])} uncalculable growth rates")
		if len(self.bad_runs['nstep']) > 0:
			print(f"Found {len(self.bad_runs['nstep'])} runs hitting the nstep limit")
		if len(self.bad_runs['nan']) > 0:
			print(f"Found {len(self.bad_runs['nan'])} runs with phi2/omega -> nan")
		if len(self.bad_runs['phi']) > 0:
			print(f"Found {len(self.bad_runs['phi'])} runs where phi does not -> 0 as n -> +/-inf")
		if len(self.bad_runs['apar']) > 0:
			print(f"Found {len(self.bad_runs['apar'])} runs where apar does not -> 0 as n -> +/-inf")
		if len(self.bad_runs['bpar']) > 0:
			print(f"Found {len(self.bad_runs['bpar'])} runs where bpar does not -> 0 as n -> +/-inf")
		if len(self.bad_runs['epar']) > 0:
			print(f"Found {len(self.bad_runs['epar'])} runs where epar does not -> 0 as n -> +/-inf")
		if len(self.runs_with_errors) > 0:
			print(f"Total: {len(self.runs_with_errors)} unique runs with errors")
		if len(self.runs_with_save_errors) > 0:
			print(f"Found {len(self.runs_with_save_errors)} runs with save errors")

	def check_all(self):
		for run in self.reader.get_all_runs():
			self.check_run(run)
		self.print_verify()
		
	def check_run(self, run):
		self.check_nan(run)
		self.check_order(run)
		self.check_convergence(run)
		self.check_nstep(run)
		self.check_phi(run)
		self.check_apar(run)
		self.check_bpar(run)
		self.check_epar(run)

	def check_nan(self, run):
		omega = self.reader('omega',run)
		phi2 = self.reader('phi2',run)
		t = self.reader('t',run)
		run_id = self.reader.get_run_id(run)
		if any([x is None for x in [omega,phi2,t]]):
			if omega is None:
				self.save_errors['omega'].add(run_id)
			if phi2 is None:
				self.save_errors['phi2'].add(run_id)
			if t is None:
				self.save_errors['t'].add(run_id)
			return
		
		if str(phi2[-1]) in ['nan','inf','-inf']:
			if str(phi2[-1]) == 'nan':
				self.bad_runs['nan'].add(run_id)
			bad_phi2_ids = [idx for idx, x in enumerate(phi2) if str(x) in ['nan','inf','-inf']]
			if len(phi2) == len(bad_phi2_ids):
				self.scan[run_id]['growth_rate'] = nan
				self.scan[run_id]['mode_frequency'] = nan
				return
			
			self.scan[run_id]['phi2'] = [x for idx, x in enumerate(phi2) if idx not in bad_phi2_ids]
			self.scan[run_id]['omega'] = [x for idx, x in enumerate(omega) if idx not in bad_phi2_ids]
			self.scan[run_id]['t'] = [x for idx, x in enumerate(t) if idx not in bad_phi2_ids]
			self.scan[run_id]['growth_rate'] = imag(omega[-1])
			self.scan[run_id]['mode_frequency'] = real(omega[-1])
			
	def check_order(self, run):
		phi2 = self.reader('phi2',run)
		t = self.reader('t',run)
		run_id = self.reader.get_run_id(run)
		if any([x is None for x in [phi2,t]]):
			self.convergence['uncalculated'].add(run_id)
			if phi2 is None:
				self.save_errors['phi2'].add(run_id)
			if t is None:
				self.save_errors['t'].add(run_id)
			return
			
		increasing = [y > x for x,y in zip(t,t[1:])]
		if not all(increasing):
			self.bad_runs['order'].add(run_id)
		positive = [x >= 0 for x in phi2]
		if not all(positive):
			self.bad_runs['negative'].add(run_id)
			self.scan[run_id]['phi2'] = [x for idx, x in enumerate(phi2) if positive[idx]]
			self.scan[run_id]['t'] = [x for idx, x in enumerate(t) if positive[idx]]
									
	def check_convergence(self, run):
		from scipy.stats import pearsonr
		omega = self.reader('omega',run)
		phi2 = self.reader('phi2',run)
		t = self.reader('t',run)
		run_id = self.reader.get_run_id(run)
		if any([x is None for x in [omega,phi2,t]]):
			if omega is None:
				self.save_errors['omega'].add(run_id)
			if phi2 is None:
				self.save_errors['phi2'].add(run_id)
			if t is None:
				self.save_errors['t'].add(run_id)
			return

		def find_gr(nt, t, phi2, gr, mf, gradgr):
			fit = polyfit(t[-nt:],log(phi2[-nt:]),1)
			fgr = fit[0]/2
			pr = abs(pearsonr(t[-nt:],log(phi2[-nt:]))[0])
			drop = log10(phi2[-1]/phi2[0])
			if gr == 0:
				diff = fgr
			else:
				diff = (fgr-gr)/gr
			if abs(diff) < 0.05 and pr > 0.99:
				return gr, mf, 'converged'
			elif pr > 0.99:
				#MARKER FOR MF SAYING POTENTIALLY WRONG
				return fgr, mf, 'converged_fit'
			elif gradgr < -0.05 and nt > 10:
				if fgr < 0:
					return fgr, nan, 'unconverged_stable'
				else:
					return gradgr, nan, 'unconverged_stable'
			elif drop < -3 and nt > 10:
				return gradgr, nan, 'unconverged_stable'
			else:
				return nan, nan, 'unconverged'
		t = [x for xi, x in enumerate(t) if phi2[xi] != 0]
		phi2 = [x for x in phi2 if x !=0]
		if len(phi2) != 0:
			gradgr = log(phi2[-1]/phi2[0])/(t[-1]-t[0])/2
			findingGR = True
			if len(phi2) > 100:
				nt = len(phi2)//10 + 10
			else:
				nt = 11
			init_nt = nt
			while findingGR:
				if nt > 11:
					nt = nt - 10
				else:
					nt = nt - 1
				try:
					new_gr, new_mf, cat = find_gr(nt = nt, t = t, phi2 = phi2, gr = self.reader('growth_rate',run), mf = self.reader('mode_frequency',run), gradgr = gradgr)
					if cat not in ['unconverged','unconverged_stable']:
						findingGR = False
					elif init_nt > 11 and nt < 110:
						findingGR = False
						nt = len(phi2)//10
					elif nt == 4:
						findingGR = False
						nt = init_nt
						
				except Exception as e:
					findingGR = False
					new_gr = nan
					new_mf = nan
					cat = 'uncalculated'
					nt = None
					print(f"ERROR ON {run_id}")
					print(e)
		else:
			new_gr = nan
			new_mf = nan
			cat = 'uncalculated'
			nt = None
		
		self.scan[run_id]['growth_rate'] = new_gr
		self.scan[run_id]['mode_frequency'] = new_mf
		self.convergence[cat].add(run_id)
		self.nts[run_id] = nt

	def check_nstep(self, run):
		omega = self.reader('omega',run)
		run_id = self.reader.get_run_id(run)
		if omega is None:
			self.save_errors['omega'].add(run_id)
			return
		lim = int(self._nstep/self._nwrite)
		if len(omega) in [lim,lim+1]:
			self.bad_runs['nstep'].add(run_id)
		
	def check_phi(self, run):
		phi = self.reader('phi',run)
		run_id = self.reader.get_run_id(run)
		if phi is None:
			self.save_errors['phi'].add(run_id)
			return
		try:
			if self.reader.inputs['grid_option'] != 'box':
				phi_max = max([abs(x) for x in phi])
				if phi_max != 0 and (abs(phi[0])/phi_max > 0.05 or abs(phi[-1])/phi_max > 0.05):
					self.bad_runs['phi'].add(run_id)
		except:
			pass

	def check_apar(self, run):
		apar = self.reader('apar',run)
		run_id = self.reader.get_run_id(run)
		if apar is None:
			self.save_errors['apar'].add(run_id)
			return
		try:
			if self.reader.inputs['grid_option'] != 'box':
				apar_max = max([abs(x) for x in apar])
				if apar_max != 0 and (abs(apar[0])/apar_max > 0.05 or abs(apar[-1])/apar_max > 0.05):
					self.bad_runs['apar'].add(run_id)
		except:
			pass
		
	def check_bpar(self, run):
		bpar = self.reader('bpar',run)
		run_id = self.reader.get_run_id(run)
		if bpar is None:
			self.save_errors['bpar'].add(run_id)
			return
		try:
			if self.reader.inputs['grid_option'] != 'box':
				bpar_max = max([abs(x) for x in bpar])
				if bpar_max != 0 and (abs(bpar[0])/bpar_max > 0.05 or abs(bpar[-1])/bpar_max > 0.05):
					self.bad_runs['bpar'].add(run_id)
		except:
			pass
	
	def check_epar(self, run):
		epar = self.reader('epar',run)
		run_id = self.reader.get_run_id(run)
		if epar is None:
			self.save_errors['epar'].add(run_id)
			return
		try:
			if self.reader.inputs['grid_option'] != 'box':
				epar_max = max([abs(x) for x in epar])
				if epar_max != 0 and (abs(epar[0])/epar_max > 0.05 or abs(epar[-1])/epar_max > 0.05):
					self.bad_runs['epar'].add(run_id)
		except:
			pass
