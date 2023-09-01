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
		self.save_errors = {'omega': set(), 'phi2': set(), 't': set(), 'phi': set(), 'apar': set(), 'bpar': set()}
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
		for run_id in self.scan.keys():
			self.check_run(run_id)
		self.print_verify()
		
	def check_run(self, run_id):
		self.check_nan(run_id)
		self.check_order(run_id)
		self.check_convergence(run_id)
		self.check_nstep(run_id)
		self.check_phi(run_id)
		self.check_apar(run_id)
		self.check_bpar(run_id)
		self.check_epar(run_id)

	def check_nan(self, run_id):
		data = self.scan[run_id]
		if any([x not in data for x in ['omega','phi2','t']]):
			return
		if any([data[x] is None for x in ['omega','phi2','t']]):
			if data['omega'] is None:
				self.save_errors['omega'].add(run_id)
			if data['phi2'] is None:
				self.save_errors['phi2'].add(run_id)
			if data['t'] is None:
				self.save_errors['t'].add(run_id)
			return
	
		if str(data['phi2'][-1]) in ['nan','inf','-inf','0.0','0']:
			if str(data['phi2'][-1]) == 'nan':
				self.bad_runs['nan'].add(run_id)
			bad_phi2_ids = [idx for idx, x in enumerate(data['phi2']) if str(x) in ['nan','inf','-inf','0.0','0']]
			if len(data['phi2']) == len(bad_phi2_ids):
				data['phi2'] = None
				data['omega'] = None
				data['t'] = None
				data['growth_rates_all'] = nan
				data['mode_frequencies_all'] = nan
				return
			
			data['phi2'] = [x for idx, x in enumerate(data['phi2']) if idx not in bad_phi2_ids]
			data['omega'] = [x for idx, x in enumerate(data['omega']) if idx not in bad_phi2_ids]
			data['t'] = [x for idx, x in enumerate(data['t']) if idx not in bad_phi2_ids]
			data['growth_rate'] = imag(data['omega'][-1])
			data['mode_frequency'] = real(data['omega'][-1])
		self.scan[run_id] = data
			
	def check_order(self, run_id):
		data = self.scan[run_id]
		if any([x not in data for x in ['phi2','t']]):
			self.convergence['uncalculated'].add(run_id)
			return
		if any([data[x] is None for x in ['phi2','t']]):
			if data['phi2'] is None:
				self.save_errors['phi2'].add(run_id)
			if data['t'] is None:
				self.save_errors['t'].add(run_id)
			self.convergence['uncalculated'].add(run_id)
			return
			
		increasing = [y > x for x,y in zip(data['t'],data['t'][1:])]
		if not all(increasing):
			self.bad_runs['order'].add(run_id)
		positive = [x >= 0 for x in data['phi2']]
		if not all(positive):
			self.bad_runs['negative'].add(run_id)
			self.scan[run_id]['phi2'] = [x for idx, x in enumerate(self.scan[run_id]['phi2']) if positive[idx]]
			self.scan[run_id]['t'] = [x for idx, x in enumerate(self.scan[run_id]['t']) if positive[idx]]
									
	def check_convergence(self, run_id):
		from scipy.stats import pearsonr
		data = self.scan[run_id]
		if any([x not in data for x in ['omega','phi2','t']]):
			return
		if any([data[x] is None for x in ['omega','phi2','t']]):
			if data['omega'] is None:
				self.save_errors['omega'].add(run_id)
			if data['phi2'] is None:
				self.save_errors['phi2'].add(run_id)
			if data['t'] is None:
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
		
		gradgr = log(data['phi2'][-1]/data['phi2'][0])/(data['t'][-1]-data['t'][0])/2
		findingGR = True
		if len(data['phi2']) > 100:
			nt = len(data['phi2'])//10 + 10
		else:
			nt = 11
		init_nt = nt
		while findingGR:
			if nt > 11:
				nt = nt - 10
			else:
				nt = nt - 1
			try:
				new_gr, new_mf, cat = find_gr(nt = nt, t = data['t'], phi2 = data['phi2'], gr = data['growth_rate'], mf = data['mode_frequency'], gradgr = gradgr)
				if cat not in ['unconverged','unconverged_stable']:
					findingGR = False
				elif init_nt > 11 and nt < 110:
					findingGR = False
					nt = len(data['phi2'])//10
				elif nt == 2:
					findingGR = False
					nt = len(data['phi2'])
					
			except Exception as e:
				findingGR = False
				new_gr = nan
				new_mf = nan
				cat = 'uncalculated'
				print(f"ERROR ON {run_id}")
				print(e)
		
		self.scan[run_id]['growth_rate'] = new_gr
		self.scan[run_id]['mode_frequency'] = new_mf
		self.convergence[cat].add(run_id)
		self.nts[run_id] = nt

	def check_nstep(self, run_id):
		data = self.scan[run_id]
		if 'omega' not in data:
			return
		if data['omega'] is None:
			self.save_errors['omega'].add(run_id)
			return
		lim = int(self._nstep/self._nwrite)
		if len(data['omega']) in [lim,lim+1]:
			self.bad_runs['nstep'].add(run_id)
		
	def check_phi(self, run_id):
		data = self.scan[run_id]
		if 'phi' not in data:
			return
		if data['phi'] is None:
			self.save_errors['phi'].add(run_id)
			return
		try:
			phi_max = max([abs(x) for x in data['phi']])
			if phi_max != 0 and (abs(data['phi'][0])/phi_max > 0.05 or abs(data['phi'][-1])/phi_max > 0.05):
				self.bad_runs['phi'].add(run_id)
		except:
			pass

	def check_apar(self, run_id):
		data = self.scan[run_id]
		if 'apar' not in data:
			return
		if data['apar'] is None:
			self.save_errors['apar'].add(run_id)
			return
		try:
			apar_max = max([abs(x) for x in data['apar']])
			if apar_max != 0 and (abs(data['apar'][0])/apar_max > 0.05 or abs(data['apar'][-1])/apar_max > 0.05):
				self.bad_runs['apar'].add(run_id)
		except:
			pass
		
	def check_bpar(self, run_id):
		data = self.scan[run_id]
		if 'bpar' not in data:
			return
		if data['bpar'] is None:
			self.save_errors['bpar'].add(run_id)
			return
		try:
			bpar_max = max([abs(x) for x in data['bpar']])
			if bpar_max != 0 and (abs(data['bpar'][0])/bpar_max > 0.05 or abs(data['bpar'][-1])/bpar_max > 0.05):
				self.bad_runs['bpar'].add(run_id)
		except:
			pass
	
	def check_epar(self, run_id):
		data = self.scan[run_id]
		if 'epar' not in data:
			return
		if data['epar'] is None:
			self.save_errors['epar'].add(run_id)
			return
		try:
			epar_max = max([abs(x) for x in data['epar']])
			if epar_max != 0 and (abs(data['epar'][0])/epar_max > 0.05 or abs(data['epar'][-1])/epar_max > 0.05):
				self.bad_runs['epar'].add(run_id)
		except:
			pass
