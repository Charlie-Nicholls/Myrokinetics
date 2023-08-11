from numpy import array, trapz, nan, log10, floor, append

def avg_kperp2_f(run,f):
	ky = run['ky']
	ky2 = ky**2
	gds2 = array(run['gds2'])
	J = array(run['jacob'])
	theta = array(run['theta'])
	
	df = abs(array(run[f]))
	#Sometimes the magnitude of df is such that squaring the value takes it out of python's size limit for floats, but since we're integrating df2 on top&bottom we can just scale it and this change cancels itself out
	df = df/10**floor(log10(max(df)))
	df2 = df**2
	
	top_grand = ky2*gds2*J*df2
	bottom_grand = df2*J
	
	top = trapz(top_grand,theta)
	bottom = trapz(bottom_grand,theta)
	
	return top/bottom

def QL(run_ids,gyro_data):
	#try:
	errors = []
	for run_id in run_ids:
		run = gyro_data[run_id]
		if run['phi'] is None or run['apar'] is None or run['bpar'] is None or str(run['growth_rate']) in ['None','nan','inf','-inf'] or str(run['ky']) in ['None','nan','inf','-inf']:
			errors.append(run_id)
	run_ids = [x for x in run_ids if x not in errors]
	
	kys = []
	lookup = {}
	for run_id in run_ids:
		ky = gyro_data[run_id]['ky']
		kys.append(ky)
		lookup[ky] = run_id
	kys.sort()

	grs = []
	for ky in kys:
		run_id = lookup[ky]
		grs.append(gyro_data[run_id]['growth_rate'] if gyro_data[run_id]['growth_rate'] > 0 else 0)
	
	kp_phis = array(())
	kp_apars = array(())
	kp_bpars = array(())
	max_phis = array(())
	max_apars = array(())
	max_bpars = array(())
	for ky in kys:
		run_id = lookup[ky]
		run = gyro_data[run_id]
		kp_phis = append(kp_phis,avg_kperp2_f(run, 'phi'))
		kp_apars = append(kp_apars,avg_kperp2_f(run, 'apar'))
		kp_bpars = append(kp_bpars,avg_kperp2_f(run, 'bpar'))

		max_phis = append(max_phis,max(abs(array(run['phi']))))
		max_apars = append(max_apars,max(abs(array(run['apar']))))
		max_bpars = append(max_bpars,max(abs(array(run['bpar']))))

	#Sometimes the values are such that squaring takes it out of python's size limit for floats, we can scale all the values and it will cancel itself out
	
	nruns = len(kys)
	for i in range(nruns):
		scale = floor(log10(max_phis[i]))
		max_phis[i] = max_phis[i]/10**scale
		max_apars[i] = max_apars[i]/10**scale
		max_bpars[i] = max_bpars[i]/10**scale

	max_phis2 = max_phis**2
	max_apars2 = max_apars**2
	max_bpars2 = max_bpars**2

	err = []
	for i in range(nruns):
		if 'nan' in [str(max_phis2[i]),str(max_apars2[i]),str(max_bpars2[i]),str(kp_phis[i]),str(kp_apars[i]),str(kp_bpars[i]),str(grs[i])]:
			err.append(i)
	if len(err) == nruns:
		return nan
			
	max_phis2 = array([max_phis2[i] for i in range(nruns) if i not in err])
	max_apars2 = array([max_apars2[i] for i in range(nruns) if i not in err])
	max_bpars2 = array([max_bpars2[i] for i in range(nruns) if i not in err])
	kp_phis = array([kp_phis[i] for i in range(nruns) if i not in err])
	kp_apars = array([kp_apars[i] for i in range(nruns) if i not in err])
	kp_bpars = array([kp_bpars[i] for i in range(nruns) if i not in err])
	grs = array([grs[i] for i in range(nruns) if i not in err])
	kys = array([kys[i] for i in range(nruns) if i not in err])

	intergrand = grs*(1/kp_phis + max_apars2/(max_phis2*kp_apars) + max_bpars2/(max_phis2*kp_bpars))
	
	QL = trapz(intergrand,kys)
	print(intergrand,kys,QL)
	return QL
	#except:
		#return nan
