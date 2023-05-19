from numpy import array, trapz, nan, log10, floor

def avg_kperp2_f(indexes,data,f,ky):
	p,i,j,k = indexes
	ky2 = ky**2
	gds2 = array(data['gds2'][p][i][j][k])
	J = array(data['jacob'][p][i][j][k])
	theta = array(data['theta'][p][i][j][k])
	
	df = abs(array(data[f][p][i][j][k]))
	#Sometimes the magnitude of df is such that squaring the value takes it out of python's size limit for floats, but since we're integrating df2 on top&bottom we can just scale it and this change cancels itself out
	df = df/10**floor(log10(max(df)))
	df2 = df**2
	
	top_grand = ky2*gds2*J*df2
	bottom_grand = df2*J
	
	top = trapz(top_grand,theta)
	bottom = trapz(bottom_grand,theta)
	
	return top/bottom

def QL(indexes,data,kys):
	p,i,j = indexes
	try:
		nruns = len(kys)
		kys = array(kys)
		grs = array([data['growth_rates_all'][p][i][j][k] for k in range(nruns)])
		kp_phis = array([avg_kperp2_f((p,i,j,k), data,'phi',ky) for k, ky in enumerate(kys)])
		kp_apars = array([avg_kperp2_f((p,i,j,k), data,'apar',ky) for k, ky in enumerate(kys)])
		kp_bpars = array([avg_kperp2_f((p,i,j,k), data,'bpar',ky) for k, ky in enumerate(kys)])

		max_phis = array([max(abs(array(data['phi'][p][i][j][k]))) for k in range(nruns)])
		max_apars = array([max(abs(array(data['apar'][p][i][j][k]))) for k in range(nruns)])
		max_bpars = array([max(abs(array(data['bpar'][p][i][j][k]))) for k in range(nruns)])

		#Sometimes the values are such that squaring takes it out of python's size limit for floats, we can scale all the values and it will cancel itself out
		for k in range(nruns):
			scale = floor(log10(max_phis[k]))
			max_phis[k] = max_phis[k]/10**scale
			max_apars[k] = max_apars[k]/10**scale
			max_bpars[k] = max_bpars[k]/10**scale

		max_phis2 = max_phis**2
		max_apars2 = max_apars**2
		max_bpars2 = max_bpars**2

		err = set()
		for k in range(nruns):
			if 'nan' in [str(max_phis2[k]),str(max_apars2[k]),str(max_bpars2[k]),str(kp_phis[k]),str(kp_apars[k]),str(kp_bpars[k]),str(grs[k])]:
				err.add(k)
		max_phis2 = array([max_phis2[k] for k in range(nruns) if k not in err])
		max_apars2 = array([max_apars2[k] for k in range(nruns) if k not in err])
		max_bpars2 = array([max_bpars2[k] for k in range(nruns) if k not in err])
		kp_phis = array([kp_phis[k] for k in range(nruns) if k not in err])
		kp_apars = array([kp_apars[k] for k in range(nruns) if k not in err])
		kp_bpars = array([kp_bpars[k] for k in range(nruns) if k not in err])
		grs = array([grs[k] for k in range(nruns) if k not in err])
		kys = array([kys[k] for k in range(nruns) if k not in err])

		intergrand = grs*(1/kp_phis + max_apars2/(max_phis2*kp_apars) + max_bpars2/(max_phis2*kp_bpars))
		QL = trapz(intergrand,kys)

		return QL
	except:
		return nan

def mset_avg_kperp2_f(data,f):
	ky2 = data['ky']**2
	gds2 = data['gds2']
	J = data['jacob']
	df2 = abs(data[f])**2
	theta = data['theta']
	
	top_grand = ky2*gds2*J*df2
	bottom_grand = df2*J
	
	top = trapz(top_grand,theta)
	bottom = trapz(bottom_grand,theta)
	
	return top/bottom
	
def mset_QL(mset):
	nruns = len(mset.runs)
	kys = array([mset[i]['data']['ky'][0] for i in range(nruns)])
	grs = array([mset[i]['data'].gr for i in range(nruns)])
	kp_phis = array([mset_avg_kperp2_f(mset[i]['data'],'phi') for i in range(nruns)])
	kp_apars = array([mset_avg_kperp2_f(mset[i]['data'],'apar') for i in range(nruns)])
	kp_bpars = array([mset_avg_kperp2_f(mset[i]['data'],'bpar') for i in range(nruns)])
	max_phis = array([max(abs(mset[i]['data']['phi'])**2) for i in range(nruns)])
	max_apars = array([max(abs(mset[i]['data']['apar'])**2) for i in range(nruns)])
	max_bpars = array([max(abs(mset[i]['data']['bpar'])**2) for i in range(nruns)])
	
	intergrand = grs*(1/kp_phis + max_apars/(max_phis*kp_apars) + max_bpars/(max_phis*kp_bpars))
	QL = trapz(intergrand,kys)
	
	return QL


def QL_nonan(indexes,data,kys):
	p,i,j = indexes
	nruns = len(kys)
	kys = array(kys)
	grs = array([data['growth_rates_all'][p][i][j][k] for k in range(nruns)])
	kp_phis = array([avg_kperp2_f((p,i,j,k), data,'phi',ky) for k, ky in enumerate(kys)])
	kp_apars = array([avg_kperp2_f((p,i,j,k), data,'apar',ky) for k, ky in enumerate(kys)])
	kp_bpars = array([avg_kperp2_f((p,i,j,k), data,'bpar',ky) for k, ky in enumerate(kys)])

	max_phis = array([max(abs(array(data['phi'][p][i][j][k]))) for k in range(nruns)])
	max_apars = array([max(abs(array(data['apar'][p][i][j][k]))) for k in range(nruns)])
	max_bpars = array([max(abs(array(data['bpar'][p][i][j][k]))) for k in range(nruns)])

	#Sometimes the values are such that squaring takes it out of python's size limit for floats, we can scale all the values and it will cancel itself out
	for k in range(nruns):
		scale = floor(log10(max_phis[k]))
		max_phis[k] = max_phis[k]/10**scale
		max_apars[k] = max_apars[k]/10**scale
		max_bpars[k] = max_bpars[k]/10**scale

	max_phis2 = max_phis**2
	max_apars2 = max_apars**2
	max_bpars2 = max_bpars**2

	err = set()
	for k in range(nruns):
		if 'nan' in [str(max_phis2[k]),str(max_apars2[k]),str(max_bpars2[k]),str(kp_phis[k]),str(kp_apars[k]),str(kp_bpars[k]),str(grs[k])]:
			err.add(k)
	max_phis2 = array([max_phis2[k] for k in range(nruns) if k not in err])
	max_apars2 = array([max_apars2[k] for k in range(nruns) if k not in err])
	max_bpars2 = array([max_bpars2[k] for k in range(nruns) if k not in err])
	kp_phis = array([kp_phis[k] for k in range(nruns) if k not in err])
	kp_apars = array([kp_apars[k] for k in range(nruns) if k not in err])
	kp_bpars = array([kp_bpars[k] for k in range(nruns) if k not in err])
	grs = array([grs[k] for k in range(nruns) if k not in err])
	kys = array([kys[k] for k in range(nruns) if k not in err])

	intergrand = grs*(1/kp_phis + max_apars2/(max_phis2*kp_apars) + max_bpars2/(max_phis2*kp_bpars))
	QL = trapz(intergrand,kys)

	return QL
