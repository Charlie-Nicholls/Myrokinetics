from numpy import array, trapz, nan

def avg_kperp2_f(indexes,data,f,kys):
	p,i,j,k = indexes
	ky2 = kys[k]**2
	gds2 = array(data['gds2'][p][i][j][k])
	J = array(data['jacob'][p][i][j][k])
	df2 = abs(array(data[f][p][i][j][k]))**2
	theta = array(data['theta'][p][i][j][k])
	
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
		kp_phis = array([avg_kperp2_f((p,i,j,k), data,'phi',kys) for k in range(len(kys))])
		kp_apars = array([avg_kperp2_f((p,i,j,k), data,'apar',kys) for k in range(len(kys))])
		kp_bpars = array([avg_kperp2_f((p,i,j,k), data,'bpar',kys) for k in range(len(kys))])
		max_phis = array([max(abs(array(data['phi'][p][i][j][k]))**2) for k in range(len(kys))])
		max_apars = array([max(abs(array(data['apar'][p][i][j][k]))**2) for k in range(len(kys))])
		max_bpars = array([max(abs(array(data['bpar'][p][i][j][k]))**2) for k in range(len(kys))])
		
		intergrand = grs*(1/kp_phis + max_apars/(max_phis*kp_apars) + max_bpars/(max_phis*kp_bpars))
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
