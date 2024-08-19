from numpy import array, trapz, nan, log10, floor, append
from copy import copy
from scipy.interpolate import InterpolatedUnivariateSpline

def avg_kperp2_f(run, reader, f):
	ky = reader('ky',run)
	ky2 = ky**2
	gds2 = array(reader('gds2',run))
	J = array(reader('jacob',run))
	theta = array(reader('theta',run))
	
	df = abs(array(reader(f,run)))
	#Sometimes the magnitude of df is such that squaring the value takes it out of python's size limit for floats, but since we're integrating df2 on top&bottom we can just scale it and this change cancels itself out
	df = df/10**floor(log10(max(df)))
	df2 = df**2
	
	top_grand = ky2*gds2*J*df2
	bottom_grand = df2*J
	
	top = trapz(top_grand,theta)
	bottom = trapz(bottom_grand,theta)
	
	return top/bottom

def QL(run_ids, reader, returnlist = False):
	try:
		errors = []
		for run_id in run_ids:
			run = reader.get_run_from_id(run_id) 
			if reader('phi',run) is None or reader('apar',run) is None or reader('bpar',run) is None or str(reader('growth_rate',run)) in ['None','nan','inf','-inf'] or str(reader('ky',run)) in ['None','nan','inf','-inf']:
				errors.append(run_id)
		run_ids = [x for x in run_ids if x not in errors]
		
		kys = []
		lookup = {}
		for run_id in run_ids:
			ky = reader('ky',reader.get_run_from_id(run_id))
			kys.append(ky)
			lookup[ky] = run_id
		kys.sort()

		grs = []
		for ky in kys:
			run_id = lookup[ky]
			run = reader.get_run_from_id(run_id) 
			grs.append(reader('growth_rate',run) if reader('growth_rate',run) > 0 else 0)
		
		kp_phis = array(())
		kp_apars = array(())
		kp_bpars = array(())
		max_phis = array(())
		max_apars = array(())
		max_bpars = array(())
		for ky in kys:
			run_id = lookup[ky]
			run = reader.get_run_from_id(run_id) 
			kp_phis = append(kp_phis,avg_kperp2_f(run, reader, 'phi'))
			kp_apars = append(kp_apars,avg_kperp2_f(run, reader, 'apar'))
			kp_bpars = append(kp_bpars,avg_kperp2_f(run, reader, 'bpar'))

			max_phis = append(max_phis,max(abs(array(reader('phi',run)))))
			max_apars = append(max_apars,max(abs(array(reader('apar',run)))))
			max_bpars = append(max_bpars,max(abs(array(reader('bpar',run)))))
			
		nruns = len(kys)
		err = []
		for i in range(nruns):
			if 'nan' in [str(max_phis[i]),str(max_apars[i]),str(max_bpars[i]),str(kp_phis[i]),str(kp_apars[i]),str(kp_bpars[i]),str(grs[i])]:
				err.append(i)
		if len(err) == nruns:
			return nan, [[nan],[nan]]
				
		max_phis = array([max_phis[i] for i in range(nruns) if i not in err])
		max_apars = array([max_apars[i] for i in range(nruns) if i not in err])
		max_bpars = array([max_bpars[i] for i in range(nruns) if i not in err])
		kp_phis = array([kp_phis[i] for i in range(nruns) if i not in err])
		kp_apars = array([kp_apars[i] for i in range(nruns) if i not in err])
		kp_bpars = array([kp_bpars[i] for i in range(nruns) if i not in err])
		grs = array([grs[i] for i in range(nruns) if i not in err])
		kys = array([kys[i] for i in range(nruns) if i not in err])

		intergrand = grs*(1/kp_phis + max_apars/(max_phis*kp_apars) + max_bpars/(max_phis*kp_bpars))
		
		QL = trapz(intergrand,kys)
		if returnlist:
			return QL, [list(intergrand),list(kys)]
		else:
			return QL
	except Exception as e:
		print(e)
		if returnlist:
			return nan, [[nan],[nan]]
		else:
			return nan

def cal_flow_shear_ql(reader, run, ge):
	kys = reader['ky']
	thets = reader['theta0']
	#Ignores -ve theta0, unsure if correct
	thets = [x for x in thets if x >=0]
	ql_ky = []
	th_maxs = []
	for ky in kys:
		#print(run, ky)
		rn = copy(run)
		rn['ky'] = ky
		if reader.inputs['scan_format'] == 'point':
			thets = set()
			rids = reader.get_run_list(rn)
			for rid in rids:
				thets.add(reader.data['gyro'][rid]['theta0'])
			thets = list(thets)
			thets.sort()
		grs = []
		ql_th0 = []
		subthets = []
		for th in thets:
			rn['theta0'] = th
			gr = reader('growth_rate',rn)
			if str(gr) not in ['nan','inf','-inf','None']:
				grs.append(gr)
				ql_th0.append(reader('ql_metric',rn))
				subthets.append(th)
		if len(grs) == 0:
			ql_ky.append(0)
			continue
		gam = max(grs)
		if gam > ge/10:
			if 'shear' in run:
				shear = run['shear']
			else:
				shear = reader['shear'][0]
			th_max = (ge / (shear * gam))
			th_maxs.append(th_max)
		else:
			ql_ky.append(0)
			th_maxs.append(0)
			continue
		spline = InterpolatedUnivariateSpline(subthets,ql_th0)
		if th_max not in subthets:
			ql_th0.append(spline(th_max))
			subthets.append(th_max)
		ql_th0 = array(ql_th0)
		subthets = array(subthets)
		ql_ky.append(trapz(ql_th0[subthets<=th_max],subthets[subthets<=th_max])/th_max)
	ql = trapz(ql_ky,kys)
	reader.data['_ql_fs'] = {}
	reader.data['_ql_fs'][ge] = {}
	reader.data['_ql_fs'][ge]['theta0_maxs'] = th_maxs
	reader.data['_ql_fs'][ge]['ql_ky'] = ql_ky
	return ql, ql_ky
