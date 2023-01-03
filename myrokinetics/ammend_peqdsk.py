def AmmendPEQDSK(peq_file = None, geq_file = None, peq = None, geq = None):
	if geq_file is None and geq is None:
		print("ERROR: No GEQDSK file or object provided")
		return
	if geq is None:
		from geqdsk_reader import geqdsk
		geq = geqdsk(geq_file)
	if peq is None:
		from .peqdsk_reader import peqdsk
		peq = peqdsk(peq_file)
		
	psi_n = peq['psinorm']
	rho = []
	for psiN in psi_n:
		if psiN < 1.0e-3:
			psi_n = psi_n[psi_n != psiN]	
		else:
			fs = geq.flux_surface(psiN = psiN)
			rho.append((max(fs['R']) - min(fs['R']))/2)
	rhonorm = rho/max(rho)
	
	f = open(peq_file,'a')
	f.write(f"{len(rho)+1} psinorm rho rhonorm")
	f.write("\n 0.0000000   0.0000000   0.0000000")
	for i in range(len(rho)):
		f.write(f"\n {psi_n[i]:.7f}   {rho[i]:.7f}   {rhonorm[i]:.7f}")  
	f.close()
	return
