import os
from numpy import full, real, imag, array, loadtxt, transpose, savez
from Myrokinetics.ncdf2dict import ncdf2dict as readnc
from Myrokinetics.inputs import scan_inputs
from uuid import uuid4

def convert_box_to_myro(filename = None, run_data = None, directory = "./", input_file = None, QuickSave = False):
	
	inputs = scan_inputs(input_file = input_file, directory = directory)
	gyro_data = {}
	group_data = {}
	only = set({'omega','ky','kx'})
	if not QuickSave:
		only = only | set({'phi','bpar','apar','phi2','t','theta', 'gds2', 'jacob','ql_metric_by_mode', 'phi2_by_mode', 'phi2_by_kx', 'phi2_by_ky'})
	if inputs['non_linear']:
		only = only | set({'heat_flux_tot'})
	data_keys = ['growth_rate','mode_frequency','omega','phi','bpar','apar','epar','phi2','parity','ql_metric']
	group_keys = ['phi2_avg','t','theta', 'gds2', 'jacob','heat_flux_tot','phi2_by_kx', 'phi2_by_ky']
	
	if run_data == None:
		run_data = readnc(f"{directory}/{filename}.out.nc",only=only)
		
	group_key = run_data['attributes']['id']
	group_data[group_key] = {}
	kxs = list(run_data['kx'])
	kys = list(run_data['ky'])
	inputs.inputs['dimension_1']['values'] = kys
	inputs.inputs['dimension_2']['values'] = kxs
	inputs.load_dimensions()
	inputs.write_scan_input()
	dimensions = inputs.dimensions
	single_parameters = inputs.single_parameters
	
	
	equilibrium = {}
	psiN = single_parameters['psin'].values[0]
	equilibrium[psiN] = {}
	equilibrium[psiN]['shear'] = single_parameters['shear'].values[0]
	equilibrium[psiN]['beta_prime'] = single_parameters['beta_prime'].values[0]
	
	gyro_keys = {}
	for dim in dimensions.values():
		gyro_keys[dim.name] = {}
		for val in dim.values:
			gyro_keys[dim.name][val] = set()

	for yi, ky in enumerate(run_data['ky']):
		for xi, kx in enumerate(run_data['kx']):
			run_key = str(uuid4())
			gyro_keys['ky'][ky].add(run_key)
			gyro_keys['kx'][kx].add(run_key)
			
			gyro_data[run_key] = {}
			gyro_data[run_key]['kx'] = kx
			gyro_data[run_key]['ky'] = ky
			gyro_data[run_key]['group_key'] = group_key
			for key in data_keys:
				gyro_data[run_key][key] = None

			for key in only:
				try:
					key_data = run_data[key]				
					if key == 'omega':
						om = key_data[-1,yi,xi]
						if type(om) != complex:
							om = key_data[-2,yi,xi]
						gyro_data[run_key]['growth_rate'] = imag(om)
						gyro_data[run_key]['mode_frequency'] = real(om)
						gyro_data[run_key]['omega'] = key_data[:,yi,xi].tolist()
					elif key in ['phi','apar','bpar']:
						gyro_data[run_key][key] = key_data[yi,xi,:].tolist()
						if key == 'phi':
							try:
								symsum = sum(abs(key_data[yi,xi,:] + key_data[yi,xi,::-1]))/sum(abs(key_data[yi,xi,:]))
							except:
								symsum = 1
							if  symsum > 1.5:
								gyro_data[run_key]['parity'] = 1
							elif symsum < 0.5:
								gyro_data[run_key]['parity'] = -1
							else:
								gyro_data[run_key]['parity'] = 0
					
					elif key in ['ql_metric_by_mode']:
						gyro_data[run_key]['ql_metric'] = key_data[-1,yi,xi]
					elif key in ['phi2_by_mode']:
						gyro_data[run_key]['phi2'] = key_data[:,yi,xi]
					elif key in ['t','theta', 'gds2', 'jacob','heat_flux_tot','phi2_by_kx', 'phi2_by_ky']:
						group_data[group_key][key] = key_data.tolist()
					elif key in ['phi2']:
						group_data[group_key]['phi2_avg'] = key_data.tolist()
					
					'''
					elif key in ['epar']:
						epar_path = f"{sub_dir}/itteration_{itt}.epar"
				
						bpar = key_data['bpar'][yi,xi,:]
						epar_data = loadtxt(epar_path)
						epar = []
						for l in range(len(epar_data[:,3])):
							epar.append(complex(epar_data[l,3],epar_data[l,4]))
						epar = array(epar)
						gyro_data[run_key]['epar'] = epar
					''' 
				except:
						pass

		data = {'gyro': gyro_data,
		'ideal': None,
		'group': group_data,
		'equilibrium': equilibrium,
		'_gyro_keys': gyro_keys,
		'_ideal_keys': None,
		}
	
	file_lines = None
	savez(f"{directory}/{filename}", inputs = inputs.inputs, data = data, files = file_lines)
