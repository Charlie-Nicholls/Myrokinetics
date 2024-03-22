from pathlib import Path
from .dimensions_gs2 import dimensions_list as dimensions_list_gs2
from .dimensions_cgyro import dimensions_list as dimensions_list_cgyro
from .systems import systems

dim_lookup_gs2 = {'_list': [], '_full_list': []}
for dim in dimensions_list_gs2:
	dim_lookup_gs2['_list'].append(dim.name_keys[0].lower())
	for dim_name in dim.name_keys:
		dim_name = dim_name.lower()
		dim_lookup_gs2[dim_name] = dim
		dim_lookup_gs2['_full_list'].append(dim_name)

dim_lookup_cgyro = {'_list': [], '_full_list': []}
for dim in dimensions_list_cgyro:
	dim_lookup_cgyro['_list'].append(dim.name_keys[0].lower())
	for dim_name in dim.name_keys:
		dim_name = dim_name.lower()
		dim_lookup_cgyro[dim_name] = dim
		dim_lookup_cgyro['_full_list'].append(dim_name)	

template_dir = Path(__file__).parent
template_dir.resolve()
template_dir = str(template_dir)

gs2_template = "template.gs2"
cgyro_template = "template.cgyro"

inputs_template = "default_input.in"

__all__ = ['dim_lookup_gs2','dim_lookup_cgyro','template_dir','gs2_template','cgyro_template','systems','inputs_template']