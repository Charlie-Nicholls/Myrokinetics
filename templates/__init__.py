from pathlib import Path
from .dimensions import dimensions_list
from .systems import systems

dim_lookup = {'_list': [], '_full_list': []}
for dim in dimensions_list:
	dim_lookup['_list'].append(dim.name)
	for dim_name in dim.name_keys:
		dim_name = dim_name.lower()
		dim_lookup[dim_name] = dim
		dim_lookup['_full_list'].append(dim_name)

template_dir = Path(__file__).parent
template_dir.resolve()
template_dir = str(template_dir)

gs2_template = "template.gs2"

inputs_template = "default_input.in"

__all__ = ['dim_lookup','template_dir','gs2_template','systems','inputs_template']
