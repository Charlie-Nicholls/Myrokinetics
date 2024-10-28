from pathlib import Path

template_dir = Path(__file__).parent
template_dir.resolve()
template_dir = str(template_dir)

gs2_template = "template.gs2"
cgyro_template = "template.cgyro"
tglf_template = "template.tglf"

inputs_template = "default_input.in"

__all__ = ['template_dir','gs2_template','cgyro_template','tglf_template','inputs_template']
