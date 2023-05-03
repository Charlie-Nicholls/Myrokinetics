from pathlib import Path

template_dir = Path(__file__).parent / "templates"
template_dir.resolve()
template_dir = str(template_dir)

gs2_template = "template.gs2"
gyro_job = "gyro.job"
ideal_job = "ideal.job"
