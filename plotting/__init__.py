from .scan_plotter import plot_scan
from .diagnostics_plotter import plot_diag
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar
from .set_plotter import plot_set
from .theta_plotter import plot_theta, plot_theta_scan, plot_theta_set
from .ql_plotter import plot_ql
from .slice_plotter import plot_slice

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Diag"] = plot_diag
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar
Plotters["Set"] = plot_set
Plotters["Theta"] = plot_theta_scan
Plotters["Theta_Single"] = plot_theta
Plotters["Theta_Set"] = plot_theta_set
Plotters["QL"] = plot_ql
Plotters["Slice"] = plot_slice

__all__ = ["Plotters"]
