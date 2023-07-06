from .scan_plotter import plot_scan
from .diagnostics_plotter import plot_diag, plot_diag_single, plot_diag_set
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar
from .set_plotter import plot_set
from .theta_plotter import plot_theta, plot_theta_scan, plot_theta_set
from .ql_plotter import plot_ql

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Diag"] = plot_diag
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar
Plotters["Diag_Single"] = plot_diag_single
Plotters["Diag_Set"] = plot_diag_set
Plotters["Set"] = plot_set
Plotters["Theta"] = plot_theta_scan
Plotters["Theta_Single"] = plot_theta
Plotters["Theta_Set"] = plot_theta_set
Plotters["QL"] = plot_ql

__all__ = ["Plotters"]
