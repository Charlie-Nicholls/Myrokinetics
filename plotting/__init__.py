from .scan_plotter import plot_scan
from .diagnostics_plotter import plot_diag
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar
from .theta_plotter import plot_theta
from .ql_plotter import plot_ql
from .slice_plotter import plot_slice
from .kxky_plotter import plot_kxky
from .box_diagnostics_plotter import plot_box_diag

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Diag"] = plot_diag
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar
Plotters["Theta"] = plot_theta
Plotters["QL"] = plot_ql
Plotters["Slice"] = plot_slice
Plotters["kxky"] = plot_kxky
Plotters["Box_Diag"] = plot_box_diag

__all__ = ["Plotters"]
