from .scan_plotter import plot_scan
from .diagnostics_plotter import plot_diag, plot_single
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Diag"] = plot_diag
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar
Plotters["Eigen_Single"] = plot_single

__all__ = ["Plotters"]
