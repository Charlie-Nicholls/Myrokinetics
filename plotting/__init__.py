from .scan_plotter import plot_scan
from .diagnostics_plotter import plot_diag, plot_diag_single, plot_diag_set
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Diag"] = plot_diag
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar
Plotters["Diag_Single"] = plot_diag_single
Plotters["Diag_Set"] = plot_diag_set

__all__ = ["Plotters"]
