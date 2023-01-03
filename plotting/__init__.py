from .scan_plotter import plot_scan
from .eigen_plotter import plot_eigen
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Eigen"] = plot_eigen
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar

__all__ = ["Plotters"]
