from .scan_plotter import plot_scan
from .eigen_plotter import plot_eigen, plot_single
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Eigen"] = plot_eigen
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar
Plotters["Eigen_Single"] = plot_single

__all__ = ["Plotters"]
