from .scan_plotter import plot_scan
from .diagnostics_plotter import plot_diag
from .ideal_plotter import plot_ideal
from .epar_plotter import plot_epar
from .theta_plotter import plot_theta
from .two_d_plotter import plot_2d
from .slice_plotter import plot_slice
from .kxky_plotter import plot_kxky
from .box_diagnostics_plotter import plot_box_diag
from .slider_ax import slider_axes

Plotters = {}
Plotters["Scan"] = plot_scan
Plotters["Diag"] = plot_diag
Plotters["Ideal"] = plot_ideal
Plotters["Epar"] = plot_epar
Plotters["Theta"] = plot_theta
Plotters["2D"] = plot_2d
Plotters["Slice"] = plot_slice
Plotters["kxky"] = plot_kxky
Plotters["Box_Diag"] = plot_box_diag
Plotters["Sliders"] = slider_axes

__all__ = ["Plotters"]
