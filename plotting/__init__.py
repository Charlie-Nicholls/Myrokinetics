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
from .nonlinear_plotter import plot_phi2_by_mode, plot_hflux, plot_nl_phi2, plot_zonality

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
Plotters["NL_Phi2"] = plot_phi2_by_mode
Plotters["NL_Hflux"] = plot_hflux
Plotters["NL_Phi2_by_k"] = plot_nl_phi2
Plotters["NL_Zonality"] = plot_zonality

__all__ = ["Plotters"]
