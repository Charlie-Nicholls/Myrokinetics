from .geqdsk_reader import geqdsk as readg
from matplotlib.pyplot import *
from scipy.interpolate import RectBivariateSpline
from matplotlib.widgets import TextBox
from numpy import linspace, zeros, diff, insert
from pyrokinetics import Pyro

class geqdsk_interpolator(object):
    def __init__(self, file1 = None, file2 = None):
        if file1 is None or file2 is None:
            print("ERROR: file1 and/or file2 is None")
        self.geq1 = readg(file1)
        self.geq1 = readg(file2)
        self.pyro1 = Pyro(eq_file=file1,eq_type="GEQDSK")
        self.pyro2 = Pyro(eq_file=file2,eq_type="GEQDSK")
        self.data = {}
        self.create_interpolation()
    
    def create_interpolation(self, psires = 100, res = 300):
        variables = ['Rmaj', 'Z0', 'beta_prime', 'delta', 'kappa', 'q', 'rho', 's_delta', 's_kappa', 'shat', 'shift', 'r_minor', 'f_psi', 'B0', 'dpsidr', 'pressure', 'dpressure_drho', 'bunit_over_b0', 'a_minor']
        #array variables = ['R', 'Z', 'theta', 'b_poloidal']
        #a_minor is constant with respect to psiN
        #'psi_n' also required
        constants = ['local_geometry','btccw','ipccw','s_zeta', 'zeta']
        dat = {}
        psiNs = linspace(1/psires,1,psires)
        self.psiN = psiNs
        times = linspace(0,1,res)
        py1 = self.pyro1
        py2 = self.pyro2
        for var in variables:
            dat[var] = zeros((res,psires))
        for p, psiN in enumerate(psiNs):
            py1.load_local_geometry(psi_n=psiN)
            py2.load_local_geometry(psi_n=psiN)
            lg1 = py1.local_geometry
            lg2 = py2.local_geometry
            if p == 0:
                for cvar in constants:
                    dat[cvar] = lg1[cvar]
            for t, time in enumerate(times):
                for var in variables:
                    dat[var][t][p] = (lg2[var]-lg1[var])*time + lg1[var]
        for var in variables:
            self.data[var] = RectBivariateSpline(times,psiNs,dat[var])
    
    def var(self, var = None, t = None):
        if var is None:
            print("ERROR: var must be specified")
            return
        if var not in self.data.keys():
            print(f"ERROR: var must be in {self.data.keys()}")
            return
        val = zeros((len(self.psiN)))
        for p, psiN in enumerate(self.psiN):
            val[p] = self.data[var](t,psiN)
        return val
    
    def time(self, t = None):
        if t is None or t < 0 or t > 1:
            print("ERROR: t must be specified between 0 and 1")
            return
        ne = self.var(var = 'ne', t = t)
        te = self.var(var = 'te', t = t)
        return self.peqdsk_time(ne = ne, te = te, psiN = self.psiN, t = t)
    
    def write_peqdsk(self, t = None, filename = None):
        if t is None or t < 0 or t > 1:
            print("ERROR: t must be specified between 0 and 1")
            return
        a = self.time(t = t)
        a.write_peqdsk(filename = filename)
    
    def plot(self, var = 'beta_prime'):
            def draw_fig(val = None):
                val = eval(val)
                k = self.var(var = var, t = val)
                ax.cla()
                ax.plot(self.psiN,k0,'r', label = "eq1")
                ax.plot(self.psiN,k,'k--', label = "new eq")
                ax.plot(self.psiN,k1,'b', label = "eq2")
                ax.legend(loc = 0)
                ax.set_title(var)
                fig.canvas.draw_idle()
            fig, ax = subplots()
            k0 = self.var(var = var, t = 0)
            k1 = self.var(var = var, t = 1)
            taxes = fig.add_axes([0.5, 0.01, 0.05, 0.04])
            tbox = TextBox(taxes, 'time:', initial = '0.5')
            tbox.on_submit(draw_fig)
            draw_fig("0.5")
            show()
    
    class peqdsk_time(object):
        def __init__(self, ne = None, te = None, psiN = None, t = None):
            self.ne = ne
            self.te = te
            self.psiN = psiN
            self.time = self.t = t 
        
        def write_peqdsk(self, filename = None):
            if filename is None:
                filename = f"p{self.t}.peqdsk"
            f = open(filename,'w')
            f.write(f"{len(self.psiN)} psinorm ne(10^20/m^3) dne/dpsiN\n")
            for p, psi in enumerate(self.psiN):
                differ = diff(self.ne)/diff(self.psiN)
                differ = insert(differ, 0,0)
                f.write(f" {psi:.7f}   {self.ne[p]:.7f}   {differ[p]:.7f}\n")
            f.write(f"{len(self.psiN)} psinorm te(KeV) dte/dpsiN\n")
            for p, psi in enumerate(self.psiN):
                differ = diff(self.ne)/diff(self.psiN)
                differ = insert(differ, 0,0)
                f.write(f" {psi:.7f}   {self.te[p]:.7f}   {differ[p]:.7f}\n")
            f.close()
        
        def plot(self):
            fig, ax = subplots(1,3,figsize=(15,5))
            ax[0].cla()
            ax[1].cla()
            ax[2].cla()
            ax[0].plot(self.psiN,self.ne,'r')
            ax[1].plot(self.psiN,self.te,'r')
            ax[2].plot(self.psiN,self.ne*self.te,'r')
            ax[0].set_title("Density")
            ax[1].set_title("Temperature")
            ax[2].set_title("Pressure")
            show()


