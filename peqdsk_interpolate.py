from .peqdsk_reader import peqdsk as readp
from matplotlib.pyplot import *
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
from matplotlib.widgets import TextBox
from numpy import linspace, array, diff, insert, zeros


class peqdsk_interpolator(object):
    def __init__(self, file1 = None, file2 = None):
        if file1 is None or file2 is None:
            print("ERROR: file1 and/or file2 is None")
        self.peq1 = readp(file1)
        self.peq2 = readp(file2)
        self.keys = self._find_keys()
        self.psiN = None
        self.data = {}
        self.create_interpolation()
    
    def _find_keys(self):
        keys = []
        for key in self.peq1.keys():
            if key in self.peq2.keys():
                keys.append(key)
                
        if not set(['psinorm','ne','te']).issubset(keys):
            print("ERROR: psinorm, ne or te keys not found")
        return keys
    
    def create_interpolation(self, res = 300):
        for var in ['ne','te']:#self.keys:
            if var != 'psinorm':
                self._create_var(var, res)
    def _create_var(self, var = None, res = 300):
        if var is None:
            print("ERROR: var must be specified")
            return
        if var not in self.keys:
            print(f"ERROR: var must be in {self.keys}")
            return
        if var == 'psinorm':
            print("ERROR: var can not be psinorm")
            return
        
        psi1 = self.peq1['psinorm']
        psi2 = self.peq2['psinorm']
        
        var1 = self.peq1[var]
        var2 = self.peq2[var]
        
        if len(psi2) < len(psi1):
            var1_IUS = InterpolatedUnivariateSpline(psi1,var1)
            var1 = var1_IUS(psi2)
            psi = psi2
        else: 
            var2_IUS = InterpolatedUnivariateSpline(psi2,var2)
            var2 = var2_IUS(psi1)
            self.psiN = psi = psi1
        
        time = linspace(0,1,res)
        var_2d = []
        for p in range(len(psi)):
            var_p = []
            for t in time:
                var_p.append((var2[p] - var1[p])*t + var1[p])
            var_2d.append(var_p)
        self.data[var] = RectBivariateSpline(array(psi),array(time),array(var_2d))
    
    def var(self, var = None, t = None):
        if var is None:
            print("ERROR: var must be specified")
            return
        if var not in self.data.keys():
            print(f"ERROR: var must be in {self.data.keys()}")
            return
        val = zeros((len(self.psiN)))
        for p, psiN in enumerate(self.psiN):
            val[p] = self.data[var](psiN,t)
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
      
    
    def plot(self):
            def draw_fig(val = None):
                val = eval(val)
                ne = self.var(var = 'ne', t = val)
                te = self.var(var = 'te', t = val)
                ax[0].cla()
                ax[1].cla()
                ax[2].cla()
                ax[0].plot(self.psiN,ne0,'r', label = "eq1")
                ax[0].plot(self.psiN,ne,'k--', label = "new eq")
                ax[0].plot(self.psiN,ne1,'b', label = "eq2")
                ax[1].plot(self.psiN,te0,'r', label = "eq1")
                ax[1].plot(self.psiN,te,'k--', label = "new eq")
                ax[1].plot(self.psiN,te1,'b', label = "eq2")
                ax[2].plot(self.psiN,ne0*te0,'r', label = "eq1")
                ax[2].plot(self.psiN,ne*te,'k--', label = "new eq")
                ax[2].plot(self.psiN,ne1*te1,'b', label = "eq2")
                ax[0].legend(loc = 0)
                ax[1].legend(loc = 0)
                ax[2].legend(loc = 0)
                ax[0].set_title("Density")
                ax[1].set_title("Temperature")
                ax[2].set_title("Pressure")
                fig.canvas.draw_idle()
            
            fig, ax = subplots(1,3,figsize=(15,5))
            ne0 = self.var(var = 'ne', t = 0)
            te0 = self.var(var = 'te', t = 0)
            ne1 = self.var(var = 'ne', t = 1)
            te1 = self.var(var = 'te', t = 1)
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
            print(f"Created {filename}")
        
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
