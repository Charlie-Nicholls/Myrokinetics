#Obtained from David Dickinson
# A set of utilities for working with geqdsk files

import re
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.interpolate import UnivariateSpline as smoothspline
from scipy.interpolate import RectBivariateSpline as spline2D

class geqdsk(object):

    def __init__(self, filename = None, directory = "./"):
        import os.path
        self.filename = os.path.join(directory, filename)
        self.data = {}
        self.lines = {}
        self.levels = []

        if self.filename is not None:
            self._open_and_parse()

    def __getitem__(self, key):
        return self.data[key]

    def _open_and_parse(self, filename = None):
        if self.filename is None and filename is None:
            raise ValueError("ERROR: No filename set so cannot open")
        elif self.filename is None:
            self.filename = filename
        elif filename is None:
            pass
        elif self.filename != filename:
            print("Warning: Inconsistent filenames, object has {a} but passed {b}. Using {b}".format(
                a=self.filename, b = filename))
            self.filename = filename

        # Now check that the file exists and if so read all lines
        try:
            with open(self.filename) as the_file:
                self.lines = the_file.readlines()

        except IOError:
            print("Couldn't read from file {a}".format(a = self.filename))

        # Now add in some derived data
        self.data = self._read_from_lines(self.lines)
        self._calculate_derived_data()


    def _read_from_lines(self, the_lines = None):
        if the_lines is None:
            raise ValueError("ERROR: Must pass the_lines.")

        # Dictionary of raw data
        raw_data = {}

        # Now parse the lines - make a copy first
        lines = list(the_lines)

        def strip_empty(the_list):
            return list(filter(None, the_list))

        # Header with fortran format 6a8, 3i4
        # we want the last two integers, the first can be ignored
        # everything else is a descriptor
        header = lines.pop(0)
        # To extract last two split on " ", should have at least
        # three entries.
        tmp = strip_empty(header.rstrip().split(" "))
        if len(tmp) < 3:
            raise SyntaxError("ERROR: Malformed header")

        nw = int(tmp[-2]) # Number of R points
        nh = int(tmp[-1]) # Number of Z points

        if len(tmp) == 3:
            descriptor = tmp[0]
        else:
            descriptor = " ".join(tmp[:-3])

        raw_data['nr'] = nw
        raw_data['nw'] = nw
        raw_data['nz'] = nh
        raw_data['nh'] = nh
        raw_data['descriptor'] = descriptor

        # Next we expect lots of lines of format 5e16.9
        # The first four lines corresponds to individual scalar
        # quanitities. As there's no guarantee of spaces between numbers
        # we have to use regexp to extract these.
        single_float = r'([ \-\+]\d\.\d+[Ee][\+\-]\d\d)'
        up_to_five_floats = r'^\s*' + 5*(single_float+r'?') + r'\s*$'
        pattern = re.compile(up_to_five_floats)

        # Define the variable names for these items, note some repetition
        # and some to ignore.
        ignored_name = "xdum"
        names = [
            ["rdim", "zdim", "rcentr", "rleft", "zmid"],
            ["rmaxis", "zmaxis", "simag", "sibry", "bcentr"],
            ["current", "simag", ignored_name, "rmaxis", ignored_name],
            ["zmaxis", ignored_name, "sibry", ignored_name, ignored_name]
            ]

        for nm in names:
            tmp = lines.pop(0)
            result = re.findall(pattern, tmp)[0]

            if len(result) != 5:
                raise SyntaxError("ERROR: Malformed scalars")

            for i, the_name in enumerate(nm):
                if the_name == ignored_name:
                    continue
                the_value = float(result[i])
                if the_name in raw_data.keys():
                    the_other_value = raw_data[the_name]
                    if the_other_value != the_value:
                        print(the_value, the_other_value)

#                        raise ValueError("ERROR: Inconsistent duplicate definitions of variable {a}".format(a=the_name))
                        the_value = the_other_value

                raw_data[the_name] = the_value


        # Now we can get onto reading in profiles, firstly floating point profiles with
        # five values per line. We have a number of 1D in R profiles, followed by a 2D
        # profile and then some integers etc.
        # Could try to read values one at a time, just enough for a profile or all at once
        def get_1d(npts, lines_in):
            nvals = 0
            output = np.zeros(npts)
            while nvals < npts:
                tmp = lines_in.pop(0)
                result = strip_empty(re.findall(pattern, tmp)[0])

                for i, val in enumerate(result):
                    output[nvals] = float(val)
                    nvals += 1
                    if nvals == npts:
                        break

            # Now check if we have any unused tokens, if so join them up and add
            # them back into the list
            if i < len(result) - 1 :
                leftover = " ".join(result[i:])
                lines_in.insert(0, leftover)

            return output

        def get_2d(nx, ny, lines_in):
            return get_1d(nx*ny, lines_in).reshape([ny, nx]).T

        raw_data['fpol'] = get_1d(raw_data['nr'], lines)
        raw_data['pres'] = get_1d(raw_data['nr'], lines)
        raw_data['ffprim'] = get_1d(raw_data['nr'], lines)
        raw_data['pprime'] = get_1d(raw_data['nr'], lines)
        raw_data['psirz'] = get_2d(raw_data['nr'], raw_data['nz'], lines)
        raw_data['qpsi'] = get_1d(raw_data['nr'], lines)

        # Now we expect to read two integers
        tmp = strip_empty(lines.pop(0).rstrip().split(" "))
        nbbbs, limitr = map(int, tmp)

        # Next we have two interleaved 1D arrays
        tmp = get_1d(nbbbs*2, lines)
        raw_data['rbbbs'] = tmp[0::2]
        raw_data['zbbbs'] = tmp[1::2]

        # Finally we have two more interleaved 1D arrays
        tmp = get_1d(limitr*2, lines)
        raw_data['rlim'] = tmp[0::2]
        raw_data['zlim'] = tmp[1::2]

        return raw_data

    def _calculate_derived_data(self):
        # These are the axes of the rectangular grid
        self.data['R'] = np.linspace(self.data['rleft'], self.data['rleft'] + self.data['rdim'], self.data['nr'])
        self.data['Z'] = np.linspace(self.data['zmid'] - 0.5*self.data['zdim'], self.data['zmid'] + 0.5*self.data['zdim'], self.data['nz'])

        # Make uniform psi grids
        self.data['psiN'] = np.linspace(0.0, 1.0, self.data['nr'])
        self.data['psi'] = np.linspace(self.data['simag'], self.data['sibry'], self.data['nr'])

        # Construct splines for later use
        self._psiN = spline(self.data['psi'], self.data['psiN'], ext = 2)
        self._psi = spline(self.data['psiN'], self.data['psi'], ext = 2)
        self._fpol = spline(self.data['psiN'], self.data['fpol'], ext = 2)
        self._pres = spline(self.data['psiN'], self.data['pres'], ext = 2)
        self._ffprim = spline(self.data['psiN'], self.data['ffprim'], ext = 2)
        self._pprime = spline(self.data['psiN'], self.data['pprime'], ext = 2)
        self._qpsi = smoothspline(self.data['psiN'], self.data['qpsi'], ext = 2, s = 0.01)

        self._psiRZ = spline2D(self.data['R'], self.data['Z'], self.data['psirz'])

    def _get_psiN(self, psi = None, psiN = None):
        if psi is None and psiN is None:
            return self.data['psiN']
        elif psi is not None and psiN is not None:
            raise ValueError("ERROR: Must specifiy exactly one of psi and psiN")
        elif psiN is None:
            return self._psiN(psi)
        else:
            return psiN

    def _get_value_or_derivative(self, the_spline, psi = None, psiN = None, derivative = None, **kwargs):
        # Note here we assume the spline's independent variable is psiN
        the_psiN = self._get_psiN(psi, psiN)
        if derivative is None:
            return the_spline(the_psiN)
        else:
            return the_spline.derivative(n = derivative)(the_psiN)

    def _get_psi(self, psi = None, psiN = None):
        if psiN is None and psi is None:
            return self.data['psi']
        elif psi is not None and psiN is not None:
            raise ValueError("ERROR: Must specifiy exactly one of psi and psiN")
        elif psi is None:
            return self._psi(psiN)
        else:
            return psi

    def psiN(self, psi = None, psiN = None, derivative = None):
        # Here we don't use the generic _get_value_or_derivative as that assumes
        # the spline is agains psiN
        the_psi = self._get_psi(psi = psi, psiN = psiN)
        if derivative is None:
            return self._psiN(the_psi)
        else:
            return self._psiN.derivative(n = derivative)(the_psi)

    def psi(self, psi = None, psiN = None, derivative = None):
        return self._get_value_or_derivative(self._psi, psi, psiN, derivative)

    def psiRZ(self, R = None, Z = None, derivativeR = 0, derivativeZ = 0):
        # Get's psi(R, Z) or derivatives for all {R[i], Z[i]} i.e. expect
        # pairwise RZ not the axes/dimensions. If we made grid True instead
        # we'd have the opposite behaviour and we'd try to form grids based
        # on the two input dimensions
        return self._psiRZ(R, Z, derivativeR, derivativeZ, grid = False)

    def fpol(self, psi = None, psiN = None, derivative = None):
        return self._get_value_or_derivative(self._fpol, psi, psiN, derivative)

    def pres(self, psi = None, psiN = None, derivative = None):
        return self._get_value_or_derivative(self._pres, psi, psiN, derivative)

    def ffprim(self, psi = None, psiN = None, derivative = None):
        return self._get_value_or_derivative(self._ffprim, psi, psiN, derivative)

    def pprime(self, psi = None, psiN = None, derivative = None):
        return self._get_value_or_derivative(self._pprime, psi, psiN, derivative)

    def qpsi(self, psi = None, psiN = None, derivative = None):
        return self._get_value_or_derivative(self._qpsi, psi, psiN, derivative)

    # Primarily just for plotting to help us visualise the surface
    # of interest etc.
    def flux_surface(self, psi = None, psiN = None, isgnBpol = 1):
        the_psi = self._get_psi(psi, psiN)
        the_psiN = self._get_psiN(the_psi, None)

        if the_psiN < 1.0e-3:
            print("Warning : Requested psiN very small - likely no contours available so making psiN = 1.0e-3")
            the_psiN = 1.0e-3
            the_psi = self._get_psi(None, the_psiN)

        if abs(isgnBpol) != 1:
            raise ValueError("ERROR: isgnBpol should have magnitude of exactly 1.")

        # We need to find the contours of constant psi as a funcion of R,Z
        # Probably getting data out of a matplotlib contour is the best bet.

        import matplotlib.pyplot as plt

        fig = plt.figure()
        cnt = plt.contour(self['R'], self['Z'], self['psirz'].T, levels=[the_psi])
        path = cnt.collections[0].get_paths()
        plt.close(fig)

        # Now we need to work out which path is correct, probably the one with the most points?
        # More strictly it's the one where all points are <= bdry/limitr
        # We can pick any point on the path an draw the line to the origin. If the
        # number of times this intersects the bdry is even then that point is outside
        # the boundary and hence that path can't be a valid one. But then we need to work out
        # how to find out if these things intersect (and if so how many times)
        maxR = self['R'].max() ; minR = self['R'].min()
        maxZ = self['Z'].max() ; minZ = self['Z'].min()

        results = []
        for i, case in enumerate(path):
            rr = case.vertices[:,0]
            zz = case.vertices[:,1]
            if (rr > maxR).any() or (rr < minR).any() or (zz > maxZ).any() or (zz < minZ).any():
                #Can't be this case as part of contour outside boundary
                pass
            else:
                results.append(case)

        if len(results) == 1:
            result = results[0]
        else:
            # Otherwise find the path with the least average distance from rmaxis, zmaxis
            distances = [ np.mean( np.sqrt((x.vertices[:,0] - self['rmaxis'])**2 +
                                           (x.vertices[:,1] - self['zmaxis'])**2)) for x in results]
            result = results[np.array(distances).argmin()]

        # Now we have the geometry let's save it.
        R, Z = result.vertices[:,0], result.vertices[:,1]

        # Now package these and other values on this surface
        data = {}
        data['psi'] = the_psi
        data['psiN'] = self.psiN( psi = the_psi)
        data['dpsiN_dpsi'] = self.psiN( psi = the_psi, derivative = 1)

        # Note we "roll" the R, Z arrays to ensure it starts with minimum minor radius
        # at the first position.
        nshift = -R.argmin()
        data['R'] = np.roll(R, nshift)
        data['Z'] = np.roll(Z, nshift)

        # Here we are finding out if there are any duplicate points and then removing them
        uniqueR = set(np.unique(data['R'], return_index = True)[1])
        uniqueZ = set(np.unique(data['Z'], return_index = True)[1])
        dupR = set(range(len(data['R']))) - uniqueR
        dupZ = set(range(len(data['Z']))) - uniqueZ
        dupBoth = np.array(list(dupR.intersection(dupZ)))

        #print("Found {n} duplicate R,Z points from {nn}.".format(n = len(dupBoth), nn = len(data['R'])))

        if len(dupBoth) > 0:
            data['R'] = np.delete(data['R'], dupBoth)
            data['Z'] = np.delete(data['Z'], dupBoth)

        # Derive a geometric theta -- should we subtract off rmaxis, or the average etc?
        # Typically in GS2 we use {rmaxis, zmaxis} as the reference point so do that here.
        data['theta'] = np.arctan2(data['Z'] - self['zmaxis'], data['R'] - self['rmaxis'])

        # Here we guard against the case where R, Z points are unique to floating precision
        # but close enough that the calculated theta has duplicate points. We identify duplicate
        # theta points and delete those indices from theta, R and Z. Instead of using unique as above
        # we use diff and where to find those locations where the theta grid is not increasing.
        badTheta = np.where(np.diff(-data['theta']) <= 0.0)

        if len(badTheta) > 0:
            #print("Found {n} duplicate theta points from {nn}.".format(n = len(badTheta), nn = len(data['theta'])))
            data['theta'] = np.delete(data['theta'], badTheta)
            data['R'] = np.delete(data['R'], badTheta)
            data['Z'] = np.delete(data['Z'], badTheta)

        # Record values of raw data on surface
        data['fpol'] = self.fpol( psi = the_psi)
        data['pres'] = self.pres( psi = the_psi)
        data['ffprim'] = self.ffprim( psi = the_psi)
        data['pprime'] = self.pprime( psi = the_psi)
        data['qpsi'] = self.qpsi( psi = the_psi)

        # Record derivatives of raw values on surface
        data['dfpol_dpsiN'] = self.fpol( psi = the_psi, derivative = 1)
        data['dpres_dpsiN'] = self.pres( psi = the_psi, derivative = 1)
        data['dffprim_dpsiN'] = self.ffprim( psi = the_psi, derivative = 1)
        data['dpprime_dpsiN'] = self.pprime( psi = the_psi, derivative = 1)
        data['dqpsi_dpsiN'] = self.qpsi( psi = the_psi, derivative = 1)

        # Find dpsi_dR, dpsi_dZ needed to construct Bpol
        data['dpsi_dR'] = self.psiRZ(data['R'], data['Z'], derivativeR = 1)
        data['dpsi_dZ'] = self.psiRZ(data['R'], data['Z'], derivativeZ = 1)

        # Construct various magnetic fields
        data['BR'] = data['dpsi_dZ']/data['R']
        data['BZ'] = -data['dpsi_dR']/data['R']
        data['Bpol'] = isgnBpol*np.sqrt(data['BR']**2 + data['BZ']**2)
        data['Btor'] = data['fpol']/data['R']
        data['Bmag'] = np.sqrt(data['Bpol']**2 + data['Btor']**2)

        # Could add current calculations here, but don't really
        # need them currently.

        # Add some other related parameters
        mu0 = 4.0e-7 * np.pi
        # Here beta is the total beta using the full on surface pressure and the magnetic field on axis
        # This is not the GS2 beta which is roughly a single species beta
        data['beta'] = 2 * mu0 * data['pres'] / self['bcentr']**2
        # Here we assume B is constant, as it is in the above definition of beta
        # Note GS2's beta_prime_input is dbeta_dpsiN / pres, i.e. beta * (1/Lp) = beta * (1/p) * dp/dpsiN
        data['dbeta_dpsiN'] = data['dpres_dpsiN'] * data['beta']
        # Here we estimate the shear assuming a radial coordinate of psiN -- this seems to disagree a bit (~10%)
        # with GS2 which calculates this in a somewhat different way.
        data['shear'] = data['dqpsi_dpsiN'] * data['psiN'] / data['qpsi']

        # Copy some other global data for convenience
        data['rmaxis'] = self['rmaxis']
        data['zmaxis'] = self['zmaxis']

        return data

    def _plot_psi(self):
        import matplotlib.pyplot as plt
        
        plt.figure()
        plt.contourf(self['R'], self['Z'], self['psirz'].T, 256, vmin = self['simag'], vmax = self['sibry'])
        plt.colorbar()
        plt.contour(self['R'], self['Z'], self['psirz'].T, levels=[self['sibry']], colors=['g'])
        if self.levels:
        	plt.contour(self['R'], self['Z'], self['psirz'].T, levels=self.levels, colors=['k'])
        plt.plot(self['rbbbs'],self['zbbbs'],'k--')
        plt.axis('equal')
        plt.show()
