#Obtained from David Dickinson
#Originally based on datafile.py from the BOUT++ tools
#-->Writing ability removed
#-->Add conversion to dictionary (with complex forming)

try:
    from numpy import complex, product
except ImportError:
    print("ERROR: NumPy module not available")
    raise

library = None # Record which library to use

try:
    from netCDF4 import Dataset
    library = "netCDF4"
except ImportError:
    try:
        from Scientific.IO.NetCDF import NetCDFFile as Dataset
        from Scientific.N import Int, Float, Float32
        library = "Scientific"
    except ImportError:
        try:
            from scipy.io.netcdf import netcdf_file as Dataset
            library = "scipy"
        except:
            print("DataFile: No supported NetCDF modules available")
            raise

class DataFile:
    handle = None

    def open(self, filename, format='NETCDF3_CLASSIC'):
        try:
            if library == "scipy":
                self.handle = Dataset(filename, "r", mmap=False)
            else:
                self.handle = Dataset(filename, "r")
        except IOError:
            print("Couldn't read file : '"+filename+"'")
            self.handle=None

        # Record if writing
        self.writeable = False

    def close(self):
        if self.handle != None:
            self.handle.close()
        self.handle = None

    def __init__(self, filename=None, format='NETCDF3_CLASSIC', **kwargs):
        if filename != None:
            self.open(filename, format=format)

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def read(self, name, ranges=None):
        """Read a variable from the file."""
        if self.handle == None: return None

        try:
            var = self.handle.variables[name]
        except KeyError:
            # Not found. Try to find using case-insensitive search
            var = None
            if self.handle is not none:
                for n in self.handle.variables.keys():
                    if n.lower() == name.lower():
                        print("WARNING: Reading '"+n+"' instead of '"+name+"'")
                        var = self.handle.variables[n]
            if var == None:
                return None
        ndims = len(var.dimensions)
        if ndims == 0:
            data = var.getValue()
            return data #[0]
        else:
            if ranges != None:
                if len(ranges) != 2*ndims:
                    print("Incorrect number of elements in ranges argument")
                    return None
                
                if library == "Scientific":
                    # Passing ranges to var[] doesn't seem to work
                    data = var[:]
                    if ndims == 1:
                        data = data[ranges[0]:ranges[1]]
                    elif ndims == 2:
                        data = data[ranges[0]:ranges[1], 
                                    ranges[2]:ranges[3]]
                    elif ndims == 3:
                        data = data[ranges[0]:ranges[1], 
                                    ranges[2]:ranges[3],
                                    ranges[4]:ranges[5]]
                    elif ndims == 4:
                        data = data[(ranges[0]):(ranges[1]),
                                    (ranges[2]):(ranges[3]),
                                    (ranges[4]):(ranges[5]),
                                    (ranges[6]):(ranges[7])]
                else:
                    if ndims == 1:
                        data = var[ranges[0]:ranges[1]]
                    elif ndims == 2:
                        data = var[ranges[0]:ranges[1], 
                                   ranges[2]:ranges[3]]
                    elif ndims == 3:
                        data = var[ranges[0]:ranges[1], 
                                   ranges[2]:ranges[3],
                                   ranges[4]:ranges[5]]
                    elif ndims == 4:
                        #print "Ranges = ", ranges
                        data = var[(ranges[0]):(ranges[1]),
                                   (ranges[2]):(ranges[3]),
                                   (ranges[4]):(ranges[5]),
                                   (ranges[6]):(ranges[7])]
                return data
            else:
                return var[:]

    def __getitem__(self, name):
        var = self.read(name)
        if var is None:
            raise KeyError("No variable found: "+name)
        return var

    def list(self):
        """List all variables in the file."""
        if self.handle == None: return []
        return self.handle.variables.keys()

    def keys(self):
        """List all variables in the file."""
        return self.list()

    def dimensions(self, varname):
        """Array of dimension names"""
        if self.handle == None: return None
        try:
            var = self.handle.variables[varname]
        except KeyError:
            return None
        return var.dimensions
        
    def ndims(self, varname):
        """Number of dimensions for a variable."""
        if self.handle == None: return None
        try:
            var = self.handle.variables[varname]
        except KeyError:
            return None
        return len(var.dimensions)
    
    def size(self, varname):
        """List of dimension sizes for a variable."""
        if self.handle == None: return []
        try:
            var = self.handle.variables[varname]
        except KeyError:
            return []
        
        def dimlen(d):
            dim = self.handle.dimensions[d]
            if dim != None:
                t = type(dim).__name__
                if t == 'int':
                    return dim
                return len(dim)
            return 0
        return map(lambda d: dimlen(d), var.dimensions)

    def __setitem__(self, key, value):
        self.write(key, value)

    def get_global_attributes(self):
        attributes = {}
        try:
            for name in self.handle.ncattrs():
                attributes[name] = self.handle.getncattr(name)
        except:
            print("ERROR: Couldn't extract attributes")
        return attributes

    def toDict(self,comDim='ri', only=None, **kwargs):
        """Convert object to a simple dictionary, optionally squashing complex dim."""

        ret={}
        if self.handle is None: return ret
        ci=complex(0.0,1.0)
        for varName in self.handle.variables.keys():
            if only is not None and varName not in only:
                continue
            dims=self.handle.variables[varName].dimensions
            #Not complex
            if not comDim in dims:
                ret[varName]=self[varName]
            #Complex
            else:
                #Original data
                data=self[varName]
                oldShape=data.shape

                #Which dim is complex
                dimInd=dims.index(comDim)
                if oldShape[dimInd] != 2:
                    print("ERROR: Expecting complex dimension to have length 2.")
                    return None

                #Make sure complex is last index
                indOrd=list(range(len(dims)))
                indOrd.append(indOrd.pop(dimInd))
                newShape=list(oldShape)
                newShape.pop(dimInd)
                dataT=data.transpose(indOrd)

                fdata=dataT.reshape([product(newShape),2])
                cdata=fdata[:,0]+ci*fdata[:,1]
                cdata=cdata.reshape(newShape)
                ret[varName]=cdata.copy()

        # Now try to add any attributes which may exist
        possible_attribute_names = ['attributes', '_attributes',
                                    'global_attributes', '_global_attributes']
        attribute_name_found = False
        for name in possible_attribute_names:
            if name not in ret.keys():
                attribute_name_found = True
                attribute_name = name
                break
        if not attribute_name_found:
            print("Warning: Couldn't find a unique name for attributes so not including.")
        else:
            ret[attribute_name] = self.get_global_attributes()
        return ret


def ncdf2dict(filename=None,*args,**kwargs):
    """Returns a dictionary representation of a netcdf file."""
    if filename is None:
        print("ERROR: Must pass filename (todo: extend to match a pattern)")
        return
    
    a=DataFile(filename=filename,*args,**kwargs)
    return a.toDict(*args,**kwargs)
