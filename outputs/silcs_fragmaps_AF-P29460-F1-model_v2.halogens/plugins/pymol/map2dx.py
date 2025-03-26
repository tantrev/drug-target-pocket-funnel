from __future__ import print_function

_map_header_tmpl = """GRID_PARAMETER_FILE {paramfile}
GRID_DATA_FILE {datafile}
MACROMOLECULE {molecule}
SPACING {spacing:4.3f}
NELEMENTS {npts[0]} {npts[1]} {npts[2]}
CENTER {center[0]:5.3f} {center[1]:5.3f} {center[2]:5.3f}
"""
_dx_header_tmpl = """object 1 class gridpositions counts {nx} {ny} {nz}
origin {ori[0]:12.5e} {ori[1]:12.5e} {ori[2]:12.5e}
delta {spacing:12.5e} 0 0
delta 0 {spacing:12.5e} 0
delta 0 0 {spacing:12.5e}
object 2 class gridconnections counts {nx} {ny} {nz}
object 3 class array type double rank 0 items {nvals} data follows
"""
_dx_footer_tmpl = """attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
"""

class ADGridMap:
    """Class for handling AutoDock Map files"""

    def __init__(self, fp = None, name = 'map'):

        self.name = ''
        self.npts = [0,0,0]
        self.n = [0,0,0]
        self.center = [0,0,0]
        self.origin = [0,0,0]
        self.nelem = 0
        self.spacing = 0.
        self.values = []
        self.datafile = ''
        self.molecule = ''
        self.paramfile = ''
        self.precision = 0.0001
        if fp is not None:
            self.read(fp,name)

    def read(self,fp,name='map'):
        self.name = name
        for i in range(6):
            line = fp.readline()
            if i == 0:
                self.paramfile = '' #line.split()[1]
            elif i == 1:
                self.datafile = '' #line.split()[1]
            elif i == 2:
                self.molecule = '' #line.split()[1]
            elif i == 3:
                self.spacing = float(line.split()[1])
            elif i == 4:
                self.npts = [int(float(x)) for x in line.split()[1:]]
            elif i == 5:
                self.center = [float(x) for x in line.split()[1:]]
        for i in range(3):
            self.n[i] = self.npts[i]+1
        self.nelem=self.n[0]*self.n[1]*self.n[2]
        i = 0
        while i < self.nelem:
            val = float(fp.readline())
            self.values.append(val)
            i+=1
        for i in range(3):
            self.origin[i] = self.center[i] - int(self.npts[i]/2)*self.spacing

    def meta(self):
        _dx_header_tmpl.format(**self.__dict__)
        return s

    def write(self,fp):
        fp.write(self.meta())
        for x in self.values:
            if abs(x) < self.precision:
                fp.write("0.\n")
            else:
                fp.write("%.3f\n" % x)

    def writeDX(self,fname):
        fp = open(fname,'w')
        nx = self.n[0]
        ny = self.n[1]
        nz = self.n[2]
        ori = self.origin
        spacing = self.spacing
        vals = self.values

        fp.write(_dx_header_tmpl.format(nx=nx, ny=ny, nz=nz, ori=ori, spacing=spacing, nvals=len(vals)))
        col=0;
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    fp.write(" %12.5E" % vals[k*nx*ny + j*nx + i])
                    col+=1;
                    if col==3:
                        fp.write("\n")
                        col=0
        if col != 0:
            fp.write("\n")
        fp.write(_dx_footer_tmpl)
        fp.close()

if __name__ == "__main__":
    import argparse
    import os
    import glob

    parser = argparse.ArgumentParser(description='Convert AutoDock Map format file to DX format file')
    parser.add_argument('mapfile', nargs='+',
                        help='AutoDock Map file name. The output will be saved at the same location where .map file is and the file extension of .dx will be added at the end.')
    args = parser.parse_args()

    for mapfile in args.mapfile:
        dxfile = '.'.join(mapfile.split('.')[:-1]) + '.dx'
        if os.path.exists(mapfile):
            ADGridMap(open(mapfile)).writeDX(dxfile)

