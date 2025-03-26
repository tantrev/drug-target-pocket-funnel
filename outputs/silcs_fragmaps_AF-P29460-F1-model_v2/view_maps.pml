# SILCS Sampler
#
# Copyright (c) 2016 SilcsBio, LLC. All Rights Reserved.
#
# Any unauthorized reproduction, use, or modification of this work
# is prohibited.
#
# Contact:
#
# SilcsBio, LLC
# 1100 Wicomico Street, Suite 323
# Baltimore, MD, 21230
# 410-401-9794
# info@silcsbio.com
# http://www.silcsbio.com

# load protein
proteinpdb = "AF-P29460-F1-model_v2.pdb"
ligdir = "./ligand"
mapdir = "maps"
prefix = "AF-P29460-F1-model_v2"
filename = 'view_maps.pml'

from pymol import plugins
pmgapp = plugins.get_pmgapp()

# change working directory in case the script is loaded by GUI
python
def change_workingdir():
    global workdir
    import tkMessageBox, tkFileDialog
    tkMessageBox.showerror('File not found', 'Select the folder where silcs_fragmap files are downloaded.')
    workdir = tkFileDialog.askdirectory(initialdir='.')
    if os.path.exists(os.path.join(workdir, filename)):
        os.chdir(workdir)
        sys.path.append(os.path.join(workdir, 'plugins', 'pymol'))
        cmd.do("run " + filename)

if not os.path.exists(proteinpdb):
    pmgapp.execute(change_workingdir)
python end

# load protein
cmd.load(proteinpdb)
hide everything
set transparency, 0.5

# display setup
#bg_color white

# load PyMOL plugin
import sys
sys.path.append(os.path.join('plugins', 'pymol'))
from fragmap_tools import *
pmgapp.execute(FragMap.create_dialog)

python
# adjust default variables for fragmap
def update_default():
    FragMap.mapdir.setvalue(mapdir)
    FragMap.prefix.setvalue(prefix)
    FragMap.liganddir.setvalue(os.path.join(ligdir, '%s.gfe.score.ligand.pdb' % prefix))
    FragMap.mcsilcsdir.setvalue(os.path.join(ligdir, '%s.ligand.*.pdb' % prefix))
python end
pmgapp.execute(update_default)

pmgapp.execute(FragMap.show_fragmap)

python
# turn on protein (noh) surface
def show_surface():
    FragMap.visualization_group.surf_entry_flag.set(1)
    FragMap.visualization_group.showSurface('surf')

# turn off first mc-silcs ligand
def turn_off_mcligand():
    molfile = FragMap.mcligandsloaded[0]
    FragMap.mcligandsmap[molfile]['visible'].set(0)
    FragMap.load_mcligand(molfile)
python end
pmgapp.execute(show_surface)
