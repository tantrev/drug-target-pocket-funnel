# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
PyMol FragMap Plugin

Copyright (c) 2016-2021 SilcsBio, LLC. All Rights Reserved.

Any unauthorized reproduction, use, or modification of this work
is prohibited.

Contact:

SilcsBio, LLC
1100 Wicomico Street, Suite 323
Baltimore, MD, 21230
410-401-9794
info@silcsbio.com
http://www.silcsbio.com

This work contains code derived from the APBS Plugin, Copyright 2009,
Michael G. Lerner (http://pymolapbsplugin.svn.sourceforge.net). Below
is the original copyright notice and permission notice.

# ----------------------------------------------------------------------
# APBS TOOLS is Copyright (C) 2009 by Michael G. Lerner
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Michael G. Lerner not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# MICHAEL G. LERNER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL MICHAEL G. LERNER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

This work contains code derived from the Pmw contrib directory, by
Rob W. W. Hooft. Below is the original permission notice.

# Filename dialogs using Pmw
#
# (C) Rob W.W. Hooft, Nonius BV, 1998
#
# Modifications:
#
# J. Willem M. Nissink, Cambridge Crystallographic Data Centre, 8/2002
#    Added optional information pane at top of dialog; if option
#    'info' is specified, the text given will be shown (in blue).
#    Modified example to show both file and directory-type dialog
#
# No Guarantees. Distribute Freely.
# Please send bug-fixes/patches/features to <r.hooft@euromail.com>
"""

from __future__ import division
from __future__ import generators
from __future__ import print_function

DEBUG = 0

import tempfile
import os,math,re
import string
import copy
import glob

try:
    import Tkinter
    from Tkinter import *
    import ttk
except ImportError:
    import tkinter as Tkinter
    from tkinter import *
try:
    import tkFileDialog
except:
    import tkinter.filedialog as tkFileDialog

import Pmw
import distutils.spawn # used for find_executable
import traceback
import pymol

#
# Global config variables
#

def __init__(self):
    """
    Init PyMOL, by adding FragMap Tools to the GUI under Plugins
    """
    self.menuBar.addmenuitem(
        'Plugin', 'command',
        'Launch FragMap Tools',
        label='FragMap Tools...',
        command = lambda s=self: FragMap.create_dialog(s)
    )

class pymolutil:
    """
    A quick collection of utility functions.
    """

    @staticmethod
    def getMacroMolecules():
        """returns all molecules that PyMOL knows about excluding the ligands loaded by the plugin"""
        pluginloaded = []
        for molfile in FragMap.ligandsmap:
            ligand = FragMap.ligandsmap[molfile]
            pluginloaded.append(ligand['object'])
        for molfile in FragMap.mcligandsmap:
            ligand = FragMap.mcligandsmap[molfile]
            pluginloaded.append(ligand['object'])
        return [i for i in pymolutil.getMolecules() if i not in pluginloaded]

    @staticmethod
    def getMolecules():
        """returns all molecules that PyMOL knows about"""
        return [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:molecule']

    @staticmethod
    def getMaps():
        """returns all maps that PyMOL knows about"""
        return [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:map']

pymolutil = pymolutil()

# Global variable

FRAGMAPS_AVAILABLE = {
    'apolar': {
        'label': 'Generic Apolar Map',
        'color': 'green',
        'visible': 1,
        'level': -1.2,
    },
    'hbdon': {
        'label': 'Generic Donor Map',
        'color': 'blue',
        'visible': 1,
        'level': -1.2,
    },
    'hbacc': {
        'label': 'Generic Acceptor Map',
        'color': 'red',
        'visible': 1,
        'level': -1.2,
    },
    'excl': {
        'label': 'Exclusion Map',
        'color': 'sand',
        'visible': 0,
        'level': 0.6,
        'use_surf': True,
    },
    'benc': {
        'label': 'Benzene Carbon Map',
        'color': 'purple',
        'visible': 0,
        'level': -1.2,
    },
    'prpc': {
        'label': 'Propane Carbon Map',
        'color': 'lime',
        'visible': 0,
        'level': -1.2,
    },
    'meoo': {
        'label': 'Methanol Oxygen Map',
        'color': 'tan',
        'visible': 1,
        'level': -1.2,
    },
    'forn': {
        'label': 'Formamide Nitrogen Map',
        'color': 'blue',
        'visible': 0,
        'level': -1.2,
    },
    'foro': {
        'label': 'Formamide Oxygen Map',
        'color': 'red',
        'visible': 0,
        'level': -1.2,
    },
    'mamn': {
        'label': 'Methylammonium Nitrogen Map',
        'color': 'cyan',
        'visible': 1,
        'level': -1.2,
    },
    'acec': {
        'label': 'Acetate Oxygen Map',
        'color': 'orange',
        'visible': 1,
        'level': -1.2,
    },
    'aalo': {
        'label': 'Acetaldehyde Oxygen Map',
        'color': 'red',
        'visible': 0,
        'level': -1.2,
    },
    'iminh': {
        'label': 'Imidazole Donor Nitrogen',
        'color': 'blue',
        'visible': 0,
        'level': -1.2,
    },
    'imin': {
        'label': 'Imidazole Acceptor Nitrogen',
        'color': 'red',
        'visible': 0,
        'level': -1.2,
    },
    'tipo': {
        'label': 'Water Oxygen Map',
        'color': 'black',
        'visible': 0,
        'level': -0.5,
    },
    'clbx': {
        'label': 'Chlorobenzene Map',
        'color': 'pink',
        'visible': 0,
        'level': -1.2,
    },
    'flbx': {
        'label': 'Fluorobenzene Map',
        'color': 'blue',
        'visible': 0,
        'level': -1.2,
    },
    'brbx': {
        'label': 'Bromobenzene Map',
        'color': 'brown',
        'visible': 0,
        'level': -1.2,
    },
    'clex': {
        'label': 'Chloroethane Map',
        'color': 'green',
        'visible': 0,
        'level': -1.2,
    },
    'fetx': {
        'label': 'Fluoroethane Map',
        'color': 'green',
        'visible': 0,
        'level': -1.2,
    },
    'tfec': {
        'label': 'Trifluoroethane Map',
        'color': 'green',
        'visible': 0,
        'level': -1.2,
    },
    'dmeo': {
        'label': 'Dimethyl-ethoxide Map',
        'color': 'red',
        'visible': 0,
        'level': -1.2,
    },
}

FRAGMAPS_ENABLED = (
    'apolar', 'hbdon', 'hbacc', 'excl',
    'benc', 'prpc', 'meoo', 'forn', 'foro', 'mamn', 'acec', 'aalo', 'iminh', 'imin', 'tipo',
    'flbx', 'clbx', 'brbx', 'clex', 'tfec','dmeo'
)


# https://stackoverflow.com/a/37861801/532799
class ScrolledWindow(Frame):
    """
    1. Master widget gets scrollbars and a canvas. Scrollbars are connected
    to canvas scrollregion.

    2. self.scrollwindow is created and inserted into canvas

    Usage Guideline:
    Assign any widgets as children of <ScrolledWindow instance>.scrollwindow
    to get them inserted into canvas

    __init__(self, parent, canv_w = 400, canv_h = 400, *args, **kwargs)
    docstring:
    Parent = master of scrolled window
    canv_w - width of canvas
    canv_h - height of canvas

    """

    def __init__(self, parent, canv_w = 400, canv_h = 400, *args, **kwargs):
        """Parent = master of scrolled window
        canv_w - width of canvas
        canv_h - height of canvas
        """
        super().__init__(parent, *args, **kwargs)

        self.parent = parent

        # creating a scrollbars
        self.xscrlbr = Scrollbar(self.parent, orient='horizontal')
        self.xscrlbr.grid(column=0, row=1, sticky='ew', columnspan=2)
        self.yscrlbr = Scrollbar(self.parent)
        self.yscrlbr.grid(column=1, row=0, sticky='ns')
        # creating a canvas
        self.canv = Canvas(self.parent)
        self.canv.config(relief='flat',
                         width=10,
                         heigh=10, bd=2)
        # placing a canvas into frame
        self.canv.grid(column=0, row=0, sticky='nsew')
        # accociating scrollbar comands to canvas scroling
        self.xscrlbr.config(command=self.canv.xview)
        self.yscrlbr.config(command=self.canv.yview)

        # creating a frame to inserto to canvas
        self.scrollwindow = Frame(self.parent)

        self.canv.create_window(0, 0, window=self.scrollwindow, anchor='nw')

        self.canv.config(xscrollcommand=self.xscrlbr.set,
                         yscrollcommand=self.yscrlbr.set,
                         scrollregion=(0, 0, 100, 100))

        self.yscrlbr.lift(self.scrollwindow)
        self.xscrlbr.lift(self.scrollwindow)
        self.scrollwindow.bind('<Configure>', self._configure_window)
        self.scrollwindow.bind('<Enter>', self._bound_to_mousewheel)
        self.scrollwindow.bind('<Leave>', self._unbound_to_mousewheel)

        return

    def _bound_to_mousewheel(self, event):
        self.canv.bind_all("<MouseWheel>", self._on_mousewheel)

    def _unbound_to_mousewheel(self, event):
        self.canv.unbind_all("<MouseWheel>")

    def _on_mousewheel(self, event):
        self.canv.yview_scroll(int(-1*(event.delta/120)), "units")

    def _configure_window(self, event):
        # update the scrollbars to match the size of the inner frame
        size = (self.scrollwindow.winfo_reqwidth(), self.scrollwindow.winfo_reqheight())
        self.canv.config(scrollregion='0 0 %s %s' % size)
        if self.scrollwindow.winfo_reqwidth() != self.canv.winfo_width():
            # update the canvas's width to fit the inner frame
            self.canv.config(width=self.scrollwindow.winfo_reqwidth())
        if self.scrollwindow.winfo_reqheight() != self.canv.winfo_height():
            # update the canvas's width to fit the inner frame
            self.canv.config(height=self.scrollwindow.winfo_reqheight())


# Plugin code

class FragMapPlugin:

    def set_map_filename(self, maptype, value):
        setattr(self, '{}_fragmap_filename'.format(maptype), value)

    def get_map_filename(self, maptype):
        return getattr(self, '{}_fragmap_filename'.format(maptype))

    def set_map_filename_value(self, maptype, fname):
        self.get_map_filename(maptype).setvalue(fname)

    def get_map_filename_value(self, maptype):
        return self.get_map_filename(maptype).getvalue()

    def __init__(self):
        self.fragmaps_available = FRAGMAPS_AVAILABLE
        self.fragmaps_enabled = FRAGMAPS_ENABLED
        # for ligands tab
        self.ligandsmap = {} # main data storage
        self.ligandsloaded = [] # keeping track of ligands that are loaded
        # for mcligands tab
        self.mcligandsmap = {} # main data storage
        self.mcligandsloaded = [] # keeping track of ligands that are loaded

    def create_dialog(self, app=None):
        if app is None:
            app = pymol.plugins.get_pmgapp()
        self.parent = app.root

        # Create the dialog.
        self.dialog = Pmw.Dialog(
            self.parent,
            buttons = (
                'Visualize FragMap',
                'Load Ligands',
                'Load SILCS-MC',
                'Exit FragMap tools'),
            title = 'PyMOL FragMap Tools',
            command = self.execute,
        )
        self.dialog.geometry("%dx%d" % (810, 820))
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = Tkinter.Label(
            self.dialog.interior(),
            text = """PyMOL FragMap Tools

SilcsBio, 2021 - http://silcsbio.com""",
            background = '#3857a5',
            foreground = 'white',
            pady = 10,
        )
        w.grid(row=0, column=0, padx=4, pady=4, sticky=E+W)
        self.dialog.interior().grid_columnconfigure(0, weight=1)

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.grid(row=1, column=0, padx=10, pady=10, sticky=E+W+N+S)
        self.dialog.interior().grid_rowconfigure(1, weight=1)

        #
        # Page: FragMap locations
        #
        page = self.notebook.add('FragMap Locations')
        frame=Frame(page, width=750, height=700)
        frame.pack(expand=True, fill=BOTH)

        group = Pmw.Group(frame, tag_text='Locations')
        group.pack(fill='both', expand=1, padx=10, pady=5)

        sw = ScrolledWindow(group.interior())
        innerFrame = Frame(sw.scrollwindow)
        innerFrame.grid(row=0, column=0, padx=10, pady=10, sticky=W+N+E+S)
        group.interior().grid_columnconfigure(0, weight=1)

        # below seemed to be needed to expand the innerFrame to fill the page
        group.interior().grid_rowconfigure(0, weight=1)

        # add mapdir and prefix
        self.selection = Pmw.EntryField(
            innerFrame,
            labelpos='w',
            label_text='Selection to use: ',
            value='(polymer)',
        )
        self.mapdir = Pmw.EntryField(
            innerFrame,
            labelpos='w',
            label_text='Map folder: ',
            value='2b_gen_maps',
        )
        if pymol.cmd.get_object_list():
            default_object_name = pymol.cmd.get_object_list()[0]
        else:
            default_object_name = ""
        self.prefix = Pmw.EntryField(
            innerFrame,
            labelpos='w',
            label_text='Map prefix: ',
            value=default_object_name,
        )
        self.liganddir = Pmw.EntryField(
            innerFrame,
            labelpos='w',
            label_text='Ligand folder: ',
            value='',
        )
        self.mcsilcsdir = Pmw.EntryField(
            innerFrame,
            labelpos='w',
            label_text='SILCS-MC folder: ',
            value='',
        )

        alignlabels = []
        components = (self.mapdir, self.prefix)
        for entry in components:
            entry.pack(fill='x', padx = 20, pady = 2) # vertical
            alignlabels.append(entry)

        # register FragMap filename forms
        for fragmap_type in FragMap.fragmaps_enabled:
            fragmap = FragMap.fragmaps_available[fragmap_type]
            setter = lambda value: self.set_map_filename_value(fragmap_type, value)
            label = FileDialogButtonClassFactory.get(setter)
            filename_placeholder = '{{prefix}}.{}.gfe.dx'.format(fragmap_type)
            if fragmap_type == 'excl':
                filename_placeholder = '{prefix}.excl.dx'
            element = Pmw.EntryField(
                innerFrame,
                labelpos = 'w',
                label_text = '{}: '.format(fragmap['label']),
                value = os.path.join('{mapdir}',filename_placeholder),
            )
            self.set_map_filename(fragmap_type, element)
            self.get_map_filename(fragmap_type).pack(
                fill = 'x', padx = 20, pady = 2
            )
            alignlabels.append(self.get_map_filename(fragmap_type))
        Pmw.alignlabels(alignlabels)

        label = Tkinter.Label(
            innerFrame,
            pady = 10,
            justify=LEFT,
            text = """ """,
        )
        label.pack()

        #
        # Page: Visualization page
        #
        page = self.notebook.add('Visualization')
        frame = Frame(page, width=800, height=800)
        frame.pack(expand=True, fill=BOTH)

        group = VisualizationGroup(
            frame,
            tag_text='Visualization',
            visgroup_num=1
        )
        self.visualization_group = group
        group.pack(fill='both', expand=1, padx=10, pady=10)

        #
        # Page: Ligand page
        #
        page = self.notebook.add('Ligands')
        group = LigandsGroup(
            page,
            tag_text='Ligands',
            visgroup_num=1
        )
        self.ligands_group = group
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)

        #
        # Page: SILCS-MC page
        #
        page = self.notebook.add('SILCS-MC')
        group = McSilcsGroup(
            page,
            tag_text='SILCS-MC',
            visgroup_num=1
        )
        self.mcsilcs_group = group
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)

        #
        # Page: About page
        #
        page = self.notebook.add('About')
        group = Pmw.Group(page, tag_text='About PyMOL FragMaps Tools')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        text = """This plugin integrates PyMOL (http://PyMOL.org/) with FragMaps from SILCS simulation.

Documentation may be found at
http://silcsbio.com

In the simplest case,

1) Load a structure into PyMOL.
2) Start this plugin.
3) Click the "Visualize FragMap" button.

Contact info@silcsbio.com for support.

Created by SilcsBio, 2016
"""
        #
        # Add this as text in a scrollable pane.
        # Code based on Caver plugin
        # http://loschmidt.chemi.muni.cz/caver/index.php
        #
        interior_frame=Frame(group.interior())
        bar=Scrollbar(interior_frame,)
        text_holder=Text(interior_frame,yscrollcommand=bar.set,background="#ddddff",font="Times 14")
        bar.config(command=text_holder.yview)
        text_holder.insert(END,text)
        text_holder.pack(side=LEFT,expand="yes",fill="both")
        bar.pack(side=LEFT,expand="yes",fill="y")
        interior_frame.pack(expand="yes",fill="both")

        self.notebook.setnaturalsize()
        self.showAppModal()

    def showAppModal(self):
        #self.dialog.activate() #geometry = 'centerscreenfirst',globalMode = 'nograb')
        self.dialog.show()

    def execute(self, result, refocus=True):
        """Defines action upon pressing the buttons in the bottom row"""

        if result == 'FragMap Help':
            import webbrowser
            webbrowser.open("http://www.silcsbio.com")

        elif result == 'Load Ligands':
            #self.show_ligands()
            self.ligand_dialog()

        elif result == 'Load SILCS-MC':
            self.mcligand_dialog()

        elif result == 'Visualize FragMap':
            self.show_fragmap()

        else:
            self.quit()

    def show_fragmap(self):
        """Load FragMaps and show the "Visualization" tab in the plugin window"""

        for fragmap_type in FragMap.fragmaps_enabled:
            fragmap = FragMap.fragmaps_available[fragmap_type]
            fname = self.get_map_filename_value(fragmap_type)
            fname = fname.format(**{
                'mapdir': self.mapdir.getvalue(),
                'prefix': self.prefix.getvalue(),
            })

            # check if autodock map file exists but no dx format
            map_fname = '.'.join(fname.split('.')[:-1] + ['map'])
            if not os.path.exists(fname) and os.path.exists(map_fname):
                try:
                    ADGridMap(open(map_fname)).writeDX(fname)
                except:
                    pass
            if os.path.exists(fname):
                pymol.cmd.load(fname, '%s_%s' % ('map', fragmap_type))
            else:
                print("%s (%s) does not exists" % (fragmap['label'], fname))
                #self.visualization_group_1.maps[fragmap_type]['visible'].set(0)

        self.visualization_group.refresh()
        self.notebook.tab('Visualization').focus_set()
        self.notebook.selectpage('Visualization')

    def ligand_dialog(self):
        """Show file dialog for loading more ligands molecules."""

        options = {}
        options['filetypes'] = [
            ('Coordinate Files', "*.pdb *.mol2 *.sdf *.sd"),
        ]
        options['title'] = 'Load Ligands'
        molfiles = tkFileDialog.askopenfilenames(**options)
        self.liganddir.setvalue(' '.join(molfiles))
        self.show_ligands()

    def show_ligands(self):
        """Load ligands that matches with pattern provided by "liganddir"
        variable and show the "Ligands" tab in the plugin window.
        """

        liganddir = self.liganddir.getvalue()
        for pattern in liganddir.split(' '):
            molfiles = glob.glob(pattern)
            for molfile in molfiles:
                if molfile in self.ligandsloaded:
                    continue

                self.load_ligand(molfile)

                # makes the only first ligands to be visible by default
                if len(self.ligandsloaded) == 1:
                    self.ligandsmap[molfile]['visible'].set(1)
                    self.load_ligand(molfile)

        self.ligands_group.refresh()
        self.notebook.tab('Ligands').focus_set()
        self.notebook.selectpage('Ligands')

    def load_ligand(self, molfile):
        """Load Ligand into PyMOL. If the ligand is already loaded,
        it will simply turn on and off."""

        if molfile in self.ligandsloaded:
            ligand = self.ligandsmap[molfile]
            if ligand['visible'].get() == 1:
                pymol.cmd.enable(ligand['object'])
            else:
                pymol.cmd.disable(ligand['object'])
            return

        # ligandsmap dictionary carries ligand related properties
        # as well as visualization state
        self.ligandsmap[molfile] = {
            'visible': IntVar(),
            'gfe_checked': IntVar(),
            'color_checked': IntVar(),
            'lgfe': None,
            'le': None,
            'object': None,
        }
        self.ligandsloaded.append(molfile)
        
        molfile_base = os.path.basename(molfile)
        molfile_pre = os.path.splitext(molfile_base)[0]
        while molfile_pre in pymol.cmd.get_names():
            molfile_pre += ".1" 
        ligand = self.ligandsmap[molfile]
        pymol.cmd.load(molfile, object=molfile_pre, zoom=0)
        ligand['object'] = pymol.cmd.get_names()[-1]
        ligand['visible'].set(0)
        pymol.cmd.hide('line', ligand['object'])
        pymol.cmd.show('stick', '%s and not elem H' % ligand['object'])
        pymol.cmd.disable(ligand['object'])

        pymol.stored.bfactor = 0
        pymol.stored.nclass = 0
        pymol.stored.nheavy = 0
        selection = ligand['object']
        if molfile.endswith('.pdb'):
            # read LGFE from B factor column
            pymol.cmd.iterate(selection, "stored.bfactor = stored.bfactor + b")
            pymol.cmd.iterate("%s and not segid NCLA" % selection, "stored.nclass = stored.nclass + 1")
            pymol.cmd.iterate("%s and not elem H" % selection, "stored.nheavy = stored.nheavy + 1")
            with open(molfile) as fp:
                for line in fp:
                    if line.startswith("REMARK"):
                        if line.strip().split()[1] == "LGFE":
                            ligand['lgfe'] = float(line.strip().split()[-1])
                        if line.strip().split()[1] == "LE":
                            ligand['le'] = float(line.strip().split()[-1])
                    if line[:3] == 'END':
                        break

        if molfile.endswith('.sdf'):
            with open(molfile) as fp:
                for line in fp:
                    if line.endswith('<atom.prop.CLASS>\n'):
                        line = next(fp)
                        silcsclass = line.split()
                        pymol.stored.nclass = len(silcsclass)

                    if line.endswith('<LGFE>\n'):
                        line = next(fp)
                        ligand['lgfe'] = float(line)
                    
                    if line.endswith('<LE>\n'):
                        line = next(fp)
                        ligand['le'] = float(line)

                    if line.endswith('<atom.dprop.GFE>\n'):
                        line = next(fp)
                        gfe = [float(_) for _ in line.split()]
                        pymol.stored.bfactor = sum(gfe)

                    if line.startswith("$$$$"):
                        break

            for i in range(len(gfe)):
                if silcsclass[i] == 'NCLA': continue
                pymol.cmd.alter("%s and rank %d" % (selection, i), 'b=%.2f' % gfe[i])
                pymol.stored.nclass += 1

    def toggle_gfe_ligand(self, molfile, label):
        """Toggle atomic GFE contribution either by label or color"""

        ligand = self.ligandsmap[molfile]
        if label == 'color':
            if ligand['color_checked'].get():
                pymol.cmd.spectrum('b', 'green_white_red', ligand['object'], minimum=-1.0, maximum=1.0)
            else:
                cid = self.ligandsloaded.index(molfile)
                colors = ['g', 'c', 'm', 'y', 's', 'w', 'b', 'o', 'p', 'k']
                color_code = colors[cid % len(colors)]
                color = getattr(pymol.util, 'cba%s' % color_code)
                color(ligand['object'])
        if label == 'label':
            if ligand['gfe_checked'].get():
                #pymol.cmd.iterate('%s and not segid NCLA' % ligand['object'], 'label ')
                selection = '%s and not elem H' % ligand['object']
                pymol.cmd.label(selection, '"%.2f" % b')
                pymol.cmd.set('label_font_id', 10)
                pymol.cmd.set('label_color', 'green', '%s and b < 0' % selection)
                pymol.cmd.set('label_color', 'red', '%s and b > 0' % selection)
            else:
                pymol.cmd.hide('label', ligand['object'])

    def ligand_zoom(self, molfile):
        """Zoom ligand"""

        ligand = self.ligandsmap[molfile]
        pymol.cmd.zoom(ligand['object'], buffer=5)

    def define_pocket(self, target, radius):
        """Show side chains of protein atoms near selected ligands"""

        objects = []
        for molfile in self.ligandsloaded:
            ligand = self.ligandsmap[molfile]
            if ligand['visible'].get() == 1:
                objects.append(ligand['object'])
        selection = ["( %s and not elem H )" % o for o in objects]
        pymol.cmd.select('sele', "( byres ( ( %s ) around %d ) and %s ) and not elem H" % (' or '.join(selection), int(radius), target))
        pymol.cmd.create('pocket', 'sele')
        pymol.cmd.show('stick', 'pocket')
        pymol.util.cbag('pocket')

    def mcligand_dialog(self):
        """Show file dialog for loading more ligands molecules."""

        options = {}
        options['filetypes'] = [
            ('Coordinate Files', '*.pdb *.mol2 *.sdf')
        ]
        options['title'] = 'Load Ligands'
        molfiles = tkFileDialog.askopenfilenames(**options)
        self.mcsilcsdir.setvalue(' '.join(molfiles))
        self.show_mcsilcs()

    def show_mcsilcs(self):
        """Load SILCS-MC ligands and show the "SILCS-MC" tab in the plugin window"""

        mcsilcsdir = self.mcsilcsdir.getvalue()
        molfiles = glob.glob(mcsilcsdir)
        for molfile in molfiles:
            if molfile in self.mcligandsloaded:
                continue

            self.load_mcligand(molfile)

            # makes the only first ligands to be visible by default
            if len(self.mcligandsloaded) == 1:
                self.mcligandsmap[molfile]['visible'].set(1)
                self.load_mcligand(molfile)

        pymol.cmd.sync(timeout=5.0)
        self.mcsilcs_group.refresh()
        self.notebook.tab('SILCS-MC').focus_set()
        self.notebook.selectpage('SILCS-MC')

    def load_mcligand(self, molfile, frame=None):
        """Load Ligand into PyMOL"""

        if molfile in self.mcligandsloaded:
            ligand = self.mcligandsmap[molfile]
            if ligand['visible'].get() == 1:
                pymol.cmd.enable(ligand['object'])
            else:
                pymol.cmd.disable(ligand['object'])

            if frame == 'prev':
                frame = ligand['current_frame'].get() - 1
                if frame < 1: frame = 0
            if frame == 'next':
                frame = ligand['current_frame'].get() + 1
                if frame > ligand['nframes']: frame = ligand['nframes']
            if frame and frame != ligand['current_frame'].get():
                pymol.cmd.set('state', frame, ligand['object'])
                ligand['current_frame'].set(frame)
                ligand['lgfe'].set("%.2f kcal/mol" % ligand['lgfe_arr'][int(float(frame))-1])
            return

        # ligandsmap dictionary carries ligand related properties
        # as well as visualization state
        self.mcligandsmap[molfile] = {
            'visible': IntVar(),
            'lgfe': StringVar(),
            'object': None,
            'nframes': 0,
            'current_frame': IntVar(),
            'lgfe_arr': [],
        }
        self.mcligandsloaded.append(molfile)
        molfile_base = os.path.basename(molfile)
        molfile_pre = os.path.splitext(molfile_base)[0]
        while molfile_pre in pymol.cmd.get_names():
            molfile_pre += ".1" 

        ligand = self.mcligandsmap[molfile]
        pymol.cmd.load(molfile, object=molfile_pre, multiplex=0, zoom=0)
        ligand['object'] = pymol.cmd.get_names()[-1]
        ligand['nframes'] = pymol.cmd.count_states(ligand['object'])
        ligand['current_frame'].set(0)
        ligand['lgfe'].set('NA')
        ligand['visible'].set(0)
        pymol.cmd.hide('line', ligand['object'])
        pymol.cmd.show('stick', '%s and not elem H' % ligand['object'])
        pymol.cmd.disable(ligand['object'])

        # compute LGFE scores of each frame
        lgfe = 0
        nclass = 0
        nheavy = 0
        factor = None

        if molfile.endswith('.pdb'):
            for line in open(molfile):
                if line.startswith("REMARK"):
                    if line.strip().split()[1] == "LGFE":
                        lgfe = float(line.strip().split()[-1])
                if line[:3] == 'END':
                    factor = float(nheavy) / nclass
                    ligand['lgfe_arr'].append(lgfe)
                    lgfe = 0
                    continue
                if line[:4] == 'ATOM':
                    atomname = line[11:17].strip()
                    atomtype = atomname[0]
                    clss = line[72:76]
                    nclass += 1 if clss != 'NCLA' else 0
                    nheavy += 1 if atomtype != 'H' else 0

        if molfile.endswith('.sdf'):
            # we stopped renormalizing LGFE scores in 2017
            factor = 1
            with open(molfile) as fp:
                for line in fp:
                    if line.startswith('$$$$'):
                        lgfe = lgfe * factor
                        ligand['lgfe_arr'].append(lgfe)
                        lgfe = 0
                        continue
                    if line.endswith('<atom.dprop.GFE>\n'):
                        line = next(fp)
                        gfe = [float(_) for _ in line.split()]
                    if line.endswith('<LGFE>\n'):
                        line = next(fp)
                        lgfe = float(line)

        # wait until each ligand is completely loaded
        pymol.cmd.sync(timeout=5.0)

        # set frame to 1
        self.load_mcligand(molfile, 1)

    def mcligand_zoom(self, molfile):
        """Zoom ligand"""
        ligand = self.mcligandsmap[molfile]
        pymol.cmd.zoom(ligand['object'], buffer=5)

    def quit(self):
        self.dialog.destroy() # stops CPU hogging, perhaps fixes Ubuntu bug MGL

    def fixColumns(self,sel):
        """
        Make sure that everything fits into the correct columns.
        This means doing some rounding. It also means getting rid of
        chain, occupancy and b-factor information.
        """
        #pymol.cmd.alter_state(1,'all','(x,y,z)=(int(x*1000)/1000.0, int(y*1000)/1000.0, int(z*1000)/1000.0)')
        #pymol.cmd.alter_state(1,'all','(x,y,z)=float("%.2f"%x),float("%.2f"%y),float("%.2f"%z)')
        pymol.cmd.alter_state(1,'all','(x,y,z)=float("%.3f"%x),float("%.3f"%y),float("%.3f"%z)')
        #pymol.cmd.alter(sel,'chain=""')
        pymol.cmd.alter(sel,'b=0')
        pymol.cmd.alter(sel,'q=0')

FragMap = FragMapPlugin()

# PmwGroups

class VisualizationGroup(Pmw.Group):
    """FragMap visualization (adjust level)"""

    def __init__(self, *args, **kwargs):
        my_options = 'visgroup_num'.split()
        for option in my_options:
            # use these options as attributes of this class
            # and remove them from the kwargs dict before
            # passing on to Pmw.Group.__init__().
            setattr(self, option, kwargs.pop(option))

        kwargs['tag_text'] = kwargs['tag_text'] + ' (%s)'%self.visgroup_num
        Pmw.Group.__init__(self, *args, **kwargs)

        FragMap.maps = copy.deepcopy(FragMap.fragmaps_available)
        for k in FragMap.fragmaps_enabled:
            visible = FragMap.maps[k]['visible']
            FragMap.maps[k]['visible'] = IntVar()
            FragMap.maps[k]['visible'].set(visible)
        self.refresh()

    def refresh(self):
        things_to_kill = 'error_label update_buttonbox'.split()
        for thing in things_to_kill:
            try:
                getattr(self, thing).destroy()
                delattr(self, thing)
                #print "destroyed",thing
            except AttributeError:
                #print "couldn't destroy",thing

                # note: this attributeerror will also hit if getattr(self,thing) misses.
                # another note: both halves of the if/else make an update_buttonbox.
                # if you rename the one in the top half to something else, you'll cause nasty Pmw errors.
                pass

        sw = ScrolledWindow(self.interior())
        self.sw = sw

        fragmapFrame = Frame(self.sw.scrollwindow)
        fragmapFrame.grid(row=0, column=0, padx=10, pady=10, sticky=W+N+E+S)
        self.interior().grid_columnconfigure(0, weight=1)

        # decoration
        self.interior().grid_rowconfigure(0, weight=1)
        separator = Frame(self.interior(), height=2, bd=1, relief=SUNKEN)
        separator.grid(row=1, column=0, sticky=W+E, padx=10, pady=4)

        # protein surface
        surfaceFrame = Frame(self.interior())
        surface_label_0 = Tkinter.Label( surfaceFrame, text = "Surface Type", )
        surface_label_1 = Tkinter.Label( surfaceFrame, text = "Molecule" )
        surface_entry_flag = IntVar()
        surface_entry_surf = Checkbutton(
            surfaceFrame,
            text = "Protein Surface",
            variable=surface_entry_flag,
            command = lambda surf_type='surf': self.showSurface(surf_type),
        )
        protein_list_surf = Pmw.ComboBox(
            surfaceFrame,
            scrolledlist_items=pymolutil.getMacroMolecules(),
            entryfield_entry_state='readonly',
        )
        cartoon_entry_flag = IntVar()
        surface_entry_cartoon = Checkbutton(
            surfaceFrame,
            text = "Protein Cartoon",
            variable=cartoon_entry_flag,
            command = lambda surf_type='cartoon': self.showSurface(surf_type),
        )
        protein_list_cartoon = Pmw.ComboBox(
            surfaceFrame,
            scrolledlist_items=pymolutil.getMacroMolecules(),
            entryfield_entry_state='readonly',
        )
        surface_entry_2 = Checkbutton(
            surfaceFrame,
            text = "Exclusion Map",
            variable = FragMap.maps['excl']['visible'],
            command = lambda name='excl': self.updateFragMap(name),
        )
        surface_label_0.grid(row=0, column=0)
        surface_label_1.grid(row=0, column=1)
        surface_entry_surf.grid(row=1, column=0, sticky=W)
        surface_entry_cartoon.grid(row=2, column=0, sticky=W)
        protein_list_surf.grid(row=1, column=1, padx=30, sticky=W)
        protein_list_cartoon.grid(row=2, column=1, padx=30, sticky=W)
        surfaceFrame.grid_columnconfigure(1, weight=1)
        surfaceFrame.grid(row=2, column=0, padx=10, ipady=5, sticky=E+W)
        if pymolutil.getMacroMolecules():
            protein_list_surf.selectitem(0)
            protein_list_cartoon.selectitem(0)
        self.surf_entry_flag = surface_entry_flag
        self.cartoon_entry_flag = cartoon_entry_flag
        self.protein_list_surf = protein_list_surf
        self.protein_list_cartoon = protein_list_cartoon

        if not [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:map'] and [i for i in pymol.cmd.get_names() if pymol.cmd.get_type(i)=='object:molecule']:
            self.error_label = Tkinter.Label(
                fragmapFrame,
                pady = 10,
                justify=LEFT,
                text = '''You must have at least a molecule and a map loaded.
If you have a molecule and a map loaded, please click "Update"''',
            )
            self.error_label.pack()
            self.update_buttonbox = Pmw.ButtonBox(fragmapFrame, padx=0)
            self.update_buttonbox.pack()
            self.update_buttonbox.add('Update',command=self.refresh)
            return

        header_label_0 = Tkinter.Label(
            fragmapFrame,
            pady = 2,
            justify=LEFT,
            text = "FragMap Type",
        )
        header_label_1 = Tkinter.Label(
            fragmapFrame,
            pady = 2,
            justify=LEFT,
            text = "GFE Level",
        )
        header_label_0.grid(row=0, column=0)
        header_label_1.grid(row=0, column=1, columnspan=3)

        offset = 1
        for rowidx, map_name in enumerate(FragMap.fragmaps_enabled):
            if map_name == 'excl':
                continue

            #map_group = Pmw.Group(self.mm_group.interior(),ring_borderwidth=0,tag_pyclass = None)
            map_entry = Checkbutton(
                fragmapFrame,
                text = FragMap.maps[map_name]['label'],
                variable = FragMap.maps[map_name]['visible'],
                command = lambda name=map_name: self.updateFragMap(name),
            )
            #map_entry.pack(side=LEFT)
            map_entry.grid(row=rowidx+offset, column=0, sticky='w')
            setattr(self, 'map_%s' % map_name, map_entry)

            map_minus = Tkinter.Button(
                fragmapFrame,
                text='-', padx=0, pady=0, width=2,
                command=lambda value='-', name=map_name: self.updateFragMap(name, value),
            )
            map_plus = Tkinter.Button(
                fragmapFrame,
                text='+', padx=0, pady=0, width=2,
                command=lambda value='+', name=map_name: self.updateFragMap(name, value),
            )
            map_bar = Tkinter.Scale(
                fragmapFrame,
                from_ = -2,
                to = 2,
                orient=HORIZONTAL,
                resolution=0.1,
                digit=2,
                showvalue=False,
                bg='#7f7f7f',
                activebackground='#7f7f7f',
                command=lambda value, name=map_name: self.updateFragMap(name, value),
            )
            map_level_label = StringVar()
            map_level = Tkinter.Label(
                fragmapFrame,
                pady = 2,
                justify=RIGHT,
                anchor=E,
                textvariable=map_level_label,
            )
            map_level_unit = Tkinter.Label(
                fragmapFrame,
                pady = 2,
                justify=LEFT,
                text='kcal/mol',
            )
            map_bar.set(FragMap.maps[map_name]['level'])
            map_bar.grid(row=rowidx+offset, column=2)
            map_minus.grid(row=rowidx+offset, column=1)
            map_plus.grid(row=rowidx+offset, column=3)
            map_level.grid(row=rowidx+offset, column=4, padx=10)
            map_level_unit.grid(row=rowidx+offset, column=5)
            setattr(self, 'map_%s_level' % map_name, map_bar)
            setattr(self, 'map_%s_level_label' % map_name, map_level_label)
            #map_group.pack(fill = 'both', expand = 1, padx = 6, pady = 2, side=TOP)

        surface_entry_2.grid(row=3, column=0, sticky=W)

    def isFragMapReady(self):
        for map_name in FragMap.fragmaps_enabled:
            if not hasattr(self, 'map_%s_level'%map_name): return False
        return True

    def showSurface(self, surface_type):
        object_name = getattr(self, 'protein_list_%s' % surface_type).get()
        flag = getattr(self, '%s_entry_flag' % surface_type).get()
        new_object_name = 'protein-%s' % surface_type
        if flag:
            pymol.cmd.select('sele', "%s not elem H" % (object_name))
            pymol.cmd.create(new_object_name, 'sele')
            pymol.cmd.hide('everything', new_object_name)
            if surface_type == 'surf':
                pymol.cmd.show('surface', new_object_name)
                pymol.cmd.color('gray60', new_object_name)
                pymol.cmd.zoom(object_name)
            if surface_type == 'cartoon':
                pymol.cmd.show('cartoon', new_object_name)
                pymol.util.cbc(new_object_name)
                pymol.cmd.zoom(object_name)
        else:
            pymol.cmd.hide("everything", new_object_name)

    def updateFragMap(self, map_name, value=None):
        map_objects_loaded = pymol.cmd.get_names()
        fragmap = FragMap.maps[map_name]
        use_surf = fragmap.get('use_surf', False)
        map_object = 'map_%s'%map_name
        if use_surf:
            fragmap_name = 'map_%s_surf' % map_name
            fragmap_fn = pymol.cmd.isosurface
        else:
            fragmap_name = 'map_%s_mesh' % map_name
            fragmap_fn = pymol.cmd.isomesh

        map_visible = FragMap.maps[map_name]['visible'].get()
        map_level = float(FragMap.maps[map_name]['level'])
        if value:
            if value in ('-', '+'):
                map_level += -0.1 if value == '-' else 0.1
            else:
                map_level = float(value)
        FragMap.maps[map_name]['level'] = map_level
        try:
            getattr(self, 'map_%s_level' % map_name).set(map_level)
            getattr(self, 'map_%s_level_label' % map_name).set(str(map_level))
        except:
            pass
        if value in ('-', '+'):
            # updateFragMap will be called again via scroll bar event
            return
        if map_object not in map_objects_loaded:
            #print("Map not found: %s" % map_object)
            return

        fragmap_loaded = fragmap_name in map_objects_loaded
        if fragmap_loaded:
            if map_visible:
                pymol.cmd.enable(fragmap_name)
                fragmap_fn(fragmap_name, map_object, map_level)
            else:
                pymol.cmd.disable(fragmap_name)
                return
        else:
            if map_visible and map_name != 'excl':
                pymol.cmd.map_double(map_object)
            if map_visible:
                fragmap_fn(fragmap_name, map_object, map_level)

        if map_visible and not fragmap_loaded:
            # first time loading the map
            pymol.cmd.color(FragMap.maps[map_name]['color'], fragmap_name)

    def _validateFragMapLevel(self, text):
        try:
            float(text)
        except:
            return -1
        if self.isFragMapReady():
            self.updateFragMap()
        return 1

# Ligands

class LigandsGroup(Pmw.Group):
    """Ligands tab"""

    def __init__(self,*args,**kwargs):
        my_options = 'visgroup_num'.split()
        for option in my_options:
            # use these options as attributes of this class
            # and remove them from the kwargs dict before
            # passing on to Pmw.Group.__init__().
            setattr(self,option,kwargs.pop(option))

        kwargs['tag_text'] = kwargs['tag_text'] + ' (%s)'%self.visgroup_num
        Pmw.Group.__init__(self,*args,**kwargs)
        self.refresh()

    def refresh(self):
        things_to_kill = 'error_label update_buttonbox'.split()
        for thing in things_to_kill:
            try:
                getattr(self,thing).destroy()
                delattr(self,thing)
                #print "destroyed",thing
            except AttributeError:
                #print("couldn't destroy",thing)
                pass
            except:
                pass

        if not FragMap.ligandsloaded:
            self.error_label = Tkinter.Label(
                self.interior(),
                pady = 10,
                justify=LEFT,
                text = '''You must have at least a molecule and a map loaded.
If you have a molecule and a map loaded, please click "Update"''',
            )
            self.error_label.pack()
            self.update_buttonbox = Pmw.ButtonBox(self.interior(), padx=0)
            self.update_buttonbox.pack()
            self.update_buttonbox.add('Update',command=self.refresh)
            return

        ligandFrame = Frame(self.interior())
        ligandFrame.grid(row=0, column=0, padx=10, pady=10, sticky=W+N+E+S)
        self.interior().grid_columnconfigure(0, weight=1)

        header_label_0 = Tkinter.Label(
            ligandFrame,
            text = "Ligand Filename",
        )
        header_label_1 = Tkinter.Label(
            ligandFrame,
            text = "LGFE score",
        )
        header_label_2 = Tkinter.Label(
            ligandFrame,
            text = "LE score",
        )
        header_label_3 = Tkinter.Label(
            ligandFrame,
            text = "GFE",
        )
        header_label_0.grid(row=0, column=0)
        header_label_1.grid(row=0, column=1)
        header_label_2.grid(row=0, column=2)
        header_label_3.grid(row=0, column=3)

        offset = 1
        filenamelengthcap = 50
        for rowidx, molfile in enumerate(FragMap.ligandsloaded):
            ligand = FragMap.ligandsmap[molfile]
            filename = molfile
            if len(filename) > filenamelengthcap:
                filename = '... %s' % molfile[-filenamelengthcap:]

            ligand_entry = Checkbutton(
                ligandFrame,
                text = filename,
                variable = ligand['visible'],
                onvalue=1,
                offvalue=0,
                command = lambda name=molfile: FragMap.load_ligand(name),
            )
            ligand_entry.grid(row=rowidx+offset, column=0, sticky='w')
            ligand['checkbox'] = ligand_entry

            ligand_lgfe_label = Tkinter.Label(
                ligandFrame,
                text = "%.2f kcal/mol" % ligand['lgfe'] if ligand['lgfe'] else "NA",
            )
            ligand_lgfe_label.grid(row=rowidx+offset, column=1, sticky='w', padx=5)

            ligand_le_label = Tkinter.Label(
                ligandFrame,
                text = "%.2f kcal/mol" % ligand['le'] if ligand['le'] else "NA",
            )
            ligand_le_label.grid(row=rowidx+offset, column=2, sticky='w', padx=5)

            gfe_group = Frame(ligandFrame)
            ligand_gfe_label_toggle = Checkbutton(
                gfe_group,
                text = 'GFE',
                variable = ligand['gfe_checked'],
                command = lambda name=molfile: FragMap.toggle_gfe_ligand(name, 'label'),
            )
            ligand_gfe_color_toggle = Checkbutton(
                gfe_group,
                text = 'Color',
                variable = ligand['color_checked'],
                command = lambda name=molfile: FragMap.toggle_gfe_ligand(name, 'color'),
            )
            ligand_gfe_label_toggle.grid(row=0, column=0, sticky='w')
            ligand_gfe_color_toggle.grid(row=0, column=1, sticky='w')
            gfe_group.grid(row=rowidx+offset, column=3, sticky='w', padx=5)

            ligand_zoom_button = Tkinter.Button(
                ligandFrame,
                text='Zoom', padx=0, pady=2, width=8,
                command=lambda name=molfile: FragMap.ligand_zoom(name),
            )
            ligand_zoom_button.grid(row=rowidx+offset, column=4, sticky='w', padx=5)

        # decoration
        self.interior().grid_rowconfigure(0, weight=1)
        separator = Frame(self.interior(), height=2, bd=1, relief=SUNKEN)
        separator.grid(row=1, column=0, sticky=W+E, padx=10, pady=4)

        # ligand pocket
        pocketFrame = Frame(self.interior())
        pocket_label_0 = Tkinter.Label( pocketFrame, text = "Show Protein Atoms Near the Selected Ligands:", )
        pocket_label_1 = Tkinter.Label( pocketFrame, text = "Protein Molecule:" )
        protein_list = Pmw.ComboBox(
            pocketFrame,
            scrolledlist_items=pymolutil.getMacroMolecules(),
            entryfield_entry_state='readonly',
        )
        protein_list.selectitem(0)
        pocket_radius = Pmw.EntryField(
            pocketFrame,
            labelpos='w',
            label_text='Radius:',
            value='5',
            entry_width=4,
        )
        pocket_label_2 = Tkinter.Label( pocketFrame, text = "Å" )
        pocket_button = Tkinter.Button(
            pocketFrame,
            text='Define Pocket',
            pady=2,
            command=lambda: FragMap.define_pocket(protein_list.get(), pocket_radius.get())
        )
        pocket_label_0.grid(row=0, column=0, columnspan=5)
        pocket_label_1.grid(row=1, column=0, padx=5)
        protein_list.grid(row=1, column=1, padx=5)
        pocket_radius.grid(row=1, column=2, padx=5)
        pocket_label_2.grid(row=1, column=3, padx=2)
        pocket_button.grid(row=1, column=4, padx=5)
        pocketFrame.grid(row=2, column=0, padx=10, pady=5)


# SILCS-MC

class McSilcsGroup(Pmw.Group):
    """SILCS-MC tab"""

    def __init__(self,*args,**kwargs):
        my_options = 'visgroup_num'.split()
        for option in my_options:
            # use these options as attributes of this class
            # and remove them from the kwargs dict before
            # passing on to Pmw.Group.__init__().
            setattr(self,option,kwargs.pop(option))

        kwargs['tag_text'] = kwargs['tag_text'] + ' (%s)'%self.visgroup_num
        Pmw.Group.__init__(self,*args,**kwargs)
        self.refresh()

    def refresh(self):
        things_to_kill = 'error_label update_buttonbox'.split()
        for thing in things_to_kill:
            try:
                getattr(self,thing).destroy()
                delattr(self,thing)
            except AttributeError:
                #print "couldn't destroy",thing
                pass

        if not FragMap.mcligandsloaded:
            self.error_label = Tkinter.Label(
                self.interior(),
                pady = 10,
                justify=LEFT,
                text = '''You must have at least a molecule and a map loaded.
If you have a molecule and a map loaded, please click "Update"''',
            )
            self.error_label.pack()
            self.update_buttonbox = Pmw.ButtonBox(self.interior(), padx=0)
            self.update_buttonbox.pack()
            self.update_buttonbox.add('Update',command=self.refresh)
            return

        header_label_0 = Tkinter.Label(
            self.interior(),
            text = "Ligand Filename",
        )
        header_label_1 = Tkinter.Label(
            self.interior(),
            text = "Conformation",
        )
        header_label_2 = Tkinter.Label(
            self.interior(),
            text = "LGFE",
        )
        header_label_0.grid(row=0, column=0)
        header_label_1.grid(row=0, column=1)
        header_label_2.grid(row=0, column=3)

        offset = 1
        filenamelengthcap = 50
        for rowidx, molfile in enumerate(FragMap.mcligandsloaded):
            ligand = FragMap.mcligandsmap[molfile]

            filename = molfile
            if len(filename) > filenamelengthcap:
                filename = '... %s' % molfile[-filenamelengthcap:]

            ligand_entry = Checkbutton(
                self.interior(),
                text = filename,
                variable = ligand['visible'],
                onvalue=1,
                offvalue=0,
                command = lambda name=molfile: FragMap.load_mcligand(name),
            )
            ligand_entry.grid(row=rowidx+offset, column=0, sticky='w')

            conf_group = Frame(self.interior())
            conf_left = Tkinter.Button(
                conf_group,
                text='<', padx=2, pady=0, width=4,
                command=lambda value='prev', name=molfile: FragMap.load_mcligand(name, value),
            )
            conf_right = Tkinter.Button(
                conf_group,
                text='>', padx=2, pady=0, width=4,
                command=lambda value='next', name=molfile: FragMap.load_mcligand(name, value),
            )
            conf_bar = Tkinter.Scale(
                conf_group,
                from_ = 1,
                to = ligand['nframes'],
                orient=HORIZONTAL,
                resolution=1,
                digit=1,
                showvalue=False,
                bg='#7f7f7f',
                activebackground='#7f7f7f',
                variable=ligand['current_frame'],
                command=lambda value, name=molfile: FragMap.load_mcligand(name, value),
            )
            conf_group.grid(row=rowidx+offset, column=1, padx=5)
            conf_bar.set(0)
            conf_bar.grid(row=0, column=1)
            conf_left.grid(row=0, column=0)
            conf_right.grid(row=0, column=2)

            conf_frame_label = Tkinter.Label(
                self.interior(),
                width=5,
                anchor=E,
                textvariable=ligand['current_frame'],
            )
            conf_frame_label.grid(row=rowidx+offset, column=2, padx=5)

            conf_lgfe_level = Tkinter.Label(
                self.interior(),
                textvariable=ligand['lgfe'],
            )
            conf_lgfe_level.grid(row=rowidx+offset, column=3, padx=5)

            ligand_zoom_button = Tkinter.Button(
                self.interior(),
                text='Zoom', padx=0, pady=2, width=8,
                command=lambda name=molfile: FragMap.mcligand_zoom(name),
            )
            ligand_zoom_button.grid(row=rowidx+offset, column=4, sticky='w', padx=5)

        self.pack(fill = 'both', expand = 1, padx = 10, pady = 10, side=TOP)


# PmwExtensions

"""
This contains all of the visualization groups that we'll use for our
PMW interface.
"""

# Generically useful PMW extensions

import os,fnmatch,time

def _errorpop(master,text):
    d=Pmw.MessageDialog(
        master,
        title="Error",
        message_text=text,
        buttons=("OK",)
    )
    d.component('message').pack(ipadx=15,ipady=15)
    d.activate()
    d.destroy()

class PmwFileDialog(Pmw.Dialog):
    """File Dialog using Pmw"""
    def __init__(self, parent = None, **kw):
        # Define the megawidget options.
        optiondefs = (
	    ('filter',    '*',              self.newfilter),
	    ('directory', os.getcwd(),      self.newdir),
	    ('filename',  '',               self.newfilename),
	    ('historylen',10,               None),
	    ('command',   None,             None),
            ('info',      None,             None),
	    )
        self.defineoptions(kw, optiondefs)
        # Initialise base class (after defining options).
        Pmw.Dialog.__init__(self, parent)

        self.withdraw()

        # Create the components.
        interior = self.interior()

        if self['info'] is not None:
            rowoffset=1
            dn = self.infotxt()
            dn.grid(row=0,column=0,columnspan=2,padx=3,pady=3)
        else:
            rowoffset=0

        dn = self.mkdn()
        dn.grid(row=0+rowoffset,column=0,columnspan=2,padx=3,pady=3)
        del dn

        # Create the directory list component.
        dnb = self.mkdnb()
        dnb.grid(row=1+rowoffset,column=0,sticky='news',padx=3,pady=3)
        del dnb

        # Create the filename list component.
        fnb = self.mkfnb()
        fnb.grid(row=1+rowoffset,column=1,sticky='news',padx=3,pady=3)
        del fnb

        # Create the filter entry
        ft = self.mkft()
        ft.grid(row=2+rowoffset,column=0,columnspan=2,padx=3,pady=3)
        del ft

        # Create the filename entry
        fn = self.mkfn()
        fn.grid(row=3+rowoffset,column=0,columnspan=2,padx=3,pady=3)
        fn.bind('<Return>',self.okbutton)
        del fn

        # Buttonbox already exists
        bb=self.component('buttonbox')
        bb.add('OK',command=self.okbutton)
        bb.add('Cancel',command=self.cancelbutton)
        del bb

        Pmw.alignlabels([self.component('filename'),
			 self.component('filter'),
			 self.component('dirname')])

    def infotxt(self):
        """ Make information block component at the top """
        return self.createcomponent(
                'infobox',
                (), None,
                Tkinter.Label, (self.interior(),),
                width=51,
                relief='groove',
                foreground='darkblue',
                justify='left',
                text=self['info']
            )

    def mkdn(self):
        """Make directory name component"""
        return self.createcomponent(
	    'dirname',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['directory'],
	    entryfield_entry_width=40,
            entryfield_validate=self.dirvalidate,
	    selectioncommand=self.setdir,
	    labelpos='w',
	    label_text='Directory:')

    def mkdnb(self):
        """Make directory name box"""
        return self.createcomponent(
	    'dirnamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='directories',
	    labelpos='n',
	    hscrollmode='none',
	    dblclickcommand=self.selectdir)

    def mkft(self):
        """Make filter"""
        return self.createcomponent(
	    'filter',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filter'],
	    entryfield_entry_width=40,
	    selectioncommand=self.setfilter,
	    labelpos='w',
	    label_text='Filter:')

    def mkfnb(self):
        """Make filename list box"""
        return self.createcomponent(
	    'filenamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='files',
	    labelpos='n',
	    hscrollmode='none',
	    selectioncommand=self.singleselectfile,
	    dblclickcommand=self.selectfile)

    def mkfn(self):
        """Make file name entry"""
        return self.createcomponent(
	    'filename',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filename'],
	    entryfield_entry_width=40,
            entryfield_validate=self.filevalidate,
	    selectioncommand=self.setfilename,
	    labelpos='w',
	    label_text='Filename:')

    def dirvalidate(self,string):
        if os.path.isdir(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL

    def filevalidate(self,string):
        if string=='':
            return Pmw.PARTIAL
        elif os.path.isfile(string):
            return Pmw.OK
        elif os.path.exists(string):
            return Pmw.PARTIAL
        else:
            return Pmw.OK

    def okbutton(self):
        """OK action: user thinks he has input valid data and wants to
           proceed. This is also called by <Return> in the filename entry"""
        fn=self.component('filename').get()
        self.setfilename(fn)
        if self.validate(fn):
            self.canceled=0
            self.deactivate()

    def cancelbutton(self):
        """Cancel the operation"""
        self.canceled=1
        self.deactivate()

    def tidy(self,w,v):
        """Insert text v into the entry and at the top of the list of
           the combobox w, remove duplicates"""
        if not v:
            return
        entry=w.component('entry')
        entry.delete(0,'end')
        entry.insert(0,v)
        list=w.component('scrolledlist')
        list.insert(0,v)
        index=1
        while index<list.index('end'):
            k=list.get(index)
            if k==v or index>self['historylen']:
                list.delete(index)
            else:
                index=index+1
        w.checkentry()

    def setfilename(self,value):
        if not value:
            return
        value=os.path.join(self['directory'],value)
        dir,fil=os.path.split(value)
        self.configure(directory=dir,filename=value)

        c=self['command']
        if callable(c):
            c()

    def newfilename(self):
        """Make sure a newly set filename makes it into the combobox list"""
        self.tidy(self.component('filename'),self['filename'])

    def setfilter(self,value):
        self.configure(filter=value)

    def newfilter(self):
        """Make sure a newly set filter makes it into the combobox list"""
        self.tidy(self.component('filter'),self['filter'])
        self.fillit()

    def setdir(self,value):
        self.configure(directory=value)

    def newdir(self):
        """Make sure a newly set dirname makes it into the combobox list"""
        self.tidy(self.component('dirname'),self['directory'])
        self.fillit()

    def singleselectfile(self):
        """Single click in file listbox. Move file to "filename" combobox"""
        cs=self.component('filenamebox').curselection()
        if cs!=():
            value=self.component('filenamebox').get(cs)
            self.setfilename(value)

    def selectfile(self):
        """Take the selected file from the filename, normalize it, and OK"""
        self.singleselectfile()
        value=self.component('filename').get()
        self.setfilename(value)
        if value:
            self.okbutton()

    def selectdir(self):
        """Take selected directory from the dirnamebox into the dirname"""
        cs=self.component('dirnamebox').curselection()
        if cs!=():
            value=self.component('dirnamebox').get(cs)
            dir=self['directory']
            if not dir:
                dir=os.getcwd()
            if value:
                if value=='..':
                    dir=os.path.split(dir)[0]
                else:
                    dir=os.path.join(dir,value)
            self.configure(directory=dir)
            self.fillit()

    def askfilename(self,directory=None,filter=None):
        """The actual client function. Activates the dialog, and
           returns only after a valid filename has been entered
           (return value is that filename) or when canceled (return
           value is None)"""
        if directory!=None:
            self.configure(directory=directory)
        if filter!=None:
            self.configure(filter=filter)
        self.fillit()
        self.canceled=1 # Needed for when user kills dialog window
        self.activate()
        if self.canceled:
            return None
        else:
            return self.component('filename').get()

    lastdir=""
    lastfilter=None
    lasttime=0
    def fillit(self):
        """Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
        if self.lastdir==self['directory'] and self.lastfilter==self['filter'] and self.lasttime>os.stat(self.lastdir)[8]:
            return
        self.lastdir=self['directory']
        self.lastfilter=self['filter']
        self.lasttime=time.time()
        dir=self['directory']
        if not dir:
            dir=os.getcwd()
        dirs=['..']
        files=[]
        try:
            fl=os.listdir(dir)
            fl.sort()
        except os.error as arg:
            if arg[0] in (2,20):
                return
            raise
        for f in fl:
            if os.path.isdir(os.path.join(dir,f)):
                dirs.append(f)
            else:
                filter=self['filter']
                if not filter:
                    filter='*'
                if fnmatch.fnmatch(f,filter):
                    files.append(f)
        self.component('filenamebox').setlist(files)
        self.component('dirnamebox').setlist(dirs)

    def validate(self,filename):
        """Validation function. Should return 1 if the filename is valid,
           0 if invalid. May pop up dialogs to tell user why. Especially
           suited to subclasses: i.e. only return 1 if the file does/doesn't
           exist"""
        return 1

class PmwExistingFileDialog(PmwFileDialog):
    def filevalidate(self,string):
        if os.path.isfile(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL

    def validate(self,filename):
        if os.path.isfile(filename):
            return 1
        elif os.path.exists(filename):
            _errorpop(self.interior(),"This is not a plain file")
            return 0
        else:
            _errorpop(self.interior(),"Please select an existing file")
            return 0

class FileDialogButtonClassFactory:
    def get(fn,filter='*'):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class FileDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.

            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                Tkinter.Button.__init__(*(self, master, cnf), **kw)
                self.configure(command=self.set)
            def set(self):
                fd = PmwFileDialog(self.master,filter=filter)
                fd.title('Please choose a file')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)

# Map to DX conversion

class ADGridMap:
    """Class for handling AutoDock Map files"""

    def __init__(self, fp=None, name='map'):

        self.name = ''
        self.npts = [0, 0, 0]
        self.n = [0, 0, 0]
        self.center = [0, 0, 0]
        self.origin = [0, 0, 0]
        self.nelem = 0
        self.spacing = 0.
        self.values = []
        self.datafile = ''
        self.molecule = ''
        self.paramfile = ''
        self.precision = 0.0001
        if fp is not None:
            self.read(fp,name)

    def read(self, fp, name='map'):
        self.name = name
        lines = [fp.readline() for _ in range(6)] # first six lines
        self.paramfile = '' #line.split()[1]
        self.datafile = '' #line.split()[1]
        self.molecule = '' #line.split()[1]
        self.spacing = float(lines[3].split()[1])
        self.npts = [int(x) for x in lines[4].split()[1:]]
        self.center = [float(x) for x in lines[5].split()[1:]]
        self.n = [self.npts[i]+1 for i in range(3)]
        self.nelem = self.n[0] * self.n[1] * self.n[2]
        i = 0
        while i < self.nelem:
            val = float(fp.readline())
            self.values.append(val)
            i += 1
        self.origin = [self.center[i] - int(self.npts[i]/2)*self.spacing for i in range(3)]

    def meta(self):
        s= """GRID_PARAMETER_FILE {paramfile}
GRID_DATA_FILE {datafile}
MACROMOLECULE {molecule}
SPACING {spacing:4.3f}
NELEMENTS {npts[0]} {npts[1]} {npts[2]}
CENTER {center[0]:5.3f} {center[1]:5.3f} {center[2]:5.3f}
""".format(**self.__dict__)
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

        fp.write("""object 1 class gridpositions counts {nx} {ny} {nz}
origin {ori[0]:12.5e} {ori[1]:12.5e} {ori[2]:12.5e}
delta {spacing:12.5e} 0 0
delta 0 {spacing:12.5e} 0
delta 0 0 {spacing:12.5e}
object 2 class gridconnections counts {nx} {ny} {nz}
object 3 class array type double rank 0 items {nvals} data follows
""".format(nx=nx, ny=ny, nz=nz, ori=ori, spacing=spacing, nvals=len(vals)))
        col = 0;
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    fp.write(" %12.5E" % vals[k*nx*ny + j*nx + i])
                    col += 1;
                    if col == 3:
                        fp.write("\n")
                        col = 0
        if col != 0:
            fp.write("\n")
        fp.write("""attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
""")
        fp.close()


