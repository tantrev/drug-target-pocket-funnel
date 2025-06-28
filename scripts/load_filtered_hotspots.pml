import os
from pymol import cmd

base_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
hotspot_path = os.path.join(base_dir, "outputs", "silcs_hotspots_report.AF-P29460.astex_library", "hotspots_sites.pdb")

cmd.load(hotspot_path, "hotspots_all")

cmd.select("siteA_sel", "hotspots_all and index 120+5")
cmd.select("siteB_sel", "hotspots_all and index 99+9+61+24")

cmd.create("siteA_hotspots", "siteA_sel")
cmd.create("siteB_hotspots", "siteB_sel")

cmd.delete("siteA_sel")
cmd.delete("siteB_sel")

cmd.disable("hotspots_all")
cmd.hide("everything", "hotspots_all")

cmd.show("spheres", "siteA_hotspots")
cmd.set("sphere_scale", 0.3, "siteA_hotspots")
cmd.color("cyan", "siteA_hotspots")

cmd.show("spheres", "siteB_hotspots")
cmd.set("sphere_scale", 0.3, "siteB_hotspots")
cmd.color("magenta", "siteB_hotspots")

cmd.zoom("siteA_hotspots or siteB_hotspots")
