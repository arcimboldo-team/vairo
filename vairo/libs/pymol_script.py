import logging
import os
import subprocess
import tempfile
from libs import utils


def create_pymol_session(a_air):
    script = """
from pymol import cmd
cmd.set("bg_rgb", "0xffffff")
cmd.set("antialias", '2')
cmd.set("ribbon_sampling", '10')
cmd.set("hash_max", '220')
cmd.set("dash_length", '0.10000')
cmd.set("dash_gap", '0.30000')
cmd.set("cartoon_sampling", '14')
cmd.set("cartoon_loop_quality", '6.00000')
cmd.set("cartoon_rect_length", '1.10000')
cmd.set("cartoon_oval_length", '0.80000')
cmd.set("cartoon_oval_quality", '10.00000')
cmd.set("cartoon_tube_quality", '9.00000')
cmd.set("dash_width", '3.00000')
cmd.set("transparency", '0.60000')
cmd.set("two_sided_lighting", '0')
cmd.set("sculpt_vdw_weight", '0.45000')
cmd.set("sculpt_field_mask", '2047')
cmd.set("ray_shadow", 'off')
cmd.set("auto_color_next", '2')
cmd.set("button_mode_name", '3-Button Viewing')
cmd.set("mouse_selection_mode", '2')
cmd.set("cartoon_nucleic_acid_mode", '2')
cmd.set("cartoon_putty_quality", '11.00000')
cmd.set("cartoon_ring_mode", '1')
cmd.set("cartoon_ladder_color", 'cyan')
cmd.set("cartoon_nucleic_acid_color", 'cyan')
cmd.set("ray_trace_mode", '1')
cmd.set("sculpt_min_weight", '2.25000')
cmd.set("mesh_negative_color", 'grey30')
cmd.set("ray_transparency_oblique_power", '1.00000')
cmd.set("movie_quality", '60')
cmd.set("use_shaders", 'on')
cmd.set("volume_bit_depth", '8')
cmd.set("mesh_as_cylinders", 'on')
cmd.set("line_as_cylinders", 'on')
cmd.set("ribbon_as_cylinders", 'on')
cmd.set("nonbonded_as_cylinders", 'on')
cmd.set("nb_spheres_quality", '3')
cmd.set("alignment_as_cylinders", 'on')
cmd.set("dot_as_spheres", 'on')
cmd.set("valence", 'off')
"""
    i = 1
    if a_air.output.ranked_list:
        pdb_list = [a_air.output.ranked_list[0]]
        pdb_list.extend([experimental for experimental in a_air.output.experimental_list if
                         experimental.name == a_air.output.best_experimental])
        pdb_list.extend([ranked for ranked in a_air.output.ranked_list[1:] if ranked.filtered])
        pdb_list.extend([experimental for experimental in a_air.output.experimental_list if
                         experimental.name != a_air.output.best_experimental])
        pdb_list.extend(a_air.output.templates_list)

        for pdb_in in pdb_list:
            script += f'cmd.load("{pdb_in.split_path}", "{pdb_in.name}")\n'
            if pdb_in in a_air.output.ranked_list:
                script += (f'cmd.spectrum(expression="b", palette="rainbow_rev", selection="{pdb_in.name}", minimum=0, '
                           f'maximum=100)\n')
            elif pdb_in in a_air.output.experimental_list:
                script += f'cmd.color("gray70", "{pdb_in.name}")\n'
            script += f'cmd.show_as("cartoon", "{pdb_in.name}")\n'

            if pdb_in.accepted_interfaces:
                for interface in pdb_in.get_interfaces_with_path():
                    script += f'cmd.load("{interface.path}", "{utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.show_as("sticks", "{utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.show("surface", "{utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.set("transparency", "0.50000")\n'
                    script += f'cmd.color("lime", "resn ALA and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("density", "resn ARG and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("deepsalmon", "resn ASN and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("warmpink", "resn ASP and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("paleyellow", "resn CYS and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("tv_red", "resn GLN and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("ruby", "resn GLU and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("slate", "resn HIS and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("forest", "resn ILE and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("smudge", "resn LEU and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("deepblue", "resn LYS and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("sand", "resn MET and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("gray40", "resn PHE and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("gray20", "resn PRO and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("tv_orange", "resn SER and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("brown", "resn THR and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("palegreen", "resn TRP and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("wheat", "resn TYR and {utils.get_file_name(interface.path)}")\n'
                    script += f'cmd.color("pink", "resn VAL and {utils.get_file_name(interface.path)}")\n'

        if a_air.output.conservation_ranked_path:
            script += f'cmd.load("{a_air.output.conservation_ranked_path}", "{utils.get_file_name(a_air.output.conservation_ranked_path)}")\n'
            script += f'cmd.spectrum("b", "blue_white_red", "{utils.get_file_name(a_air.output.conservation_ranked_path)}", "0", "100")\n'

        
        script += f'cmd.reset()\n'

        if a_air.output.ranked_list[0].get_interfaces_with_path():
            for interface in a_air.output.ranked_list[0].get_interfaces_with_path():
                script += f'cmd.disable("*")\n'
                script += f'cmd.enable("{a_air.output.ranked_list[0].name}")\n'
                script += f'cmd.enable("{utils.get_file_name(interface.path)}")\n'
                script += f'cmd.orient("{a_air.output.ranked_list[0].name}")\n'
                script += f'cmd.zoom("{utils.get_file_name(interface.path)}")\n'
                script += (f'cmd.scene(key="{i}: {interface.name} interface", action="store", message="Zoom '
                           f'into interface {interface.name}")\n')
                i += 1

        script += f'cmd.disable("*")\n'
        script += f'cmd.enable("{a_air.output.ranked_list[0].name}")\n'
        script += f'cmd.enable("{a_air.output.best_experimental}")\n'
        script += f'cmd.orient("{a_air.output.ranked_list[0].name}")\n'
        script += (f'cmd.scene(key="{i}: Reset best prediction (and experimental)", action="store", message="Reset best '
                f'prediction (and experimental)")\n')

        script += f'cmd.disable("*")\n'
        script += f'cmd.enable("{a_air.output.ranked_list[0].name}")\n'
        script += f'cmd.enable("{a_air.output.best_experimental}")\n'
        script += f'cmd.orient("{a_air.output.ranked_list[0].name}")\n'

        for zoom in a_air.pymol_show_list:
            script += f'cmd.zoom(("{zoom}"), 10, complete=1)\n'
            script += f'cmd.scene(key="{i}: Residues {zoom}", action="store", message="Zoom into residues {zoom}")\n'
            i += 1

        script += f'cmd.save("{a_air.output.pymol_session_path}")\n'
        script += 'cmd.quit()\n'

        try:
            with tempfile.TemporaryDirectory() as tmpdirname:
                pymol_script = os.path.join(tmpdirname, 'script_pymol.py')
                with open(pymol_script, 'w+') as f_out:
                    f_out.write(script)
                cmd = 'which pymol'
                subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                cmd = f'pymol -ckq {pymol_script}'
                out, err = subprocess.Popen(cmd, shell=True, env={}, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
                                            stderr=subprocess.STDOUT).communicate()
        except Exception as e:
            logging.error('Error creating a PyMOL session. PyMOL might not be in the path. Skipping.')
            pass
