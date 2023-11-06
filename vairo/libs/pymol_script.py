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
cmd.set("surface_quality", '4')
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
cmd.set("max_threads", '4')
cmd.set("button_mode_name", '3-Button Viewing')
cmd.set("mouse_selection_mode", '2')
cmd.set("cartoon_nucleic_acid_mode", '2')
cmd.set("cartoon_putty_quality", '11.00000')
cmd.set("cartoon_ring_mode", '1')
cmd.set("cartoon_ladder_color", 'cyan')
cmd.set("cartoon_nucleic_acid_color", 'cyan')
cmd.set("ray_trace_mode", '1')
cmd.set("sculpt_min_weight", '2.25000')
cmd.set("surface_negative_color", 'grey50')
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
"""
    i = 1
    if a_air.output.ranked_list:
        pdb_list = [ranked for ranked in a_air.output.ranked_list if ranked.filtered]
        pdb_list.extend(a_air.output.experimental_list)
        pdb_list.extend(a_air.output.templates_list)
        for j, pdb_in in enumerate(pdb_list):
            script += f'cmd.load("{pdb_in.split_path}", "{pdb_in.name}")\n'
            if j != 0:
                script += f'cmd.disable("{pdb_in.name}")\n'
            for interface in pdb_in.interfaces:
                script += f'cmd.load("{interface.path}", "{utils.get_file_name(interface.path)}")\n'
                script += f'cmd.disable("{utils.get_file_name(interface.path)}")\n'

        for zoom in a_air.pymol_show_list:
            key = f'F{i}'
            script += f'cmd.zoom("center", {zoom})\n'
            script += f'cmd.scene(key="{key}", action="store", message="Zoom into residues {zoom}")\n'
            script += f'cmd.scene(key="{key}", action="rename", new_key="{key}: User zoom")\n'
            i += 1

    script += f'cmd.reset()\n'
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
