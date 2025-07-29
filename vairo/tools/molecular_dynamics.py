#!/usr/bin/env python3

import argparse
import shlex
import shutil
import subprocess
import os
import sys
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import matplotlib.ticker as ticker
from Bio.PDB import PDBParser, PDBIO
import numpy as np

current_directory = os.path.dirname(os.path.abspath(__file__))
target_directory = os.path.abspath(os.path.join(current_directory, '..', '..'))
sys.path.append(target_directory)
from vairo.libs import features, bioutils, plots, utils, structures, template_modifications, global_variables

def generate_rmds_plots():
    def process_and_plot_results(reference_results, plot_file, plot_title):
        colors = ['blue', 'orange', 'm', 'g']
        plt.figure(figsize=(12, 6))
        for i, (model, pdbs) in enumerate(reference_results.items()):
            sorted_items = sorted(pdbs.items(), key=lambda x: os.path.basename(x[0]))
            y_vals = [float(rmsd) for _, rmsd in sorted_items]
            x_vals = list(range(1, len(y_vals) + 1))
            plt.plot(x_vals, y_vals, marker='o', linestyle='-', color=colors[i], label=model)

        ax = plt.gca()
        ax.set_xlabel('Model Number')
        ax.set_ylabel('RMSD')
        ax.set_title(plot_title)
        ax.grid(True)
        ax.set_ylim(0, 10)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_locator(ticker.NullLocator())

        plt.legend()
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')

    def split_models_in_pdb(pdb_path):
        name = utils.get_file_name(pdb_path)
        output_dir = os.path.join(os.getcwd(), name)
        utils.create_dir(name, True)
        structure = bioutils.get_structure(pdb_path)
        io = PDBIO()
        for model in structure:
            io.set_structure(model)
            out_file = os.path.join(output_dir, f"model_{model.id:03d}.pdb")
            io.save(out_file)
        return output_dir

    def extract_residue_ranges(res_list):
        results_list = []
        for seq, start_resnum, chain_id in res_list:
            end_resnum = start_resnum + len(seq) - 1
            results_list.append((chain_id, start_resnum, end_resnum))
        return results_list

    def prepare_lsqkab_input(reference_ranges, target_ranges):
        lines = []
        for reference, target in zip(reference_ranges, target_ranges):
            ref_chain, ref_start, ref_end = reference
            target_chain, target_start, target_end = target
            lines.append(f'FIT RESIDUE CA {ref_start} to {ref_end} CHAIN {ref_chain}')
            lines.append(f'MATCH {target_start} {target_end} CHAIN {target_chain}')
        return lines

    def write_lsqkab_input_file(lsqkab_input_lines, pdb1, pdb2):
        with open('lsqkab_input.cards', 'w') as f:
            f.write('lsqkab ')
            f.write(f'xyzinf {pdb1} ')
            f.write(f'xyzinm {pdb2} ')
            f.write('DELTAS deltas.out ')
            f.write('xyzout output.pdb << END-lsqkab\n')
            f.write('output deltas \n')
            f.write('output XYZ \n')
            for line in lsqkab_input_lines:
                f.write(line + '\n')
            f.write(f'end \n')
            f.write(f'END-lsqkab')
        stdout, stderr = subprocess.Popen(['bash', 'lsqkab_input.cards'], stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE).communicate()
        text = stdout.decode('utf-8', errors='ignore')
        lines = text.splitlines()
        pattern = r"RMS\s+XYZ\s+DISPLACEMENT\s*=\s*([\d\.]+)"
        for line in lines:
            m = re.search(pattern, line)
            if m:
                rmsd = m.group(1)
                return rmsd
        return None

    def rmsd_ca(atoms1, atoms2):
        coords1 = np.array([a.get_coord() for a in atoms1])
        coords2 = np.array([a.get_coord() for a in atoms2])
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff * diff, axis=1)))

    def calculate_translation_list(reference_numbering, reference_pdb, target_pdb):
        reference_struct = bioutils.get_structure(reference_pdb)
        target_struct = bioutils.get_structure(target_pdb)
        reference_sequence = bioutils.extract_sequence_msa_from_pdb(reference_pdb)
        target_sequence = bioutils.extract_sequence_msa_from_pdb(target_pdb)
        translation_list = []

        for reference_chain, query_start, query_end in reference_numbering:
            query_ca = [reference_struct[0][reference_chain][i]['CA'] for i in range(query_start - 1, query_end)]
            query_seq = reference_sequence[reference_chain][query_start - 1:query_end]
            best_hit = None
            for target_chain, target_chain_seq in target_sequence.items():
                target_start = target_chain_seq.find(query_seq)
                while target_start != -1:
                    target_end = target_start + len(query_seq)
                    target_ca = [
                        target_struct[0][target_chain][res_id]['CA']
                        for res_id in range(target_start + 1, target_end + 1)
                    ]
                    score = rmsd_ca(query_ca, target_ca)
                    if best_hit is None or score < best_hit[3]:
                        best_hit = (target_chain, target_start + 1, target_end, score)
                    target_start = target_chain_seq.find(query_seq, target_start + 1)
            if best_hit is not None:
                translation_list.append((best_hit[0], best_hit[1], best_hit[2]))
            else:
                sys.exit(1)
        return translation_list

    reference = 'Crystal'
    input_dict = {
        'greengreen': {
            'title': 'RMSD Interfaces; D1-D1 % Dimer Reference',
            'references': {
                'Crystal': 'GreenTetramerCrystal_clean.pdb',
                'VAIRO': 'greeboundDimer_super.pdb'
            },
            'frames': {
                'VAIRO': [('Crystal D1-D1 Dimer trajectory', 'frames_GGxtal_200.pdb'),
                          ('VAIRO D1-D1 Dimer trajectory', 'frames_GGvairo_200.pdb')]
            },
            'superpose_list': [('YNGKTYTANLKAD', 84, 'A'),
                               ('YNGKTYTANLKAD', 84, 'B'),
                               ('DVSFNFGSEN', 128, 'A'),
                               ('DVSFNFGSEN', 128, 'B'),
                               ('LDQNGVASLTN', 173, 'A'),
                               ('LDQNGVASLTN', 173, 'B'),
                               ]
        },
        'blueblue': {
            'title': 'RMSD Interfaces; D2-D2 % Dimer Reference',
            'references': {
                'Crystal': 'BlueTetramerCrystal.pdb',
                'VAIRO': 'blueboundDimer_clean_cut4md_renumbered_super.pdb'
            },
            'frames': {
                'VAIRO': [('Crystal D2-D2 Dimer trajectory', 'frames_BBxtal_200.pdb'),
                          ('VAIRO D2-D2 Dimer trajectory', 'frames_BBvairo_200.pdb')]
            },
            'superpose_list': [('NVNFYDVTSGATVTNG', 199, 'B'),
                               ('NVNFYDVTSGATVTNG', 199, 'A'),
                               ('AAQYADKKLNTRTANT', 241, 'B'),
                               ('AAQYADKKLNTRTANT', 241, 'A'),
                               ]
        },
        'greenblue': {
            'title': 'RMSD Interfaces; D1-D2 % Dimer Reference',
            'references': {
                'Crystal': 'BlueTetramerCrystal.pdb',
                'VAIRO': 'greenbluefromgbdMaskBlueNoNaiveT2_cut4md_clean_renumbered_super.pdb'
            },
            'frames': {
                'VAIRO': [('Crystal D1-D2 Dimer trajectory', 'frames_GBxtal1_200.pdb'),
                          ('VAIRO D1-D2 Dimer trajectory', 'frames_GBvairo1_200.pdb')]
            },
            'superpose_list': [('SAVAANTANNTPAIAGNL', 59, 'A'),
                               ('YAINTTDNSN', 190, 'A'),
                               ('SVNADNQGQVNVANVVAAINSKYF', 217, 'D'),
                               ('LKDQKIDVNSVGYFKAPHTFTV', 264, 'D'),
                               ]
        },
        'greentetramer': {
            'title': 'RMSD Interfaces; D1-D1/D1-D2 % Tetramer Reference',
            'references': {
                'Crystal': 'GreenTetramerCrystal_clean.pdb',
                'VAIRO': 'greentetramerTilefromPiecesgreenSeqMSAnomaskR0_clean_cut4md_renumbered_super.pdb'
            },
            'frames': {
                'VAIRO': [('Crystal D1-D1/D1-D2 Tetramer trajectory', 'frames_GTxtal_200.pdb'),
                          ('VAIRO D1-D1/D1-D2 Tetramer trajectory', 'frames_GTvairo_200.pdb')]
            },
            'superpose_list': [('SAVAANTANNTPAIAGNL', 59, 'A'),
                               ('SAVAANTANNTPAIAGNL', 59, 'E'),
                               ('YAINTTDNSN', 190, 'A'),
                               ('YAINTTDNSN', 190, 'E'),
                               ('SVNADNQGQVNVANVVAAINSKYF', 217, 'B'),
                               ('SVNADNQGQVNVANVVAAINSKYF', 217, 'D'),
                               ('LKDQKIDVNSVGYFKAPHTFTV', 264, 'B'),
                               ('LKDQKIDVNSVGYFKAPHTFTV', 264, 'D'),
                               ('YNGKTYTANLKAD', 84, 'A'),
                               ('YNGKTYTANLKAD', 84, 'B'),
                               ('YNGKTYTANLKAD', 84, 'D'),
                               ('YNGKTYTANLKAD', 84, 'E'),
                               ('DVSFNFGSEN', 128, 'A'),
                               ('DVSFNFGSEN', 128, 'B'),
                               ('DVSFNFGSEN', 128, 'D'),
                               ('DVSFNFGSEN', 128, 'E'),
                               ('LDQNGVASLTN', 173, 'A'),
                               ('LDQNGVASLTN', 173, 'B'),
                               ('LDQNGVASLTN', 173, 'D'),
                               ('LDQNGVASLTN', 173, 'E'),
                               ]
        },
        'bluetetramer': {
            'title': 'RMSD Interfaces; D2-D2/D1-D2 % Tetramer Reference',
            'references': {
                'Crystal': 'BlueTetramerCrystal.pdb',
                'VAIRO': 'blueboundTetramerfromPiecesSeqgbTR1_clean_cut4md_renumbered_super.pdb'
            },
            'frames': {
                'VAIRO': [('Crystal D2-D2/D1-D2 Tetramer trajectory', 'frames_BTxtal_200.pdb'),
                          ('VAIRO D2-D2/D1-D2 Tetramer trajectory', 'frames_BTvairo_200.pdb')],
            },
            'superpose_list': [('SAVAANTANNTPAIAGNL', 59, 'E'),
                               ('SAVAANTANNTPAIAGNL', 59, 'A'),
                               ('YAINTTDNSN', 190, 'E'),
                               ('YAINTTDNSN', 190, 'A'),
                               ('SVNADNQGQVNVANVVAAINSKYF', 217, 'B'),
                               ('SVNADNQGQVNVANVVAAINSKYF', 217, 'D'),
                               ('LKDQKIDVNSVGYFKAPHTFTV', 264, 'B'),
                               ('LKDQKIDVNSVGYFKAPHTFTV', 264, 'D'),
                               ('NVNFYDVTSGATVTNG', 199, 'B'),
                               ('NVNFYDVTSGATVTNG', 199, 'D'),
                               ('NVNFYDVTSGATVTNG', 199, 'A'),
                               ('NVNFYDVTSGATVTNG', 199, 'E'),
                               ('AAQYADKKLNTRTANT', 241, 'B'),
                               ('AAQYADKKLNTRTANT', 241, 'D'),
                               ('AAQYADKKLNTRTANT', 241, 'A'),
                               ('AAQYADKKLNTRTANT', 241, 'E'),
                               ]
        },
    }
    generate_list = ['bluetetramer', 'greentetramer', 'greenblue', 'greengreen', 'blueblue']
    # generate_list = ['greenblue']
    for generate in generate_list:
        old_path = os.getcwd()
        path = os.path.join(os.getcwd(), generate)
        os.chdir(path)
        superpose_list = input_dict[generate]['superpose_list']
        references = input_dict[generate]['references']
        frames = input_dict[generate]['frames']

        reference_superpose_list = extract_residue_ranges(superpose_list)
        results_dict = {ref: {} for ref in references.keys()}
        superpose_translation_dict = {reference: reference_superpose_list}
        superpose_translation_dict.update({
            key: calculate_translation_list(reference_superpose_list, references[reference], ref)
            for key, ref in references.items()
            if key != reference
        })

        for frame_type, frame_values in frames.items():
            for frame_name, frame_value in frame_values:
                dir_path = split_models_in_pdb(frame_value)
                for reference_type, reference_values in references.items():
                    results_dict[reference_type][frame_name] = {}
                    lsqkab_input = prepare_lsqkab_input(superpose_translation_dict[frame_type],
                                                        superpose_translation_dict[reference_type])
                    for pdb in os.listdir(dir_path):
                        pdb_path = os.path.join(dir_path, pdb)
                        rmsd = write_lsqkab_input_file(lsqkab_input, reference_values, pdb_path)
                        results_dict[reference_type][frame_name][pdb_path] = rmsd

        for ref in references:
            title = input_dict[generate]['title'].replace('%', ref)
            file = re.sub(r'[^\w\-]', '_', title)
            process_and_plot_results(results_dict[ref], file, title)
        os.chdir(old_path)

def postprocess_run(input_path: str, suffix: str):
    def gmx_run(cmd, stdin=""):
        result = subprocess.run(
            shlex.split(cmd),
            input=stdin.encode() if stdin else None,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode, cmd, result.stdout, result.stderr
            )

    def energy_jobs():
        return [
            ("em1.edr", f"potential_em1_{suffix}.xvg", "10 0\n"),
            ("em2.edr", f"potential_em2_{suffix}.xvg", "10 0\n"),
            ("nvt.edr", f"temperature_{suffix}.xvg", "15 0\n"),
            ("npt.edr", f"pressure_{suffix}.xvg", "17 0\n"),
            ("npt.edr", f"density_{suffix}.xvg", "23 0\n"),
        ]

    def return_path(name):
        return os.path.join(input_path, name)

    output_path = os.path.join(os.getcwd(), 'processed_run')
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)
    os.chdir(output_path)

    for edr, xvg, sel in energy_jobs():
        gmx_run(f"gmx energy -f {return_path(edr)} -o {xvg}", stdin=sel)

    gmx_run(
        f"gmx trjconv -s {return_path('em1.tpr')} -f {return_path('md.xtc')} "
        f"-o whole_{suffix}.xtc -pbc whole",
        stdin="1\n",
    )
    gmx_run(
        f"gmx trjconv -s {return_path('em1.tpr')} -f whole_{suffix}.xtc "
        f"-o cluster_{suffix}.xtc -pbc cluster -center",
        stdin="1\n1\n1\n",
    )
    gmx_run(
        f"gmx rms -s {return_path('md.tpr')} -f cluster_{suffix}.xtc "
        f"-o rmsd_{suffix}.xvg -tu ns", stdin="4\n4\n"
    )
    gmx_run(
        f"gmx rms -s {return_path('em1.tpr')} -f cluster_{suffix}.xtc "
        f"-o rmsd_origin_{suffix}.xvg -tu ns", stdin="4\n4\n"
    )
    gmx_run(
        f"gmx gyrate -s {return_path('md.tpr')} -f cluster_{suffix}.xtc "
        f"-o gyrate_{suffix}.xvg", stdin="1\n"
    )
    gmx_run(
        f"gmx trjconv -f cluster_{suffix}.xtc -b 0 -e 100000 "
        f"-o traj_1000_{suffix}.xtc -skip 10"
    )
    gmx_run(
        f"gmx trjconv -f traj_1000_{suffix}.xtc -s {return_path('md.tpr')} "
        f"-o frames_{suffix}.pdb -skip 5",
        stdin="1\n",
    )

def generate_pca_plots():
    print('test')

def preprocess_run(input_pdb: str, mdp_folder: str):
    def run_cmd(command, input_text=None):
        print(f"Running: {' '.join(command)}")
        result = subprocess.run(
            command,
            input=input_text,
            text=True,
            capture_output=True
        )
        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            raise RuntimeError(f"Command failed: {' '.join(command)}")


    base_name = utils.get_file_name(input_pdb)
    required_files = ["em1.mdp", "em2.mdp", "md.mdp", "npt.mdp", "nvt.mdp"]
    output_path = os.path.join(os.getcwd(), f'{base_name}_prepared_gmx_files')
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)
    os.chdir(output_path)
    shutil.copy2(input_pdb, output_path)
    shutil.copytree(mdp_folder, output_path, dirs_exist_ok=True)
    missing_files = [f for f in required_files if not os.path.exists(os.path.join(output_path, f))]
    if missing_files:
        print('Missing .mdp files')
        sys.exit(1)

    # Step 1: pdb2gmx
    run_cmd(
        ["gmx", "pdb2gmx", "-f", input_pdb, "-o", f"{base_name}_processed.gro", "-water", "spce"],
        input_text="15\n"
    )
    print("done define force field and create topology")

    # Step 2: editconf
    run_cmd(
        ["gmx", "editconf", "-f", f"{base_name}_processed.gro", "-o", f"{base_name}_box.gro", "-c", "-d", "2.0", "-bt",
         "dodecahedron"])
    print("done box generation")

    # Step 3: solvate
    run_cmd(["gmx", "solvate", "-cp", f"{base_name}_box.gro", "-cs", "spc216.gro", "-o", f"{base_name}_solv.gro", "-p",
             "topol.top"])
    print("done solvating")

    # Step 4: grompp (ions)
    run_cmd(["gmx", "grompp", "-f", "em1.mdp", "-c", f"{base_name}_solv.gro", "-p", "topol.top", "-o", "ions.tpr"])
    print("done generate ions.tpr")

    # Step 5: genion
    run_cmd(["gmx", "genion", "-s", "ions.tpr", "-o", f"{base_name}_solv_ions.gro", "-p", "topol.top", "-pname", "NA",
             "-nname", "CL", "-neutral"])
    print("done")

    # Step 6: Energy minimization (em1)
    run_cmd(["gmx", "grompp", "-f", "em1.mdp", "-c", f"{base_name}_solv_ions.gro", "-p", "topol.top", "-o", "em1.tpr"])
    run_cmd(["gmx", "mdrun", "-v", "-deffnm", "em1"])

    # Step 7: Energy minimization (em2)
    run_cmd(["gmx", "grompp", "-f", "em2.mdp", "-c", "em1.gro", "-p", "topol.top", "-o", "em2.tpr"])
    run_cmd(["gmx", "mdrun", "-v", "-deffnm", "em2"])

    # Step 8: NVT equilibration
    run_cmd(["gmx", "grompp", "-f", "nvt.mdp", "-c", "em2.gro", "-r", "em2.gro", "-p", "topol.top", "-o", "nvt.tpr"])
    run_cmd(["gmx", "mdrun", "-deffnm", "nvt"])

    # Step 9: NPT equilibration
    run_cmd(
        ["gmx", "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro", "-t", "nvt.cpt", "-p", "topol.top", "-o",
         "npt.tpr"])
    run_cmd(["gmx", "mdrun", "-deffnm", "npt"])

    # Step 10: Production MD
    run_cmd(
        ["gmx", "grompp", "-f", "md.mdp", "-c", "npt.gro", "-r", "npt.gro", "-t", "npt.cpt", "-p", "topol.top", "-o",
         "md.tpr"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="GROMACS automation.")
    subparsers = ap.add_subparsers(dest="command", required=True)

    parser_rmsd = subparsers.add_parser("rmsd", help="Generate RMSD plots")
    parser_pca = subparsers.add_parser("pca", help="Generate PCA plots")
    parser_pre = subparsers.add_parser("pre", help="Prepare files for GROMACS")
    parser_pre.add_argument("--input_pdb", required=True, help="Pdb file (.pdb)")
    parser_pre.add_argument("--mdp_folder", required=True, help=" MDPs folder (.mdps)")
    parser_post = subparsers.add_parser("post", help="Analyze the trajectory")
    parser_post.add_argument("--input_path", required=True, help="Trajectory file (.xtc)")
    parser_post.add_argument("--suffix", required=True, help="Suffix name")



    # Parse args
    args = ap.parse_args()

    # Command routing
    if args.command == "post":
        postprocess_run(args.input_path, args.suffix)
    elif args.command == 'pre':
        preprocess_run(args.input_pdb, args.mdp_folder)
    elif args.command == "rmsd":
        generate_rmds_plots()
    elif args.command == 'pca':
        generate_pca_plots()

