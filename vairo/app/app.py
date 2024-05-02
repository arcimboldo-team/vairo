#! /usr/bin/env python3

import io
import os
import pickle
import shutil
import sys
from pathlib import Path

from flask import Flask, render_template, request, jsonify, send_file
from werkzeug.utils import secure_filename
import tempfile

target_directory = os.path.dirname(Path(__file__).absolute().parent.parent)
sys.path.append(target_directory)
from libs import bioutils, features
from tools import utilities

app = Flask(__name__)


def transform_dict(inputDict: dict):
    result = {}
    repeatedKeys = {}
    for key in inputDict.keys():
        for value in inputDict.getlist(key):
            if key in repeatedKeys:
                if isinstance(repeatedKeys[key], list):
                    repeatedKeys[key].append(value)
                else:
                    repeatedKeys[key] = [repeatedKeys[key], value]
            else:
                repeatedKeys[key] = value

    for key, value in repeatedKeys.items():
        if isinstance(value, str) and not value.strip():
            continue
        parts = key.split('-')
        labels = [part for part in parts if not part.isdigit()]
        numbers = [part for part in parts if part.isdigit()]
        current_level = result
        for label, number in zip(labels, numbers):
            if label not in current_level:
                current_level[label] = {}
            if number not in current_level[label]:
                current_level[label][number] = {}
            current_level = current_level[label][number]
        current_level[labels[-1]] = value
    return result


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/features')
def show_modfeatures():
    return render_template('modfeatures.html')


@app.route('/modfeaturesinfo')
def show_modfeaturesinfo():
    return render_template('modfeaturesinfo.html')

@app.route('/output')
def show_output():
    return render_template('output.html')

@app.route('/', methods=["POST"])
def form_vairo():
    input_dict = transform_dict(request.form)
    files_dict = transform_dict(request.files)
    print(input_dict)
    print(files_dict)
    output_path = input_dict.get('output')
    config_file = os.path.join(output_path, 'config.bor')
    files_path = os.path.join(output_path, 'input_files')
    if not os.path.exists(files_path):
        os.makedirs(files_path)
    seq_str = ''
    if 'sequence' in input_dict:
        seq_str = 'sequences:\n'
        for seq_id, seq_info in input_dict['sequence'].items():
            if seq_info.get('input') == 'file':
                file = files_dict['sequence'][seq_id].get('fasta')
                filename = secure_filename(file.filename)
                seq_path = os.path.join(files_path, f'{seq_id}_{filename}.fasta')
                file.save(seq_path)
                seq_str += f'  name: {filename}\n'
            else:
                seq_path = os.path.join(files_path, f'sequence_{seq_id}.fasta')
                fasta_info = seq_info.get('text')
                with open(seq_path, 'w') as f_out:
                    f_out.write(f'>seq\n{fasta_info}\n')
            seq_str += f'- fasta_path: {seq_path}\n'
            copies = seq_info.get('copies')
            positions = seq_info.get('positions')
            mutations = seq_info.get('mutations')
            if copies is not None:
                seq_str += f"  num_of_copies: {copies}\n"
            if positions is not None:
                seq_str += f"  positions: {positions}\n"
            if mutations is not None:
                seq_str += f"  mutations:\n"
                for mutationID, mutationInfo in seq_info['mutations'].items():
                    res = mutationInfo.get('res')
                    pos = mutationInfo.get('pos')
                    if res is not None and pos is not None:
                        seq_str += f"    -'{res}': {pos}\n"

    feat_str = ''
    if 'feature' in input_dict:
        feat_str = 'features:\n'
        for feat_id, feat_info in input_dict['feature'].items():
            file = files_dict['feature'][feat_id].get('pkl')
            filename = secure_filename(sequence.filename)
            pkl_path = os.path.join(files_path, f'{feat_id}_{filename}')
            file.save(pkl_path)
            feat_str += f'- path: {pkl_path}\n'
            feat_str += f"  keep_msa: {'True' if feat_info.get('addmsa') is not None else 'False'}\n"
            feat_str += f"  keep_templates: {'True' if feat_info.get('addtemplates') is not None else 'False'}\n"
            pos = feat_info.get('pos')
            regionfeat = feat_info.get('regionfeat')
            regionquery = feat_info.get('regionquery')
            msa_mask = feat_info.get('mask')
            sequence = files_dict['feature'][feat_id].get('sequence')
            if pos is not None:
                feat_str += f"  positions: {pos}\n"
            if regionfeat is not None:
                feat_str += f"  positions_feature: {regionfeat}\n"
            if regionquery is not None:
                feat_str += f"  positions_query: {regionquery}\n"
            if msa_mask is not None:
                feat_str += f"  msa_mask: {msa_mask}\n"
            if sequence is not None:
                filename = secure_filename(sequence.filename)
                feat_fasta_path = os.path.join(files_path, f'{feat_id}_{filename}')
                sequence.save(feat_fasta_path)
                feat_str += f"  sequence: {feat_fasta_path}\n"

    library_str = ''
    if 'library' in input_dict:
        feat_str = 'append_library:\n'
        for library_id, library_info in input_dict['library'].items():
            lib_path = os.path.join(files_path, f'lib_{library_id}')
            if os.path.exists(lib_path):
                shutil.rmtree(lib_path)
            os.makedirs(lib_path)
            if library_info.get('input') == 'folder':
                lib_folder = files_dict['library'][library_id].get('folder')
                for file in lib_folder:
                    filename = secure_filename(file.filename)
                    file.save(os.path.join(lib_path, filename))
            else:
                file = files_dict['library'][library_id].get('fasta')
                filename = secure_filename(file.filename)
                file.save(os.path.join(lib_path, filename))
            feat_str += f"- path: {lib_path}\n"
            feat_str += f"  add_to_msa: {'True' if library_info.get('addmsa') is not None else 'False'}\n"
            feat_str += f"  add_to_templates: {'True' if library_info.get('addtemplates') is not None else 'False'}\n"
            regionlib = library_info.get('lib')
            regionquery = library_info.get('query')
            if regionlib is not None:
                feat_str += f"  positions_library: {regionlib}\n"
            if regionquery is not None:
                feat_str += f"  positions_query: {regionquery}\n"

    template_str = ''
    if 'template' in input_dict:
        template_str = 'templates:\n'
        for template_id, template_info in input_dict['template'].items():
            if template_info.get('input') == 'code':
                pdb_path = template_info.get('code')
            else:
                file = files_dict['template'][template_id].get('file')
                filename = secure_filename(file.filename)
                pdb_path = os.path.join(files_path, f'{template_id}_{filename}')
                file.save(pdb_path)
            template_str += f"- pdb: {pdb_path}\n"
            template_str += f"  add_to_msa: {'True' if template_info.get('addmsa') is not None else 'False'}\n"
            template_str += f"  add_to_templates: {'True' if template_info.get('addtemplates') is not None else 'False'}\n"
            template_str += f"  generate_multimer: {'True' if template_info.get('multimer') is not None else 'False'}\n"
            template_str += f"  aligned: {'True' if template_info.get('aligned') is not None else 'False'}\n"

            modify = template_info.get('modify')
            if modify is not None:
                template_str += f"  modifications:\n"
                for modify_id, modify_info in template_info['modify'].items():
                    chain = modify_info.get('where')
                    delete = modify_info.get('delete')
                    pos = modify_info.get('pos')
                    modify_residues = True if modify_info.get('residues') is not None else False
                    template_str += f"    - chain: {chain}\n"
                    if delete is not None:
                        template_str += f"      delete_residues: {delete}\n"
                    if pos is not None:
                        template_str += f"      position: {pos}\n"
                    aminos = modify_info.get('amino')
                    if aminos is not None and modify_residues:
                        template_str += f"  replace:\n"
                        for amino_id, amino_info in aminos.items():
                            pos = amino_info.get('pos')
                            select = amino_info.get('select')
                            if pos is not None:
                                template_str += f"      - residues: {pos}\n"
                            if select == 'residue':
                                template_str += f"        by: {select}\n"
                            else:
                                file = files_dict['template'][template_id]['modify'][modify_id]['amino'][amino_id].get('fasta')
                                filename = secure_filename(file.filename)
                                fasta_path = os.path.join(files_path, f'{template_id}_{modify_id}_{amino_id}_{filename}')
                                file.save(fasta_path)
                                template_str += f"        by: {fasta_path}\n"

    with open(config_file, 'w') as f_out:
        f_out.write(f'mode: {input_dict.get("mode")}\n')
        f_out.write(f'output_dir: {output_path}\n')
        f_out.write(f'af2_dbs_path: {input_dict.get("databases")}\n')
        f_out.write(f'run_af2: True\n')
        f_out.write(seq_str)
        f_out.write(feat_str)
        f_out.write(library_str)
        f_out.write(template_str)
    return jsonify({})


@app.route('/generate-multimer', methods=["POST"])
def generate_multimer():
    try:
        pdb_data = request.form.get('templateData')
        results_dict = {}
        with tempfile.NamedTemporaryFile(mode='w+') as pdb_input:
            pdb_input.write(pdb_data)
            pdb_input.flush()
            chain_dict = bioutils.split_pdb_in_chains(pdb_path=pdb_input.name)
            multimer_chain_dict = dict(sorted(bioutils.generate_multimer_chains(pdb_input.name, chain_dict).items()))
            for key, values in multimer_chain_dict.items():
                results_dict[key] = []
                for value in values:
                    results_dict[key].append(bioutils.extract_sequence_msa_from_pdb(value)[key])

            return jsonify(results_dict)
    except Exception as e:
        print(e)
        return jsonify({}), 500


@app.route('/read-pkl', methods=["POST"])
def read_pkl():
    try:
        pkl_file = request.files.get('featuresFile')
        pkl_data = pickle.load(pkl_file.stream)
        length = pkl_data['seq_length'][0]
        return {
            'num_msa': int(pkl_data.get('num_msa', 0)),
            'num_templates': int(pkl_data.get('num_templates', 0)),
            'msa_coverage': pkl_data.get('msa_coverage', [1] * length),
            'templates_coverage': pkl_data.get('templates_coverage', [1] * length)
        }
    except Exception as e:
        print(e)
        return jsonify({}), 500


@app.route('/read-features-info', methods=["POST"])
def read_features():
    try:
        pkl_file = request.files.get('fileFeatures')
        region = request.form.get('rangeFeatures')
        ini_identity = int(request.form.get('iniIdentity'))
        end_identity = int(request.form.get('endIdentity'))
        run_uniprot = True if request.form.get('runUniprot') == 'true' else False
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_file.write(pkl_file.stream.read())
            return_dict = utilities.extract_features_info(temp_file.name, region, ini_identity, end_identity,
                                                          run_uniprot)
            return return_dict
    except Exception as e:
        print(e)
        return jsonify({}), 500


@app.route('/modify-pkl', methods=["POST"])
def modify_pkl():
    try:
        pkl_file = request.files.get('fileFeatures')
        min_identity = float(request.form.get('iniIdentity'))
        max_identity = float(request.form.get('endIdentity'))
        delete_list = request.form.get('deleteList').replace(" ", "").split(',')
        with tempfile.NamedTemporaryFile(delete=False) as features_file:
            features_file.write(pkl_file.stream.read())
            new_feature = features.create_features_from_file(pkl_in_path=features_file.name)

        if new_feature:
            new_feature.delete_by_id(delete_list)
            new_feature.delete_by_range(min_identity, max_identity)
            with tempfile.NamedTemporaryFile(delete=False) as new_file:
                new_feature.write_pkl(new_file.name)
                new_file.flush()
                return send_file(
                    io.BytesIO(new_file.read()),
                    mimetype='image/plain',
                    as_attachment=True,
                    download_name='modified_features.pkl')
    except Exception as e:
        print(e)
        return jsonify({}), 500


if __name__ == '__main__':
    app.json.sort_keys = False
    app.run()
