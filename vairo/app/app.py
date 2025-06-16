#! /usr/bin/env python3

import io
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import re
from pathlib import Path

from flask import Flask, render_template, request, jsonify, send_file, session
from werkzeug.utils import secure_filename

target_directory = os.path.dirname(Path(__file__).absolute().parent.parent)
sys.path.append(target_directory)
from vairo import libs
from libs import bioutils, features, utils, alphafold_classes
from tools import utilities

app = Flask(__name__)
app.secret_key = 'secret'

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
def show_index():
    return render_template('index.html', active_page="index", sub_page=None)

@app.route('/parameterization')
def show_parameterization():
    return render_template('parameterization.html', active_page="run", sub_page="parameterization")

@app.route('/input')
def show_input():
    return render_template('input.html', active_page="run", sub_page="input")

@app.route('/output')
def show_output():
    session['lastmodified'] = 0.0
    return render_template('output.html', active_page="run", sub_page="output")

@app.route('/features')
def show_modfeatures():
    return render_template('modfeatures.html', active_page="tools", sub_page=None)

@app.route('/modfeaturesinfo')
def show_modfeaturesinfo():
    return render_template('modfeaturesinfo.html', active_page="tools", sub_page=None)

@app.route('/check-output', methods=['POST'])
def check_output():
    folder = request.form.get('folder')
    html_path = os.path.join(folder, 'output', 'output.html')
    yml_path = os.path.join(folder, 'config.yml')
    session['html_path'] = html_path
    session['yml_path'] = yml_path
    if folder and os.path.isdir(folder):
        html_exists = True if os.path.exists(html_path) else False
        yml_exists = True if os.path.exists(yml_path) else False
        return jsonify({"status": "success", "yml_exists": yml_exists, "html_exists": html_exists})
    return jsonify({"status": "success", "yml_exists": False, "html_exists": False})

@app.route('/check-databases', methods=['POST'])
def check_databases():
    folder = request.form.get('folder')
    try:
        af2_db = alphafold_classes.AlphaFoldPaths(af2_dbs_path=folder)
        exist = af2_db.validate_db_paths()
        if exist:
            return jsonify({"status": "success", "exists": True})
    except Exception as e:
        print(e)
        pass
    return jsonify({"status": "error", "exists": False})

@app.route('/read-yml', methods=["GET"])
def read_yml():
    config_file = os.path.join(session['yml_path'])
    if os.path.exists(config_file):
        read_data = open(config_file, 'r').read()
        return jsonify({'status': 'error', "data": read_data})
    else:
        return jsonify({'status': 'error', "data": None})


@app.route('/run-vairo', methods=["POST"])
def run_vairo():
    run_vairo_path = os.path.join(target_directory, 'run_vairo.py')
    param_dict = transform_dict(request.form)
    write_yml(output_path=session['yml_path'], output_str=param_dict.get('text'))
    #subprocess.Popen(["nq", run_vairo_path, session['yml_path']])
    return jsonify({'status': 'success'})


@app.route('/form-vairo', methods=["POST"])
def form_vairo():
    try:
        save_session(request.form)
        param_dict = transform_dict(request.form)
        files_dict = transform_dict(request.files)
        output_path = param_dict.get('output')
        files_path = os.path.join(output_path, 'param_files')
        if not os.path.exists(files_path):
            os.makedirs(files_path)
    
        runaf2 = True if param_dict.get('runaf2') is not None else False
        config_str = f'mode: {param_dict.get("mode")}\n'
        config_str += f'output_dir: {output_path}\n'
        config_str += f'af2_dbs_path: {param_dict.get("databases")}\n'
        config_str += f"run_af2: {runaf2}\n"

        if 'sequence' in param_dict:
            config_str += 'sequences:\n'
            for seq_id, seq_info in param_dict['sequence'].items():
                if seq_info.get('input') == 'file':
                    file = files_dict['sequence'][seq_id].get('fasta')
                    filename = secure_filename(file.filename)
                    seq_path = os.path.join(files_path, f'{seq_id}_{filename}')
                    file.save(seq_path)
                    config_str += f'  name: {filename}\n'
                else:
                    seq_path = os.path.join(files_path, f'sequence_{seq_id}.fasta')
                    fasta_info = seq_info.get('text')
                    with open(seq_path, 'w') as f_out:
                        f_out.write(f'>seq\n{fasta_info}\n')
                config_str += f'- fasta_path: {seq_path}\n'
                copies = seq_info.get('copies')
                positions = seq_info.get('positions')
                mutations = seq_info.get('mutations')
                if copies is not None:
                    config_str += f"  num_of_copies: {copies}\n"
                if positions is not None:
                    config_str += f"  positions: {positions}\n"
                if mutations is not None:
                    config_str += f"  mutations:\n"
                    for mutation_id, mutation_info in seq_info['mutations'].items():
                        res = mutation_info.get('res')
                        pos = mutation_info.get('pos')
                        if res is not None and pos is not None:
                            config_str += f"    -'{res}': {pos}\n"

        if param_dict.get("mode") == 'guided':
            if 'template' in param_dict:
                config_str += 'templates:\n'
                for template_id, template_info in param_dict['template'].items():
                    if template_info.get('radio') == 'code':
                        pdb_path = template_info.get('code')
                    else:
                        file = files_dict['template'][template_id].get('file')
                        filename = secure_filename(file.filename)
                        pdb_path = os.path.join(files_path, f'{template_id}_{filename}')
                        file.save(pdb_path)
                    config_str += f"- pdb: {pdb_path}\n"
                    config_str += f"  add_to_msa: {'True' if template_info.get('addmsa') is not None else 'False'}\n"
                    config_str += f"  add_to_templates: {'True' if template_info.get('addtemplates') is not None else 'False'}\n"
                    config_str += f"  generate_multimer: {'True' if template_info.get('multimer') is not None else 'False'}\n"
                    config_str += f"  aligned: {'True' if template_info.get('aligned') is not None else 'False'}\n"

                    modify = template_info.get('modify')
                    if modify is not None:
                        config_str += f"  modifications:\n"
                        for modify_id, modify_info in template_info['modify'].items():
                            chain = modify_info.get('where')
                            delete = modify_info.get('delete')
                            pos = modify_info.get('pos')
                            modify_residues = True if modify_info.get('residues') is not None else False
                            config_str += f"    - chain: {chain}\n"
                            if delete is not None:
                                config_str += f"      delete_residues: {delete}\n"
                            if pos is not None:
                                config_str += f"      position: {pos}\n"
                            aminos = modify_info.get('amino')
                            if aminos is not None and modify_residues:
                                config_str += f"  mutations:\n"
                                for amino_id, amino_info in aminos.items():
                                    pos = amino_info.get('pos')
                                    select = amino_info.get('select')
                                    if pos is not None:
                                        config_str += f"      - numbering_residues: {pos}\n"
                                    if select == 'residue':
                                        config_str += f"        mutate_with: {select}\n"
                                    else:
                                        file = files_dict['template'][template_id]['modify'][modify_id]['amino'][
                                            amino_id].get('fasta')
                                        filename = secure_filename(file.filename)
                                        fasta_path = os.path.join(files_path,
                                                                  f'{template_id}_{modify_id}_{amino_id}_{filename}')
                                        file.save(fasta_path)
                                        config_str += f"        mutate_with: {fasta_path}\n"


            if 'feature' in param_dict:
                config_str += 'features:\n'
                for feat_id, feat_info in param_dict['feature'].items():
                    file = files_dict['feature'][feat_id].get('pkl')
                    filename = secure_filename(file.filename)
                    pkl_path = os.path.join(files_path, f'{feat_id}_{filename}')
                    file.save(pkl_path)
                    config_str += f'- path: {pkl_path}\n'
                    config_str += f"  keep_msa: {'True' if feat_info.get('addmsa') is not None else 'False'}\n"
                    config_str += f"  keep_templates: {'True' if feat_info.get('addtemplates') is not None else 'False'}\n"
                    pos = feat_info.get('pos')
                    regionfeat = feat_info.get('regionfeat')
                    regionquery = feat_info.get('regionquery')
                    msa_mask = feat_info.get('mask')
                    sequence = files_dict['feature'][feat_id].get('sequence')
                    if pos is not None:
                        config_str += f"  positions: {pos}\n"
                    if regionfeat is not None:
                        config_str += f"  numbering_features: {regionfeat}\n"
                    if regionquery is not None:
                        config_str += f"  numbering_query: {regionquery}\n"
                    if msa_mask is not None:
                        config_str += f"  msa_mask: {msa_mask}\n"
                    if sequence is not None:
                        filename = secure_filename(sequence.filename)
                        feat_fasta_path = os.path.join(files_path, f'{feat_id}_{filename}')
                        sequence.save(feat_fasta_path)
                        config_str += f"  sequence: {feat_fasta_path}\n"

            if 'library' in param_dict:
                config_str += 'append_library:\n'
                for library_id, library_info in param_dict['library'].items():
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
                    config_str += f"- path: {lib_path}\n"
                    config_str += f"  add_to_msa: {'True' if library_info.get('addmsa') is not None else 'False'}\n"
                    config_str += f"  add_to_templates: {'True' if library_info.get('addtemplates') is not None else 'False'}\n"
                    regionlib = library_info.get('lib')
                    regionquery = library_info.get('query')
                    if regionlib is not None:
                        config_str += f"  numbering_library: {regionlib}\n"
                    if regionquery is not None:
                        config_str += f"  numbering_query: {regionquery}\n"

        write_yml(output_path=session["yml_path"], output_str=config_str)
        return jsonify({'status': 'success'})
    except Exception as e:
        print(e)
        return jsonify({'status': 'error'}), 500

def write_yml(output_path: str, output_str: str):
    with open(output_path, 'w') as f_out:
        f_out.write(output_str)

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
        return jsonify({'status': 'error'}), 500

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

@app.route('/check-update', methods=["GET"])
def check_update():
    if not os.path.exists(session['html_path']):
        return jsonify({'changed': False, 'error': 'File not found'})
    last_modified = os.path.getmtime(session['html_path'])
    changed = last_modified > session['lastmodified']
    session['lastmodified'] = last_modified
    if changed:
        with open(session['html_path'], 'r', encoding='utf-8') as f:
            html_content = f.read()
        body_content = re.sub(r'<footer\b[^>]*>.*?</footer>', '', html_content, flags=re.DOTALL | re.IGNORECASE)
    else:
        body_content = None
    return jsonify({
        'changed': changed,
        'content': body_content
    })

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
            return jsonify(return_dict)
    except Exception as e:
        print(e)
        return jsonify({'status': 'error'}), 500

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
        return jsonify({'status': 'error'}), 500
    except Exception as e:
        print(e)
        return jsonify({'status': 'error'}), 500


@app.route('/save-form-data', methods=['POST'])
def save_form_data():
    try:
        save_session(request.form)
        return jsonify({'status': 'success'})
    except Exception as e:
        print(e)
        return jsonify({'status': 'error'}), 500

def save_session(input_dict: dict):
    session['form_data'] = input_dict

@app.route('/load-form-data', methods=['GET'])
def load_form_data():
    form_data = session.get('form_data', {})
    return jsonify({'status': 'success', 'data': form_data})

if __name__ == '__main__':
    app.json.sort_keys = False
    app.run()
