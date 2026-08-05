#! /usr/bin/env python3
import os
import sys
import csv
import copy
import glob
import shutil
import logging
import argparse
from datetime import datetime
import yaml
from libs import utils, bioutils, analysis_shift_vairo
from tools import utilities

SHIFT_ONLY_KEYS = ('models', 'analysis', 'template', 'use_product')

FORCED_BASE_SETTINGS = {
    'run_af2': True,
    'mode': 'guided',
}

DEFAULT_TEMPLATE_FLAGS = {
    'add_to_msa': True,
    'add_to_templates': True,
    'aligned': False,
    'generate_multimer': False
}

CONFIG_KEY_ORDER = ('output_dir', 'run_af2', 'mode', 'af2_dbs_path', 'sequences')

def _order_config(cfg: dict) -> dict:
    ordered = {k: cfg[k] for k in CONFIG_KEY_ORDER if k in cfg}
    for k in cfg:
        if k not in ordered and k != 'templates':
            ordered[k] = cfg[k]
    if 'templates' in cfg:
        ordered['templates'] = cfg['templates']
    return ordered


def _aggregate_shift_summaries(shift_summaries, out_csv: str):
    rows = []
    for model_name, path in shift_summaries:
        try:
            with open(path) as fh:
                reader = csv.reader(fh)
                next(reader, None)
                for row in reader:
                    rows.append([model_name] + row)
        except Exception:
            continue
    if rows:
        with open(out_csv, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(["model", "combo_name", "sequence"])
            writer.writerows(rows)
        logging.error(f'Wrote aggregated shift summary: {out_csv} ({len(rows)} models)')
    return out_csv


def _build_combined_fasta(fasta_out: str, chain_to_fasta: dict) -> str:
    with open(fasta_out, 'w') as fh:
        for chain, fasta in chain_to_fasta.items():
            sequence = bioutils.extract_sequence(fasta)
            fh.write(f'>{chain}\n{sequence}\n')
    return fasta_out


def _link(src: str, dst: str) -> str:
    src = os.path.abspath(src)
    if os.path.islink(dst) or os.path.exists(dst):
        try:
            os.remove(dst)
        except OSError:
            pass
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)
    return dst


def _setup_run(bor_path: str, gen_dir: str, run_dir: str):
    tag = os.path.splitext(os.path.basename(bor_path))[0][len('model_'):]
    os.makedirs(run_dir, exist_ok=True)
    run_output = os.path.join(run_dir, 'output')

    model_pdb = os.path.join(gen_dir, f'model_{tag}.pdb')
    model_fasta = os.path.join(gen_dir, f'model_{tag}.fasta')
    link_pdb = (_link(model_pdb, os.path.join(run_dir, f'{tag}.pdb'))
                if os.path.exists(model_pdb) else model_pdb)
    if os.path.exists(model_fasta):
        _link(model_fasta, os.path.join(run_dir, f'{tag}.fasta'))

    with open(bor_path) as fh:
        cfg = yaml.load(fh, Loader=yaml.SafeLoader) or {}
    cfg['output_dir'] = run_output
    if cfg.get('templates'):
        cfg['templates'][0]['pdb'] = link_pdb
    run_config = os.path.join(run_dir, f'{tag}.bor')
    with open(run_config, 'w') as fh:
        yaml.dump(cfg, fh, default_flow_style=False, sort_keys=False)
    return run_config, run_output, tag


def _parse_model(model: dict, model_idx: int):
    if not isinstance(model, dict):
        raise Exception(f"models[{model_idx}] must be a mapping with 'pdb' and 'modifications'")
    pdb_path = model.get('pdb')
    if not pdb_path or not os.path.exists(pdb_path):
        raise Exception(f"models[{model_idx}]: 'pdb' is missing or the file does not exist: {pdb_path}")
    entries = model.get('modifications')
    if not entries:
        raise Exception(
            f"models[{model_idx}]: 'modifications' must list at least one "
            f"{{chain, fasta, shift}} entry")

    chain_to_fasta = {}
    chain_steps_dict = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise Exception(
                f"models[{model_idx}]: each 'modifications' entry must be a mapping with "
                f"chain/fasta/shift, got: {entry!r}")
        chain = entry.get('chain')
        if chain is None:
            raise Exception(f"models[{model_idx}]: a modifications entry is missing 'chain'")
        chain = str(chain)
        fasta = entry.get('fasta')
        if not fasta or not os.path.exists(fasta):
            raise Exception(f"models[{model_idx}] chain {chain}: fasta not found: {fasta}")
        shift = entry.get('shift', [0, 0])
        if isinstance(shift, (list, tuple)):
            if len(shift) != 2:
                raise Exception(f"models[{model_idx}] chain {chain}: shift must be [min, max]")
            shift_min, shift_max = sorted((int(shift[0]), int(shift[1])))
        else:
            shift_min = shift_max = int(shift)
        chain_to_fasta[chain] = fasta
        chain_steps_dict[chain] = [shift_min, shift_max]

    return pdb_path, chain_to_fasta, chain_steps_dict


def _auto_split_resids(input_load):
    glycines = int(input_load.get('glycines', 50))
    lengths = []
    for seq in (input_load.get('sequences') or []):
        fasta = seq.get('fasta_path')
        if not fasta or not os.path.exists(fasta):
            return None
        length = len(bioutils.extract_sequence(fasta))
        lengths.extend([length] * int(seq.get('num_of_copies', 1)))
    if len(lengths) < 2:
        return None
    resids, pos = [], 0
    for length in lengths[:-1]:
        pos += length
        resids.append(pos) 
        pos += glycines 
    return resids


def _query_sequence(input_load):
    glycines = int(input_load.get('glycines', 50))
    parts = []
    for s in (input_load.get('sequences') or []):
        fasta = s.get('fasta_path')
        if not fasta:
            continue
        seq = bioutils.extract_sequence(fasta)
        parts.extend([seq] * int(s.get('num_of_copies', 1)))
    return ('G' * glycines).join(parts)


def _chain_query_offsets(chain_to_fasta, query_seq):
    offsets, seen = {}, {}
    for chain, fasta in chain_to_fasta.items():
        seq = bioutils.extract_sequence(fasta)
        k = seen.get(seq, 0)
        pos, start = -1, 0
        for _ in range(k + 1):
            pos = query_seq.find(seq, start) if query_seq else -1
            if pos < 0:
                break
            start = pos + 1
        offsets[chain] = pos
        seen[seq] = k + 1
    return offsets


def _chain_query_positions(chain_to_fasta, input_load):
    expanded = []
    for s in (input_load.get('sequences') or []):
        fp = s.get('fasta_path')
        seq = bioutils.extract_sequence(fp) if fp else None
        expanded.extend([seq] * int(s.get('num_of_copies', 1)))
    positions, used = {}, set()
    for chain, fasta in chain_to_fasta.items():
        seq = bioutils.extract_sequence(fasta)
        pos = None
        for idx, exp_seq in enumerate(expanded):
            if idx not in used and exp_seq == seq:
                pos = idx + 1
                used.add(idx)
                break
        positions[chain] = pos
    return positions


def _parse_use_product(raw):
    if isinstance(raw, (list, tuple)):
        if raw and all(isinstance(g, (list, tuple)) for g in raw):
            groups = [list(g) for g in raw]
        else:
            groups = [list(raw)]                # a single flat group of chain ids
        return True, groups
    return bool(raw), None


def _run_analysis(input_load, workdir, input_path=None):
    # The analysis always runs; the optional 'analysis:' block only adds configuration
    # (references, surfaces, split residues, filter thresholds, ...).
    acfg = input_load.get('analysis') or {}
    split_resids = acfg.get('split_residues')
    if split_resids is None:
        split_resids = _auto_split_resids(input_load)
        if split_resids:
            logging.error(f'analysis: auto split residues from sequences = {split_resids}')
    elif not isinstance(split_resids, (list, tuple)):
        split_resids = [split_resids]
    split_resids = [int(r) for r in split_resids] if split_resids else None

    logging.error('Running post-VAIRO shift analysis...')
    analysis_shift_vairo.run_analysis(
        workdir=workdir,
        reference=acfg.get('reference'),
        split_resids=split_resids,
        ranks=tuple(acfg.get('ranks', (0, 1, 2, 3, 4))),
        cluster_cutoff=float(acfg.get('cluster_cutoff', 5.0)),
        positive_refs=acfg.get('similar_to'),
        negative_refs=acfg.get('different_from'),
        surfaces=acfg.get('surfaces'),
        distance_filters=acfg.get('distance_filters'),
        disulfides=acfg.get('disulfides'),
        assembly_template=acfg.get('assembly_template'),
        config_yml=input_path,          # show the CURRENT input .yml in the report (not a stale copy)
        **{k: acfg[k] for k in acfg
           if k.startswith(('min_', 'max_'))})


def main():
    try:
        utils.create_logger()
        logging.error('')
        logging.error('SHIFT VAIRO')
        logging.error('--------------')
        logging.error('')
        parser = argparse.ArgumentParser(
            description='Generate shifted VAIRO models (by threading + trimming templates) and launch VAIRO on each.')
        parser.add_argument('input_path', nargs='?', help='Path to the input file (.yml)')
        parser.add_argument('-check', '--check', action='store_true',
                            help='Validate the input file and exit')
        parser.add_argument('--generate-only', action='store_true',
                            help='Generate the shifted models and VAIRO configs but do NOT launch VAIRO')
        parser.add_argument('--analyse-only', action='store_true',
                            help='Skip generation/launch; run the results analysis on an existing output dir')
        args = parser.parse_args()
        if not args.input_path:
            raise SystemExit

        input_path = os.path.abspath(args.input_path)
        logging.error(f'Timestamp: {datetime.now()}')
        if not os.path.exists(input_path):
            raise Exception('The given configuration file path does not exist or you do not have read permissions')
        logging.error(f'Reading the SHIFT VAIRO configuration file at {input_path}')
        with open(input_path) as f:
            input_load = yaml.load(f, Loader=yaml.SafeLoader)
        if not isinstance(input_load, dict):
            raise Exception('The input file did not parse into a mapping')

        models = input_load.get('models')
        if not models:
            raise Exception("The input file must contain a non-empty 'models' list")

        base_config = copy.deepcopy(input_load)
        for key in SHIFT_ONLY_KEYS:
            base_config.pop(key, None)
        base_config.update(FORCED_BASE_SETTINGS)
        overrides = dict(input_load.get('template') or {})
        template_flags = dict(DEFAULT_TEMPLATE_FLAGS)
        template_flags.update(overrides)
        base_config['templates'] = [dict(template_flags, pdb='PLACEHOLDER')]
        base_config = _order_config(base_config)
        utils.check_input(base_config)

        if args.check:
            for model_idx, model in enumerate(models):
                _parse_model(model, model_idx) 
            logging.error('Input file is CORRECT')
            raise SystemExit

        root = os.path.abspath(input_load.get('output_dir') or os.path.join(os.path.dirname(input_path), 'output'))
        if args.analyse_only:
            _run_analysis(input_load, root, input_path)
            return
        config_dir = os.path.join(root, 'config')
        models_dir = os.path.join(root, 'models')
        runs_dir = os.path.join(root, 'runs')
        for d in (root, config_dir, models_dir, runs_dir):
            utils.create_dir(d)

        shutil.copy2(input_path, os.path.join(root, 'shift_input.yml'))

        base_config_path = os.path.join(config_dir, 'base_config.yml')
        with open(base_config_path, 'w') as f:
            yaml.dump(base_config, f, default_flow_style=False, sort_keys=False)
        try:
            query_seq = _query_sequence(input_load)
        except Exception:
            query_seq = ''

        single_model = len(models) == 1
        run_specs = []
        shift_summaries = []
        for model_idx, model in enumerate(models):
            pdb_path, chain_to_fasta, chain_steps_dict = _parse_model(model, model_idx)
            model_name = utils.get_file_name(pdb_path)
            model_tag = f'model_{model_idx}_{model_name}'
            gen_dir = os.path.join(models_dir, model_tag)
            utils.create_dir(gen_dir)
            model_base_config = copy.deepcopy(base_config)
            template_path = os.path.abspath(pdb_path)
            combined_fasta = _build_combined_fasta(os.path.join(config_dir, f'combined_{model_tag}.fasta'), chain_to_fasta)
            use_product, chain_groups = _parse_use_product(
                model.get('use_product', input_load.get('use_product', True)))
            one_chain = bool(model.get('one_chain', False))
            chain_query_offsets = None
            if one_chain:
                chain_query_offsets = _chain_query_offsets(chain_to_fasta, query_seq)
                for chain, off in chain_query_offsets.items():
                    if off < 0:
                        logging.error(f'  chain {chain}: fasta not found in the query sequence; '
                                      f'cannot place it by query position.')
            else:
                positions = _chain_query_positions(chain_to_fasta, input_load)
                mods = [{'chain': c, 'position': p} for c, p in positions.items() if p]
                if mods:
                    model_base_config['templates'][0]['modifications'] = mods
                for c, p in positions.items():
                    if not p:
                        logging.error(f'  chain {c}: fasta not among the query sequences; '
                                      f'no template position assigned.')
            model_base_path = os.path.join(config_dir, f'base_config_{model_tag}.yml')
            with open(model_base_path, 'w') as f:
                yaml.dump(model_base_config, f, default_flow_style=False, sort_keys=False)

            mode = f'grouped {chain_groups}' if chain_groups else ('product' if use_product else 'lineal')

            logging.error(f'Model {model_idx} ({model_name}): modifying chains {chain_steps_dict} '
                          f'(one_chain={one_chain}, combos={mode})')

            model_configs = utilities.generate_shift_models(
                fasta_path=combined_fasta,
                template_path=template_path,
                chain_steps_dict=chain_steps_dict,
                bor_file=model_base_path,
                use_product=use_product,
                chain_groups=chain_groups,
                one_chain=one_chain,
                chain_query_offsets=chain_query_offsets,
                output_dir=gen_dir) or []


            logging.error(f'Model {model_idx} ({model_name}): generated {len(model_configs)} VAIRO configs')
            for bor in model_configs:
                tag = os.path.splitext(os.path.basename(bor))[0][len('model_'):]
                run_id = tag if single_model else f'{model_tag}_{tag}'
                run_config, run_output, tag = _setup_run(bor, gen_dir, os.path.join(runs_dir, run_id))
                run_specs.append((run_config, run_output, tag))
            summary_csv = os.path.join(gen_dir, 'shift_summary.csv')
            if os.path.exists(summary_csv):
                shift_summaries.append((model_name, summary_csv))

        if not run_specs:
            raise Exception('No VAIRO configs were generated. Check the generate_shift_models output above.')

        _aggregate_shift_summaries(shift_summaries, os.path.join(root, 'shift_summary.csv'))

        if args.generate_only:
            logging.error(f'--generate-only: {len(run_specs)} runs prepared under {runs_dir} (not launching VAIRO):')
            for run_config, _, _ in run_specs:
                logging.error(f'  {run_config}')
            logging.error('SHIFT VAIRO finished (generate-only).')
            logging.error(f'Timestamp: {datetime.now()}')
            return

        logging.error(f'Launching VAIRO on {len(run_specs)} runs...')
        launched = skipped = 0
        for run_config, run_output, tag in run_specs:
            try:
                already_done = os.path.isdir(run_output) and bool(utils.read_rankeds(input_path=run_output))
            except Exception:
                already_done = False
            if already_done:
                logging.error(f'  -> SKIP {tag} (already finished: rankeds in {run_output})')
                skipped += 1
                continue
            logging.error(f'  -> VAIRO {run_config}')
            try:
                bioutils.run_vairo(yml_path=run_config, input_path=run_output)
                launched += 1
            except Exception:
                logging.error(f'VAIRO run failed or produced no rankeds for {tag}; continuing.', exc_info=True)
                continue
        logging.error(f'VAIRO launches: {launched} run, {skipped} skipped (already finished).')

        _run_analysis(input_load, root, input_path)

        logging.error(f'SHIFT VAIRO finished. {len(run_specs)} runs processed.')
        logging.error(f'Timestamp: {datetime.now()}')

    except SystemExit:
        sys.exit()
    except Exception:
        logging.error(f'Timestamp: {datetime.now()}')
        logging.error('ERROR:', exc_info=True)


if __name__ == "__main__":
    main()
