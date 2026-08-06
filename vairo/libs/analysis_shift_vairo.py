import os
import re
import io
import copy
import glob
import base64
import shutil
import logging
import itertools

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from jinja2 import Environment, FileSystemLoader
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from Bio.PDB import PDBParser, PDBIO, Superimposer, Chain, Model, Structure

from libs import bioutils


METRIC_DIRECTION = {
    "rmsd_similar": True,
    "rmsd_different": False,
    "surface_accuracy": False,
    "plddt": False,
    "rmsd_template": True,
    "rmsd_ref": True,
    "interface_area": False,
    "interface_energy": True,
}

METRIC_ORDER = ["rmsd_similar", "rmsd_different", "rmsd_template", "surface_accuracy",
                "plddt", "rmsd_ref", "interface_area", "interface_energy"]


def _lower_is_better(col):
    return METRIC_DIRECTION.get(col, True)


def _read_structure(path):
    return PDBParser(QUIET=True).get_structure("pdb", path)


def _write_structure(structure, path):
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(path)


def _ca_map(structure):
    """{(chain_id, resnum): CA atom} for every standard residue with a CA."""
    return {(ch.id, res.id[1]): res["CA"]
            for ch in structure[0] for res in ch
            if res.id[0] == " " and "CA" in res}


def _ca_atoms(path):
    return _ca_map(_read_structure(path))


def _chain_ca_maps(structure):
    """{chain_id: {resnum: CA atom}} for every chain that has at least one CA."""
    out = {}
    for ch in structure[0]:
        cas = {res.id[1]: res["CA"] for res in ch if res.id[0] == " " and "CA" in res}
        if cas:
            out[ch.id] = cas
    return out


def _residue(structure, chain_id, resnum):
    for ch in structure[0]:
        if ch.id == chain_id:
            for res in ch:
                if res.id[1] == int(resnum):
                    return res
    return None


def _chain_count(pdb):
    return sum(1 for _ in _read_structure(pdb)[0])


def _min_atom_distance(res_a, res_b, atom_names=None):
    """Closest approach between two residues. With ``atom_names=(a, b)`` only those
    named atoms are compared (e.g. ("SG", "SG") for a disulfide); otherwise every
    non-hydrogen atom is. NaN when either residue/atom is missing."""
    if res_a is None or res_b is None:
        return float("nan")

    def _select(res, name):
        if name:
            return [res[name]] if name in res else []
        return [a for a in res if a.element != "H"]

    atoms_a = _select(res_a, atom_names[0] if atom_names else None)
    atoms_b = _select(res_b, atom_names[1] if atom_names else None)
    if not atoms_a or not atoms_b:
        return float("nan")
    coords_a = np.array([a.get_coord() for a in atoms_a])
    coords_b = np.array([a.get_coord() for a in atoms_b])
    return float(np.sqrt(((coords_a[:, None, :] - coords_b[None, :, :]) ** 2).sum(-1)).min())



# A shift tag is chain-letter(s) + signed integer, joined by "_": e.g. "A2_B-1".
_TAG_TOKEN_RE = re.compile(r"^([A-Za-z]+)(-?\d+)$")


def parse_tag(tag):
    """'A2_B-1' -> {'A': 2, 'B': -1}. Unparsable tokens are ignored ({} if none match)."""
    out = {}
    for token in str(tag).split("_"):
        match = _TAG_TOKEN_RE.match(token)
        if match:
            out[match.group(1)] = int(match.group(2))
    return out


def _as_list(x):
    """None/'' -> []; a single string -> [string]; anything else -> list(x)."""
    if not x:
        return []
    return [x] if isinstance(x, str) else list(x)


def _ref_column(kind, path):
    """Column name for an RMSD-to-reference, e.g. ('similar', '/p/6R17.pdb') -> 'rmsd_similar_6R17'."""
    base = re.sub(r"[^0-9A-Za-z]+", "_", os.path.splitext(os.path.basename(str(path)))[0])
    return f"rmsd_{kind}_{base}"



def _run_dirs(workdir):
    """Yield (run_dir, tag, template_pdb) for each VAIRO run under workdir/runs/*."""
    for run_dir in sorted(glob.glob(os.path.join(workdir, "runs", "*"))):
        if not os.path.isdir(run_dir):
            continue
        templates = glob.glob(os.path.join(run_dir, "*.pdb"))
        template = templates[0] if templates else None
        tag = (os.path.splitext(os.path.basename(template))[0] if template
               else os.path.basename(run_dir))
        if parse_tag(tag):
            yield run_dir, tag, template


def collect_models(workdir, ranks=(0, 1, 2, 3, 4), dest=None):
    """One row per (run, rank) prediction: model name, tag, rank, pdb path, template,
    and a shift_<chain> column per shifted chain. Copies pdbs into `dest` if given."""
    rows = []
    for run_dir, tag, template in _run_dirs(workdir):
        shifts = parse_tag(tag)
        out_dir = os.path.join(run_dir, "output")
        for rank in ranks:
            src = os.path.join(out_dir, f"ranked_{rank}.pdb")
            if not os.path.exists(src):
                continue
            if dest:
                os.makedirs(dest, exist_ok=True)
                dst = os.path.join(dest, f"ranked_{rank}_{tag}.pdb")
                shutil.copy2(src, dst)
                src = dst
            row = {"model": os.path.basename(src), "tag": tag, "rank": rank,
                   "pdb": src, "template": template}
            row.update({f"shift_{chain}": shift for chain, shift in shifts.items()})
            rows.append(row)

    df = pd.DataFrame(rows)
    logging.error(f"collect_models: {len(df)} models, "
                  f"{df['tag'].nunique() if len(df) else 0} shift combinations")
    return df


def collect_runs(workdir):
    """One dict per run: {'tag', 'html'} where html is the run's output.html (or None)."""
    runs = []
    for run_dir, tag, _ in _run_dirs(workdir):
        html = os.path.join(run_dir, "output", "output.html")
        runs.append({"tag": tag, "html": html if os.path.exists(html) else None})
    logging.error(f"collect_runs: {len(runs)} runs, "
                  f"{sum(1 for r in runs if r['html'])} with an output.html")
    return runs



def superpose(fixed_pdb, moving_pdb, output_pdb=None):
    """Global least-squares (Kabsch) CA RMSD of moving_pdb onto fixed_pdb over their
    common (chain, resnum) atoms. Optionally writes the superposed moving structure."""
    fixed_struct = _read_structure(fixed_pdb)
    moving_struct = _read_structure(moving_pdb)
    fixed, moving = _ca_map(fixed_struct), _ca_map(moving_struct)
    common = sorted(set(fixed) & set(moving))
    if not common:
        logging.error(f"  no common CA atoms: {os.path.basename(fixed_pdb)} vs "
                      f"{os.path.basename(moving_pdb)}")
        return float("nan")
    sup = Superimposer()
    sup.set_atoms([fixed[k] for k in common], [moving[k] for k in common])
    if output_pdb:
        sup.apply(list(moving_struct.get_atoms()))
        _write_structure(moving_struct, output_pdb)
    return round(float(sup.rms), 4)


def add_template_rmsd(df, superposed_dir=None):
    """Global least-squares CA RMSD of each model to its own template."""
    if superposed_dir:
        os.makedirs(superposed_dir, exist_ok=True)
    rmsds = []
    for _, row in df.iterrows():
        out = os.path.join(superposed_dir, row["model"]) if superposed_dir else None
        try:
            rmsd = superpose(row["template"], row["pdb"], output_pdb=out)
        except Exception as exc:
            logging.error(f"  template rmsd failed for {row['model']}: {exc}")
            rmsd = float("nan")
        rmsds.append(rmsd)
    df = df.copy()
    df["rmsd_template"] = rmsds
    logging.error(f"add_template_rmsd: global least-squares CA RMSD to template; "
                  f"median {np.nanmedian(rmsds):.2f} A")
    return df


def add_perchain_template_rmsd(df):
    """Per-chain CA RMSD to the template (one column per chain, plus the worst-chain
    value), fitting each chain on its own so a single bad chain is not averaged away."""
    per_row, chain_ids = [], set()
    for _, row in df.iterrows():
        vals = {}
        template_chains = _chain_ca_maps(_read_structure(row["template"]))
        model_chains = _chain_ca_maps(_read_structure(row["pdb"]))
        for cid in set(template_chains) & set(model_chains):
            common = sorted(set(template_chains[cid]) & set(model_chains[cid]))
            sup = Superimposer()
            sup.set_atoms([template_chains[cid][r] for r in common],
                          [model_chains[cid][r] for r in common])
            vals[cid] = round(sup.rms, 4)
            chain_ids.add(cid)
        per_row.append(vals)

    df = df.copy()
    for cid in sorted(chain_ids):
        df[f"rmsd_template_{cid}"] = [v.get(cid, float("nan")) for v in per_row]
    df["rmsd_template_chain_max"] = [max(v.values()) if v else float("nan") for v in per_row]
    if chain_ids:
        logging.error(f"add_perchain_template_rmsd: per-chain LS RMSD for chains "
                      f"{sorted(chain_ids)}; median worst-chain "
                      f"{np.nanmedian(df['rmsd_template_chain_max']):.2f} A")
    return df


def add_rmsd(df, against, column, superposed_dir=None):
    """CA RMSD of each model against a reference. `against` is either a fixed pdb path
    or the name of a df column holding a per-row reference path."""
    per_row = against in df.columns
    if not per_row and not (against and os.path.exists(str(against))):
        raise ValueError(f"reference not found: {against}")
    if superposed_dir:
        os.makedirs(superposed_dir, exist_ok=True)

    values = []
    for _, row in df.iterrows():
        target = row.get(against) if per_row else against
        if not target or not os.path.exists(str(target)):
            values.append(float("nan"))
            continue
        out = os.path.join(superposed_dir, row["model"]) if superposed_dir else None
        try:
            rmsd = superpose(str(target), row["pdb"], output_pdb=out)
        except Exception as exc:
            logging.error(f"  {column} failed for {row['model']}: {exc}")
            rmsd = float("nan")
        values.append(rmsd)

    df = df.copy()
    df[column] = values
    return df


def add_plddt(df):
    """Mean CA pLDDT (AlphaFold confidence, stored in the B-factor) per model."""
    vals = []
    for _, row in df.iterrows():
        try:
            bfactors = [a.get_bfactor() for a in _read_structure(row["pdb"]).get_atoms()
                        if a.get_name() == "CA"]
            vals.append(round(float(np.mean(bfactors)), 2) if bfactors else float("nan"))
        except Exception as exc:
            logging.error(f"  plddt failed for {row['model']}: {exc}")
            vals.append(float("nan"))
    df = df.copy()
    df["plddt"] = vals
    if np.isfinite(vals).any():
        logging.error(f"add_plddt: median {np.nanmedian(vals):.1f}")
    return df


def add_geometry_checks(df, distance_filters=None, disulfides=None, pdb_col="pdb", chain_copies=None):
    specs = []
    for chain_a, res_a, chain_b, res_b, max_dist in (distance_filters or []):
        chain_a, res_a, chain_b, res_b, max_dist = chain_a, int(res_a), chain_b, int(res_b), float(max_dist)
        specs.append((f"dist_{chain_a}{res_a}_{chain_b}{res_b}",
                      chain_a, res_a, chain_b, res_b, None, max_dist,
                      f"{chain_a}{res_a}–{chain_b}{res_b} within {max_dist:g} Å"))
    for chain_a, res_a, chain_b, res_b in (disulfides or []):
        chain_a, res_a, chain_b, res_b = chain_a, int(res_a), chain_b, int(res_b)
        specs.append((f"ss_{chain_a}{res_a}_{chain_b}{res_b}",
                      chain_a, res_a, chain_b, res_b, ("SG", "SG"), 2.5,
                      f"disulfide {chain_a}{res_a}–{chain_b}{res_b} (SG–SG ≤ 2.5 Å)"))
    if not specs:
        return df, []

    def copy_pairs(chain_a, chain_b):
        a_positions = chain_copies.get(chain_a, [chain_a]) if chain_copies else [chain_a]
        b_positions = chain_copies.get(chain_b, [chain_b]) if chain_copies else [chain_b]
        return list(zip(a_positions, b_positions))

    columns = {spec[0]: [] for spec in specs}
    for _, row in df.iterrows():
        structure = _read_structure(row[pdb_col])
        for col, chain_a, res_a, chain_b, res_b, atoms, _, _ in specs:
            dists = [_min_atom_distance(_residue(structure, pa, res_a),
                                        _residue(structure, pb, res_b), atom_names=atoms)
                     for pa, pb in copy_pairs(chain_a, chain_b)]
            dists = [d for d in dists if np.isfinite(d)]
            columns[col].append(round(max(dists), 2) if dists else float("nan"))

    df = df.copy()
    for col, vals in columns.items():
        df[col] = vals
    filters = [(col, max_allowed, label) for col, *_, max_allowed, label in specs]
    logging.error(f"add_geometry_checks: {[f[0] for f in filters]} on {pdb_col}"
                  f"{' across all copies' if chain_copies else ''}")
    return df, filters


def _parse_surfaces(surfaces):
    """-> {'buried': {(chain, res)}, 'exposed': {(chain, res)}}. Each item is either a
    [chain, res, 'b'|'e'] list (like distance_filters/disulfides) or a compact 'A38e' string."""
    expected = {"buried": set(), "exposed": set()}
    for item in surfaces or []:
        try:
            if isinstance(item, (list, tuple)):
                chain, res, status = str(item[0]), int(item[1]), str(item[2]).lower()
            else:
                chain, status, res = item[0], item[-1].lower(), int(item[1:-1])
            expected["buried" if status == "b" else "exposed"].add((chain, res))
        except (IndexError, ValueError):
            logging.error(f"  bad surface item, skipping: {item!r}")
    return expected


def add_surface(df, surfaces, workdir=None, pdb_col="pdb", chain_copies=None):
    """Fraction of expected buried/exposed residues that come out correctly (SASA-based).
    With `chain_copies`, the expected residues are replicated to every copy of each chain."""
    expected = _parse_surfaces(surfaces)
    if chain_copies:
        expected = {kind: {(cid, res) for (chain, res) in residues
                           for cid in chain_copies.get(chain, [chain])}
                    for kind, residues in expected.items()}
    if not (expected["buried"] or expected["exposed"]):
        return df

    workdir = workdir or os.path.join(os.path.dirname(df.iloc[0]["pdb"]), "surfaces")
    os.makedirs(workdir, exist_ok=True)
    accuracies = []
    for _, row in df.iterrows():
        out = os.path.join(workdir, f"surface_{row['model']}")
        data = bioutils.analyse_surface_residues(row[pdb_col], out, expected)
        accuracies.append(float(data.get("accuracy", float("nan"))))
    df = df.copy()
    df["surface_accuracy"] = accuracies
    logging.error(f"add_surface: {int(np.sum(np.isfinite(accuracies)))}/{len(df)} scored")
    return df


def split_chains(pdb_in, pdb_out, split_resids):
    """Cut a single-chain pdb into chains A, B, ... at the given residue boundaries
    (so PISA can find interfaces). Returns False if fewer than 2 chains result."""
    splits = sorted(int(s) for s in split_resids)
    model = _read_structure(pdb_in)[0]
    new_model = Model.Model(0)
    chains = {}
    for res in list(model.get_residues()):
        idx = sum(1 for s in splits if res.id[1] > s)
        cid = chr(ord("A") + idx)
        if cid not in chains:
            chains[cid] = Chain.Chain(cid)
            new_model.add(chains[cid])
        new_res = copy.copy(res)
        new_res.detach_parent()
        chains[cid].add(new_res)

    if len(chains) < 2:
        return False
    new_structure = Structure.Structure("split")
    new_structure.add(new_model)
    _write_structure(new_structure, pdb_out)
    return True


def add_interface_metrics(df, split_resids, workdir=None, pdb_col="pdb"):
    """Total PISA interface area / binding energy per model. Multi-chain models are
    used directly; a single-chain model is first cut at `split_resids`."""
    workdir = workdir or os.path.join(os.path.dirname(df.iloc[0]["pdb"]), "interfaces")
    os.makedirs(workdir, exist_ok=True)
    areas, energies, counts = [], [], []
    for _, row in df.iterrows():
        interfaces = []
        if _chain_count(row[pdb_col]) >= 2:
            pdb_for_pisa = row[pdb_col]
        elif split_resids and split_chains(row[pdb_col], os.path.join(workdir, row["model"]),
                                           split_resids):
            pdb_for_pisa = os.path.join(workdir, row["model"])   # single chain -> cut to domains
        else:
            pdb_for_pisa = None
            logging.error(f"  {row['model']}: single chain, no usable split, skipping")
        if pdb_for_pisa:
            try:
                interfaces = bioutils.find_interface_from_pisa(pdb_for_pisa, workdir) or []
            except Exception as exc:
                logging.error(f"  PISA failed for {row['model']}: {exc}")

        if interfaces:
            areas.append(round(sum(i.area for i in interfaces), 3))
            energies.append(round(sum(i.deltaG for i in interfaces), 3))
            counts.append(len(interfaces))
        else:
            areas.append(float("nan"))
            energies.append(float("nan"))
            counts.append(0)
    df = df.copy()
    df["interface_area"] = areas
    df["interface_energy"] = energies
    df["n_interfaces"] = counts
    logging.error(f"add_interface_metrics: {int(np.sum(np.isfinite(areas)))}/{len(df)} scored")
    return df


def add_clusters(df, cutoff=5.0):
    """Group models that are structurally the same (average-linkage on pairwise CA RMSD)."""
    df = df.copy()
    paths = list(df["pdb"])
    if len(paths) < 2:
        df["cluster"] = 1 if paths else []
        return df

    atom_maps = [_ca_atoms(p) for p in paths]
    common = sorted(set.intersection(*(set(a) for a in atom_maps)))
    ref = [atom_maps[0][k] for k in common]
    coords = []
    for atoms in atom_maps:
        moving = [atoms[k] for k in common]
        sup = Superimposer()
        sup.set_atoms(ref, moving)
        rot, tran = sup.rotran
        coords.append((np.array([a.coord for a in moving]) @ rot + tran).ravel())

    dist = pdist(np.vstack(coords)) / np.sqrt(len(common))
    labels = fcluster(linkage(dist, method="average"), t=cutoff, criterion="distance")
    df["cluster"] = labels
    logging.error(f"add_clusters: {len(set(labels))} clusters at {cutoff} A "
                  f"({len(common)} common CA atoms)")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Assembly building
#
# Replicate each prediction (an n-chain unit) into a symmetric template made of K
# such units, so packing-dependent metrics can be scored on the full assembly.
# Each copy gets one rigid transform placing the whole prediction onto it.
# ─────────────────────────────────────────────────────────────────────────────

def _place_chain(chain, rotation, translation):
    """A detached deep copy of `chain` with the rigid transform applied to every atom.
    The parent is severed first so deepcopy does not clone the whole parent structure."""
    saved_parent = chain.parent
    chain.parent = None
    placed = copy.deepcopy(chain)
    chain.parent = saved_parent
    for atom in placed.get_atoms():
        atom.set_coord(atom.get_coord() @ rotation + translation)
    return placed


def _superpose_unit(prediction_cas, copy_cas, pairing):
    """Fit the whole prediction onto one template copy in a single rigid move, using the
    CA atoms of the paired chains together. Returns (rms, rotation, translation) or None."""
    fixed, moving = [], []
    for pred_chain, copy_chain in pairing:
        for resnum, ca in prediction_cas[pred_chain].items():
            if resnum in copy_cas[copy_chain]:
                moving.append(ca)
                fixed.append(copy_cas[copy_chain][resnum])
    if len(moving) < 3:
        return None
    superposer = Superimposer()
    superposer.set_atoms(fixed, moving)
    rotation, translation = superposer.rotran
    return superposer.rms, rotation, translation


def _best_pairing(prediction_cas, copy_cas):
    """Which predicted chain maps to which chain of this template copy. Tries every
    same-length assignment and keeps the one with the lowest whole-unit RMSD. This is
    computed per copy AND per model because equal-length chains can fit different
    positions depending on the shift."""
    pred_chains = list(prediction_cas)
    candidates = [[cc for cc in copy_cas if len(copy_cas[cc]) == len(prediction_cas[pc])]
                  for pc in pred_chains]
    best = None
    for choice in itertools.product(*candidates):
        if len(set(choice)) != len(pred_chains):
            continue
        pairing = list(zip(pred_chains, choice))
        fit = _superpose_unit(prediction_cas, copy_cas, pairing)
        if fit and (best is None or fit[0] < best[0]):
            best = (fit[0], pairing)
    return best[1] if best else None


def build_assembly(df, template, outdir, pdb_col="pdb"):
    """Replicate each prediction into the assembly `template`, keeping its structure intact.
    Returns (df with an 'assembly' path column, chain_copies) where chain_copies maps each
    predicted chain to the list of template positions it fills, one per copy."""
    os.makedirs(outdir, exist_ok=True)
    if df.empty:
        return df.assign(assembly=[]), {}

    # Split the template into consecutive n-chain blocks, one per copy of the prediction.
    template_cas = _chain_ca_maps(_read_structure(template))
    template_chain_ids = list(template_cas)
    n = len(_chain_ca_maps(_read_structure(df.iloc[0][pdb_col])))
    copy_cas_list = [{c: template_cas[c] for c in template_chain_ids[i:i + n]}
                     for i in range(0, len(template_chain_ids) - n + 1, n)]

    assembly_paths, chain_copies = [], None
    for _, row in df.iterrows():
        pred_struct = _read_structure(row[pdb_col])
        prediction, prediction_cas = pred_struct[0], _chain_ca_maps(pred_struct)
        assembly = Model.Model(0)
        placed_at = {chain: [] for chain in prediction_cas}
        for copy_cas in copy_cas_list:
            pairing = _best_pairing(prediction_cas, copy_cas)
            if pairing is None:
                continue
            _, rotation, translation = _superpose_unit(prediction_cas, copy_cas, pairing)
            for pred_chain, copy_chain in pairing:
                placed = _place_chain(prediction[pred_chain], rotation, translation)
                placed.id = copy_chain
                assembly.add(placed)
                placed_at[pred_chain].append(copy_chain)
        structure = Structure.Structure("assembly")
        structure.add(assembly)
        out_path = os.path.join(outdir, row["model"])
        _write_structure(structure, out_path)
        assembly_paths.append(out_path)
        chain_copies = chain_copies or placed_at

    df = df.copy()
    df["assembly"] = assembly_paths
    logging.error(f"build_assembly: {len(copy_cas_list)} copies of {n} chains -> "
                  f"{n * len(copy_cas_list)}-chain assemblies; chain copies {chain_copies}")
    return df, chain_copies


# ─────────────────────────────────────────────────────────────────────────────
# Filtering, scoring and summaries
# ─────────────────────────────────────────────────────────────────────────────

def apply_filter(df, col, thr, kind, label):
    """Keep rows where col <= thr (kind='max') or col >= thr (kind='min'); NaNs always
    pass. Returns (filtered_df, funnel_stage). If nothing passes, keeps everything so the
    report is never empty. A no-op when the threshold is absent/infinite/zero."""
    def _stage(frame):
        return {"stage": label, "models": len(frame), "tags": int(frame["tag"].nunique())}

    active = thr is not None and col in df.columns and (
        (kind == "max" and np.isfinite(thr)) or (kind == "min" and thr > 0))
    if not active:
        return df, _stage(df)

    n_before = len(df)
    values = df[col].astype(float)
    keep = values.isna() | (values <= thr if kind == "max" else values >= thr)
    kept = df[keep].copy()
    logging.error(f"apply_filter[{col}]: kept {len(kept)}/{n_before} "
                  f"({'<=' if kind == 'max' else '>='} {thr:g}); {n_before - len(kept)} discarded")
    if kept.empty:
        logging.error(f"apply_filter[{col}]: nothing passed; keeping all so the report is not "
                      f"empty -- relax the {col} cutoff if unexpected")
        return df, _stage(df)
    return kept, _stage(kept)


def score_models(df, **thresholds):
    """Flag which models pass optional min_<col>/max_<col> thresholds (e.g.
    max_rmsd_template=5, min_surface_accuracy=80) in a `passed` column. Nothing
    is dropped; no aggregate score is invented."""
    df = df.copy()
    keep = pd.Series(True, index=df.index)
    for name, thr in thresholds.items():
        if thr is None:
            continue
        kind, _, col = name.partition("_")
        if col not in df.columns or kind not in ("min", "max"):
            continue
        values = df[col].astype(float)
        keep &= values.notna() & ((values >= thr) if kind == "min" else (values <= thr))
    df["passed"] = keep
    logging.error(f"score_models: {int(keep.sum())}/{len(df)} models pass the thresholds")
    return df


def choose_focus(df):
    """Pick the metric the ranking should focus on, by priority:
      1. a "should resemble" reference    -> rmsd_similar   (lower is better)
      2. only a "should differ" reference -> rmsd_different (higher is better)
      3. otherwise                        -> rmsd_template  (global RMSD to the template, lower better)
    Returns (column, ascending, human_description). surface_accuracy is intentionally NOT a
    candidate: the expected surface is a filter (min_surface_accuracy discards models below
    it), not a selection. rmsd_template is both a ranking metric and a filter.
    """
    def has(col):
        return col in df.columns and df[col].notna().any()

    for col, desc in (("rmsd_similar", "similarity to the reference the prediction should resemble"),
                      ("rmsd_different", "difference from the reference the prediction should not resemble"),
                      ("rmsd_template", "global least-squares RMSD to the template")):
        if has(col):
            return col, _lower_is_better(col), desc
    return None, True, "no ranking metric available"


def summarise_shifts(df, focus_col=None, ascending=True):
    """One row per shift combination - the 'which shift is best' table, sorted by
    the focus metric (best value of that shift)."""
    if df.empty:
        return pd.DataFrame()
    cols = [c for c in METRIC_ORDER if c in df.columns and df[c].notna().any()]
    grouped = df.groupby("tag")
    out = grouped.agg(models=("model", "count"))
    for col in cols:
        out[f"{col}_mean"] = grouped[col].mean()
        out[f"{col}_best"] = grouped[col].min() if _lower_is_better(col) else grouped[col].max()
    for extra, how, name in (("passed", "sum", "passed"), ("cluster", "nunique", "clusters")):
        if extra in df.columns:
            out[name] = grouped[extra].agg(how)
    out = out.reset_index()
    key = f"{focus_col}_best" if focus_col else None
    if key and key in out.columns:
        out = out.sort_values(key, ascending=ascending)
    return out.reset_index(drop=True)


def summarise_chains(df, metric=None):
    """One row per (chain, shift): the mean and best `metric` over the models that
    share that shift, so you can see which shift each chain prefers."""
    if not metric:
        return pd.DataFrame()
    best = min if _lower_is_better(metric) else max
    rows = []
    for col in df.columns:
        if not col.startswith("shift_"):
            continue
        chain = col[len("shift_"):]
        for shift, group in df[[col, metric]].dropna().groupby(col)[metric]:
            rows.append({"chain": chain, "shift": int(shift), "metric": metric,
                         "models": len(group), "metric_mean": group.mean(),
                         "metric_best": best(group)})
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Report artefacts: alignment GIF and HTML report
# ─────────────────────────────────────────────────────────────────────────────

def build_shift_gif(runs, outdir, filename="shift_alignment.gif", duration=150, plot_index=0):
    """Stitch each run's alignment plot (extracted from its output.html) into one GIF,
    ordered by shift and captioned with the shift tag."""
    png_re = re.compile(r'base64,([A-Za-z0-9+/=]+)"')

    def _plot(html):
        matches = png_re.findall(open(html).read())
        if len(matches) <= plot_index:
            return None
        return Image.open(io.BytesIO(base64.b64decode(matches[plot_index]))).convert("RGB")

    def _key(run):
        shifts = parse_tag(run["tag"])
        return tuple(shifts[c] for c in sorted(shifts)), run["tag"]

    frames = []
    for run in sorted((r for r in (runs or []) if r.get("html")), key=_key):
        try:
            im = _plot(run["html"])
        except Exception:
            im = None
        if im is not None:
            frames.append((run["tag"], im))
    if not frames:
        logging.error("build_shift_gif: no alignment plots found; skipping")
        return None

    band = 46
    font = ImageFont.truetype("DejaVuSans-Bold.ttf", 26)
    width = max(im.width for _, im in frames)
    height = max(im.height for _, im in frames)

    canvases = []
    for tag, im in frames:
        canvas = Image.new("RGB", (width, height + band), "white")
        canvas.paste(im, ((width - im.width) // 2, band + (height - im.height) // 2))
        ImageDraw.Draw(canvas).text((14, band // 2), f"shift {tag}", fill="black",
                                    font=font, anchor="lm")
        canvases.append(canvas)

    path = os.path.join(outdir, filename)
    canvases[0].save(path, save_all=True, append_images=canvases[1:],
                     duration=duration, loop=0, optimize=True, disposal=2)
    logging.error(f"build_shift_gif: {len(canvases)} frames -> {path}")
    return path


def _template_name(path):
    """Human name for a template from its superposed path: the run-dir name with the
    leading 'model_<n>_' stripped, falling back to the pdb basename."""
    real = os.path.realpath(str(path))
    return re.sub(r"^model_\d+_", "", os.path.basename(os.path.dirname(real))) \
        or os.path.splitext(os.path.basename(real))[0]


def _ref_name(ref):
    return os.path.splitext(os.path.basename(str(ref)))[0]


def _read_vairo_assets(tmpl_dir):
    """Pull the <style> and <script> blocks out of VAIRO's output.html so the shift
    report can reuse the same styling. Returns ('', '') if unavailable."""
    try:
        src = open(os.path.join(tmpl_dir, "output.html")).read()
    except Exception:
        return "", ""
    css = "".join(re.findall(r"<style type=.text/css.>.*?</style>", src, re.S))
    js = "".join(re.findall(r"<script.*?</script>", src, re.S))
    return css, js


def write_html_report(df, shifts, chains, outdir, reference="", thresholds=None,
                      runs=None, comparisons=None, focus=None, config_yml=None,
                      collected=None, funnel=None):
    """Render the shift-analysis HTML report from the scored tables. `data` is the single
    context dict handed to the Jinja template; each block below fills one section of it."""
    focus_col, focus_desc = focus or (None, None)
    focus_lower = _lower_is_better(focus_col) if focus_col else True
    output_html = os.path.join(outdir, "shift_analysis.html")
    tmpl_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "templates")
    css, js = _read_vairo_assets(tmpl_dir)

    def fmt(v):
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return ""
        return f"{v:.2f}" if isinstance(v, float) else str(v)

    def cellfmt(v):
        return fmt(float(v)) if v is not None and np.isfinite(v) else ""

    best = shifts.iloc[0] if len(shifts) else None

    template_names = {_template_name(p) for p in df["template"].dropna().unique()}
    template_group_name = next(iter(template_names)) if len(template_names) == 1 else "Template"

    # --- per-model table columns (friendly labels; some cells pair two metrics) ---
    model_cols, col_labels, col_paired = [], {}, {}

    def _add_col(key, label, pair=None):
        if key in df.columns and df[key].notna().any():
            model_cols.append(key)
            col_labels[key] = label
            if pair and pair in df.columns:
                col_paired[key] = pair

    for col, kind, ref in (comparisons or []):
        _add_col(col, f"rmsd {_ref_name(ref)} ({kind})")
    _add_col("rmsd_template", f"rmsd {template_group_name} (template)")
    for col in ("surface_accuracy", "plddt", "rmsd_ref"):
        _add_col(col, col)
    _add_col("interface_energy", "pisa_energies", "interface_area")   # cell = "energy (area)"

    similar_names = [_ref_name(r) for c, k, r in (comparisons or []) if k == "similar"]
    different_names = [_ref_name(r) for c, k, r in (comparisons or []) if k == "different"]
    pretty_base = {
        "rmsd_similar": f"rmsd {similar_names[0]} (similar)" if len(similar_names) == 1 else "rmsd (best similar)",
        "rmsd_different": f"rmsd {different_names[0]} (different)" if len(different_names) == 1 else "rmsd (best different)",
        "rmsd_template": f"rmsd {template_group_name} (template)",
    }

    def _tag_label(col):
        for base, pretty in pretty_base.items():
            if col == f"{base}_mean":
                return f"{pretty} mean"
            if col == f"{base}_best":
                return f"{pretty} best"
        return col

    metric_cols = [c for c in shifts.columns if c not in ("tag", "models", "passed", "clusters")]
    tag_cols = ["models"] + metric_cols
    tag_labels = {"models": "models", **{c: _tag_label(c) for c in metric_cols}}

    input_yaml = input_yaml_name = ""
    if config_yml and os.path.exists(str(config_yml)):
        try:
            input_yaml = open(config_yml).read()
            input_yaml_name = os.path.basename(config_yml)
        except Exception as exc:
            logging.error(f"could not read config yml {config_yml}: {exc}")

    data = {
        "reference": os.path.basename(str(reference)), "analysis_dir": outdir,
        "thresholds": thresholds or {}, "focus_desc": focus_desc,
        "focus_metric": focus_col, "focus_lower": focus_lower,
        "input_yaml": input_yaml, "input_yaml_name": input_yaml_name,
        "best_tag": best["tag"] if best is not None else None,
        "best_tag_shifts": ", ".join(f"chain {c} shift {v}" for c, v in
                                     sorted(parse_tag(best["tag"]).items())) if best is not None else None,
        "tag_cols": tag_cols, "tag_labels": tag_labels,
        "tag_rows": [{"tag": r["tag"], "cells": {c: fmt(r.get(c)) for c in tag_cols}}
                     for _, r in shifts.iterrows()],
        "model_cols": model_cols, "col_labels": col_labels,
        "model_rows": [], "chain_stats": {}, "funnel": [],
        "surface_rows": [], "plots": [],
        "cluster_rows": [], "cluster_cols": [], "cluster_plots": [],
        "files": sorted(os.path.basename(p) for p in glob.glob(os.path.join(outdir, "*.csv"))),
    }

    def _cell(r, col):
        val = fmt(r.get(col))
        if col in col_paired and val:
            return f"{val} ({fmt(r.get(col_paired[col]))})"
        return val

    top = df.sort_values(focus_col, ascending=focus_lower) if focus_col else df
    tmin = df["rmsd_template"].min() if "rmsd_template" in df.columns else None
    for _, r in top.iterrows():
        rowclass = "template-row" if tmin is not None and r.get("rmsd_template") == tmin else ""
        data["model_rows"].append({"model": r["model"], "tag": r["tag"], "rowclass": rowclass,
                                   "cells": {c: _cell(r, c) for c in model_cols}})
    if len(top):
        data["best_model"] = top.iloc[0]["model"]

    if "surface_accuracy" in df.columns and df["surface_accuracy"].notna().any():
        sub = df[df["surface_accuracy"].notna()].sort_values("surface_accuracy", ascending=False)
        data["surface_rows"] = [{"tag": r["tag"], "model": r["model"],
                                 "accuracy": fmt(float(r["surface_accuracy"]))}
                                for _, r in sub.iterrows()]

    if not chains.empty:
        for chain, sub in chains.groupby("chain"):
            metric_name = sub["metric"].iloc[0]
            top_mean = sub["metric_mean"].min() if _lower_is_better(metric_name) else sub["metric_mean"].max()
            data["chain_stats"][chain] = {
                "metric": metric_name,
                "rows": [{"shift": int(r["shift"]), "n": int(r["models"]),
                          "mean": fmt(r["metric_mean"]), "best": fmt(r["metric_best"]),
                          "best_of_chain": bool(r["metric_mean"] == top_mean)}
                         for _, r in sub.sort_values("shift").iterrows()]}

    collected_n, collected_tags = collected if collected else (len(df), df["tag"].nunique())
    data["funnel"].append({"stage": "collected", "models": collected_n,
                           "tags": collected_tags, "pct": 100.0})
    for stage in (funnel or []):
        data["funnel"].append({**stage,
                               "pct": round(100.0 * stage["models"] / max(collected_n, 1), 1)})
    if "passed" in df.columns:
        sub = df[df["passed"]]
        if len(sub) != len(df):
            data["funnel"].append({"stage": "passed thresholds", "models": len(sub),
                                   "tags": int(sub["tag"].nunique()),
                                   "pct": round(100.0 * len(sub) / max(collected_n, 1), 1)})

    if "cluster" in df.columns:
        data["cluster_cols"] = ["cluster", "members", "shifts", "models"]
        data["cluster_rows"] = [{"cluster": str(cid), "members": str(len(sub)),
                                 "shifts": ", ".join(sorted(sub["tag"].unique())),
                                 "models": ", ".join(sorted(sub["model"]))}
                                for cid, sub in df.groupby("cluster")]

    data["runs"] = []
    for run in (runs or []):
        href = os.path.relpath(run["html"], outdir) if run.get("html") else None
        data["runs"].append({"tag": run["tag"], "html": href})

    gif_path = build_shift_gif(runs, outdir)
    data["gif"] = os.path.basename(gif_path) if gif_path else None

    refs = []
    for col, kind, ref in (comparisons or []):
        if col in df.columns and np.isfinite(df[col]).any():
            refs.append({"name": os.path.basename(ref), "col": col, "kind": kind})
    data["comparison_refs"] = [{"name": r["name"]} for r in refs]
    data["comparison_rows"] = []
    if refs:
        first = refs[0]
        sub = df[np.isfinite(df[first["col"]])].sort_values(
            first["col"], ascending=(first["kind"] == "similar"))
        for _, r in sub.iterrows():
            cells = [{"rmsd": cellfmt(r.get(rf["col"]))} for rf in refs]
            data["comparison_rows"].append({"model": r["model"], "cells": cells})

    sub = df[df["rmsd_template"].notna()].sort_values("rmsd_template")
    data["template_group_name"] = template_group_name
    data["template_rows"] = [{"model": r["model"], "rmsd": cellfmt(r["rmsd_template"])}
                             for _, r in sub.iterrows()]

    data["best_detail"] = ""
    if focus_col:
        idx = df[focus_col].idxmin() if focus_lower else df[focus_col].idxmax()
        v = fmt(float(df.loc[idx, focus_col]))
        data["best_detail"] = {
            "rmsd_similar": f"CA RMSD {v} Å to the reference it should resemble",
            "rmsd_different": f"CA RMSD {v} Å to the reference it should differ from",
            "surface_accuracy": f"{v}% surface accuracy",
            "rmsd_template": f"global least-squares CA RMSD {v} Å to its template",
        }.get(focus_col, v)

    data["template_bests"] = []
    tmp = df[df["rmsd_template"].notna()].copy()
    tmp["template_name"] = tmp["template"].map(_template_name)
    for name, group in tmp.groupby("template_name"):
        b = group.loc[group["rmsd_template"].idxmin()]
        data["template_bests"].append({
            "name": name, "model": b["model"], "rmsd": fmt(float(b["rmsd_template"]))})

    env = Environment(loader=FileSystemLoader(tmpl_dir))
    html = env.from_string(open(os.path.join(tmpl_dir, "shift_output.html")).read()).render(
        data=data, vairo_css=css, vairo_js=js)
    with open(output_html, "w") as fh:
        fh.write(html)
    logging.error(f"HTML report: {output_html}")
    return output_html



def run_analysis(workdir, reference=None, outdir=None, split_resids=None,
                 ranks=(0, 1, 2, 3, 4), cluster_cutoff=5.0,
                 positive_refs=None, negative_refs=None, surfaces=None,
                 max_rmsd_template=None, min_surface_accuracy=0.0,
                 distance_filters=None, disulfides=None, assembly_template=None,
                 config_yml=None, html=True, **thresholds):

    max_rmsd_template = float("inf") if max_rmsd_template is None else float(max_rmsd_template)
    min_surface_accuracy = float(min_surface_accuracy or 0.0)

    outdir = outdir or os.path.join(workdir, "analysis")
    os.makedirs(outdir, exist_ok=True)

    df = collect_models(workdir, ranks=ranks, dest=outdir)
    if df.empty:
        raise Exception(f"no models found under {workdir}")
    collected = (len(df), df["tag"].nunique())

    tmpl_sup_dir = None if _as_list(positive_refs) else os.path.join(outdir, "superposed", "template")
    df = add_template_rmsd(df, superposed_dir=tmpl_sup_dir)
    df = add_perchain_template_rmsd(df)
    rlabel = "∞" if not np.isfinite(max_rmsd_template) else f"{max_rmsd_template:g}"
    df, stage = apply_filter(df, "rmsd_template_chain_max", max_rmsd_template, "max",
                             f"per-chain template CA RMSD ≤ {rlabel} Å")
    funnel = [stage]

    if reference:
        df = add_rmsd(df, reference, "rmsd_ref")

    comparisons = []
    for kind, refs in (("similar", _as_list(positive_refs)),
                       ("different", _as_list(negative_refs))):
        for ref in refs:
            if not os.path.exists(str(ref)):
                logging.error(f"comparison reference not found, skipping: {ref}")
                continue
            col = _ref_column(kind, ref)
            sup_dir = os.path.join(outdir, "superposed", col)
            df = add_rmsd(df, str(ref), col, superposed_dir=sup_dir)
            comparisons.append((col, kind, str(ref)))

    similar_cols = [c for c, k, _ in comparisons if k == "similar"]
    different_cols = [c for c, k, _ in comparisons if k == "different"]
    if similar_cols:
        df["rmsd_similar"] = df[similar_cols].min(axis=1)
    if different_cols:
        df["rmsd_different"] = df[different_cols].min(axis=1)

    df = add_plddt(df)

    pack_col, chain_copies = "pdb", None
    if assembly_template and os.path.exists(str(assembly_template)):
        df, chain_copies = build_assembly(df, str(assembly_template),
                                          os.path.join(outdir, "assemblies"))
        pack_col = "assembly"

    df, geom_filters = add_geometry_checks(df, distance_filters, disulfides,
                                           pdb_col=pack_col, chain_copies=chain_copies)
    for col, mx, label in geom_filters:
        df, stage = apply_filter(df, col, mx, "max", label)
        funnel.append(stage)

    if surfaces:
        df = add_surface(df, surfaces, workdir=os.path.join(outdir, "surfaces"),
                         pdb_col=pack_col, chain_copies=chain_copies)
        df, stage = apply_filter(df, "surface_accuracy", min_surface_accuracy, "min",
                                 f"surface accuracy ≥ {min_surface_accuracy:g}%")
        funnel.append(stage)
    if split_resids:
        df = add_interface_metrics(df, split_resids, workdir=os.path.join(outdir, "interfaces"),
                                   pdb_col=pack_col)
    if cluster_cutoff:
        df = add_clusters(df, cutoff=cluster_cutoff)

    df = score_models(df, **thresholds)
    focus_col, ascending, focus_desc = choose_focus(df)
    logging.error(f"Ranking focus: {focus_desc} ({focus_col}, "
                  f"{'lower' if ascending else 'higher'} is better)")
    shifts = summarise_shifts(df, focus_col=focus_col, ascending=ascending)
    chains = summarise_chains(df, metric=focus_col)

    for name, table in (("models", df), ("shifts", shifts), ("chains", chains)):
        table.to_csv(os.path.join(outdir, f"{name}.csv"), index=False)
    logging.error(f"Wrote models.csv / shifts.csv / chains.csv to {outdir}")

    if not (config_yml and os.path.exists(str(config_yml))):
        config_yml = next((p for p in (os.path.join(workdir, "shift_input.yml"),
                                       os.path.join(workdir, "config", "base_config.yml"))
                           if os.path.exists(p)), None)
    if html:
        write_html_report(df, shifts, chains, outdir, reference=reference or "",
                          thresholds=dict(thresholds, split_resids=split_resids),
                          runs=collect_runs(workdir), comparisons=comparisons,
                          focus=(focus_col, focus_desc), config_yml=config_yml,
                          collected=collected, funnel=funnel)

    if len(shifts):
        logging.error(f"Best shift: {shifts.iloc[0]['tag']}")
    return df, shifts, chains
