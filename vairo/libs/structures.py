import dataclasses
import os
import statistics
import sys
import shutil
from typing import Dict, List
from libs import utils, bioutils


@dataclasses.dataclass(frozen=True)
class CCAnalysisOutput:
    coord: List[float]
    module: float
    angle: float


class BinariesPath:
    def __init__(self, binaries_path):
        pd2cc_path: str
        cc_analysis_path: str
        hinges_path: str
        spong_path: str

        if sys.platform == "darwin":
            self.cc_analysis_path = os.path.join(binaries_path, 'cc_analysis_mac')
            self.pd2cc_path = os.path.join(binaries_path, 'pdb2cc_mac')
            self.hinges_path = os.path.join(binaries_path, 'hinges_mac')
            self.spong_path = os.path.join(binaries_path, 'spong_mac')
        else:
            self.cc_analysis_path = os.path.join(binaries_path, 'cc_analysis_linux')
            self.pd2cc_path = os.path.join(binaries_path, 'pdb2cc_linux')
            self.hinges_path = os.path.join(binaries_path, 'hinges_linux')
            self.spong_path = os.path.join(binaries_path, 'spong_linux')


@dataclasses.dataclass(frozen=True)
class Cluster:
    name: str
    label: str
    path: str
    relative_path: str
    encoded_path: bytes
    rankeds: Dict
    templates: Dict


@dataclasses.dataclass(frozen=True)
class Hinges:
    decreasing_rmsd_middle: float
    decreasing_rmsd_total: float
    one_rmsd: float
    middle_rmsd: float
    min_rmsd: float
    overlap: float
    groups: List


@dataclasses.dataclass
class FeaturesInput:
    path: str
    keep_msa: int
    keep_templates: int
    msa_mask: List[int]
    replace_sequence: str
    numbering_features: List[int]
    numbering_query: List[int]
    mutate_residues: dict
    num_msa: int = dataclasses.field(default=0)
    num_templates: int = dataclasses.field(default=0)

    def add_information(self, num_msa: int = 0, num_templates: int = 0):
        self.num_msa = num_msa
        self.num_templates = num_templates


@dataclasses.dataclass
class Library:
    path: str
    aligned: str
    add_to_msa: bool
    add_to_templates: bool
    numbering_query: List[int]
    numbering_library: List[int]
    num_msa: int = dataclasses.field(default=0)
    num_templates: int = dataclasses.field(default=0)

    def add_information(self, num_msa: int = 0, num_templates: int = 0):
        self.num_msa = num_msa
        self.num_templates = num_templates


@dataclasses.dataclass(frozen=True)
class Dendogram:
    dendogram_list: List[str]
    dendogram_plot: str
    encoded_dendogram_plot: bytes


@dataclasses.dataclass(frozen=True)
class Alignment:
    aligned_columns: int
    total_columns: int
    evalue: str
    identities: int
    hhr_path: str
    mapping: Dict
    chain: str


@dataclasses.dataclass(frozen=True)
class GanttPlot:
    plot_both: bytes
    legend_both: str
    plot_template: bytes
    legend_template: str
    plot_msa: bytes
    legend_msa: str


@dataclasses.dataclass
class Interface:
    name: str
    res_chain1: List[int]
    res_chain2: List[int]
    chain1: str
    chain2: str
    se_gain1: float
    se_gain2: float
    solvation1: float
    solvation2: float
    area: float
    deltaG: float
    nhb: float
    path: str = dataclasses.field(default=None)

    def set_structure(self, path: str):
        self.path = path


@dataclasses.dataclass(frozen=True)
class Frobenius:
    template: str
    dist_coverage: float
    encoded_dist_plot: bytes
    dist_plot: str
    ang_coverage: float
    ang_plot: str
    encoded_ang_plot: bytes
    core: int


@dataclasses.dataclass(frozen=True)
class PdbRanked:
    pdb: str
    rmsd: float
    aligned_residues: int
    total_residues: int
    qscore: float


class Pdb:
    def __init__(self, path: str):
        self.path: str
        self.name: str
        self.split_path: str = None
        self.compactness: float
        self.ramachandran: float
        self.ah: int
        self.bs: int
        self.total_residues: int
        self.interfaces: List[Interface] = []
        self.accepted_interfaces: bool = False

        self.path = path
        self.name = utils.get_file_name(path)

    def set_path(self, path: str):
        self.path = path

    def set_split_path(self, path: str):
        self.split_path = path

    def set_compactness(self, compactness: float):
        self.compactness = compactness

    def set_ramachandran(self, ramachandran: float):
        self.ramachandran = ramachandran

    def set_secondary_structure(self, ah: int, bs: int, total_residues: int):
        self.ah = ah
        self.bs = bs
        self.total_residues = total_residues

    def set_interfaces(self, interfaces: List[Interface]):
        self.interfaces = interfaces

    def set_accepted_interfaces(self, value: bool):
        self.accepted_interfaces = value

    def get_interfaces_with_path(self) -> List[Interface]:
        return [interface for interface in self.interfaces if interface.path is not None]


class ExperimentalPdb(Pdb):
    def __init__(self, path: str):
        super().__init__(path=path)
        self.split_path = path


class TemplateExtracted(Pdb):
    def __init__(self, path: str):
        super().__init__(path=path)
        self.percentage_list: List[float]
        self.identity: float
        self.sequence_msa: str
        self.template: str = ''
        self.originalseq_path: str = ''

    def set_template(self, template, originalseq_path: str):
        self.template = template
        self.originalseq_path = originalseq_path
        if self.template is not None:
            shutil.copy2(self.template.template_originalseq_path, self.originalseq_path)
        else:
            shutil.copy2(self.split_path, self.originalseq_path)

    def add_percentage(self, percentage_list: List[float]):
        self.percentage_list = percentage_list

    def set_identity(self, identity: float):
        self.identity = identity

    def set_sequence_msa(self, sequence_msa: str):
        self.sequence_msa = sequence_msa


class Ranked(Pdb):
    def __init__(self, path: str):
        super().__init__(path=path)
        self.minimized_path: str
        self.plddt: int
        self.superposition_templates: List[PdbRanked] = []
        self.superposition_experimental: List[PdbRanked] = []
        self.potential_energy: float = None
        self.frobenius_plots: List[Frobenius] = []
        self.filtered: bool = False
        self.best: bool = False
        self.qscore: float
        self.qscore_dict: Dict = {}
        self.encoded: bytes

    def set_plddt(self):
        if self.split_path is not None:
            split = bioutils.read_bfactors_from_residues(pdb_path=self.split_path)
        else:
            split = bioutils.read_bfactors_from_residues(pdb_path=self.path)
        split = [item for sublist in split.values() for item in sublist]
        plddt_list = [value for value in split if value is not None]
        self.plddt = round(statistics.mean(map(float, plddt_list)), 2)

    def set_qscore(self, qscore: float):
        self.qscore = round(qscore, 3) if qscore is not None else qscore

    def set_ranked_to_qscore_dict(self, qscore: float, ranked_name: str):
        self.qscore_dict[ranked_name] = round(qscore, 3) if qscore is not None else qscore

    def set_filtered(self, filtered: bool):
        self.filtered = filtered

    def set_best(self, best: bool):
        self.best = best

    def set_minimized_path(self, path: str):
        self.minimized_path = path

    def add_template(self, template: PdbRanked):
        self.superposition_templates.append(template)

    def add_experimental(self, experimental: PdbRanked):
        self.superposition_experimental.append(experimental)

    def set_potential_energy(self, potential_energy: float):
        self.potential_energy = potential_energy

    def set_encoded(self, path: str):
        self.encoded = utils.encode_data(path)

    def sort_template_rankeds(self):
        self.superposition_templates.sort(key=lambda x: (x.qscore is None, x.qscore), reverse=True)

    def sort_experimental_rankeds(self):
        self.superposition_experimental.sort(key=lambda x: (x.aligned_residues is None, x.aligned_residues),
                                             reverse=True)
