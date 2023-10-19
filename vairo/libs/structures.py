import dataclasses
import os
import sys
from typing import Dict, List
from libs import utils


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
    msa_delete: List[int]
    sequence: str
    positions: List[int]
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
    positions: Dict
    positions_list: List[str]
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


@dataclasses.dataclass(frozen=True)
class InterfaceTemplate:
    template: str
    dist_coverage: float
    core: int
    dist_plot: str
    encoded_dist_plot: bytes


@dataclasses.dataclass
class Interface:
    name: str
    res_list: List[int]
    chain1: str
    chain2: str
    se_gain1: float
    se_gain2: float
    solvation1: float
    solvation2: float
    interface_template: List[InterfaceTemplate] = dataclasses.field(default_factory=list)

    def add_frobenius_information(self, template: str, dist_coverage: float, core: int, dist_plot: str):
        interface = InterfaceTemplate(
            template=template,
            dist_coverage=dist_coverage,
            dist_plot=dist_plot,
            encoded_dist_plot=utils.encode_data(dist_plot),
            core=core
        )
        self.interface_template.append(interface)


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
        self.split_path: str
        self.compactness: float
        self.ramachandran: float
        
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

class Template (Pdb):
    def __init__(self, path: str):
        super().__init__(path=path)
        self.percentage_list: List[float]
        self.identity: float
        self.sequence_msa: str

    def add_percentage(self, percentage_list: List[float]):
        self.percentage_list = percentage_list

    def set_identity(self, identity: float):
        self.identity = identity

    def set_sequence_msa(self, sequence_msa: str):
        self.sequence_msa = sequence_msa

class Ranked (Pdb):
    def __init__(self, path: str):
        super().__init__(path=path)
        self.minimized_path: str
        self.plddt: int
        self.ah: int
        self.bs: int
        self.total_residues: int
        self.superposition_templates: List[PdbRanked] = []
        self.superposition_experimental: List[PdbRanked] = []
        self.mapping: Dict = {}
        self.potential_energy: float = None
        self.interfaces: List[Interface] = []
        self.frobenius_plots: List[Frobenius] = []
        self.filtered: bool = False
        self.best: bool = False
        self.rmsd: float
        self.rmsd_dict: Dict = {}
        self.encoded: bytes

    def set_plddt(self, plddt: float):
        self.plddt = round(plddt)

    def set_rmsd(self, rmsd: float):
        self.rmsd = round(rmsd, 2) if rmsd is not None else rmsd

    def set_ranked_to_rmsd_dict(self, rmsd: float, ranked_name: str):
        self.rmsd_dict[ranked_name] = round(rmsd, 2) if rmsd is not None else rmsd

    def set_filtered(self, filtered: bool):
        self.filtered = filtered

    def set_best(self, best: bool):
        self.best = best

    def set_mapping(self, mapping: Dict):
        self.mapping = mapping

    def set_minimized_path(self, path: str):
        self.minimized_path = path

    def add_template(self, template: PdbRanked):
        self.superposition_templates.append(template)

    def add_experimental(self, experimental: PdbRanked):
        self.superposition_experimental.append(experimental)

    def add_interface(self, interface: Interface):
        self.interfaces.append(interface)

    def set_secondary_structure(self, ah: int, bs: int, total_residues: int):
        self.ah = ah
        self.bs = bs
        self.total_residues = total_residues

    def set_potential_energy(self, potential_energy: float):
        self.potential_energy = potential_energy

    def add_interfaces_frobenius_plot(self, plot: str):
        self.interfaces_frobenius_plot.append(plot)

    def set_encoded(self, path: str):
        self.encoded = utils.encode_data(path)

    def add_frobenius_plot(self, template: str, dist_plot: str, ang_plot: str, dist_coverage: float,
                           ang_coverage: float, core: float):
        frobenius = Frobenius(
            template=template,
            dist_plot=dist_plot,
            encoded_dist_plot=utils.encode_data(dist_plot),
            ang_plot=ang_plot,
            encoded_ang_plot=utils.encode_data(ang_plot),
            dist_coverage=dist_coverage,
            ang_coverage=ang_coverage,
            core=core
        )
        self.frobenius_plots.append(frobenius)

    def sort_template_rankeds(self):
        self.superposition_templates.sort(key=lambda x: (x.qscore is None, x.qscore), reverse=True)

    def sort_experimental_rankeds(self):
        self.superposition_experimental.sort(key=lambda x: (x.qscore is None, x.qscore), reverse=True)
