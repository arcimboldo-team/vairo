import dataclasses
from typing import Dict, List

from libs import utils

@dataclasses.dataclass(frozen=True)
class AlignmentDatabase:
    chain: str
    fasta_path: str
    database_path: str


@dataclasses.dataclass(frozen=True)
class Alignment:
    aligned_columns: int
    total_columns: int
    evalue: str
    identities: int
    hhr_path: str
    extracted_path: str
    database: AlignmentDatabase


@dataclasses.dataclass
class Interface:
    name: str
    res_list: List[int]
    dist_coverage: float = dataclasses.field(init=False, repr=False, default=None)
    core: int = dataclasses.field(init=False, repr=False, default=None)
    dist_plot: str = dataclasses.field(init=False, repr=False, default=None)
    encoded_dist_plot: bytes = dataclasses.field(init=False, repr=False, default=None)

    def add_frobenius_information(self, dist_coverage: float, core: int, dist_plot: str):
        self.dist_coverage = dist_coverage
        self.core = core
        self.dist_plot = dist_plot
        self.encoded_dist_plot = utils.encode_data(dist_plot)


@dataclasses.dataclass(frozen=True)
class Frobenius:
    template: str
    dist_coverage: float
    encoded_dist_plot: bytes
    dist_plot: str
    ang_coverage: float
    ang_plot: str
    encoded_ang_plot: bytes


@dataclasses.dataclass(frozen=True)
class OpenmmEnergies:
    kinetic: str
    potential: List[int]


@dataclasses.dataclass(frozen=True)
class TemplateRanked:
    template: str
    rmsd: float
    aligned_residues: int
    total_residues: int


class Ranked:
    def __init__(self, ranked_path):
        self.path: str
        self.name: str
        self.split_path: str
        self.minimized_path: str
        self.plddt: float
        self.ah: int
        self.bs: int
        self.total_residues: int
        self.superposition_templates: List[TemplateRanked] = []
        self.mapping: Dict = {}
        self.energies: OpenmmEnergies = None
        self.interfaces: List[Interface] = []
        self.frobenius_plots: List[Frobenius] = []
        self.filtered: bool = False
        self.best: bool = False

        self.path = ranked_path
        self.name = utils.get_file_name(ranked_path)
    
    def set_path(self, path: str):
        self.path = path

    def set_plddt(self, plddt: float):
        self.plddt = plddt

    def set_filtered(self, filtered: bool):
        self.filtered = filtered

    def set_best(self, best: bool):
        self.best = best

    def set_mapping(self, mapping: Dict):
        self.mapping = mapping

    def set_split_path(self, path: str):
        self.split_path = path

    def set_minimized_path(self, path: str):
        self.minimized_path = path

    def add_template(self, template: TemplateRanked):
        self.superposition_templates.append(template)

    def add_interface(self, interface: Interface):
        self.interfaces.append(interface)

    def set_secondary_structure(self, ah: int, bs: int, total_residues: int):
        self.ah = ah
        self.bs = bs
        self.total_residues = total_residues

    def set_energies(self, energies: OpenmmEnergies):
        self.energies = energies

    def add_interfaces_frobenius_plot(self, plot: str):
        self.interfaces_frobenius_plot.append(plot)

    def add_frobenius_plot(self, template: str, dist_plot: str, ang_plot: str, dist_coverage: float, ang_coverage: float):
        frobenius = Frobenius(
                        template=template, 
                        dist_plot=dist_plot, 
                        encoded_dist_plot=utils.encode_data(dist_plot),
                        ang_plot=ang_plot,
                        encoded_ang_plot=utils.encode_data(ang_plot),
                        dist_coverage=dist_coverage, 
                        ang_coverage=ang_coverage
                        )
        self.frobenius_plots.append(frobenius)

    def sort_template_rankeds(self):
        self.superposition_templates.sort(key=lambda x: (x.rmsd is None, x.rmsd))