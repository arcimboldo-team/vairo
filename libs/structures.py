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


@dataclasses.dataclass(frozen=True)
class InterfaceName:
    name: str
    res_list: List[int]


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
        self.interfaces: List[InterfaceName] = []
        self.frobenius_plot: List[str] = []
        self.filtered: bool = False

        self.path = ranked_path
        self.name = utils.get_file_name(ranked_path)
    
    def set_path(self, path: str):
        self.path = path

    def set_plddt(self, plddt: float):
        self.plddt = plddt

    def set_filtered(self, filtered: bool):
        self.filtered = filtered

    def set_mapping(self, mapping: Dict):
        self.mapping = mapping

    def set_split_path(self, path: str):
        self.split_path = path

    def set_minimized_path(self, path: str):
        self.minimized_path = path

    def add_template(self, template: TemplateRanked):
        self.superposition_templates.append(template)

    def add_interface(self, interface: InterfaceName):
        self.interfaces.append(interface)

    def set_secondary_structure(self, ah: int, bs: int, total_residues: int):
        self.ah = ah
        self.bs = bs
        self.total_residues = total_residues

    def set_energies(self, energies: OpenmmEnergies):
        self.energies = energies

    def add_frobenius_plot(self, plot):
        self.frobenius_plot.append(plot)