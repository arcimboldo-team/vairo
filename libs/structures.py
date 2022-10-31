import dataclasses
from typing import List

@dataclasses.dataclass(frozen=True)
class InterfaceName:
    name: str
    res_list: List[int]

@dataclasses.dataclass(frozen=True)
class AlignmentSequence:
    fasta_path: str
    database_path: str
    hhr_path: str

@dataclasses.dataclass(frozen=True)
class AlignmentDatabase:
    chain: str
    fasta_path: str
    database_path: str