import dataclasses
from typing import List

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
