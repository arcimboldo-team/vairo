import dataclasses
from typing import List

@dataclasses.dataclass(frozen=True)
class InterfaceName:
    name: str
    res_list: List[int]
