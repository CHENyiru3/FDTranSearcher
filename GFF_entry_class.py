from typing import Tuple
from dataclasses import dataclass
# Define the GFFEntry dataclass
@dataclass
class GFFEntry:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float
    strand: str
    phase: str
    attributes: dict

    def get_position_key(self) -> Tuple[str, int, int, str]:
        return (self.seqid, self.start, self.end, self.strand)