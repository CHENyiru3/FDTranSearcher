from dataclasses import dataclass
from typing import Tuple
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

# analyze gff format attributes
def parse_gff_attributes(attr_string: str) -> dict:
    attrs = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = value
    return attrs