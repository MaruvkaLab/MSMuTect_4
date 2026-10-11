from dataclasses import dataclass


@dataclass
class Mutation:
    seq: str
    insertion: bool = False
    deletion: bool = False
    substitution: bool = False
    position: int = 0
    enters_or_exits_locus: bool = False # deletion begins inside locus, and exits it, or begins outside and enters it


    def __eq__(self, other):
        return (self.seq == other.seq and self.insertion == other.insertion and self.deletion == other.deletion
                and self.substitution == other.substitution and self.position==other.position)
