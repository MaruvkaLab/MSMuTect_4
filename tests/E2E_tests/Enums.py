from enum import Enum, auto

class STAGE(Enum):  # stages of analysis
    BAMS = auto()
    HISTOGRAMS = auto()
    ALLELES = auto()
    MUTATIONS_FULL = auto()
    MUTATIONS_EFFICIENT = auto() # doesnt calculate all alleles


class EXECUTION_TYPE(Enum): # how to run msmutect
    BASH = auto()
    EXECUTABLE = auto()
    DOCKER = auto()

