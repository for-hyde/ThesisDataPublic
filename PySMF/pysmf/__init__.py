from .core.bamobject import BamObject
#from .core.sparsematrix import SparseMatrix
from .analysis.singlesite import call_SMF
from .analysis.multisite import call_SMF_averages  

__all__ = [
    call_SMF,
    call_SMF_averages,
    BamObject
]

