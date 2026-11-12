"""PETfold / PETcofold: RNA structure prediction of multiple alignments.

Python interface to the C++ edition of PETfold 2.2 and PETcofold 3.4
(Seemann et al., NAR 2008; Bioinformatics 2011).

    >>> import petfold
    >>> r = petfold.fold("example.fasta")
    >>> r.structure, r.score
    >>> r.paired_pet          # numpy (L x L) PET base pair reliabilities
    >>> c = petfold.cofold("rna1.fasta", "rna2.fasta")
"""
from importlib import resources as _resources

from . import _core
from ._api import CofoldResult, FoldResult, cofold, fold, read_fasta
from ._core import EngineError, PetfoldError

__version__ = _core.__version__
__all__ = ["fold", "cofold", "read_fasta", "FoldResult", "CofoldResult", "PetfoldError", "EngineError",
           "data_dir", "__version__"]


def data_dir() -> str:
    """Directory holding the bundled Pfold grammar files (scfg.rate, article.grm)."""
    return str(_resources.files(__name__) / "data")


_core.set_default_data_dir(data_dir())
