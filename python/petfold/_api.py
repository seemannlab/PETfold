"""High-level Python API around petfold._core."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from . import _core

AlignmentLike = Union[str, os.PathLike, Mapping[str, str], Sequence[Tuple[str, str]]]


def read_fasta(path: Union[str, os.PathLike]) -> List[Tuple[str, str]]:
    """Reads a FASTA file into a list of (identifier, sequence) pairs."""
    records: List[Tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                records.append((line[1:], ""))
            elif records:
                records[-1] = (records[-1][0], records[-1][1] + line.strip())
    return records


def _records(alignment: AlignmentLike) -> Tuple[List[str], List[str]]:
    if isinstance(alignment, (str, os.PathLike)):
        recs = read_fasta(alignment)
    elif isinstance(alignment, Mapping):
        recs = list(alignment.items())
    else:
        recs = [(str(i), str(s)) for i, s in alignment]
    if not recs:
        raise ValueError("empty alignment")
    return [r[0] for r in recs], [r[1] for r in recs]


def _read_arg(value: Optional[str]) -> str:
    """Like the command-line tools: a file name is replaced by its first line."""
    if not value:
        return ""
    if os.path.isfile(value):
        with open(value) as fh:
            return fh.readline().rstrip("\r\n")
    return value


@dataclass
class FoldResult:
    """Result of :func:`fold`. Arrays are in alignment coordinates (L columns);
    removed gap columns are zero. Pair matrices are symmetric."""

    structure: str
    score: float
    ensemble_diversity: float
    pfold_structure: str
    constraints: str
    tree: str
    gap_columns: List[int]
    ids: List[str]
    sequences: List[str]
    paired_pet: np.ndarray
    single_pet: np.ndarray
    paired_evo: np.ndarray
    single_evo: np.ndarray
    paired_thermo: np.ndarray
    single_thermo: np.ndarray
    suboptimals: List[Tuple[str, float]] = field(default_factory=list)
    log: str = ""

    @property
    def base_pairs(self) -> List[Tuple[int, int]]:
        """Base pairs (0-based alignment columns) of the PETfold structure."""
        return _pairs(self.structure)

    def pair_table(self, min_reliability: float = 0.0):
        """(i, j, PET reliability) for all i < j above a threshold."""
        i, j = np.triu_indices(self.paired_pet.shape[0], 1)
        v = self.paired_pet[i, j]
        keep = v > min_reliability
        return list(zip(i[keep].tolist(), j[keep].tolist(), v[keep].tolist()))


@dataclass
class CofoldResult:
    """Result of :func:`cofold`. Arrays cover both alignments (L1 + L2 columns)."""

    structure: str
    structure1: str
    structure2: str
    intermolecular: str
    partial_structure1: str
    partial_structure2: str
    petfold_structure1: str
    petfold_structure2: str
    delta1: float
    delta2: float
    partial_prob1: Tuple[float, float]
    partial_prob2: Tuple[float, float]
    score: float
    common_ids: List[str]
    gap_columns1: List[int]
    gap_columns2: List[int]
    paired_pet: np.ndarray
    single_pet: np.ndarray
    log: str = ""


def _pairs(db: str) -> List[Tuple[int, int]]:
    stacks = {"(": [], "[": [], "{": [], "<": []}
    close = {")": "(", "]": "[", "}": "{", ">": "<"}
    out = []
    for k, c in enumerate(db):
        if c in stacks:
            stacks[c].append(k)
        elif c in close and stacks[close[c]]:
            out.append((stacks[close[c]].pop(), k))
    return sorted(out)


def fold(
    alignment: AlignmentLike,
    *,
    tree: Optional[str] = None,
    structure: Optional[str] = None,
    suboptimal: int = 0,
    seed: Optional[int] = None,
    gap: float = 0.25,
    evocon_bp: float = 0.9,
    evocon_ss: float = 1.0,
    alpha: float = 0.2,
    beta: float = 1.0,
    use_evo: bool = True,
    use_thermo: bool = True,
    verbose: bool = False,
    grammar_dir: Optional[str] = None,
) -> FoldResult:
    """Predicts the consensus structure of an RNA alignment with PETfold.

    ``alignment`` is a FASTA path, a ``{id: sequence}`` mapping or a sequence of
    ``(id, sequence)`` pairs. ``tree`` and ``structure`` accept a string or a file
    name (first line is used), exactly like ``-t``/``-s`` of the command line.
    """
    ids, seqs = _records(alignment)
    d = _core.fold(ids, seqs, tree=_read_arg(tree), structure=_read_arg(structure), suboptimal=suboptimal,
                   seed=seed, gap=gap, evocon_bp=evocon_bp, evocon_ss=evocon_ss, alpha=alpha, beta=beta,
                   use_evo=use_evo, use_thermo=use_thermo, verbose=verbose, grammar_dir=grammar_dir or "")
    return FoldResult(**d)


def cofold(
    alignment1: AlignmentLike,
    alignment2: AlignmentLike,
    *,
    tree: Optional[str] = None,
    structure: Optional[str] = None,
    structure1: Optional[str] = None,
    structure2: Optional[str] = None,
    delta: float = 0.9,
    gamma: float = 0.1,
    extstem: bool = False,
    legacy_cutpoint: bool = False,
    gap: float = 0.25,
    evocon_bp: float = 0.9,
    evocon_ss: float = 1.0,
    alpha: float = 0.2,
    beta: float = 1.0,
    use_evo: bool = True,
    use_thermo: bool = True,
    verbose: bool = False,
    grammar_dir: Optional[str] = None,
) -> CofoldResult:
    """Predicts the joint structure of two RNA alignments with PETcofold."""
    i1, s1 = _records(alignment1)
    i2, s2 = _records(alignment2)
    d = _core.cofold(i1, s1, i2, s2, tree=_read_arg(tree), structure=_read_arg(structure),
                     structure1=_read_arg(structure1), structure2=_read_arg(structure2), delta=delta, gamma=gamma,
                     extstem=extstem, legacy_cutpoint=legacy_cutpoint, gap=gap, evocon_bp=evocon_bp,
                     evocon_ss=evocon_ss, alpha=alpha, beta=beta, use_evo=use_evo, use_thermo=use_thermo,
                     verbose=verbose, grammar_dir=grammar_dir or "")
    return CofoldResult(**d)
