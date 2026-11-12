# PETfold / PETcofold — C++ edition with Python bindings


## Installation

Both routes need the **ViennaRNA C library** 2.x (tested with 2.6.3) with its `RNAlib2.pc` pkg-config file:

```bash
conda install -c conda-forge -c bioconda viennarna     # or build ViennaRNA from source
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig      # if pkg-config does not find RNAlib2
```

### Python package (includes the command-line tools)

```bash
pip install .                 # from this directory
```

This installs the `petfold` Python module and the commands `PETfold`, `PETcofold`, `drawdot` and `drawphyl`. The Pfold grammar files are bundled; `PETFOLDBIN` is no longer needed.

### Command-line tools only (CMake)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /usr/local     # optional; installs grammar files to share/petfold
```

Options: `-DCMAKE_CXX_STANDARD=23`, `-DPETFOLD_SANITIZE=ON` (AddressSanitizer + UBSan), `-DPETFOLD_BUILD_PYTHON=ON`.

Grammar files are searched in this order: `--grammar-dir`, `$PETFOLDBIN`, the current directory (as in 2.2), the installed data directory. The build directory contains a copy.

## Command line

Options and output are the same as in PETfold 2.2 and PETcofold 3.4 (see `PETfold --help`, `PETcofold --help`):

```bash
PETfold  -f example/example.fasta
PETcofold -f example/example1.fasta -f example/example2.fasta
```

New options: `--seed` (reproducible `--suboptimal` sampling), `--grammar-dir`, `--version`, and for PETcofold `--legacy-cutpoint` (reproduces the 3.4 cut point handling, for comparisons). `--noevo` and `--nothermo` existed before but were undocumented.

## Python

```python
import petfold

r = petfold.fold("example/example.fasta")          # path, {id: seq} or [(id, seq), ...]
print(r.structure, r.score)
r.paired_pet        # (L x L) numpy array, PET base pair reliabilities (alignment columns)
r.single_pet        # (L,) unpaired reliabilities
r.paired_evo, r.paired_thermo     # the two model components
r.base_pairs        # [(i, j), ...] of the predicted structure
r.pair_table(0.5)   # [(i, j, reliability), ...]

r = petfold.fold(aln, tree="(a,b,(c,d))", structure="((((....))))", beta=0.5)
r = petfold.fold(aln, suboptimal=1000, seed=1)       # r.suboptimals = [(structure, score), ...]

c = petfold.cofold("rna1.fasta", "rna2.fasta", extstem=True)
print(c.structure)          # '{' '}' partial structures, '(' ')' step-2 pairs, '[' ']' interactions
c.paired_pet                # (L1+L2) x (L1+L2) duplex reliabilities
```

All parameters of the command line are keyword arguments (`gap`, `evocon_bp`, `evocon_ss`, `alpha`, `beta`, `use_evo`, `use_thermo`, and for `cofold` `delta`, `gamma`, `extstem`, `legacy_cutpoint`). Errors raise `petfold.PetfoldError`. Verbose output is returned in `result.log`. The computation releases the GIL, so alignments can be folded in parallel threads.

Arrays are in alignment coordinates: columns removed because of too many gaps are zero, as in the pp-file written by `-r`.


## Layout

| Path | Content |
|---|---|
| `src/pfold/` | Pfold engine (B. Knudsen), ported to C++ in namespace `pfold` with scoped memory management |
| `src/petfold/` | PETfold/PETcofold models, pipelines and command-line front ends (namespace `petfold`) |
| `src/apps/` | `main()` of the four programs; `drawdot`, `drawphyl` |
| `python/` | pybind11 module `petfold._core` and the Python package |
| `data/` | Pfold grammar files `scfg.rate`, `article.grm` |
| `example/` | test data |


