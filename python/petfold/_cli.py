"""Console entry points: PETfold, PETcofold, drawdot, drawphyl (same options as the C++ tools)."""
import sys

from . import _core


def _run(fn):
    sys.stdout.flush()
    sys.exit(fn(sys.argv[1:]))


def petfold():
    _run(_core.petfold_main)


def petcofold():
    _run(_core.petcofold_main)


def drawdot():
    _run(_core.drawdot_main)


def drawphyl():
    _run(_core.drawphyl_main)
