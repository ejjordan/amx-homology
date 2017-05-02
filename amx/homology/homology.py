#!/usr/bin/env python

"""
PROTEIN HOMOLOGY
Fix or mutate a protein structure
"""

from amx import *

init()
make_step(settings.step)
write_mdp()
if state.pdb_source: get_pdb(state.pdb_source)
else: get_start_structure(state.start_structure)
