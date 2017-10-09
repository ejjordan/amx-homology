#!/usr/bin/env python

"""
PROTEIN STRUCTURE REPAIR
Fix a protein structure from the pdb
"""

from amx import *
import os

init()
make_step(settings.step)
if state.pdb_source: get_pdb(state.pdb_source)
else: get_start_structure(state.start_structure)
state['pdbinfo']=get_pdb_sequence()
state['chains_info']=get_all_chains(**state)
write_ali_file(**state)
export_modeller_settings(**state)
bash('python ../amx/homology/modeller_script.py',cwd=state['step']+'/', log='modeller')
get_best_structure(**state)

