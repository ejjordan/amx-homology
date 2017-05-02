#!/usr/bin/env python

"""
modeller homology modeling script
"""

from modeller import *
from modeller.automodel import *
import sys

execfile('settings-homology.py')
class mymodel(automodel):
        def special_patches(self, aln):
                chains=[template_struct_chain]
                start_res=[starting_residue]
                if other_chains:
                        for chain in other_chains.keys(): chains.append(chain)
                        for start in other_chains.values(): start_res.append(start)
                self.rename_segments(segment_ids=chains,renumber_residues=start_res)
doalign2d = True
if doalign2d:
	env = environ()
        env.io.hetatm = True
        env.io.water = False
	aln = alignment(env)
        model_string=('FIRST:'+template_struct_chain, 'LAST:'+template_struct_chain)
        if other_chains:model_string=('FIRST:@ ', 'END:')
	mdl = model(env,
		file=template_struct,
		model_segment=(model_string))
	aln.append_model(mdl,
		align_codes=template_struct+template_struct_chain,
		atom_files=template_struct+'.pdb')
	aln.append(file=target_seq+'.ali', align_codes=target_seq)
	aln.align2d()
	aln.write(file='align2d.ali', alignment_format='PIR')
	aln.write(file='align2d.pap', alignment_format='PAP')
	afile = 'align2d.ali'
else:
	env = environ()
	aln = alignment(env)
	afile = 'align2d-custom.ali'
a = mymodel(env,
	alnfile=afile,
	knowns=template_struct+template_struct_chain,
	sequence=target_seq,
	assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = n_models
a.make()
