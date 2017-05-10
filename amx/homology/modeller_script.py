#!/usr/bin/env python

"""
modeller homology modeling script
"""

from modeller import *
from modeller.automodel import *
import sys,json

settings_fn = 'settings.json'
with open(settings_fn,'r') as fp:
        settings = json.load(fp)

chains=[chain.encode('ascii','ignore') for chain in settings['chains_info']]
start_resids=[settings['chains_info'][chain]['startres'] for chain in chains]
template_struct='start-structure'
target_name=settings['target_name'].encode('ascii','ignore')
class mymodel(automodel):
	def special_patches(self, aln):
		self.rename_segments(segment_ids=chains,renumber_residues=start_resids)

doalign2d = True
env = environ()
env.io.hetatm = True
env.io.water = False
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters
aln = alignment(env)
model_string=('FIRST:@ ', 'END:')
mdl = model(env,
            file=template_struct,
            model_segment=(model_string))
aln.append_model(mdl,
                 align_codes=template_struct,
                 atom_files=template_struct+'.pdb')
aln.append(file=target_name+'.ali', align_codes=target_name)
aln.align2d()
aln.write(file='align2d.ali', alignment_format='PIR')
aln.write(file='align2d.pap', alignment_format='PAP')
afile = 'align2d.ali'

a = mymodel(env,
	alnfile=afile,
	knowns=template_struct,
	sequence=target_name,
	assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = settings['many_models']
a.make()
