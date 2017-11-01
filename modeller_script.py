#!/usr/bin/env python

"""
modeller homology modeling script
"""

from modeller import *
from modeller.automodel import *
import sys,json

settings_fn = 'settings.json'
with open(settings_fn,'r') as fp: settings = json.load(fp)

refinement=settings['refinement'] #fast, slow
chains=[chain.encode('ascii','ignore') for chain in settings['chains_info']]
start_resids=[settings['chains_info'][chain]['startres'] for chain in chains]
template_struct='start-structure'
target_name=settings['target_name'].encode('ascii','ignore')
chains_info=settings['chains_info']
res_to_model=[str(miss) for chain in chains_info for miss in
			  chains_info[chain]['missing']]

class mymodel(automodel):
	def special_patches(self, aln):
		self.rename_segments(segment_ids=chains,renumber_residues=start_resids)

		
env = environ()
env.io.hetatm = False
env.io.water = False
env.io.atom_files_directory = ['./']
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

afile = '%s.ali'%target_name
a = mymodel(env,
			alnfile=afile,
			knowns='start-structure',
			assess_methods=(assess.DOPE),
			sequence='%s'%target_name)

if refinement=='fast':
	a.library_schedule = autosched.fast
	a.max_var_iterations = 300
	a.md_level = refine.fast
if refinement=='slow':
	a.library_schedule = autosched.slow
	a.max_var_iterations = 300
	a.md_level = refine.slow


a.starting_model = 1
a.ending_model = settings['many_models']

a.make()
