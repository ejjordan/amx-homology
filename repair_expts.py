{
	'repair':{
		'tags':['aamd','protein','homology'],
		'script':'repair.py',
		'params':None,
		'extensions':['homology_tools.py'],
		'settings':"""

		step: repair
		pdb source: none
		start structure: ./inputs/2gs7.pdb
		target name: egfr
		many models: 1
		number HETATMs: {}
		numbering: 'uniprot'

		#---which chains and residues to keep
		template chain:| {
		'A':{'startres':0,'stopres':0},
		'B':{'startres':0,'stopres':0},
		}


"""},
}
