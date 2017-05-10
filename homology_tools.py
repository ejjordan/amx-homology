#!/usr/bin/env python

import os,re

"""
set of tools to facilitate homology modeling
"""

dnacodemap_3to1 = {'UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S','UCC':'s','UCA':'S','UCG':'S',
                   'UAU':'Y','UAC':'Y','UAA':'STOP','UAG':'STOP','UGU':'C','UGC':'C','UGA':'STOP',
                   'UGG':'W','CUU':'L','CUC':'L','CUA':'L','CUG':'L','CCU':'P','CCC':'P','CCA':'P',
                   'CCG':'P','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','CGU':'R','CGC':'R','CGA':'R',
                   'CGG':'R','AUU':'I','AUC':'I','AUA':'I','AUG':'M','ACU':'T','ACC':'T','ACA':'T',
                   'ACG':'T','AAU':'N','AAC':'N','AAA':'K','AAG':'K','AGU':'S','AGC':'S','AGA':'R',
                   'AGG':'R','GUU':'V','GUC':'V','GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A',
                   'GCG':'A','GAU':'D','GAC':'D','GAA':'E','GAG':'E','GGU':'G','GGC':'G','GGA':'G',
                   'GGG':'G',} 
	
dnacodemap_1to3 = dict([(val,key) for key,val in dnacodemap_3to1.items()])

aacodemap_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
                  'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X'}

aacodemap_1to3 = dict([(val,key) for key,val in aacodemap_3to1.items()])

def export_modeller_settings(**state):

	"""
	Write a settings file for the modeller script.
	"""
    
	import json
	settings=state['settings']
	settings['chains_info']=state['chains_info']
	settings_fn = os.path.join(state['step'],'settings.json')
	with open(settings_fn,'w') as fp:
		json.dump(settings, fp, ensure_ascii=False)


def write_ali_file(fasta_linelen=50,**state):

	"""
	Write an ALI file for MODELLER.
	"""

	settings=state['settings']
	step=state['step']
	target_name=settings['target_name']
	target_seqs=state['chains_info']
	with open(os.path.join(step,'%s.ali'%target_name),'w') as fp:
		fp.write('>P1;'+target_name+'\n')
		fp.write('sequence:'+target_name+':::::::0.00:0.00\n')
		for chain_id in target_seqs:
			seq = target_seqs[chain_id]['sequence']
			chopped = [seq[j*fasta_linelen:(j+1)*fasta_linelen] 
					   for j in range(len(seq)/fasta_linelen+1)]
			chopped = [i for i in chopped if len(i) > 0]
			for i,seg in enumerate(chopped): fp.write(seg+'\n')
			fp.write('/\n')
		fp.write('*\n')


def sequence_from_pdb():

	"""
	get the sequence, and missing residues, from a properly formatted pdb file
	"""

	with open(filename) as fp: lines = fp.readlines()
	regex_seqres = '^SEQRES\s+[0-9]+\s+([A-Z])\s+[0-9]+\s+(.+)'
	regex_remark = '^REMARK\s300\s([A-Z]+)\s+'
	#---if SEQRES is present we get the sequence from it
	#---note that the seqres protocol below should handle missing residues even if they exist
	#---...at the beginning of the target sequence
	if any([re.match(regex_seqres,line) for line in lines]):
		seqresli = [li for li,l in enumerate(lines) if re.match(regex_seqres,l)]
		seqraw = [re.findall(regex_seqres,lines[li])[0] for li in seqresli]
		sequence = ''.join([''.join([aacodemap_3to1[j] for j in i[1].split()]) 
							for i in seqraw if i[0] == chain])
		missingli = [re.findall('^REMARK\s+([0-9]+)\sMISSING RESIDUES',l)[0]
					 for li,l in enumerate(lines) 
					 if re.match('^REMARK\s+([0-9]+)\sMISSING RESIDUES',l)]
		if missingli != []:
			if len(missingli)>1: raise Exception('cannot parse multiple MISSING RESIDUE notes')
			missingli = str(missingli[0])
			startres = int([re.findall('^REMARK\s+'+missingli+
									   '\s+[A-Z]{3}\s+[A-Z]\s+([0-9]+)',l)[0]
							for li,l in enumerate(lines)
							if re.match('^REMARK\s+'+missingli+
										'\s+[A-Z]{3}\s+[A-Z]\s+[0-9]+',l)][0])
		else: startres = int([line for line in lines if re.match('^ATOM',line)][0][22:25+1])
	elif any([re.match(regex_remark,line) for line in lines]):
		seqresli = [li for li,l in enumerate(lines) if re.match(regex_remark,l)]
		seqraw = [re.findall(regex_remark,lines[li])[0] for li in seqresli]
		sequence = ''.join(seqraw)
		startres = int([line for line in lines if re.match('^ATOM',line)][0][22:25+1])
	else: return False #file must contain either SEQRES or REMARK 300

def sequence_from_atoms():

	"""
	gets the chain ids and resnums from 'start-structure.pdb'
	"""

	import Bio		
	import Bio.PDB		
	import warnings

	parser = Bio.PDB.PDBParser()		
	with warnings.catch_warnings():
		warnings.simplefilter(action="ignore")
		structure = parser.get_structure('this_pdb', os.path.join(
			state['step'],'start-structure.pdb'))
        #---extract residue ID and name for non HETATM elements of all chains in the PDB
        seqs = {c.id:[(i.id[1],i.resname) for i in c.get_residues() if i.id[0]==' ']
                for c in structure.get_chains()}		
        sequence_info={}
        for chain in seqs.keys():
			sequence=''.join([aacodemap_3to1[i] for i in zip(*seqs[chain])[1]])	
			startres = int(seqs[chain][0][0])
			stopres = int(seqs[chain][-1][0])
			sequence_info[chain]={'startres':startres,'stopres':stopres,'sequence':sequence}
	return sequence_info

def get_all_chains(**state):
  		  
	"""
 	Extract the sequence of a PDB file using biopython
	
	Note that this may fail if there are missing residues - test
	"""
 	
	settings=state['settings']
	if 'start_structure' in settings and settings['start_structure']:
		sequence_info=sequence_from_atoms()
	else: raise Exception('\n[ERROR] reading a pdb is not yet implemented')
	
	template_seqs={}
	if type(settings['template_chain'])==list:
		for chain_info in settings['template_chain']:                        
			if not len(chain_info)==3:
				raise Exception('\n[ERROR] you must supply a list of lists of' +
								' chains you want to include in format ' +
								'["chain_id", start_resnum, stop_resnum] ' +
								'e.g. [["A",0,-1],["B",0,-1]]')
			temp_chain=chain_info[0]
			startres=int(chain_info[1])
			stopres=int(chain_info[2])
			start_residue = sequence_info[temp_chain]['startres']
			stop_residue = sequence_info[temp_chain]['stopres']
			sequence = sequence_info[temp_chain]['sequence']
			if startres==0 and stopres==0:
				template_seqs[temp_chain]={'startres':start_residue,'sequence':sequence}
			else:
				template_sequence = sequence[startres-start_residue:stopres-stop_residue]
				template_seqs[temp_chain]={'startres':startres,'sequence':template_sequence}
	else: 
		raise Exception('\n[ERROR] you must supply a list of lists of' +
						' chains you want to include in format ' +
						'["chain_id", start_resnum, stop_resnum] ' +
						'e.g. [["A",100,300],["B",0,0]] but ' +
						'"%s" was supplied'%settings['template_chain'])

	for chain_id in template_seqs:
		if 'number_HETATMs' in settings and chain_id in settings['number_HETATMs']:
			HETATMs=settings['number_HETATMs'][chain_id]
		else: HETATMs=0
		target_sequence = list(template_seqs[chain_id]['sequence']+'.'*HETATMs)
		if 'point_mutation' in settings and chain_id in settings['point_mutation']:
			if not type(settings['point_mutation']==dict):
				raise Exception('\n[ERROR] you must supply a dict of form ' +
								'{"chain_id":[mut1,mut2,..],..} but ' +
								'"%s" was supplied'%settings['point_mutation'])
			for mutation in settings['point_mutation'][chain_id]:
				from_mutation,mutation_location,to_mutation = re.findall(
					'([A-Z])([0-9]+)([A-Z])',mutation)[0]
				mutation_location = int(mutation_location)
				startres=template_seqs[chain_id]['startres']
				if not target_sequence[mutation_location-startres]==from_mutation:
					raise Exception('\n[ERROR] you have asked to change %s at residue ' +
									'%d to %s but this residue is actually ' +
									'%s!'%(from_mutation,mutation_location,to_mutation,
										   target_sequence[mutation_location-startres]))
				target_sequence[mutation_location-startres] = to_mutation
		template_seqs[chain_id]['sequence']=''.join(target_sequence)                                
	return 	template_seqs

def get_best_structure(**state):

	"""
	Select the structure with the lowest DOPE score
	"""

	settings=state['settings']
	with open(os.path.join(state['step'],'log-modeller'),'r') as fp: lines = fp.readlines()
	regex_log = '^('+settings['target_name']+'\.[A-Z0-9]+\.pdb)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)'
	results = [re.findall(regex_log,l)[0] for l in lines if re.match(regex_log,l)]
	results = [(i[0],float(i[1]),float(i[2])) for i in results]
	best = sorted(results,key=lambda x:x[1])[0][0]
	with open(os.path.join(state['step'],best)) as fp: lines = fp.readlines()
	atom_record_start = [ll for ll,l in enumerate(lines) if l[:4]=='ATOM'][0]-1
	hetatm_records = [','.join([line[17:20],line[21],line[22:26]]) for line in lines if line[:6]=='HETATM']
	hetatms={}
	for hetatm in hetatm_records:
		if hetatm in hetatms: hetatms[hetatm]+=1
		else: hetatms[hetatm]=1
	seqres = ""
	for chain,details in state['chains_info'].items():
		sequence=details['sequence'].strip('.')
		seq = [aacodemap_1to3[i] for i in sequence]
		seqlen = len(seq)
		nrows = seqlen/13+(0 if seqlen%13==0 else 1)
		chunks = [seq[i*13:(i+1)*13] for i in range(nrows)]
		additions = ""
		for cnum,chunk in enumerate(chunks):
			additions += 'SEQRES  %-2d  %s %-4d  '%(cnum+1,chain,len(sequence))+' '.join(chunk)+'\n'
		seqres += additions
	for het,num in hetatms.items():
		resname,chain_id,resnum=het.split(',')
		seqres += 'HET   %-3s  %s %-4d %5d\n'%(resname.strip(),chain_id,int(resnum),num)
	lines.insert(atom_record_start,seqres)
	with open(os.path.join(state['step'],settings['target_name']+'.pdb'),'w') as fp:
		for line in lines: fp.write(line)
	with open(os.path.join(state['step'],'best_structure_path'),'w') as fp: 
		fp.write(settings['target_name']+'.pdb'+'\n')
