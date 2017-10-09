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


def get_pdb_sequence():
	"""
	Get the sequence info if it is a properly constructed PDB file
	"""
	with open(state.here+'start-structure.pdb','r') as fp: pdb=fp.read()
	pdblines=pdb.split('\n')
	#return the one letter code sequence and a dictionary of missing residues
	#---! note this may fail for pdbs from NMR data
	regex_seqres = '^SEQRES\s+[0-9]+\s+([A-Z])\s+[0-9]+\s+(.+)'
	regex_remark = '^REMARK\s465\s{3}\d?\s{1,2}([A-Z]{3})\s([A-Z])\s+(\d+)'
	regex_dbref = '^DBREF\s{2}(\w{4})\s([A-Z])\s+(\d+)\s+(\d+)\s+UNP\s+\w{6}\s+(\w+)\s+(\d+)\s+(\d+)'
	regex_seqadv = '^SEQADV\s(\w{4})\s([A-Z]{3})\s([A-Z])\s+(\d+)\s+UNP\s+(\w{6})\s+(\w{3}\s*\d+\s*)?(\w+)'
	seqinfo={}
	#use 2-step regex to avoid errors from no matchs
	seqresli = [li for li,l in enumerate(pdblines) if re.match(regex_seqres,l)]
	seqraw = [re.findall(regex_seqres,pdblines[li])[0] for li in seqresli]
	missingli = [li for li,l in enumerate(pdblines) if re.match(regex_remark,l)]
	missing_res = [re.findall(regex_remark,pdblines[li])[0] for li in missingli]
	dbrefli = [li for li,l in enumerate(pdblines) if re.match(regex_dbref,l)]
	dbref = [re.findall(regex_dbref,pdblines[li])[0] for li in dbrefli]
	seqadvli = [li for li,l in enumerate(pdblines) if re.match(regex_seqadv,l)]
	seqadv = [re.findall(regex_seqadv,pdblines[li])[0] for li in seqadvli]
	chains=list(set([elt[0] for elt in seqraw]))
	for chain in chains:
		#generate a sequence for each chain from SEQRES remark
		sequence = ''.join([''.join([aacodemap_3to1[j] for j in i[1].split()]) 
							for i in seqraw if i[0] == chain])
		seqinfo[chain]={'missing':{}};seqinfo[chain]['sequence']=sequence
		#make note of any missing residues based on REMARK 465
		missing_res_chain=[res for res in missing_res if res[1]==chain]
		for res in missing_res_chain:
			seqinfo[chain]['missing'][int(res[2])]=res[0]
		#record residue numbering info from DBREF remark
		for ref in dbref:
			if ref[1]==chain:
				seqinfo[chain]['indexinfo']={'cryst_start':int(ref[2]),'cryst_end':int(ref[3]),
											 'uniprot_start':int(ref[5]),'uniprot_end':int(ref[6]),}
		#record mutations or cloning artifacts based on SEQADV remark
		for conflict in seqadv:
			if conflict[2]==chain and 'EXPRESSION' in conflict or 'CLONING' in conflict:
				#these residues will need to be removed to get the numbering right
				if 'artifact' in seqinfo[chain]:
					seqinfo[chain]['artifact'][int(conflict[3])]=conflict[1]
				else:
					seqinfo[chain]['artifact']={}
					seqinfo[chain]['artifact'][int(conflict[3])]=conflict[1]
			if 'ENGINEERED' in conflict and conflict[2]==chain:
				#these residues may need to be mutated to have the 'correct' protein sequence
				#it is left to the user to ensure that the correct sequence is used for any modeling
				uniprot_res,uniprot_idx=conflict[5].split()
				if 'mutated' in seqinfo[chain]:
					seqinfo[chain]['mutated'][int(conflict[3])]={
						'uniprot_id':conflict[4],
						'cryst':{int(conflict[3]):conflict[1]},
						'uniprot':{int(uniprot_idx):uniprot_res}}
				else:
					seqinfo[chain]['mutated']={}
					seqinfo[chain]['mutated'][int(conflict[3])]={
						'uniprot_id':conflict[4],
						'cryst':{int(conflict[3]):conflict[1]},
						'uniprot':{int(uniprot_idx):uniprot_res}}
	return seqinfo


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
	chains=[chain for chain in target_seqs]
	num_template_res=sum([len(target_seqs[chain]['template_seq']) for chain in chains])
	num_target_res=sum([len(target_seqs[chain]['target_seq']) for chain in chains])
	first_chain=sorted(chains)[0];last_chain=sorted(chains)[-1]
	first_cryst_residx=last_cryst_residx=''
	with open(os.path.join(step,'%s.ali'%target_name),'w') as fp:
		fp.write('>P1;start-structure\n')
		fp.write('structureX:start-structure:%s:%s:%s:%s:::0.00:0.00\n'%(first_cryst_residx,
																		 first_chain,
																		 last_cryst_residx,
																		 last_chain))
		for chain_id in target_seqs:
			seq = target_seqs[chain_id]['template_seq']
			chopped = [seq[j*fasta_linelen:(j+1)*fasta_linelen] 
					   for j in range(len(seq)/fasta_linelen+1)]
			chopped = [i for i in chopped if len(i) > 0]
			for i,seg in enumerate(chopped): fp.write(seg+'\n')
			fp.write('/\n')
		fp.write('*\n')
		fp.write('>P1;'+target_name+'\n')
		fp.write('sequence:%s:1:%s:%s:%s:::0.00:0.00\n'%(target_name,first_chain,num_target_res,last_chain))
		for chain_id in target_seqs:
			seq = target_seqs[chain_id]['target_seq']
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
			sequence_info[chain]={'sequence':sequence,
								  'indexinfo':{'cryst_start':startres,'cryst_end':stopres}}
	return sequence_info

def get_all_chains(**state):
  		  
	"""
	prepare the chains from a pbd for writing an ali file
	
	still needs to print a note on mutations
	"""
 	
	settings=state['settings']
	if all(['indexinfo' and 'sequence' in state['pdbinfo'][chain] for chain in state['pdbinfo']]):
		sequence_info=state['pdbinfo']
	else:
		try: sequence_info=sequence_from_atoms()
		except:	raise Exception('\n[ERROR] you must either a complete structure with no ' +
								'gaps or a properly formatted pdb file')

	except_str='\n[ERROR] you must supply a dict of chains you want to include' \
		'in format {"chain_id": {"startres":int, "stopres":int}} but ' \
		'%s was supplied'%settings['template_chain']
	template_seqs={}
	if type(settings['template_chain'])==dict:
		for chain,locs in settings['template_chain'].items():
			# keep only requested residues and use requested numbering for each chain
			if not all(['startres' and 'stopres' in locs]):
				raise Exception(except_str)
			startres=locs['startres']
			stopres=locs['stopres']

			sequence = sequence_info[chain]['sequence']
			#first we need to clip off any residues which are crystalization artifacts
			if 'artifact' in sequence_info[chain]:
				if not 'cryst_start' in sequence_info[chain]['indexinfo']:
					raise Exception('[ERROR] The pdb file you supplied has SEQADV but no DBREF remarks' +
									'please ensure that your pdb file is properly formatted')
				artifacts=sequence_info[chain]['artifact'].keys()
				popped=''.join(['\t%d-%s\n'%(art,sequence_info[chain]['missing'].pop(art))
								for art in artifacts])
				print('[NOTE] removed residues from missing residues list of chain %s because ' \
					  'they are artifacts:\n%s'%(chain,popped))
				min_artifact=min(artifacts)
				max_artifact=max(artifacts)
				if min_artifact<sequence_info[chain]['indexinfo']['cryst_start']:
					before_start=sequence_info[chain]['indexinfo']['cryst_start']-min_artifact
					print('[NOTE] removing residues %s from the begining of chain %s'%(before_start,chain))
				else: before_start=0
				if max_artifact>sequence_info[chain]['indexinfo']['cryst_end']:
					after_end=max_artifact-sequence_info[chain]['indexinfo']['cryst_end']
					print('[NOTE] removing %s residues from the end of chain %s'%(after_end,chain))
				else: after_end=0
				sequence=sequence[before_start:len(sequence)-after_end]

			
			#make sequnce with dashes for use in modeller .ali file
			if 'missing' in sequence_info[chain]:
				missing=sequence_info[chain]['missing']
				if not 'cryst_start' in sequence_info[chain]['indexinfo']:
					raise Exception('[ERROR] The pdb file you supplied has REMARK 465 but no DBREF remark' +
									'please ensure that your pdb file is properly formatted')			
				begin=int(sequence_info[chain]['indexinfo']['cryst_start'])
				lost=sorted(sequence_info[chain]['missing'].keys())
				lost=[gone-begin for gone in lost]
				template_seq=''.join([AA if i not in lost else '-' for i,AA in enumerate(sequence)])
			else: template_seq=sequence

			#we need to find out what numbering scheme the user wants
			if settings['numbering']=='uniprot': first_resnum='uniprot_start'
			elif settings['numbering']=='cryst': first_resnum='cryst_start'
			elif type(settings['numbering'])==dict:	first_resnum=settings['numbering'][chain]
			else: raise Exception('\n[ERROR] "numbering" setting must be either "uniprot", ' +
								  '"cryst", or a dict in format {"chain_id":int}')
			start_idx = int(sequence_info[chain]['indexinfo'][first_resnum])
			#if startres and endres are both zero use the whole sequence
			if startres==0 and stopres==0:
				#need to supply missing residues to modeller in this format
				res_to_model=['%d:%s'%(res+1,chain) for res in lost]
				template_seqs[chain]={'startres':start_idx,'target_seq':sequence,
									  'template_seq':template_seq,'missing':res_to_model}
			elif startres-start_idx>=0 and stopres-start_idx<=len(sequence):
				pre_idx=range(startres-start_idx)
				post_idx=range(stopres-start_idx,len(sequence))
				discard_res=pre_idx+post_idx
				#replace resname with "-" if user requests not to include it
				target_sequence = ''.join([AA if i not in discard_res else '-' for i,AA
										   in enumerate(sequence)])
				#need to supply missing residues to modeller in this format
				res_to_model=['%d:%s'%(res+1,chain) for res in lost if res not in discard_res]
				template_seqs[chain]={'startres':startres,'target_seq':target_sequence,
									  'template_seq':template_seq,'missing':res_to_model}
			else: raise Exception(except_str)

	else: 
		raise Exception(except_str)
	return template_seqs

def get_all_chains_DEP(**state):
  		  
	"""
 	Extract the sequence of a PDB file using biopython
	
	Note that this may fail if there are missing residues - test
	"""
 	
	settings=state['settings']
	try: sequence_info=sequence_from_atoms()
	except:	raise Exception('\n[ERROR] you must point to a complete structure with no ' +
							'gaps [WARN:UPDATE NEEDED]')
	
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
	atom_record_start = [ll for ll,l in enumerate(lines) if l[:4]=='ATOM'][0]
	"""
	hetatm_records = [','.join([line[17:20],line[21],line[22:26]]) for line in lines if line[:6]=='HETATM']
	hetatms={}
	for hetatm in hetatm_records:
		if hetatm in hetatms: hetatms[hetatm]+=1
		else: hetatms[hetatm]=1
	"""
	seqres = ""
	if 'numbering' in settings:
		seqres+='REMARK 6   using %s residue indexing\n'%settings['numbering']
	for chain in state['chains_info']:
		additions = ""
		if 'mutated' in state['pdbinfo'][chain]:
			for mut_info in state['pdbinfo'][chain]['mutated'].values():
				uniprot_id=mut_info['uniprot_id']
				cryst_resname=mut_info['cryst'].values()[0]
				cryst_residx=mut_info['cryst'].keys()[0]
				db_name='UNP'
				db_resname=mut_info['uniprot'].values()[0]
				db_residx=mut_info['uniprot'].keys()[0]
				additions+='SEQADV auto %3s %1s %-4d  %-4s %-9s %-3s %-5d ENGINEERED\n'%(
					cryst_resname,chain,cryst_residx,db_name,uniprot_id,db_resname,db_residx)
			seqres += additions
	for chain,details in state['chains_info'].items():
		additions = ""
		sequence=details['target_seq'].strip('.')
		seq = [aacodemap_1to3[i] for i in sequence]
		seqlen = len(seq)
		nrows = seqlen/13+(0 if seqlen%13==0 else 1)
		chunks = [seq[i*13:(i+1)*13] for i in range(nrows)]
		for cnum,chunk in enumerate(chunks):
			additions += 'SEQRES  %-2d  %s %-4d  '%(cnum+1,chain,len(sequence))+' '.join(chunk)+'\n'
		seqres += additions
	"""
	for het,num in hetatms.items():
		resname,chain_id,resnum=het.split(',')
		seqres += 'HET   %-3s  %s %-4d %5d\n'%(resname.strip(),chain_id,int(resnum),num)
	"""
	lines.insert(atom_record_start,seqres)
	with open(os.path.join(state['step'],settings['target_name']+'.pdb'),'w') as fp:
		for line in lines: fp.write(line)
	with open(os.path.join(state['step'],'best_structure_path'),'w') as fp: 
		fp.write(settings['target_name']+'.pdb'+'\n')
