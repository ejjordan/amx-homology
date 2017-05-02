#!/usr/bin/env python

"""
set of tools to facilitate homology modeling
"""

@narrate
def write_ali_file(fasta_linelen=50):

	"""
	Write an ALI file for MODELLER.
	"""

	with open(wordspace.step+wordspace.target_name+'.ali','w') as fp:
		fp.write('>P1;'+wordspace.target_name+'\n')
		fp.write('sequence:'+wordspace.target_name+':::::::0.00:0.00\n')
		seq = wordspace.target_sequence
		chopped = [seq[j*fasta_linelen:(j+1)*fasta_linelen] 
			for j in range(len(seq)/fasta_linelen+1)]
		chopped = [i for i in chopped if len(i) > 0]
		for i,seg in enumerate(chopped): fp.write(seg+'\n')
                if wordspace['other_chains']:
                        fp.write('/\n')
                        for chain in wordspace['other_chains']:
                                seq = wordspace['other_chains_info'][chain]['sequence']
                                chopped = [seq[j*fasta_linelen:(j+1)*fasta_linelen] 
                                           for j in range(len(seq)/fasta_linelen+1)]
                                chopped = [i for i in chopped if len(i) > 0]
                                for i,seg in enumerate(chopped): 
                                        fp.write(seg+'\n')
                fp.write('*\n')

def extract_sequence_backup(filename,chain):
  		  
        """
 	Extract the sequence of a PDB file using biopython
        
        Note that this may fail if there are missing residues - test
        """
 		
 	import Bio		
 	import Bio.PDB		
        import warnings
                

 	parser = Bio.PDB.PDBParser()		
        with warnings.catch_warnings():
                warnings.simplefilter(action="ignore")
                structure = parser.get_structure('this_pdb',filename)		
 	#---extract residue ID and name for non HETATM elements of all chains in the PDB		
 	seqs = {c.id:[(i.id[1],i.resname) 		
 		for i in c.get_residues() if i.id[0]==' '] 		
 		for c in structure.get_chains()}		

        sequence=''.join([aacodemap_3to1[i] for i in zip(*seqs[chain])[1]])	
        startres = int(seqs[chain][0][0])
	start_residue = startres

	if 'start_residue' in wordspace and wordspace['start_residue']: 
		start_residue = int(wordspace['start_residue'])
	stop_residue = len(sequence)+startres
	if 'stop_residue' in wordspace and wordspace['stop_residue']: 
		stop_residue = wordspace['stop_residue']
        if wordspace['sequence_info']:
                wordspace['sequence_info'][chain] = zip(
                    range(start_residue,stop_residue),
                    sequence[start_residue-startres:stop_residue-startres])
        else: 
                wordspace['sequence_info'] = {
                    chain:zip(range(start_residue,stop_residue),
                              sequence[start_residue-startres:stop_residue-startres])}

	return {
		'starting_residue':start_residue,
		'sequence':sequence[start_residue-startres:stop_residue-startres],
		'filename':os.path.basename(filename).split('.pdb')[0]}
 

@narrate
def extract_sequence_pdb(filename,chain):

	"""
	Extract the sequence and staring residue from a PDB file.
	This is a holdover from the original automacs.
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
		missingli = [re.findall('^REMARK\s+([0-9]+)\sMISSING RESIDUES',l)[0] for li,l in enumerate(lines) 
			if re.match('^REMARK\s+([0-9]+)\sMISSING RESIDUES',l)]
		if missingli != []:
			if len(missingli)>1: raise Exception('cannot parse multiple MISSING RESIDUE notes')
			missingli = str(missingli[0])
			startres = int([
				re.findall('^REMARK\s+'+missingli+'\s+[A-Z]{3}\s+[A-Z]\s+([0-9]+)',l)[0] 
				for li,l in enumerate(lines)
				if re.match('^REMARK\s+'+missingli+'\s+[A-Z]{3}\s+[A-Z]\s+[0-9]+',l)][0])
		else: startres = int([line for line in lines if re.match('^ATOM',line)][0][22:25+1])
	elif any([re.match(regex_remark,line) for line in lines]):
		seqresli = [li for li,l in enumerate(lines) if re.match(regex_remark,l)]
		seqraw = [re.findall(regex_remark,lines[li])[0] for li in seqresli]
		sequence = ''.join(seqraw)
		startres = int([line for line in lines if re.match('^ATOM',line)][0][22:25+1])
	else: return False #file must contain either SEQRES or REMARK 300
	start_residue = startres
	if 'start_residue' in wordspace and wordspace['start_residue']: 
		start_residue = int(wordspace['start_residue'])
	stop_residue = len(sequence)+startres
	if 'stop_residue' in wordspace and wordspace['stop_residue']: 
		stop_residue = wordspace['stop_residue']
	wordspace['sequence_info'] = {chain:zip(range(start_residue,stop_residue),
		sequence[start_residue-startres:stop_residue-startres])}
	return {
		'starting_residue':start_residue,
		'sequence':sequence[start_residue-startres:stop_residue-startres],
		'filename':os.path.basename(filename).split('.pdb')[0]}


@narrate
def get_best_structure():

	"""
	Select the structure with the lowest DOPE.
	"""

	with open(wordspace.step+'script-single.log') as fp: lines = fp.readlines()
	regex_log = '^('+wordspace.target_name+'\.[A-Z0-9]+\.pdb)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)'
	results = [re.findall(regex_log,l)[0] for l in lines if re.match(regex_log,l)]
	results = [(i[0],float(i[1]),float(i[2])) for i in results]
	best = sorted(results,key=lambda x:x[1])[0][0]
	with open(wordspace.step+best) as fp: lines = fp.readlines()
	atom_record_start = [ll for ll,l in enumerate(lines) if l[:4]=='ATOM'][0]-1
	seqres = ""
	for chain,details in wordspace['sequence_info'].items():
		seq = [aacodemap_1to3[i] for i in zip(*details)[1]]
		seqlen = len(seq)
		nrows = seqlen/13+(0 if seqlen%13==0 else 1)
		chunks = [seq[i*13:(i+1)*13] for i in range(nrows)]
                additions = ""
		for cnum,chunk in enumerate(chunks):
			additions += 'SEQRES  %-2d  %s %-4d  '%(cnum+1,chain,len(details))+' '.join(chunk)+'\n'
		seqres += additions
	lines.insert(atom_record_start,seqres)
	with open(wordspace.step+wordspace.target_name+'.pdb','w') as fp:
		for line in lines: fp.write(line)
	with open(wordspace.step+'best_structure_path','w') as fp: 
		fp.write(wordspace.target_name+'.pdb'+'\n')
	
@narrate
def get_other_chains(send_atoms=False):
        
        """
        Add back any chains removed during modeling.
        May need to update to deal with more than one chain.
        """
        other_chain_atoms=[]
        wordspace['other_chains_info']={}
        with open(wordspace['start_structure']) as fp: other_chain_lines = fp.readlines()
        for this_chain in wordspace['other_chains']:
                sequence=extract_sequence_pdb(filename=wordspace['start_structure'],chain=this_chain)
                if not sequence:
                        sequence=extract_sequence_backup(filename=wordspace['start_structure'],chain=this_chain)
                this_chain_atoms = [ll for ll,l in enumerate(other_chain_lines) if
                                    l[:4]=='ATOM' and l[21]==this_chain]
                other_chain_atoms.append([other_chain_lines[line] for line in this_chain_atoms])
                wordspace['other_chains_info'][this_chain]={'chain_name':this_chain,
                                                            'starting_residue':
                                                            sequence['starting_residue'],
                                                            'sequence':sequence['sequence']}
        if send_atoms:
                return ''.join([atom for chain in other_chain_atoms for atom in chain])
