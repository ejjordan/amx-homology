#!/usr/bin/env python

import os,sys,re,subprocess,shutil,glob

def make_step(name):
	"""Make a new step folder. See automacs.py docstring for defining variables in the state."""
	#---note any previous states
	if state.before:
		#---TRANSMIT DATA FROM ONE STATE TO THE NEXT
		#---get step information from the last state
		prev_state = state.before[-1]
		state.stepno = prev_state['stepno'] + 1
		state.steps = list(prev_state['steps'])
	else:
		if 'stepno' not in state: state.stepno = 1
		if 'steps' not in state: state.steps = []
		else: state.stepno += 1
	state.step = 's%02d-%s'%(state.stepno,name)
	os.mkdir(state.step)
	state.here = os.path.join(state.step,'')
	#---register a log file for this step.
	#---log files are in the step folder, and use the step name with ".log"
	state.step_log_file = os.path.join(state.here,state.step+'.log')
	state.steps.append(state.step)
	#---! note that files are not recopied on re-run because make-step only runs once
	#---copy all files
	for fn in state.q('files',[]):
		if not os.path.isfile(fn): raise Exception('cannot find an item requested in files: %s'%fn)
		shutil.copyfile(fn,os.path.join(state.here,os.path.basename(fn)))
	#---copy all sources
	for dn in state.q('sources',[]):
		if os.path.basename(dn)=='.': raise Exception('not a valid source: %s'%dn)
		if not os.path.isdir(dn): raise Exception('source %s is not a directory'%dn)
		shutil.copytree(dn,os.path.join(state.here,os.path.basename(dn)))
	#---retrieve files from the file registry
	if state.before and state.file_registry:
		for fn in state.file_registry: 
			shutil.copyfile(prev_state['here']+fn,state.here+fn)

def make_sidestep(name):
	"""Make a new step folder which is off-pathway."""
	if os.path.isfile(name): raise Exception('sidestep %s already exists'%name)
	os.mkdir(name)
	state.sidestep = name
	#---! formalize this 

def copy_file(src,dst,**kwargs):
	"""Wrapper for copying files."""
	cwd = os.path.join(kwargs.get('cwd',state.here),'')
	shutil.copyfile(cwd+src,cwd+dst)
	
def copy_files(src,dst):
	"""Wrapper for copying files with a glob."""
	if not os.path.isdir(dst): raise Exception('destination for copy_files must be a directory: %s'%dst)
	for fn in glob.glob(src): shutil.copy(fn,dst)

def move_file(src,dest,**kwargs):
	"""Wrapper for moving files."""
	cwd = os.path.join(kwargs.get('cwd',state.here),'')
	shutil.move(cwd+src,cwd+dest)

###---GROMACS INTERFACE

def convert_gmx_template_to_call(spec,kwargs,strict=False):
	"""
	Use gromacs call instructions along with kwargs to make a new command.
	"""
	kwargs_copy = dict(kwargs)
	#---check that we have the right incoming arguments
	missing_kwargs = [key for key in spec['required'] if key not in kwargs]
	#---remove required kwargs from the kwargs list
	subs_required = {}
	for key in spec['required']: 
		if key in kwargs: subs_required[key] = kwargs.pop(key)
	#---flags comes in as a list of tuples, but we require a dictionary
	flags = dict(spec['flags'])
	#---remaining kwargs are overrides
	if kwargs and not strict: 
		#---loop over overrides and apply them
		for key,val in kwargs.items(): 
			flags[key if key in flags else '-'+key] = val
	elif kwargs and strict: raise Exception('unprocessed kwargs with strictly no overrides: %s'%str(kwargs))
	#---replace booleans with gromacs -flag vs -noflag syntax
	#---! note that sending e.g. noflag=True is not allowed!
	for key in list(flags.keys()):
		if type(flags[key])==bool: 
			if re.match('^no',key): 
				raise Exception('refusing to accept a kwarg to a gmx flag that starts with "no": '%key)
			if flags[key]: flags[key] = ''
			else: 
				flags[re.sub('^-(.+)$',r'-no\1',key)] = '' 
				del flags[key]
	#---use the flags to construct the call
	call = state.gmxpaths[spec['command']]+' '
	try: call += (' '.join(['%s %s'%(key,val) for key,val in flags.items()]))%subs_required
	except: raise Exception('[ERROR] failed to construct the gromacs call. '+
		'NOTE missing keywords: %s'%missing_kwargs+' NOTE spec: %s'%spec+' NOTE incoming kwargs: %s'%kwargs_copy)
	#---use the explicit specs to make a record of this call
	recorded = {'call':spec['command'],'flags':dict([(flag,
		value%subs_required if type(value)==str else value) for flag,value in flags.items()])}
	return {'call':call,'recorded':recorded}

def get_last_gmx_call(name):
	"""
	The gmx call history is loaded by convert_gmx_template_to_call. We can retrieve the last call to a 
	particular gromacs utility using this function.
	"""
	recents = [ii for ii,i in enumerate(state.history_gmx) if i['call']==name]
	if not recents: raise Exception('no record of a gmx call to %s recently'%name)
	return state.history_gmx[recents[-1]]

def register_gmx_call(command,flag,value):
	"""
	Register an automatic rule for adding flags to gromacs calls.
	"""
	if 'gmx_call_rules' not in state: state.gmx_call_rules = []
	conflicting_rules = [ii for ii,i in enumerate(state.gmx_call_rules) 
		if i['command']==command and i['flag']==flag]
	if any(conflicting_rules):
		print('[NOTE] the rules list is: %s'%conflicting_rules)
		raise Exception('incoming item in gmx_call_rules conflicts with the list (see rules list above): '+
			'command="%s",flag="%s",value="%s"'%(command,flag,value))
	state.gmx_call_rules.append(dict(command=command,flag=flag,value=value))
	#---! somehow this gets propagated to the next step, which is pretty cool. explain how this happens

def register_file(fn):
	"""
	Maintain a list of new, essential files. T
	hese have not specific categories (in contrast to e.g. ITP files)
	"""
	if not state.file_registry: state.file_registry = []
	if not os.path.isfile(state.here+fn):
		raise Exception('cannot register file because it does not exist: %s'%fn)
	if fn in state.file_registry:
		raise Exception('file %s is already in the file registry'%fn)
	state.file_registry.append(fn)

def gmx(program,**kwargs):
	"""
	Construct a GROMACS command and either run it or queue it in a script.
	"""
	#---pop off protected kwargs
	protected,_protected_kwargs = {},['log','inpipe','nonessential','custom']
	for key in _protected_kwargs: protected[key] = kwargs.pop(key,None)
	#---we always require a log
	if not protected['log']: raise Exception(
		'all calls to gmx require a log keyword (recall that we will prepend "log-" to it)')
	#---check that our program template is in the library
	if program not in state.gmxcalls and not protected['custom']:
		raise Exception('gmxcalls has no template for %s'%program)
	elif protected['custom']: program_spec = protected['custom']
	else: program_spec = state.gmxcalls[program]
	#---construct the command from the template
	call_spec = convert_gmx_template_to_call(kwargs=kwargs,spec=program_spec)
	cmd,recorded = call_spec['call'],call_spec['recorded']
	#---check for automatic overrides
	if 'gmx_call_rules' in state:
		#---check and apply all rules with the same command name
		for rule in [i for i in state.gmx_call_rules if i['command']==program]:
			#---ensure no conflicts. the gmx_call_rules will not override flags to gmx
			if '-'+rule['flag'] in recorded['flags']:
				raise Exception('gmx call rule "%s" conflicts with a flag in the call: "%s"'%(
					str(rule),str(recorded)))
			else: 
				#---update the record
				recorded['flags']['-'+rule['flag']] = rule['value']
				#---append to the command
				cmd += ' -%s %s'%(rule['flag'],str(rule['value']))
	#---use the decorator on select functions inside this (amx) module
	#---decorator is only available at run time because it comes from __init__.py
	try: gmx_run_decorated = call_reporter(gmx_run,state)
	except: gmx_run_decorated = gmx_run
	#---! wtf does skip do?
	gmx_run_decorated(cmd,log=protected['log'],
		inpipe=protected['inpipe'],nonessential=protected['nonessential'])
	#---if the run works, we log the completed command 
	if 'history_gmx' not in state: state.history_gmx = []
	state.history_gmx.append(recorded)

def gmx_run(cmd,log,nonessential=False,inpipe=None):

	"""
	Run a GROMACS command instantly and log the results to a file.
	"""

	gmx_error_strings = [
		'File input/output error:',
		'command not found',
		'Fatal error:',
		'Fatal Error:',
		'Can not open file:',
		'Invalid command line argument:',
		'Software inconsistency error']

	if log == None: raise Exception('[ERROR] gmx_run needs a log file to route output')
	#---if the log is an absolute path we drop the log there without prepending "log-"
	elif log and not os.path.basename(log)==log:
		if not os.path.isdir(os.path.dirname(log)): 
			raise Exception('cannot find directory for log %s'%log)
		log_fn = log
	#---local logs get "log-" prepended and drop in the here directory
	else: log_fn = state.here+'log-'+log
	#---previously wrote a bash-only log but it makes more sense to have a comprehensive automacs log
	output = open(log_fn,'w')
	os.chmod(log_fn,0o664)
	if inpipe == None:
		proc = subprocess.Popen(cmd,cwd=state.here,shell=True,executable='/bin/bash',
			stdout=output,stderr=output)
		proc.communicate()
	else:
		proc = subprocess.Popen(cmd,cwd=state.here,shell=True,executable='/bin/bash',
			stdout=output,stderr=output,stdin=subprocess.PIPE)
		proc.communicate(input=str(inpipe).encode())
	#---check for errors
	with open(log_fn,'r') as logfile:
		for line in logfile:
			for msg in gmx_error_strings:
				if re.search(msg,line)!=None: 
					if nonessential: print('[NOTE] command failed but it is nonessential')
					else: raise Exception('%s in %s'%(msg.strip(':'),log_fn))

def get_machine_config(hostname=None):

	"""
	!!!
	"""

	machine_config = {}
	#---!
	with open(os.path.expanduser('~/.automacs.py')) as fp: exec(fp.read(),machine_config)
	#---most of the machine configuration file are headers that are loaded into the main dictionary
	machine_config = machine_config['machine_configuration']
	this_machine = 'LOCAL'
	if not hostname:
		hostnames = [key for key in machine_config 
			if any([varname in os.environ and (
			re.search(key,os.environ[varname])!=None or re.match(key,os.environ[varname]))
			for varname in ['HOST','HOSTNAME']])]
	else: hostnames = [key for key in machine_config if re.search(key,hostname)]
	#---select a machine configuration according to the hostname
	if len(hostnames)>1: raise Exception('[ERROR] multiple machine hostnames %s'%str(hostnames))
	elif len(hostnames)==1: this_machine = hostnames[0]
	else: this_machine = 'LOCAL'
	print('[STATUS] setting gmxpaths for machine: %s'%this_machine)
	machine_config = machine_config[this_machine]
	#---! previously did some ppn calculations here
	return machine_config

def get_gmx_paths(override=False,gmx_series=False):

	"""
	!!!
	"""

	gmx4paths = {'grompp':'grompp','mdrun':'mdrun','pdb2gmx':'pdb2gmx','editconf':'editconf',
		'genbox':'genbox','make_ndx':'make_ndx','genion':'genion','genconf':'genconf',
		'trjconv':'trjconv','tpbconv':'tpbconv','vmd':'vmd','gmxcheck':'gmxcheck'}

	gmx5paths = {'grompp':'gmx grompp','mdrun':'gmx mdrun','pdb2gmx':'gmx pdb2gmx',
		'editconf':'gmx editconf','genbox':'gmx solvate','make_ndx':'gmx make_ndx',
		'genion':'gmx genion','trjconv':'gmx trjconv','genconf':'gmx genconf',
		'tpbconv':'gmx convert-tpr','gmxcheck':'gmx check','vmd':'vmd','solvate':'gmx solvate'}

	#---figure out all the paths shit here. @@@!!!!!
	machine_config = get_machine_config()
	#---basic check for gromacs version series
	suffix = '' if 'suffix' not in machine_config else machine_config['suffix']
	check_gmx = subprocess.Popen('gmx%s'%suffix,shell=True,executable='/bin/bash',
		stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	if override and 'gmx_series' in machine_config: gmx_series = machine_config['gmx_series']
	elif not gmx_series:
		#---! is this the best way to search?
		if not re.search('command not found',str(check_gmx[1])): gmx_series = 5
		else:
			check_mdrun = ' '.join(subprocess.Popen('mdrun%s -g /tmp/md.log'%suffix,shell=True,
				executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate())
			if re.search('VERSION 4',check_mdrun): gmx_series = 4
			elif not override: raise Exception('gromacs is absent')
			else: print('[NOTE] preparing gmxpaths with override')

	if gmx_series == 4: gmxpaths = dict(gmx4paths)
	elif gmx_series == 5: gmxpaths = dict(gmx5paths)
	else: raise Exception('gmx_series must be either 4 or 5')

	#---! need more consistent path behavior here

	#---modify gmxpaths according to hardware configuration
	config = machine_config
	if suffix != '': 
		if gmx_series == 5:
			for key,val in gmxpaths.items():
				gmxpaths[key] = re.sub('gmx ','gmx%s '%suffix,val)
		else: gmxpaths = dict([(key,val+suffix) for key,val in gmxpaths.items()])
	if 'nprocs' in machine_config and machine_config['nprocs'] != None: 
		gmxpaths['mdrun'] += ' -nt %d'%machine_config['nprocs']
	#---use mdrun_command for quirky mpi-type mdrun calls on clusters
	if 'mdrun_command' in machine_config: gmxpaths['mdrun'] = machine_config['mdrun_command']
	#---if any utilities are keys in config we override and then perform uppercase substitutions from config
	utility_keys = [key for key in gmxpaths if key in machine_config]
	if any(utility_keys):
		for name in utility_keys:
			gmxpaths[name] = machine_config[name]
			for key,val in machine_config.items(): 
				gmxpaths[name] = re.sub(key.upper(),str(val),gmxpaths[name])
		del name
	#---even if mdrun is customized in config we treat the gpu flag separately
	if 'gpu_flag' in machine_config: gmxpaths['mdrun'] += ' -nb %s'%machine_config['gpu_flag']	
	#---export the gmxpaths to the state
	state.gmxpaths = gmxpaths