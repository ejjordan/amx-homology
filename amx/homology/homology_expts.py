{

'homology':{
#####
####
###
##
#
'tags':['aamd','protein','homology'],
'script':'homology.py',
'params':None,
'extensions':[],
'settings':"""

step : homology
modeller path : mod9.14
homology method : point
template : inputs/STRUCTURE.pdb
template chain : A
other chains : None
target name : alk_active_F1174L
point mutation : F1174L
target sequence : None
many models : 1
number HETATMs : 0
files : ['inputs/alk_active.pdb']

"""}
}
