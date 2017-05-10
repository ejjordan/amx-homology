{

    'homology':{

        'tags':['aamd','protein','homology'],
        'script':'homology.py',
        'params':None,
        'extensions':['homology_tools.py'],
        'settings':"""

        step : homology
        pdb source: None
        start structure: 'inputs/STRUCT.pdb'
        template chain : [['A',0,0]]
        point mutation : {}
        number HETATMs : {}
        target name : alk_active_F1174L
        many models : 1
        """}
}
