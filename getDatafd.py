import os
import django
import pandas as pd
import numpy as np
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true"
django.setup()
from xchem_db.models import *

x = Refinement.objects.filter(outcome__gte=4).filter(outcome__lte=5)
xtal = [i.crystal_name for i in x]

new = {
    'Xtal Name' : [j.crystal_name for j in [i.crystal_name for i in x]],
    'Protein' : [j.target_name for j in [i.target for i in xtal]],
    'Smiles' : [c.smiles for c in [i.compound for i in xtal]],
    'Resolution' : [c.res for c in x],
    'Rfree' :[c.r_free for c in x],
    'lig_confidence': [c.lig_confidence_string for c in x],
    'RMSD_Angles' : [c.rmsd_angles for c in x],  
    'RMSD_bonds' : [c.rmsd_bonds for c in x],     
    'Ramachandran Outliers' : [c.ramachandran_outliers for c in x],
    'CIF' : [c.cif for c in x],
    'Bound Conf' : [c.bound_conf for c in x],  
    'Ligand Bound Conf' : [c.lig_bound_conf for c in x],  
    'Latest PDB' : [c.pdb_latest for c in x],
    'Latest MTZ' : [c.mtz_latest for c in x]
    }
df = pd.DataFrame(new)
df.to_csv('./Data/mock.csv')
