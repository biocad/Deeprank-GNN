import pandas as pd
import json
import requests
#from flatten_json import flatten
from pathlib import Path
import shutil
import re
from tempfile import TemporaryDirectory
from typing import Any, List, Set, Tuple
from pathlib import Path

import numpy as np
from Bio import PDB
from Bio.PDB.Polypeptide import standard_aa_names, three_to_one
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser

def parse_chains(comp_name:str):
    chains= comp_name.split('_')[-1]
    ab_chains, ag_chains = chains.split('-')
    ab_chains_list = ab_chains.split('+')
    ag_chains_list = ag_chains.split('+')
    return ab_chains_list, ag_chains_list

# def get_pssm_dict_by_id(pdb_id):
#     url=f'https://3dcons.cnb.csic.es/pssm_json/{pdb_id}'
#     page = requests.get(url,verify=False)
#     d = json.loads(page.text)
#     df_chain={}
#     for k,v in d.items():
#         df_chain[k]=pd.DataFrame()
#         for row in v:
#             df_chain[k]=pd.concat([df_chain[k],pd.DataFrame.from_dict(flatten(row),orient='index').T],ignore_index=True)
#     return df_chain
    

def post_process_pssm(df_chain,iter_number=2,pss_type='pssm'):
    columns_dict={'pdbresi':'res_id','pdbresn':'aa','seqresi':'index','seqresn':'aa',
    'A': f'iter_{iter_number}_{pss_type}_0', ' R': f'iter_{iter_number}_{pss_type}_1', ' N': f'iter_{iter_number}_{pss_type}_2', ' D': f'iter_{iter_number}_{pss_type}_3', ' C': f'iter_{iter_number}_{pss_type}_4', ' Q':f'iter_{iter_number}_{pss_type}_5', ' E':f'iter_{iter_number}_{pss_type}_6', ' G':f'iter_{iter_number}_{pss_type}_7', ' H':f'iter_{iter_number}_{pss_type}_8', ' I':f'iter_{iter_number}_{pss_type}_9', ' L':f'iter_{iter_number}_{pss_type}_10', ' K':f'iter_{iter_number}_{pss_type}_11', ' M':f'iter_{iter_number}_{pss_type}_12', ' F':f'iter_{iter_number}_{pss_type}_13', ' P':f'iter_{iter_number}_{pss_type}_14', ' S':f'iter_{iter_number}_{pss_type}_15', ' T':f'iter_{iter_number}_{pss_type}_16', ' W':f'iter_{iter_number}_{pss_type}_17', ' Y':f'iter_{iter_number}_{pss_type}_18', ' V':f'iter_{iter_number}_{pss_type}_19',
    'IC':f'iter_{iter_number}_a'}
    df_chain_new={}
    for k_df,df in df_chain.items():
        df_chain_new[k_df]=pd.DataFrame()
        for k,v in columns_dict.items():
            df_chain_new[k_df][k]=df_chain[k_df][v]
    return df_chain_new

def get_ab_chains_name(df,pdb_id):
    return df.loc[df.pdb_id_b==pdb_id.upper(),'ab_chain_ids_b'].values[0].split(':')
    
def get_pssm_matrix(pdb_id,ab_chains,path=None):
    df_chains=get_pssm_dict_by_id(pdb_id)
    df_chains_new=post_process_pssm(df_chains)
    #ab_chains=get_ab_chains_name(df_anbase,pdb_id)
    print(df_chains_new)
    df_antibody=pd.concat([df_chains_new[ab_chains[0]],df_chains_new[ab_chains[0]]])
    df_antigen=df_chains_new[list(set(df_chains_new.keys())-set(ab_chains))[0]]
    print(123)
    print(df_antibody)
    if path is not None:
        if not df_antibody.empty:
            path.mkdir(parents=True,exist_ok=True)
            df_antibody.to_csv(path / f'{pdb_id}.A.pdb.pssm',sep='\t',index=False)
        if not df_antigen.empty:
            path.mkdir(parents=True,exist_ok=True)
            df_antigen.to_csv(path / f'{pdb_id}.B.pdb.pssm',sep='\t',index=False)
    return df_antibody, df_antigen

def reindex_chain(chain):
    for index,residue in enumerate(chain):
        residue._id = (" ", index, " ")
    return chain

def join_chains(chains: List[PDB.Chain.Chain], 
                new_name: str) -> List[PDB.Chain.Chain]:
    if len(chains)==1:
        chains[0]._id = new_name
        return reindex_chain(chains[0])
    else:
        first_chain_len = 1000
        for _ in chains[0]:
            first_chain_len += 1    
        residue_counter = 0
        for chn in chains[1:]:
            for residue in chn:
                residue_counter += 1
                old_id = list(residue.id)
                old_id[1] = first_chain_len + residue_counter
                residue.id = tuple(old_id)
                chains[0].add(residue)
        chains[0]._id = new_name
        return reindex_chain(chains[0])


def init_structure(chains: List[PDB.Chain.Chain]) -> Structure:
    new_struct = PDB.StructureBuilder.StructureBuilder()
    new_struct.init_structure("")
    new_struct.init_model("")
    new_struct = new_struct.get_structure()
    for chain in chains:
        new_struct[""].add(chain)
    
    return new_struct

def get_chains_to_merge(path_to_structure:Path, ab_chain_names:list[str]):
    chains_to_merge_1 = []
    chains_to_merge_2 = []
    parser = PDBParser()
    struct = parser.get_structure("structure", path_to_structure)
    for chain in struct[0].get_chains():
        if chain.id in ab_chain_names:
            chains_to_merge_1.append(chain)
        else:
            chains_to_merge_2.append(chain)
    if chains_to_merge_1:
        return chains_to_merge_1, chains_to_merge_2
    ab_chain_names=['A','B']
    chains_to_merge_2 = []
    for chain in struct[0].get_chains():
        if chain.id in ab_chain_names:
            chains_to_merge_1.append(chain)
        else:
            chains_to_merge_2.append(chain)
    return chains_to_merge_1, chains_to_merge_2


def merge_chains(chains_to_merge_1: List[PDB.Chain.Chain], 
                 chains_to_merge_2: List[PDB.Chain.Chain], 
                 path_to_output: Path):
    """
    Here we merge chains of antibody and write it to temporary PDB file
    """
    merged_1 = join_chains(chains_to_merge_1, new_name="A")
    merged_2 = join_chains(chains_to_merge_2, new_name="B")
    struct = init_structure([merged_1, merged_2])
    io=PDBIO()
    io.set_structure(struct)
    io.save(str(path_to_output))

def merge_ab_chains(source,dest,ab_chains,ag_chains):
    ab_chains,ag_chains=reindex_chain(ab_chains[0]),reindex_chain(ag_chains[0])
    struct = init_structure([ab_chains, ag_chains])
    io=PDBIO()
    io.set_structure(struct)
    io.save(str(dest))

def create_joined(source,dest,ab_chains=None):
    if ab_chains is None:
        path_to_complex=source
        ab_chains,ag_chains=parse_chains(source.stem)
    if len(ab_chains)==1:
        ab_chains_answ=ab_chains
        ab_chains,ag_chains=get_chains_to_merge(source,['A'])
        ab_chains,ag_chains=reindex_chain(ab_chains[0]),reindex_chain(ag_chains[0])
        struct = init_structure([ab_chains, ag_chains])
        io=PDBIO()
        io.set_structure(struct)
        io.save(str(dest))

        return ab_chains_answ
    ab_chains,ag_chains=get_chains_to_merge(source,ab_chains)
    merge_chains(ab_chains,ag_chains,dest)
    return ab_chains

def copy_joined(path,dest):
    for d in path.iterdir():
        pdb_id=d.stem
        print(d.stem)
        (dest/'pdb'/pdb_id).mkdir(parents=True,exist_ok=True)
        #(dest/'pssm'/pdb_id).mkdir(parents=True,exist_ok=True)
        (dest/'ref'/pdb_id).mkdir(parents=True,exist_ok=True)
        for f in d.iterdir():
            if f.stem!='real':
                ab_chains=create_joined(f , dest/ 'pdb'/d.stem / f"{f.stem}.pdb")
                ab_chains_str,_=parse_chains(f.stem)
            else:
                ab_chains=create_joined(f ,  dest/'pdb'/d.stem / f"{pdb_id}_{f.stem.split('_')[-1]}.pdb",ab_chains)
                ab_chains=create_joined(f ,  dest/'ref'/d.stem / f"{pdb_id}.pdb",ab_chains)
        try:
            get_pssm_matrix(pdb_id,ab_chains_str,dest/'pssm'/ pdb_id)
        except Exception as e:
            print(f'Exception on {pdb_id}:{e}')


def get_chains_to_merge(path_to_structure:Path, ab_chain_names:list[str]):
    chains_to_merge_1 = []
    chains_to_merge_2 = []
    parser = PDBParser()
    struct = parser.get_structure("structure", path_to_structure)
    for chain in struct[0].get_chains():
        if chain.id in ab_chain_names:
            chains_to_merge_1.append(chain)
        else:
            chains_to_merge_2.append(chain)
    if chains_to_merge_1:
        return chains_to_merge_1, chains_to_merge_2
    ab_chain_names=['A','B']
    chains_to_merge_2 = []
    for chain in struct[0].get_chains():
        if chain.id in ab_chain_names:
            chains_to_merge_1.append(chain)
        else:
            chains_to_merge_2.append(chain)
    return chains_to_merge_1, chains_to_merge_2




#df_anbase=pd.read_csv('anbase_summary.csv')
#copy_joined(Path('/data/user/shapoval/benchmark'),Path('/data/user/shapoval/joined_benchmark'))
# df_anbase=pd.read_csv('anbase_summary.csv')
# get_pssm_matrix('1bql',df_anbase,True)
