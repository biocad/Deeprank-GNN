import pandas as pd
import json
import requests
from flatten_json import flatten

def get_pssm_dict_by_id(pdb_id):
    url=f'https://3dcons.cnb.csic.es/pssm_json/{pdb_id}'
    page = requests.get(url,verify=False)
    d = json.loads(page.text)
    df_chain={}
    for k,v in d.items():
        df_chain[k]=pd.DataFrame()
        for row in v:
            df_chain[k]=pd.concat([df_chain[k],pd.DataFrame.from_dict(flatten(row),orient='index').T],ignore_index=True)
    return df_chain
    

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
    
def get_pssm_matrix(pdb_id,df_anbase,save=False):
    df_chains=get_pssm_dict_by_id(pdb_id)
    df_chains_new=post_process_pssm(df_chains)
    ab_chains=get_ab_chains_name(df_anbase,pdb_id)
    df_antibody=pd.concat([df_chains_new[ab_chains[0]],df_chains_new[ab_chains[0]]])
    df_antigen=df_chains_new[list(set(df_chains_new.keys())-set(ab_chains))[0]]
    if save:
        df_antibody.to_csv(f'{pdb_id}.A.pdb.pssm',sep='\t',index=False)
        df_antigen.to_csv(f'{pdb_id}.B.pdb.pssm',sep='\t',index=False)
    return df_antibody, df_antigen

df_anbase=pd.read_csv('anbase_summary.csv')
get_pssm_matrix('1bql',df_anbase,True)
