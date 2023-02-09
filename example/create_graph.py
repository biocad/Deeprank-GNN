from functools import partial
from pathlib import Path
from deeprank_gnn.ResidueGraph import ResidueGraph
import pandas as pd
import multiprocessing as mp
from functools import partial
import h5py

def create_and_save_graph(pdb_path,graph_dir,biopython,tmp_dir):
    g=ResidueGraph(pdb_path,biopython=biopython,tmp_dir=tmp_dir)
    (graph_dir / g.name.split('_')[-2]).mkdir(parents=True,exist_ok=True)
    with h5py.File(graph_dir / g.name.split('_')[-2] / f'{g.name}.hdf5', 'w') as f5:
        g.nx2h5(f5)

def get_comp_list(comp_ids,metrics_path):
    comp_list=[]
    for comp in comp_ids:
        comp_list+=pd.read_csv(metrics_path / f'{comp}.csv',index_col=0).index.values.tolist()
    return comp_list

if __name__=='__main__':
    # paths and constants
    train_ids=Path('../train_ids.txt')
    metrics_path=Path('/mnt/volume/metrics')
    complex_path=Path('/mnt/volume/consistent_alpha_hedge')
    graph_dir=Path('./graphs')
    tmp_dir=Path('./tmp_dir')
    nproc=70
    ################
    with open(train_ids) as f:
        comp_ids = [line.rstrip('\n') for line in f]
    comp_list=get_comp_list(comp_ids,metrics_path)
    pdb_paths=[complex_path / comp /'docked.pdb' for comp in comp_list]
    pool = mp.Pool(nproc)
    part_process = partial(create_and_save_graph, graph_dir=graph_dir,biopython=True,tmp_dir=tmp_dir)
    pool.map(part_process, pdb_paths)
