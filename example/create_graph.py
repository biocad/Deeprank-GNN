from deeprank_gnn.GraphGenMP import GraphHDF5

pdb_path = './data/pdb/1ATN/'
pssm_path = './data/pssm/1ATN/'
ref = './data/ref/1ATN/'

pdb_path = './data/pdb/1bql/'
pssm_path = './data/pssm/1bql/'
ref = './data/ref/1bql/'

GraphHDF5(pdb_path=pdb_path, ref_path=ref, pssm_path=pssm_path, biopython=True,
          graph_type='residue', outfile='1bql_residue.hdf5',
          nproc=1, tmpdir='./tmpdir')
