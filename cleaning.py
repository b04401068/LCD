'''data cleaning
extract metabolic genes
cutting na threshold
filling na'''

import pickle
import pandas as pd
import numpy as np

def save_obj( data, name ):
    with open(name+'.pkl','wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL);
    #np.save( name+'.npy', np.array(data));


def load_obj( name ):
    with open(name+'.pkl','rb') as f:
        return pickle.load(f);

meta_gene = np.load( 'meta_gene.npy');
meta_gene = [ i for i in meta_gene];

dep = pd.read_csv('demeter2.csv', index_col = 0);
meta_frame = dep.reindex([meta_gene]);
print(meta_frame)
save_obj( meta_frame, 'kegg_meta_demeter');


'''
meta_dep_cellless_frame = meta_dep_cellless_frame.dropna(thresh=1000, axis = 1);
meta_dep_cellless_frame = meta_dep_cellless_frame.dropna(thresh=500);

meta_dep_cellless_frame = meta_dep_cellless_frame.fillna(value=0);

meta_nn_dep = meta_dep_cellless_frame.values.tolist();
meta_nn_gene = meta_dep_cellless_frame.index.tolist();

save_obj(meta_nn_dep,'meta_nn_cellless_dep');
save_obj(meta_nn_gene,'meta_nn_cellless_gene');


cell = [];
cell_primary = [];
primary_dict = {};
cell_line = meta_dep_cellless_frame.columns.tolist();
for i in cell_line:
    start = 0;
    start = i.find('_',0);
    cell.append(i[0:start]);
    cell_primary.append(i[start+1:]);
    primary_dict[i[start+1:]] = 0;
primary = [ i for i in primary_dict.keys()];
cell_y = np.zeros((len(primary), len(cell)));

for i in range(len(cell)):
    for j in range(len(primary)):
        if cell_primary[i] == primary[j]:
            break;
    cell_y[j][i] = 1;
save_obj(cell, 'cell_name');

save_obj(cell_y, 'cellless_primary_binary');
save_obj(primary, 'cellless_primary_name');

print(cell_y);
print(cell_y[:,4]);
print(primary);

'''
