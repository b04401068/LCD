import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon as poly
from matplotlib.collections import PatchCollection
import math

np.random.seed(1);

def draw_polygon( weight, meta_map):
    
    patches = [];
    colors = [];
    
    for i in meta_map.keys():
        
        #decide the weight for each gene on the map
        if i in weight:
            w = weight[i];
        else:
            w = 0;
        
        #the index on the map
        place = meta_map[i];
        
        #add the index to the polygon
        for index in place:
            x = (np.array(index[0::2]))/1321.0;
            y = (788-np.array(index[1::2]))/788.0;
            coord = np.array([x,y]);
            coord = coord.T;
            polygon = poly( coord, True);
            patches.append(polygon);
            #modify the range of weight
            colors.append(100.0*w);

    #add control to show the map
    x = [0.001,0.002,0.003];
    y = [0.003,0.003,0.002];
    coord = np.array([x,y]);
    coord = coord.T;
    colors.append(0.0);
    polygon = poly( coord, True);
    patches.append(polygon);
    colors = np.array(colors);

    #if pos is True:
    #    p = PatchCollection(patches, cmap=plt.cm.Greys_r, alpha = 0.5);
    #else:
    #    p = PatchCollection(patches, cmap=plt.cm.Greys, alpha = 0.5);
    
    #add the plygon, set red blue map, alpha is the parameter for transparency of polygon in the same area
    p = PatchCollection(patches, cmap=plt.cm.RdBu_r, alpha = 0.5);
    
    p.set_array(np.array(colors));
    low = np.min(colors);
    high = np.max(colors);
    if np.abs(low) < np.abs(high):
        bar = np.abs(high);
    else:
        bar = np.abs(low);
    
    #set the clim
    p.set_clim(vmin = -bar, vmax=bar);
    
    #print(colors.shape);
    
    return p;


def draw_activation( cell_weight, name, meta_map, title):
    #total = len(cell_weight);
    #total = int(math.sqrt(total));
    fig, axes = plt.subplots(1,1,figsize=(13.21,7.88));
    if total > 1:
        ax = axes.flatten();
    else:
        ax = [axes];
    p = [];
    
    for i in cell_weight:
        #i[2] is weight
        p.append( draw_polygon(i[2], meta_map));
    for i in range(len(ax)):
        ax[i].add_collection(p[i]);
        ax[i].axis('off');
        ax[i].set_title(cell_weight[i][0]+' '+cell_weight[i][1]);

    #if pos is True:
    #    fig.patch.set_facecolor('black');
    fig.suptitle(title);
    #if pos is True:
    #    fig.savefig(name+'_pos.png', facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight');
    #else:
    #    fig.savefig(name+'_neg.png', bbox_inches='tight');
    fig.savefig(name+'.png', facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight');
    plt.close('all')

def pick_n( a2, y, n):
    Y = np.zeros(a2.shape);
    Y = Y + y;
    penalty = np.ones(a2.shape);
    penalty = penalty + y * 10;
    loss_2 = np.multiply((a2 - Y), penalty);
    loss_2 = np.sum(np.multiply(loss_2,loss_2),axis = 0, keepdims=True);
    index = np.argsort(loss_2).flatten()[0:n];
    #index = index.flatten();
    #index = np.argsort(loss_2[index]);
    return np.int64(index)

def pick_n_hidden( a1, y, n):
    #penalty = penalty + y * 10;
    a1 = a1[y,:];
    index = np.argsort(a1).flatten()[-n:];
    #index = index.flatten();
    #index = np.argsort(loss_2[index]);
    return np.int64(index)

def make_map( act, cell, primary, gene):
    cell_weight = [];
    for i in range(len(cell)):
        cell_act = act[:,i].flatten();
        cell_dic = dict( zip( gene,cell_act ));
        cell_weight.append( ( cell[i], primary[i], cell_dic ) )
    return cell_weight

meta_map = np.load('./data/meta_map.npy');
meta_map = meta_map.item();
gene = np.load("./data/meta_nn_cellless_gene.npy");

#primary_index = np.argmax(primary_binary, axis = 0);
#cell_primary = np.array([primary[i] for i in primary_index]);

primary_parameters = primary_parameters.item();
tp53_parameters = tp53_parameters.item();
gene = gene.flatten();

for i in []:
    #y = np.zeros((len(primary),1));
    #y[i] = 1;
    
    #index_ori = pick_n(a2, i, 9);
    cell_weight = make_map( dep[:,first], cell[first], cell_primary[first], gene);
    draw_activation( cell_weight, fd+nam+'_'+str(i)+'_top1_dep', meta_map, nam+'_'+str(i)+'_dep');
