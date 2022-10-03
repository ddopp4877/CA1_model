from bmtk.builder import NetworkBuilder
from bmtk.builder.auxi.node_params import positions_list, xiter_random
from bmtk.utils.sim_setup import build_env_bionet
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
from math import exp
import numpy as np
import pandas as pd
import random
import h5py
import synapses

synapses.load()
syn = synapses.syn_params_dicts()

seed = 999

rng = np.random.default_rng(seed)

net = NetworkBuilder("biophysical")
# amount of cells
numAAC = 15#147
numOLM =  16#164
numPV =  55#553
numPyr =  3115#31150

# arrays for cell location csv
cell_name = []
cell_x = []
cell_y = []
cell_z = []
# amount of cells per layer
numAAC_inSO = int(round(numAAC*0.238))
numAAC_inSP = int(round(numAAC*0.7))
numAAC_inSR = int(round(numAAC*0.062))
numPV_inSO = int(round(numPV*0.238))
numPV_inSP = int(round(numPV*0.701))
numPV_inSR = int(round(numPV*0.0596))
#not commented out previously
#totalCellNum = numAAC_inSO + numAAC_inSP + numAAC_inSR + numCCK_inSO + numCCK_inSP + numCCK_inSR + numCCK_inSLM + numNGF_inSR + numNGF_inSLM + numPV_inSO + numPV_inSP + numPV_inSR
totalCellNum = numAAC_inSO + numAAC_inSP + numAAC_inSR +  numPV_inSO + numPV_inSP + numPV_inSR + numPyr + numOLM


def AAC_to_PYR(src, trg, a, x0, sigma, max_dist):
    if src.node_id == trg.node_id:
        return 0

    sid = src.node_id
    tid = trg.node_id

    src_pos = src['positions']
    trg_pos = trg['positions']

    dist = np.sqrt((src_pos[0] - trg_pos[0]) ** 2 + (src_pos[1] - trg_pos[1]) ** 2 + (src_pos[2] - trg_pos[2]) ** 2)
    prob = a * exp(-((dist - x0) ** 2) / (2 * sigma ** 2))
    #print(dist)
    #prob = (prob/100)
    #print(prob)

    if dist <= max_dist:
        global count
        count = count + 1
    if dist <= max_dist and np.random.uniform() < prob:
        connection = 1
        #print("creating {} synapse(s) between cell {} and {}".format(1,sid,tid))
    else:
        connection = 0
    return connection

def n_connections(src, trg, max_dist, prob=0.1):
    """Referenced by add_edges() and called by build() for every source/target pair. For every given target/source
    pair will connect the two with a probability prob (excludes self-connections)"""
    if src.node_id == trg.node_id:
        return 0

    src_pos = src['positions']
    trg_pos = trg['positions']
    dist = np.sqrt((src_pos[0] - trg_pos[0]) ** 2 + (src_pos[1] - trg_pos[1]) ** 2 + (src_pos[2] - trg_pos[2]) ** 2)
    if dist <= max_dist:
        if np.random.uniform() > prob:
            return 0
        else:
            return 1
# total 400x1000x450
# Order from top to bottom is SO,SP,SR,SLM total


def make_layer_grid(xstart,ystart,zstart,x_length,y_length,z_length,min_dist):
    x_grid = np.arange(xstart, x_length+min_dist, min_dist)
    y_grid = np.arange(ystart, y_length+min_dist, min_dist)
    z_grid = np.arange(zstart, z_length+min_dist, min_dist)
    xx, yy, zz = np.meshgrid(x_grid, y_grid, z_grid)
    return np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T


#wrapper function for adding nodes to net, since there are not many params to change each time
def setNodes(netObj,Number,posList,popName,mTemplate):
    inds = rng.choice(np.arange(0, posList.shape[0]), Number, replace=False)

    # Place cell
    netObj.add_nodes(N=Number, pop_name=popName,
                positions=positions_list(positions=posList[inds, :]),
                mem_potential='e',
                model_type='biophysical',
                model_template=mTemplate,
                rotation_angle_zaxis=(np.pi/2), #90 degrees
                morphology=None)

    return np.delete(posList, inds, 0)



def setEdges(netObj,src,dest,conParams,dynamics_name,dist_range,secName,secID,secx):
    conn = netObj.add_edges(source={'pop_name': src}, target={'pop_name': dest},
                     
                     connection_rule=n_connections,
                     connection_params={'prob': conParams[0], 'max_dist': conParams[1]},  
                     delay=0.8,
                     syn_weight = 18,#syn[dynamics_name]['initW_lognormal_mean'],#syn[dynamics_name]['initW_lognormal_mean'],
                     dynamics_params=dynamics_name,
                     model_template=syn[dynamics_name]['level_of_detail'],
                     distance_range=dist_range,
                     target_sections=[secName],
                     sec_id = secID,  # check and is working putting syn on right location
                     sec_x = secx)
    return conn

############################## x,y,z,     xlen,ylen,zlen, space
pos_list_SO = make_layer_grid( 0,0,320,   400,1000,450,  20)
pos_list_SP = make_layer_grid( 0,0,290,   400,1000,320,   8)
pos_list_SR = make_layer_grid( 0,0, 80,   400,1000,290,  20)
pos_list_SLM = make_layer_grid(0,0,  0,   400,1000, 79,  20)

########## add cells and connections to the network

pos_list_SO = setNodes(net,numAAC_inSO,pos_list_SO,'AAC','hoc:axoaxoniccell')
pos_list_SO = setNodes(net,numOLM,pos_list_SO,'OLM','hoc:olmcell')
pos_list_SO = setNodes(net,numPV_inSO,pos_list_SO,'PV','hoc:pvbasketcell')

pos_list_SP = setNodes(net,numPyr,pos_list_SP,'Pyr','hoc:pyramidalcell')
pos_list_SP = setNodes(net,numAAC_inSP,pos_list_SP,'AAC','hoc:axoaxoniccell')
pos_list_SP = setNodes(net,numPV_inSP,pos_list_SP,'PV','hoc:pvbasketcell')

pos_list_SR = setNodes(net,numAAC_inSR,pos_list_SR,'AAC','hoc:axoaxoniccell')
pos_list_SR = setNodes(net,numPV_inSR,pos_list_SR,'PV','hoc:pvbasketcell')

#############        src     dest   prob      maxd     json              dist_range    secName sid secloc
conn1 = setEdges(net,'AAC','Pyr',[ 0.072,     400],'CHN2PN.json', [1,            1.1 ],'axon', 6, 0.5)
conn2 = setEdges(net,'Pyr','AAC',[ 0.009635,  400],'PN2CHN.json', [0.0,        400.0],'apical', 6, 0.5)
conn3 = setEdges(net,'PV','Pyr', [ 0.05366,   400],'PV2PN.json',  [0.0,        400.0],'somatic',0, 0.5)
conn4 = setEdges(net,'Pyr','PV', [ 0.0238,    400],'PV2PN.json',  [0.0,        400.0],'apical', 6, 0.5)
conn5 = setEdges(net,'PV','AAC', [ 0.135,     400],'PV2CHN.json', [0.0,        400.0],'somatic',0, 0.5)
conn6 = setEdges(net,'PV','PV',  [ 0.135,     400],'PV2PV.json',  [0.0,        400.0],'somatic',0, 0.5)
conn7 = setEdges(net,'OLM','Pyr',[ 0.08300,   400],'OLM2PN.json', [0.0,        400.0],'apical', 4, 0.5)
conn8 = setEdges(net,'OLM','AAC',[ 0.0800,    400],'OLM2CHN.json',[0.0,        400.0],'apical', 4, 0.5)
conn9 = setEdges(net,'OLM','PV', [ 0.0800,    400],'OLM2PV.json', [0.0,        400.0],'apical', 4, 0.5)
conn10 = setEdges(net,'OLM','OLM',[ 0.0800,    400],'OLM2OLM.json',[0.0,        400.0],'basal',  0, 0.9)
conn11 = setEdges(net,'Pyr','OLM',[  0.1320,   400],'OLM2OLM.json',[0.0,        400.0],'basal',  2, 0.5)

def lognormal(source, target,m,s):

    mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
    std = np.sqrt(np.log((s / m) ** 2 + 1))
    synaptic_weight = np.random.lognormal(mean, std, 1)
    return synaptic_weight


#connList = [conn1,conn2,conn3,conn4,conn5,conn6,conn7,conn8,conn9,conn10,conn11]
#for conn in connList:
#    conn.add_properties('syn_weight', rule=lognormal,rule_params={'m': syn[conn._edge_type_properties['dynamics_params']]["initW_lognormal_mean"],
#    's': syn[conn._edge_type_properties['dynamics_params']]["initW_lognormal_std"]},  dtypes=np.float64)

net.build()
net.save(output_dir='network')

"""


#use the below on the first build, then add mechanisms directory to the circuit config manually and comment out what is below

build_env_bionet(base_dir='./',
                network_dir='./network',
                config_file='config.json',
                tstop=t_stim, dt=0.1,
                report_vars=['v'],
                components_dir='biophys_components',
                
                spikes_inputs=[('bgpn', 'CA1_inputs/bg_pn_spikes.h5')],

                v_init=-70,
                compile_mechanisms=False)
"""



#print('Number of background spikes to PN cells: {}'.format(psg.n_spikes()))




