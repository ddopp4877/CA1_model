import glob
import json
import os
from bmtk.simulator.bionet.pyfunction_cache import add_synapse_model
from neuron import h
import numpy as np
import random

#Give same syn weight everytime with set seed
seed = 999
random.seed(seed)
np.random.seed(seed)

def Bg2Int(syn_params, sec_x, sec_id):
    """Create a bg2int synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.bg2int(sec_x, sec=sec_id)

    if syn_params.get('initW'):
        m = 0.4
        s = 0.02
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        lsyn.initW = float(np.random.lognormal(mean,std)) # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick() 
    if syn_params.get('taun1'):
        lsyn.taun1 = float(syn_params['taun1'])
    if syn_params.get('taun2'):
        lsyn.taun2 = float(syn_params['taun2'])
    if syn_params.get('gNMDAmax'):
        lsyn.gNMDAmax = float(syn_params['gNMDAmax'])
    if syn_params.get('enmda'):
        lsyn.enmda = float(syn_params['enmda'])
    if syn_params.get('taua1'):
        lsyn.taua1 = float(syn_params['taua1'])
    if syn_params.get('taua2'):
        lsyn.taua2 = float(syn_params['taua2'])
    if syn_params.get('gAMPAmax'):
        lsyn.gAMPAmax = float(syn_params['gAMPAmax'])
    if syn_params.get('eampa'):
        lsyn.eampa = float(syn_params['eampa'])
    return lsyn

def bg2int(syn_params, xs, secs):
    """Create a list of bg2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Bg2Int(syn_params, x, sec)
        syns.append(syn)
    return syns

def Bg2Pyr(syn_params, sec_x, sec_id):
    """Create a bg2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.bg2pyr(sec_x, sec=sec_id)

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    if syn_params.get('taun1'):
        lsyn.taun1 = float(syn_params['taun1'])
    if syn_params.get('taun2'):
        lsyn.taun2 = float(syn_params['taun2'])
    if syn_params.get('gNMDAmax'):
        lsyn.gNMDAmax = float(syn_params['gNMDAmax'])
    if syn_params.get('enmda'):
        lsyn.enmda = float(syn_params['enmda'])
    if syn_params.get('taua1'):
        lsyn.taua1 = float(syn_params['taua1'])
    if syn_params.get('taua2'):
        lsyn.taua2 = float(syn_params['taua2'])
    if syn_params.get('gAMPAmax'):
        lsyn.gAMPAmax = float(syn_params['gAMPAmax'])
    if syn_params.get('eampa'):
        lsyn.eampa = float(syn_params['eampa'])
    return lsyn

def bg2pyr(syn_params, xs, secs):
    """Create a list of bg2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Bg2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns

def Pyr2Chn(syn_params, sec_x, sec_id):
    """Create a pyr2int synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """
    lsyn = h.pyr2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_ampa'):
        lsyn.AlphaTmax_ampa = float(syn_params['AlphaTmax_ampa']) # par.x(21)
    if syn_params.get('Beta_ampa'):
        lsyn.Beta_ampa = float(syn_params['Beta_ampa']) # par.x(22)
    if syn_params.get('Cdur_ampa'):
        lsyn.Cdur_ampa = float(syn_params['Cdur_ampa']) # par.x(23)
    if syn_params.get('gbar_ampa'):
        lsyn.gbar_ampa = float(syn_params['gbar_ampa']) # par.x(24)
    if syn_params.get('Erev_ampa'):
        lsyn.Erev_ampa = float(syn_params['Erev_ampa']) # par.x(16)

    if syn_params.get('AlphaTmax_nmda'):
        lsyn.AlphaTmax_nmda = float(syn_params['AlphaTmax_nmda']) # par.x(25)
    if syn_params.get('Beta_nmda'):
        lsyn.Beta_nmda = float(syn_params['Beta_nmda']) # par.x(26)
    if syn_params.get('Cdur_nmda'):
        lsyn.Cdur_nmda = float(syn_params['Cdur_nmda']) # par.x(27)
    if syn_params.get('gbar_nmda'):
        lsyn.gbar_nmda = float(syn_params['gbar_nmda']) # par.x(28)
    if syn_params.get('Erev_nmda'):
        lsyn.Erev_nmda = float(syn_params['Erev_nmda']) # par.x(16)
        
    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax']) # par.x(16)
    
    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        log_weight = float(np.random.lognormal(mean,std, 1))
        if log_weight >= float(5*m):
            log_weight = float(5*m)
        lsyn.initW = float(log_weight) # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

    if syn_params.get('bACH'):
        lsyn.bACH = float(syn_params['bACH']) # par.x(17)
    if syn_params.get('aDA'):
        lsyn.aDA = float(syn_params['aDA']) # par.x(18)
    if syn_params.get('bDA'):
        lsyn.bDA = float(syn_params['bDA']) # par.x(19)
    if syn_params.get('wACH'):
        lsyn.wACH = float(syn_params['wACH']) # par.x(20)
    
    return lsyn

def pyr2chn(syn_params, xs, secs):
    """Create a list of pyr2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Pyr2Chn(syn_params, x, sec)
        syns.append(syn)
    return syns

def Pyr2Pv(syn_params, sec_x, sec_id):
    """Create a pyr2int synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.pyr2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_ampa'):
        lsyn.AlphaTmax_ampa = float(syn_params['AlphaTmax_ampa'])  # par.x(21)
    if syn_params.get('Beta_ampa'):
        lsyn.Beta_ampa = float(syn_params['Beta_ampa'])  # par.x(22)
    if syn_params.get('Cdur_ampa'):
        lsyn.Cdur_ampa = float(syn_params['Cdur_ampa'])  # par.x(23)
    if syn_params.get('gbar_ampa'):
        lsyn.gbar_ampa = float(syn_params['gbar_ampa'])  # par.x(24)
    if syn_params.get('Erev_ampa'):
        lsyn.Erev_ampa = float(syn_params['Erev_ampa'])  # par.x(16)

    if syn_params.get('AlphaTmax_nmda'):
        lsyn.AlphaTmax_nmda = float(syn_params['AlphaTmax_nmda'])  # par.x(25)
    if syn_params.get('Beta_nmda'):
        lsyn.Beta_nmda = float(syn_params['Beta_nmda'])  # par.x(26)
    if syn_params.get('Cdur_nmda'):
        lsyn.Cdur_nmda = float(syn_params['Cdur_nmda'])  # par.x(27)
    if syn_params.get('gbar_nmda'):
        lsyn.gbar_nmda = float(syn_params['gbar_nmda'])  # par.x(28)
    if syn_params.get('Erev_nmda'):
        lsyn.Erev_nmda = float(syn_params['Erev_nmda'])  # par.x(16)

    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax'])  # par.x(16)
        
    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight)  # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW  # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW  # par.x(2) * lsyn.initW
    # delay = float(syn_params['initW']) # par.x(3) + delayDistance
    # lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1'])  # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2'])  # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1'])  # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2'])  # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1'])  # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1'])  # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2'])  # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2'])  # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF'])  # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f'])  # par.x(15)

    if syn_params.get('bACH'):
        lsyn.bACH = float(syn_params['bACH'])  # par.x(17)
    if syn_params.get('aDA'):
        lsyn.aDA = float(syn_params['aDA'])  # par.x(18)
    if syn_params.get('bDA'):
        lsyn.bDA = float(syn_params['bDA'])  # par.x(19)
    if syn_params.get('wACH'):
        lsyn.wACH = float(syn_params['wACH'])  # par.x(20)

    return lsyn

def pyr2pv(syn_params, xs, secs):
    """Create a list of pyr2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Pyr2Pv(syn_params, x, sec)
        syns.append(syn)
    return syns

def PV2CHN(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.int2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba']) # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba']) # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba']) # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba']) # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba']) # par.x(16)
        
    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax']) # par.x(16)

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight)  # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

    
    return lsyn

def pv2chn(syn_params, xs, secs):
    """Create a list of int2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = PV2CHN(syn_params, x, sec)
        syns.append(syn)
    return syns

def PV2PV(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.int2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba'])  # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba'])  # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba'])  # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba'])  # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba'])  # par.x(16)

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight)  # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW  # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW  # par.x(2) * lsyn.initW
    # delay = float(syn_params['initW']) # par.x(3) + delayDistance
    # lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1'])  # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2'])  # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1'])  # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2'])  # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1'])  # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1'])  # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2'])  # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2'])  # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF'])  # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f'])  # par.x(15)

    return lsyn

def pv2pv(syn_params, xs, secs):
    """Create a list of int2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = PV2PV(syn_params, x, sec)
        syns.append(syn)
    return syns

def Chn2Pyr(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.chn2pyr(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba']) # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba']) # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba']) # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba']) # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba']) # par.x(16)
    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax'])
        

    
    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = m# np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = s#np.sqrt(np.log((s/m)**2 + 1))
        log_weight = float(np.random.lognormal(mean,std, 1)*0.001)
        if log_weight >= float(3*m):
            
            log_weight = float(3*m)
        
            
        lsyn.initW = float(log_weight) # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()
        

    if syn_params.get('Wmax'):
        lsyn.Wmax = lsyn.initW*3#float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = 0.000001#float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

 
    return lsyn

def chn2pyr(syn_params, xs, secs):
    """Create a list of int2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Chn2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns

def PV2Pyr(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """
    lsyn = h.int2pyr(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba']) # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba']) # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba']) # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba']) # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba']) # par.x(16)

    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax']) # par.x(16)

    if  syn_params.get('tau_d_GABAA'):
        lsyn.tau_d_GABAA = float(syn_params['tau_d_GABAA'])

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight) #(lsyn.initW)# par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()
    
    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

    
    return lsyn

def pv2pyr(syn_params, xs, secs):
    """Create a list of int2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = PV2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns

def Pyr2Pyr(syn_params, sec_x, sec_id):
    """Create a pyr2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.pyr2pyr(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_ampa'):
        lsyn.AlphaTmax_ampa = float(syn_params['AlphaTmax_ampa']) # par.x(21)
    if syn_params.get('Beta_ampa'):
        lsyn.Beta_ampa = float(syn_params['Beta_ampa']) # par.x(22)
    if syn_params.get('Cdur_ampa'):
        lsyn.Cdur_ampa = float(syn_params['Cdur_ampa']) # par.x(23)
    if syn_params.get('gbar_ampa'):
        lsyn.gbar_ampa = float(syn_params['gbar_ampa']) # par.x(24)
    if syn_params.get('Erev_ampa'):
        lsyn.Erev_ampa = float(syn_params['Erev_ampa']) # par.x(16)

    if syn_params.get('AlphaTmax_nmda'):
        lsyn.AlphaTmax_nmda = float(syn_params['AlphaTmax_nmda']) # par.x(25)
    if syn_params.get('Beta_nmda'):
        lsyn.Beta_nmda = float(syn_params['Beta_nmda']) # par.x(26)
    if syn_params.get('Cdur_nmda'):
        lsyn.Cdur_nmda = float(syn_params['Cdur_nmda']) # par.x(27)
    if syn_params.get('gbar_nmda'):
        lsyn.gbar_nmda = float(syn_params['gbar_nmda']) # par.x(28)
    if syn_params.get('Erev_nmda'):
        lsyn.Erev_nmda = float(syn_params['Erev_nmda']) # par.x(16)
    
    if syn_params.get('initW'):
        m = 2
        s = 1
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        lsyn.initW = float(np.random.lognormal(mean,std)) # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick() 

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

    if syn_params.get('bACH'):
        lsyn.bACH = float(syn_params['bACH']) # par.x(17)
    if syn_params.get('aDA'):
        lsyn.aDA = float(syn_params['aDA']) # par.x(18)
    if syn_params.get('bDA'):
        lsyn.bDA = float(syn_params['bDA']) # par.x(19)
    if syn_params.get('wACH'):
        lsyn.wACH = float(syn_params['wACH']) # par.x(20)
    
    return lsyn

def pyr2pyr(syn_params, xs, secs):
    """Create a list of pyr2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Pyr2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns

def OLM2PV(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.int2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba'])  # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba'])  # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba'])  # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba'])  # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba'])  # par.x(16)

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight)  # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW  # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW  # par.x(2) * lsyn.initW
    # delay = float(syn_params['initW']) # par.x(3) + delayDistance
    # lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1'])  # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2'])  # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1'])  # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2'])  # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1'])  # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1'])  # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2'])  # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2'])  # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF'])  # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f'])  # par.x(15)

    return lsyn

def olm2pv(syn_params, xs, secs):
    """Create a list of int2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = OLM2PV(syn_params, x, sec)
        syns.append(syn)
    return syns

def OLM2CHN(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.int2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba'])  # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba'])  # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba'])  # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba'])  # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba'])  # par.x(16)

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight)  # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW  # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW  # par.x(2) * lsyn.initW
    # delay = float(syn_params['initW']) # par.x(3) + delayDistance
    # lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1'])  # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2'])  # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1'])  # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2'])  # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1'])  # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1'])  # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2'])  # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2'])  # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF'])  # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f'])  # par.x(15)

    return lsyn

def olm2chn(syn_params, xs, secs):
    """Create a list of int2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = OLM2CHN(syn_params, x, sec)
        syns.append(syn)
    return syns

def OLM2OLM(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.int2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba'])  # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba'])  # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba'])  # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba'])  # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba'])  # par.x(16)
    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax'])  # par.x(16)

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight)  # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW  # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW  # par.x(2) * lsyn.initW
    # delay = float(syn_params['initW']) # par.x(3) + delayDistance
    # lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1'])  # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2'])  # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1'])  # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2'])  # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1'])  # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1'])  # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2'])  # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2'])  # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF'])  # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f'])  # par.x(15)

    return lsyn

def olm2olm(syn_params, xs, secs):
    """Create a list of int2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = OLM2OLM(syn_params, x, sec)
        syns.append(syn)
    return syns

def OLM2Pyr(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """
    lsyn = h.int2pyr(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_gaba'):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba'])  # par.x(21)
    if syn_params.get('Beta_gaba'):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba'])  # par.x(22)
    if syn_params.get('Cdur_gaba'):
        lsyn.Cdur_gaba = float(syn_params['Cdur_gaba'])  # par.x(23)
    if syn_params.get('gbar_gaba'):
        lsyn.gbar_gaba = float(syn_params['gbar_gaba'])  # par.x(24)
    if syn_params.get('Erev_gaba'):
        lsyn.Erev_gaba = float(syn_params['Erev_gaba'])  # par.x(16)
        
    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax'])  # par.x(16)

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(
            log_weight)  # (lsyn.initW)# par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW  # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW  # par.x(2) * lsyn.initW
    # delay = float(syn_params['initW']) # par.x(3) + delayDistance
    # lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1'])  # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2'])  # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1'])  # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2'])  # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1'])  # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1'])  # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2'])  # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2'])  # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF'])  # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f'])  # par.x(15)

    return lsyn

def olm2pyr(syn_params, xs, secs):
    """Create a list of int2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = OLM2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns

def Pyr2OLM(syn_params, sec_x, sec_id):
    """Create a pyr2int synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """
    lsyn = h.pyr2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_ampa'):
        lsyn.AlphaTmax_ampa = float(syn_params['AlphaTmax_ampa'])  # par.x(21)
    if syn_params.get('Beta_ampa'):
        lsyn.Beta_ampa = float(syn_params['Beta_ampa'])  # par.x(22)
    if syn_params.get('Cdur_ampa'):
        lsyn.Cdur_ampa = float(syn_params['Cdur_ampa'])  # par.x(23)
    if syn_params.get('gbar_ampa'):
        lsyn.gbar_ampa = float(syn_params['gbar_ampa'])  # par.x(24)
    if syn_params.get('Erev_ampa'):
        lsyn.Erev_ampa = float(syn_params['Erev_ampa'])  # par.x(16)

    if syn_params.get('AlphaTmax_nmda'):
        lsyn.AlphaTmax_nmda = float(syn_params['AlphaTmax_nmda'])  # par.x(25)
    if syn_params.get('Beta_nmda'):
        lsyn.Beta_nmda = float(syn_params['Beta_nmda'])  # par.x(26)
    if syn_params.get('Cdur_nmda'):
        lsyn.Cdur_nmda = float(syn_params['Cdur_nmda'])  # par.x(27)
    if syn_params.get('gbar_nmda'):
        lsyn.gbar_nmda = float(syn_params['gbar_nmda'])  # par.x(28)
    if syn_params.get('Erev_nmda'):
        lsyn.Erev_nmda = float(syn_params['Erev_nmda'])  # par.x(16)
        
    if syn_params.get('gmax'):
        lsyn.gmax = float(syn_params['gmax'])  # par.x(16)
    if syn_params.get('tau_d_NMDA'):
        lsyn.tau_d_NMDA = float(syn_params['tau_d_NMDA'])  # par.x(16)    
        

    if syn_params.get('initW'):
        m = syn_params.get('initW_lognormal_mean')
        s = syn_params.get('initW_lognormal_std')
        mean = np.log(m) - 0.5 * np.log((s / m) ** 2 + 1)
        std = np.sqrt(np.log((s / m) ** 2 + 1))
        log_weight = float(np.random.lognormal(mean, std, 1))
        if log_weight >= float(5 * m):
            log_weight = float(5 * m)
        lsyn.initW = float(log_weight)  # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick()

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW  # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW  # par.x(2) * lsyn.initW
    # delay = float(syn_params['initW']) # par.x(3) + delayDistance
    # lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1'])  # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2'])  # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1'])  # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2'])  # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1'])  # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1'])  # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2'])  # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2'])  # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF'])  # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f'])  # par.x(15)

    if syn_params.get('bACH'):
        lsyn.bACH = float(syn_params['bACH'])  # par.x(17)
    if syn_params.get('aDA'):
        lsyn.aDA = float(syn_params['aDA'])  # par.x(18)
    if syn_params.get('bDA'):
        lsyn.bDA = float(syn_params['bDA'])  # par.x(19)
    if syn_params.get('wACH'):
        lsyn.wACH = float(syn_params['wACH'])  # par.x(20)

    return lsyn

def pyr2olm(syn_params, xs, secs):
    """Create a list of pyr2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Pyr2OLM(syn_params, x, sec)
        syns.append(syn)
    return syns

def load():
    add_synapse_model(Bg2Pyr, 'bg2pyr', overwrite=False)
    add_synapse_model(Bg2Pyr, overwrite=False)
    add_synapse_model(Bg2Int, 'bg2int', overwrite=False)
    add_synapse_model(Bg2Int, overwrite=False)
    add_synapse_model(Pyr2Pyr, 'pyr2pyr', overwrite=False)
    add_synapse_model(Pyr2Pyr, overwrite=False)
    add_synapse_model(Pyr2Chn, 'pyr2chn', overwrite=False)
    add_synapse_model(Pyr2Chn, overwrite=False)
    add_synapse_model(Pyr2Pv, 'pyr2pv', overwrite=False)
    add_synapse_model(Pyr2Pv, overwrite=False)
    add_synapse_model(PV2Pyr, 'pv2pyr', overwrite=False)
    add_synapse_model(PV2Pyr, overwrite=False)
    add_synapse_model(PV2PV, 'pv2pv', overwrite=False)
    add_synapse_model(PV2PV, overwrite=False)
    add_synapse_model(PV2CHN, 'pv2chn', overwrite=False)
    add_synapse_model(PV2CHN, overwrite=False)
    add_synapse_model(Chn2Pyr, 'chn2pyr', overwrite=False)
    add_synapse_model(Chn2Pyr, overwrite=False)
    add_synapse_model(OLM2PV, 'olm2pv', overwrite=False)
    add_synapse_model(OLM2PV, overwrite=False)
    add_synapse_model(OLM2CHN, 'olm2chn', overwrite=False)
    add_synapse_model(OLM2CHN, overwrite=False)
    add_synapse_model(OLM2OLM, 'olm2olm', overwrite=False)
    add_synapse_model(OLM2OLM, overwrite=False)
    add_synapse_model(OLM2Pyr, 'olm2pyr', overwrite=False)
    add_synapse_model(OLM2Pyr, overwrite=False)
    add_synapse_model(Pyr2OLM, 'pyr2olm', overwrite=False)
    add_synapse_model(Pyr2OLM, overwrite=False)
    return

def syn_params_dicts(syn_dir='biophys_components/synaptic_models'):
    """
    returns: A dictionary of dictionaries containing all
    properties in the synapse json files
    """
    files = glob.glob(os.path.join(syn_dir,'*.json'))
    data = {}
    for fh in files:
        with open(fh) as f:
            data[os.path.basename(fh)] = json.load(f) #data["filename.json"] = {"prop1":"val1",...}
    return data
