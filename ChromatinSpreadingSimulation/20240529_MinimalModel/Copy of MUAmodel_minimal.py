import matplotlib        as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import stats
import seaborn as sns
import pandas as pd
import six
import time
import math

try:
    import cupy as xp
    import numpy as np
except ImportError:
    print('cupy was not imported. numpy instead')
    import numpy as xp
    import numpy as np

def SSA(params,TSTEP,LOOPNUM,DOX,init,loopmat):
    # ------------------------------------------------------ #
    #  Stochastic Simulation Algorithm
    #  import parameters, system size and step size
    #  run Gillespie algorithm on one-dimensional lattice
    #  output is 3D matrix where the dimension is NCELL x NNUC x LOOPNUM
    #  note that NCELL will be divided into 4 different conditions (no spacer, 1.2kb, 3kb or 5kb spacer)
    #  e.g. in case one simulates 1000 cells, 1~250 are no spacer,  251~500 are 1.2 kb,  501~750 are 3 kb,  751~1000 are 5 kb
    # ------------------------------------------------------ #

    # ---- system size ---- #  
    NCELL = init.shape[0]
    NNUC = init.shape[1]

    # ---- define looping probability ---- #  
    loop = loopmat

    # ---- define variables and constants ---- #  
    nucleosome = xp.zeros((NCELL,NNUC),dtype='f')
    nucleosome = init
    wholedynamics = xp.zeros((NCELL,NNUC,LOOPNUM+1),dtype='f')
    wholedynamics[:,:,0] = init
    tau = xp.zeros((NCELL,1),dtype='f')
    delta_n = xp.zeros((NCELL,NNUC),dtype='f')
    reactionrates = xp.zeros((NCELL,NNUC*2),dtype='f')
    hits = xp.zeros((NCELL,NNUC*2),dtype='f')
    cumrrates = xp.zeros((NCELL,NNUC*2),dtype='f')
    cumrrates_prv = xp.zeros((NCELL,NNUC*2),dtype='f')
    cumrrates_prv[:,0] = -1;

    # ---- define reaction rate ---- #  
    kfAA = xp.float32(params[0][0])
    kfAM = xp.float32(params[0][1])
    kfMM = xp.float32(params[0][2])
    kfMA = xp.float32(params[0][3])
    res = ElementSpecificReactionRate(params[1],NCELL,NNUC)
    kbA = res[0]; kdA = res[1]; kbM = res[2]; kdM = res[3];

    # ---- recruitment ---- #  
    if (DOX == 1):
        krec = xp.float32(params[2][0])
        tetpos = xp.zeros((NCELL,NNUC),dtype='f')
        tetpos[:,tet_s] = 1
        tloop = krec * xp.dot(tetpos,loop)
        tloop[:,tet_s] = krec
        if params[2][1] == 1.0:
            kdA = kdA + tloop
        else:
            kbM = kbM + tloop
            kdA = kdA + tloop
    
    # ----- main loop for SSA ----- #
    cntr = 0
    for every in range(1,LOOPNUM+1):
        unfinished = xp.ones((NCELL,1),dtype='f')
        tau = xp.zeros((NCELL,1),dtype='f')
        while (xp.sum(unfinished,axis=None) != 0):

            U = (nucleosome ==  0).astype('f')
            A = (nucleosome == -1).astype('f')
            M = (nucleosome ==  1).astype('f')
            
            Aloop = xp.dot(A,loop)
            Mloop = xp.dot(M,loop)

            syna = U * (kbA + kfAA * Aloop)
            dega = A * (kdA + kfMA * Mloop)
            synm = U * (kbM + kfMM * Mloop)
            degm = M * (kdM + kfAM * Aloop)

            reactionrates[:,0:NNUC] = dega+synm # +1 reaction
            reactionrates[:,NNUC:(NNUC*2)] = syna+degm # -1 reaction
            cumrrates = xp.cumsum(reactionrates,axis=1,dtype='f')
            rn0 = xp.tile(cumrrates[:,-1:] * (1.0 - xp.random.random((NCELL,1))),(1,NNUC*2))
            cumrrates = xp.sign(cumrrates - rn0).astype('f')
            cumrrates_prv[:,1:(NNUC*2)] = cumrrates[:,0:(NNUC*2-1)]
            hits = xp.sign(cumrrates - cumrrates_prv).astype('f')*xp.sign(reactionrates).astype('f')*xp.sign(abs(cumrrates_prv)).astype('f')
            delta_n = hits[:,0:NNUC] - hits[:,NNUC:(NNUC*2)]

            rn = xp.random.random((NCELL,1))
            dtau = xp.sum(abs(delta_n),axis=1).reshape(NCELL,1)*unfinished/xp.sum(reactionrates,axis=1).astype('f').reshape(NCELL,1)*xp.log(1.0/(1.0 - rn))

            tau = tau + dtau
            unfinished = xp.sign(xp.sign(TSTEP - tau).astype('f') + 1.0).astype('f')*xp.sign(xp.sum(nucleosome,axis=1).reshape((NCELL,1)))
            nucleosome = nucleosome + xp.tile(unfinished,(1,NNUC))*delta_n
            cntr = cntr + 1

        wholedynamics[:,:,every] = nucleosome

    return wholedynamics



def DefineElementCoordinate(NNUC):
    # ---- coordinate for each element ---- #
    global center,puro_s,puro_e,tet_s,tet_e,pEF_s,pEF_e,ppp1r12c_s,ppp1r12c_e,Cit_s,Cit_e
    offset = -5 # to set -200bp from the 3' end of pEF as 0 since it is position 0 in CUT&RUN
    center = round(NNUC/2) + offset
    ppp1r12c_s = center-16 + offset
    ppp1r12c_e = center-7 + offset
    tet_s  = center + offset;   tet_e  = center + offset;   pEF_s  = center+1 + offset;   pEF_e  = center+6 + offset;   Cit_s = center+7 + offset;   Cit_e = center+15 + offset;


def ElementSpecificReactionRate(params,NCELL,NNUC):
    # ------------------------------------------------------ #
    #  create element specific reaction rate
    #  params should be in the order of:  
    #  kdA (no reporter region), kdM (no reporter region), kbA (no reporter region), kbM (no reporter region)
    #  kbA (at mCitirne expression cassette), kbA (at mCitirne expression cassette),
    #  kbM (at mCherry expression cassette), kbM (at mCherry expression cassette)
    # ------------------------------------------------------ #

    kbA = xp.ones((NCELL,NNUC)).astype('f')
    kdA = xp.ones((NCELL,NNUC)).astype('f')
    kbM = xp.ones((NCELL,NNUC)).astype('f')
    kdM = xp.ones((NCELL,NNUC)).astype('f')

    # ---- element dependent background modification rate ---- #

    # -- non-reporter region
    kdA[:,:] = params[0]
    kdM[:,:] = params[1]
    kbA[:,:] = params[2]
    kbM[:,:] = params[3]

    # -- reporter region
    kbA[:,pEF_s:pEF_e+1] = params[4]
    kbA[:,Cit_s:Cit_e+1] = params[4]
    kbA[:,ppp1r12c_s:ppp1r12c_e+1] = params[5]

    return kbA,kdA,kbM,kdM



def ModifcationToPromoterActivity(out,ku,ka,km,citidx,mchidx):

    thresh = (ku/ka,km/ka)

    NCELL = out.shape[0]
    TNUM = out.shape[2]
    Quad = xp.zeros((4,out.shape[2])).astype('f')
    outA = (out == -1).astype('f')
    outM = (out ==  1).astype('f')
    outA1 = xp.mean(outA[:,citidx[0]:citidx[1]+1,:],axis=1)
    outM1 = xp.mean(outM[:,citidx[0]:citidx[1]+1,:],axis=1)
    C_Pact = outA1/(thresh[0] + (1.0-thresh[0])*outA1 + (thresh[1]-thresh[0])*outM1)
    outA1 = xp.mean(outA[:,mchidx[0]:mchidx[1]+1,:],axis=1)
    outM1 = xp.mean(outM[:,mchidx[0]:mchidx[1]+1,:],axis=1)
    M_Pact = outA1/(thresh[0] + (1.0-thresh[0])*outA1 + (thresh[1]-thresh[0])*outM1)

    return C_Pact.reshape((NCELL,TNUM,1)), M_Pact.reshape((NCELL,TNUM,1))



def CalcModifcationToPromoterActivityForEachCondition(out,ku,ka,km):
    NCELL = out.shape[0]
    TNUM = out.shape[2]
    citidx = (pEF_s,pEF_e)

    mchidx = (pRSV,mCh_e)
    C1,M1 = ModifcationToPromoterActivity(out[0:xp.int64(NCELL/2),:,:],ku,ka,km,citidx,mchidx)
    mchidx = (pRSV5,mCh_e5)
    C2,M2 = ModifcationToPromoterActivity(out[xp.int64(NCELL/2):NCELL,:,:],ku,ka,km,citidx,mchidx)

    Cout = xp.concatenate((C1,C2), axis=2)
    Mout = xp.concatenate((M1,M2), axis=2)

    return Cout,Mout



def CalcQuad(C_Pact,M_Pact):
    Quad = xp.zeros((4,C_Pact.shape[1],C_Pact.shape[2])).astype('f')
    Con = (C_Pact >= 0.5)
    Mon = (M_Pact >= 0.5)
    Quad[0,:,:] = xp.mean(((Con == 0.0)*(Mon == 1.0)).astype('f'),axis=0)
    Quad[1,:,:] = xp.mean(((Con == 1.0)*(Mon == 1.0)).astype('f'),axis=0)
    Quad[2,:,:] = xp.mean(((Con == 1.0)*(Mon == 0.0)).astype('f'),axis=0)
    Quad[3,:,:] = xp.mean(((Con == 0.0)*(Mon == 0.0)).astype('f'),axis=0)
    return Quad[1,:,:]+Quad[2,:,:], Quad[0,:,:]+Quad[1,:,:], Quad



