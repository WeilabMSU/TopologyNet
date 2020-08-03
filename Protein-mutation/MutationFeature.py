import os
import math
import pickle
import numpy as np
import sys
from    dionysus        import Filtration, StaticPersistence, data_dim_cmp, vertex_cmp, \
                               fill_alpha3D_complex, fill_alpha2D_complex, points_file
from    sys             import argv, exit
from    math            import sqrt

class Protein:

    def __init__(self, Name):
        self.Name = Name
        self.ResID = []
        self.AtomPos = []
        self.AtomType = []
    def setup(self):
        infile = open(self.Name+'.mutinfo')
        lines = infile.read().splitlines()
        for line in lines:
            a,b,c,d,e = line.split()
            self.ResID.append(int(a))
            self.AtomType.append(b)
            self.AtomPos.append(np.array([float(c), float(d), float(e)]))


def EleID(x, E):
    for i in range(0, len(E)):
        if x == E[i]:
            y = i
    return y

ElementList = ['C','N','O']
rs = 0.25
thr = 12.0
lth = int(np.rint(thr/rs))
small = 0.001
Cut = 12.0
lc = 45
X = np.zeros([lth,lc], float)

res_num = int(sys.argv[3])

ProtWildName = sys.argv[1]
ProtMutName = sys.argv[2]

ProtWild = Protein(ProtWildName)
ProtWild.setup()
ProtMut = Protein(ProtMutName)
ProtMut.setup()

# 0 dimensional features (pairwise interactions)
for j in range(0,2):
    ProtDigit0 = np.zeros([lth,9], float)
    if j == 0: Prot = ProtWild
    if j == 1: Prot = ProtMut
    MutSite = []
    for k in range(0, len(Prot.AtomPos)):
        if Prot.ResID[k] == res_num and Prot.AtomType[k] != 'H':
            MutSite.append((Prot.AtomPos[k], Prot.AtomType[k]))
    for l in range(0, len(MutSite)):
        for k in range(0, len(Prot.AtomPos)):
            if Prot.ResID[k] != res_num and Prot.AtomType[k] in ElementList:
                Dis = np.linalg.norm(Prot.AtomPos[k][:] - MutSite[l][0][:])
                if Dis > Cut: continue
                St = MutSite[l][1]; Pt = Prot.AtomType[k];
                Is = EleID(St, ElementList); Ip = EleID(Pt, ElementList); Ib = int(np.floor(Dis/rs))
                Ie = Is*3+Ip
                ProtDigit0[Ib,Ie] += 1.0
    X[:,j*9:j*9+9] = ProtDigit0[:,:]
X[:,18:27] = X[:,9:18] - X[:,0:9]

# Computing persistent homology
for j in range(0,2):
    if j == 0:
        Prot = ProtWild
        ProtName = ProtWildName
    else:
        Prot = ProtMut
        ProtName = ProtMutName
    pt = []
    for i in range(len(Prot.AtomPos)):
        if Prot.AtomType[i] != 'H':
            pt.append([Prot.AtomPos[i][0], Prot.AtomPos[i][1], Prot.AtomPos[i][2]])
    out_file = open(ProtName+'_All.PH', 'w')
    f = Filtration()
    fill_alpha3D_complex(pt, f)
    f.sort(data_dim_cmp)
    p = StaticPersistence(f)
    p.pair_simplices()
    smap = p.make_simplex_map(f)
    for i in p:
        if i.sign():
            b = smap[i]
            if i.unpaired():
                out_file.write(str(b.dimension())+' '+str(sqrt(b.data[0]))+' '+'inf\n')
                continue
            d = smap[i.pair()]
            out_file.write(str(b.dimension())+' '+str(sqrt(b.data[0]))+' '+str(sqrt(d.data[0]))+'\n')
    out_file.close()

# Feature for Betti-1 and Betti-2
for j in range(0,2):
    ProtDigit12 = np.zeros([lth,6], float)
    if j == 0: PHFile = open(ProtWildName+'_All.PH')
    if j == 1: PHFile = open(ProtMutName+'_All.PH')
    lines = PHFile.read().splitlines()
    for line in lines:
        dim, b, d = line.split()
        dim = int(dim); b = float(b); d = float(d);
        if d - b < small or dim == 0 or b > thr: continue
        # Birth
        intid = int(math.floor(b/rs))
        fid = (dim-1)*3
        ProtDigit12[intid, fid] += 1.0
        # Death
        intid = min(int(math.floor(d/rs)), lth-1)
        fid = (dim-1)*3+1
        ProtDigit12[intid, fid] += 1.0
        # Bar Count
        bintid = int(math.floor(b/rs)); dintid = min(int(math.floor(d/rs)), lth-1);
        fid = (dim-1)*3+2
        ProtDigit12[bintid:dintid+1, fid] += 1.0
    PHFile.close()
    X[:,27+j*6:27+j*6+6] = ProtDigit12[:,:]
X[:,39:45] = X[:,33:39] - X[:,27:33]

OutFile = open('X.npy', 'w')
np.save(OutFile, X)
OutFile.close()
