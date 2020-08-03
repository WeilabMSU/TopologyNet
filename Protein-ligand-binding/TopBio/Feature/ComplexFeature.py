import numpy as np
import pickle
import os

def GenerateFeature_interaction_ML(protein_name, working_dir,typ):

    if typ == 'pair':
        ProEleList = ['C','N','O','S']
        LigEleList = ['C','N','O','S','P','F','Cl','Br','I']
    elif typ == 'all':
        ProEleList = ['C','N','O','S','CN','CO','NO','CNO']
        LigEleList = ['C','N','O','S','P','F','Cl','Br','I','CN','CO','CS','NO','NS','OS','CNO','CNS','COS','NOS','CNOS']
    Bins = [0., 2.5, 3.0, 3.5, 4.5, 6.0, 12.0]
    def BinID(x, B):
        for iii in range(0,len(B)-1):
            if x >= B[iii] and x <= B[iii+1]:
                y = iii
        return y
    cut = 12.0

    Feature = np.zeros([len(ProEleList), len(LigEleList), len(Bins)-1], float)

    for j in range(len(ProEleList)):
        ep = ProEleList[j]
        for k in range(len(LigEleList)):
            el = LigEleList[k]
            if not os.path.exists(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_50.0_interaction.PH'):
                continue
            InFile = open(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_50.0_interaction.PH')
            lines = InFile.read().splitlines()
            for line in lines:
                dim,b,d = line.split()
                if float(d) >= cut:
                    continue
                did = BinID(float(d),Bins)
                Feature[j,k,did] += 1.0

    FeatureFlat = np.zeros([len(ProEleList)*len(LigEleList)*(len(Bins)-1)], float)
    for i in range(len(ProEleList)):
        for j in range(len(LigEleList)):
            for k in range(len(Bins)-1):
                fid = i*len(LigEleList)*(len(Bins)-1) + j*(len(Bins)-1) + k
                FeatureFlat[fid] = Feature[i,j,k]

    OutFile = open(working_dir+'/'+protein_name+'_feature_complex_interaction_ML.npy', 'w')
    np.save(OutFile, FeatureFlat)
    OutFile.close()

def GenerateFeature_interaction_1DCNN(protein_name, working_dir):
    ProEleList = ['C','N','O','S']
    LigEleList = ['C','N','O','S','P','F','Cl','Br','I']

    Bins = []
    for i in range(201):
        Bins.append(float(i)*0.25)

    def BinID(x, B):
        for iii in range(0,len(B)-1):
            if x >= B[iii] and x <= B[iii+1]:
                y = iii
        return y
    cut = 50.0

    Feature = np.zeros([len(ProEleList), len(LigEleList), len(Bins)-1], float)

    for j in range(len(ProEleList)):
        ep = ProEleList[j]
        for k in range(len(LigEleList)):
            el = LigEleList[k]
            if not os.path.exists(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_50.0_interaction.PH'):
                continue
            InFile = open(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_50.0_interaction.PH')
            lines = InFile.read().splitlines()
            for line in lines:
                dim,b,d = line.split()
                if float(d) >= cut:
                    continue
                did = BinID(float(d),Bins)
                Feature[j,k,did] += 1.0

    FeatureFlat = np.zeros([len(Bins)-1, len(ProEleList)*len(LigEleList)], float)
    for i in range(len(ProEleList)):
        for j in range(len(LigEleList)):
            fid = i*len(LigEleList)+j
            for k in range(len(Bins)-1):
                FeatureFlat[k,fid] = Feature[i,j,k]

    OutFile = open(working_dir+'/'+protein_name+'_feature_complex_interaction_1DCNN.npy', 'w')
    np.save(OutFile, FeatureFlat)
    OutFile.close()

def GenerateFeature_electrostatics_ML(protein_name, working_dir):
    Cut = 2.0

    LIGELE = ['C','N','O','S','P','F','Cl','Br','I','H']
    PROELE = ['C','N','O','S','H']

    Feature_i = []
    for ep in PROELE:
        for el in LIGELE:
            if os.path.exists(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_16.0_chg.PH'):
                InFile = open(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_16.0_chg.PH')
                Bars = []
                lines = InFile.read().splitlines()
                for line in lines:
                    a,b,c = line.split()
                    Bars.append([int(a), float(b), float(c)])
                InFile.close()
                Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
                for Bar in Bars:
                    if Bar[2] < Bar[1]:
                        LogFile.write(pdb+' '+str(Bar[0])+' '+str(Bar[1])+' '+str(Bar[2])+'\n')
                        continue
                    if Bar[2] > Cut and Bar[0] == 0: continue
                    if Bar[0] == 0:
                        Bar0Birth.append(Bar[1])
                        Bar0Death.append(Bar[2])
                if len(Bar0Birth) > 0:
                    Bar0Birth = np.asarray(Bar0Birth, float)
                    Bar0Death = np.asarray(Bar0Death, float)
                if len(Bar0Death) > 0:
                    Feature_i.append(np.mean(Bar0Death[:]))
                    Feature_i.append(np.std(Bar0Death[:]))
                    Feature_i.append(np.max(Bar0Death[:]))
                    Feature_i.append(np.min(Bar0Death[:]))
                    Feature_i.append(np.sum(Bar0Death[:]))
                    Feature_i.append(len(Bar0Death))
                else:
                    Feature_i.extend([0.]*6)
            else:
                Feature_i.extend([0.]*6)
    Fearture_i = np.asarray(Feature_i, float)

    OutFile = open(working_dir+'/'+protein_name+'_feature_complex_electrostatics_ML.npy', 'w')
    np.save(OutFile, Feature_i)

def GenerateFeature_electrostatics_1DCNN(protein_name, working_dir):
    ProEleList = ['C','N','O','S','H']
    LigEleList = ['C','N','O','S','P','F','Cl','Br','I','H']

    Bins = []
    for i in range(101):
        Bins.append(float(i)*0.01)

    def BinID(x, B):
        for iii in range(0,len(B)-1):
            if x >= B[iii] and x <= B[iii+1]:
                y = iii
        return y
    cut = 1.0

    Feature = np.zeros([len(ProEleList), len(LigEleList), len(Bins)-1], float)

    for j in range(len(ProEleList)):
        ep = ProEleList[j]
        for k in range(len(LigEleList)):
            el = LigEleList[k]
            if not os.path.exists(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_16.0_chg.PH'):
                continue
            InFile = open(working_dir+'/'+protein_name+'_'+ep+'_'+el+'_16.0_chg.PH')
            lines = InFile.read().splitlines()
            for line in lines:
                dim,b,d = line.split()
                if float(d) >= cut:
                    continue
                did = BinID(float(d),Bins)
                Feature[j,k,did] += 1.0

    FeatureFlat = np.zeros([len(Bins)-1, len(ProEleList)*len(LigEleList)], float)

    for i in range(len(ProEleList)):
        for j in range(len(LigEleList)):
            fid = i*len(LigEleList)+j
            for k in range(len(Bins)-1):
                FeatureFlat[k,fid] = Feature[i,j,k]

    OutFile = open(working_dir+'/'+protein_name+'_feature_complex_electrostatics_1DCNN.npy', 'w')
    np.save(OutFile, FeatureFlat)
    OutFile.close()

def GenerateFeature_alpha_ML(protein_name, working_dir, typ):
    Cut = 12.0

    if typ == 'carbon':
        ProEleName = ['C']
        LigEleName = ['C']
    elif typ == 'heavy':
        ProEleName = ['CNOS']
        LigEleName = ['CNOSPFClBrI']

    small = 0.01

    Feature_i = []
    InFile = open(working_dir+'/'+protein_name+'_alpha.pkl')
    BarCollection = pickle.load(InFile)
    for ep in ProEleName:
        if 'pro_'+ep in BarCollection.keys():
            Bars = BarCollection['pro_'+ep]
            Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
            for Bar in Bars:
                if Bar[2] > 12.0: Bar[2] = 12.0
                if Bar[2] - Bar[1] < small: continue
                if Bar[0] == 0:
                    Bar0Birth.append(Bar[1])
                    Bar0Death.append(Bar[2])
                if Bar[0] == 1:
                    Bar1Birth.append(Bar[1])
                    Bar1Death.append(Bar[2])
                if Bar[0] == 2:
                    Bar2Birth.append(Bar[1])
                    Bar2Death.append(Bar[2])
            if len(Bar0Birth) > 0:
                Bar0Birth = np.asarray(Bar0Birth, float)
                Bar0Death = np.asarray(Bar0Death, float)
            if len(Bar1Birth) > 0:
                Bar1Birth = np.asarray(Bar1Birth, float)
                Bar1Death = np.asarray(Bar1Death, float)
            if len(Bar2Birth) > 0:
                Bar2Birth = np.asarray(Bar2Birth, float)
                Bar2Death = np.asarray(Bar2Death, float)
            if len(Bar0Death) > 0:
                Feature_i.append(np.mean(Bar0Death[:]))
                Feature_i.append(np.std(Bar0Death[:]))
                Feature_i.append(np.max(Bar0Death[:]))
                Feature_i.append(np.min(Bar0Death[:]))
                Feature_i.append(np.sum(Bar0Death[:]))
                Feature_i.append(len(Bar0Death))
            else:
                Feature_i.extend([0.]*6)
            if len(Bar1Death) > 0:
                Feature_i.append(np.mean(Bar1Death[:] - Bar1Birth[:]))
                Feature_i.append(np.std(Bar1Death[:] - Bar1Birth[:]))
                Feature_i.append(np.max(Bar1Death[:] - Bar1Birth[:]))
                Feature_i.append(np.min(Bar1Death[:] - Bar1Birth[:]))
                Feature_i.append(np.sum(Bar1Death[:] - Bar1Birth[:]))
                Feature_i.append(Bar1Birth[np.argmax(Bar1Death[:] - Bar1Birth[:])])
                Feature_i.append(Bar1Death[np.argmax(Bar1Death[:] - Bar1Birth[:])])
                Feature_i.append(np.mean(Bar1Birth[:]))
                Feature_i.append(np.std(Bar1Birth[:]))
                Feature_i.append(np.max(Bar1Birth[:]))
                Feature_i.append(np.min(Bar1Birth[:]))
                Feature_i.append(np.sum(Bar1Birth[:]))
                Feature_i.append(np.mean(Bar1Death[:]))
                Feature_i.append(np.std(Bar1Death[:]))
                Feature_i.append(np.max(Bar1Death[:]))
                Feature_i.append(np.min(Bar1Death[:]))
                Feature_i.append(np.sum(Bar1Death[:]))
                Feature_i.append(len(Bar1Death))
            else:
                Feature_i.extend([0.]*18)
            if len(Bar2Death) > 0:
                Feature_i.append(np.mean(Bar2Death[:] - Bar2Birth[:]))
                Feature_i.append(np.std(Bar2Death[:] - Bar2Birth[:]))
                Feature_i.append(np.max(Bar2Death[:] - Bar2Birth[:]))
                Feature_i.append(np.min(Bar2Death[:] - Bar2Birth[:]))
                Feature_i.append(np.sum(Bar2Death[:] - Bar2Birth[:]))
                Feature_i.append(Bar2Birth[np.argmax(Bar2Death[:] - Bar2Birth[:])])
                Feature_i.append(Bar2Death[np.argmax(Bar2Death[:] - Bar2Birth[:])])
                Feature_i.append(np.mean(Bar2Birth[:]))
                Feature_i.append(np.std(Bar2Birth[:]))
                Feature_i.append(np.max(Bar2Birth[:]))
                Feature_i.append(np.min(Bar2Birth[:]))
                Feature_i.append(np.sum(Bar2Birth[:]))
                Feature_i.append(np.mean(Bar2Death[:]))
                Feature_i.append(np.std(Bar2Death[:]))
                Feature_i.append(np.max(Bar2Death[:]))
                Feature_i.append(np.min(Bar2Death[:]))
                Feature_i.append(np.sum(Bar2Death[:]))
                Feature_i.append(len(Bar2Death))
            else:
                Feature_i.extend([0.]*18)
            # Complex and difference
            Feature_i_p = Feature_i[-42:]
            for el in LigEleName:
                if 'com_'+ep+'_'+el in BarCollection.keys():
                    Bars = BarCollection['com_'+ep+'_'+el]
                    Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
                    for Bar in Bars:
                        if Bar[2] > 12.0: Bar[2] = 12.0
                        if Bar[2] - Bar[1] < small: continue
                        if Bar[0] == 0:
                            Bar0Birth.append(Bar[1])
                            Bar0Death.append(Bar[2])
                        if Bar[0] == 1:
                            Bar1Birth.append(Bar[1])
                            Bar1Death.append(Bar[2])
                        if Bar[0] == 2:
                            Bar2Birth.append(Bar[1])
                            Bar2Death.append(Bar[2])
                    if len(Bar0Birth) > 0:
                        Bar0Birth = np.asarray(Bar0Birth, float)
                        Bar0Death = np.asarray(Bar0Death, float)
                    if len(Bar1Birth) > 0:
                        Bar1Birth = np.asarray(Bar1Birth, float)
                        Bar1Death = np.asarray(Bar1Death, float)
                    if len(Bar2Birth) > 0:
                        Bar2Birth = np.asarray(Bar2Birth, float)
                        Bar2Death = np.asarray(Bar2Death, float)
                    if len(Bar0Death) > 0:
                        Feature_i.append(np.mean(Bar0Death[:]))
                        Feature_i.append(np.std(Bar0Death[:]))
                        Feature_i.append(np.max(Bar0Death[:]))
                        Feature_i.append(np.min(Bar0Death[:]))
                        Feature_i.append(np.sum(Bar0Death[:]))
                        Feature_i.append(len(Bar0Death))
                    else:
                        Feature_i.extend([0.]*6)
                    if len(Bar1Death) > 0:
                        Feature_i.append(np.mean(Bar1Death[:] - Bar1Birth[:]))
                        Feature_i.append(np.std(Bar1Death[:] - Bar1Birth[:]))
                        Feature_i.append(np.max(Bar1Death[:] - Bar1Birth[:]))
                        Feature_i.append(np.min(Bar1Death[:] - Bar1Birth[:]))
                        Feature_i.append(np.sum(Bar1Death[:] - Bar1Birth[:]))
                        Feature_i.append(Bar1Birth[np.argmax(Bar1Death[:] - Bar1Birth[:])])
                        Feature_i.append(Bar1Death[np.argmax(Bar1Death[:] - Bar1Birth[:])])
                        Feature_i.append(np.mean(Bar1Birth[:]))
                        Feature_i.append(np.std(Bar1Birth[:]))
                        Feature_i.append(np.max(Bar1Birth[:]))
                        Feature_i.append(np.min(Bar1Birth[:]))
                        Feature_i.append(np.sum(Bar1Birth[:]))
                        Feature_i.append(np.mean(Bar1Death[:]))
                        Feature_i.append(np.std(Bar1Death[:]))
                        Feature_i.append(np.max(Bar1Death[:]))
                        Feature_i.append(np.min(Bar1Death[:]))
                        Feature_i.append(np.sum(Bar1Death[:]))
                        Feature_i.append(len(Bar1Death))
                    else:
                        Feature_i.extend([0.]*18)
                    if len(Bar2Death) > 0:
                        Feature_i.append(np.mean(Bar2Death[:] - Bar2Birth[:]))
                        Feature_i.append(np.std(Bar2Death[:] - Bar2Birth[:]))
                        Feature_i.append(np.max(Bar2Death[:] - Bar2Birth[:]))
                        Feature_i.append(np.min(Bar2Death[:] - Bar2Birth[:]))
                        Feature_i.append(np.sum(Bar2Death[:] - Bar2Birth[:]))
                        Feature_i.append(Bar2Birth[np.argmax(Bar2Death[:] - Bar2Birth[:])])
                        Feature_i.append(Bar2Death[np.argmax(Bar2Death[:] - Bar2Birth[:])])
                        Feature_i.append(np.mean(Bar2Birth[:]))
                        Feature_i.append(np.std(Bar2Birth[:]))
                        Feature_i.append(np.max(Bar2Birth[:]))
                        Feature_i.append(np.min(Bar2Birth[:]))
                        Feature_i.append(np.sum(Bar2Birth[:]))
                        Feature_i.append(np.mean(Bar2Death[:]))
                        Feature_i.append(np.std(Bar2Death[:]))
                        Feature_i.append(np.max(Bar2Death[:]))
                        Feature_i.append(np.min(Bar2Death[:]))
                        Feature_i.append(np.sum(Bar2Death[:]))
                        Feature_i.append(len(Bar2Death))
                    else:
                        Feature_i.extend([0.]*18)
                    Feature_i.extend([i-j for i,j in zip(Feature_i[-42:], Feature_i_p[:])])
                else:
                    Feature_i.extend(Feature_i[-42:])
                    Feature_i.extend([0.]*42)
        else:
            Feature_i.extend([0.]*(126))

    Feature_i = np.asarray(Feature_i, float)
    OutFile = open(working_dir+'/'+protein_name+'_feature_complex_alpha_'+typ+'_ML.npy', 'w')
    np.save(OutFile, Feature_i)
    OutFile.close()

def GenerateFeature_alpha_1DCNN(protein_name, working_dir):
    Cut = 12.0

    small = 0.01

    rs = 0.125
    lth = int(np.round(Cut/rs))

    Feature_h_p = np.zeros([lth,8], float)
    Feature_c_p = np.zeros([lth,8], float)
    Feature_h_pl = np.zeros([lth,8], float)
    Feature_c_pl = np.zeros([lth,8], float)
    InFile = open(working_dir+'/'+protein_name+'_alpha.pkl')
    BarCollection = pickle.load(InFile)


    # h_p------------------------------------------------------
    ProEleName = ['CNOS']
    LigEleName = ['CNOSPFClBrI']
    ep = 'CNOS'
    el = 'CNOSPFClBrI'
    if 'pro_'+ep in BarCollection.keys():
        Bars = BarCollection['pro_'+ep]
        Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
        for Bar in Bars:
            if Bar[2] > 12.0: Bar[2] = 12.0
            if Bar[2] - Bar[1] < small: continue
            if Bar[0] == 0:
                Bar0Birth.append(Bar[1])
                Bar0Death.append(Bar[2])
            if Bar[0] == 1:
                Bar1Birth.append(Bar[1])
                Bar1Death.append(Bar[2])
            if Bar[0] == 2:
                Bar2Birth.append(Bar[1])
                Bar2Death.append(Bar[2])
        if len(Bar0Death) > 0:
            for i in range(len(Bar0Death)):
                did = int(np.floor(Bar0Death[i]/rs))
                Feature_h_p[did, 0] += 1.
                Feature_h_p[:did+1, 1] += 1.
        if len(Bar1Death) > 0:
            for i in range(len(Bar1Death)):
                did = int(np.floor(Bar1Death[i]/rs))
                bid = int(np.floor(Bar1Birth[i]/rs))
                Feature_h_p[bid, 2] += 1.
                Feature_h_p[did, 3] += 1.
                Feature_h_p[bid:did+1, 4] += 1.
        if len(Bar2Death) > 0:
            for i in range(len(Bar2Death)):
                did = int(np.floor(Bar2Death[i]/rs))
                bid = int(np.floor(Bar2Birth[i]/rs))
                Feature_h_p[bid, 5] += 1.
                Feature_h_p[did, 6] += 1.
                Feature_h_p[bid:did+1, 7] += 1.

    #h_pl------------------------------------------------------
    if 'pro_'+ep in BarCollection.keys() and 'com_'+ep+'_'+el in BarCollection.keys():
        Bars = BarCollection['com_'+ep+'_'+el]
        Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
        for Bar in Bars:
            if Bar[2] > 12.0: Bar[2] = 12.0
            if Bar[2] - Bar[1] < small: continue
            if Bar[0] == 0:
                Bar0Birth.append(Bar[1])
                Bar0Death.append(Bar[2])
            if Bar[0] == 1:
                Bar1Birth.append(Bar[1])
                Bar1Death.append(Bar[2])
            if Bar[0] == 2:
                Bar2Birth.append(Bar[1])
                Bar2Death.append(Bar[2])
        if len(Bar0Death) > 0:
            for i in range(len(Bar0Death)):
                did = int(np.floor(Bar0Death[i]/rs))
                Feature_h_pl[did, 0] += 1.
                Feature_h_pl[:did+1, 1] += 1.
        if len(Bar1Death) > 0:
            for i in range(len(Bar1Death)):
                did = int(np.floor(Bar1Death[i]/rs))
                bid = int(np.floor(Bar1Birth[i]/rs))
                Feature_h_pl[bid, 2] += 1.
                Feature_h_pl[did, 3] += 1.
                Feature_h_pl[bid:did+1, 4] += 1.
        if len(Bar2Death) > 0:
            for i in range(len(Bar2Death)):
                did = int(np.floor(Bar2Death[i]/rs))
                bid = int(np.floor(Bar2Birth[i]/rs))
                Feature_h_pl[bid, 5] += 1.
                Feature_h_pl[did, 6] += 1.
                Feature_h_pl[bid:did+1, 7] += 1.

    # c_p------------------------------------------------------
    ep = 'C'
    el = 'C'
    if 'pro_'+ep in BarCollection.keys():
        Bars = BarCollection['pro_'+ep]
        Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
        for Bar in Bars:
            if Bar[2] > 12.0: Bar[2] = 12.0
            if Bar[2] - Bar[1] < small: continue
            if Bar[0] == 0:
                Bar0Birth.append(Bar[1])
                Bar0Death.append(Bar[2])
            if Bar[0] == 1:
                Bar1Birth.append(Bar[1])
                Bar1Death.append(Bar[2])
            if Bar[0] == 2:
                Bar2Birth.append(Bar[1])
                Bar2Death.append(Bar[2])
        if len(Bar0Death) > 0:
            for i in range(len(Bar0Death)):
                did = int(np.floor(Bar0Death[i]/rs))
                Feature_c_p[did, 0] += 1.
                Feature_c_p[:did+1, 1] += 1.
        if len(Bar1Death) > 0:
            for i in range(len(Bar1Death)):
                did = int(np.floor(Bar1Death[i]/rs))
                bid = int(np.floor(Bar1Birth[i]/rs))
                Feature_c_p[bid, 2] += 1.
                Feature_c_p[did, 3] += 1.
                Feature_c_p[bid:did+1, 4] += 1.
        if len(Bar2Death) > 0:
            for i in range(len(Bar2Death)):
                did = int(np.floor(Bar2Death[i]/rs))
                bid = int(np.floor(Bar2Birth[i]/rs))
                Feature_c_p[bid, 5] += 1.
                Feature_c_p[did, 6] += 1.
                Feature_c_p[bid:did+1, 7] += 1.

    #h_pl------------------------------------------------------
    if 'pro_'+ep in BarCollection.keys() and 'com_'+ep+'_'+el in BarCollection.keys():
        Bars = BarCollection['com_'+ep+'_'+el]
        Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
        for Bar in Bars:
            if Bar[2] > 12.0: Bar[2] = 12.0
            if Bar[2] - Bar[1] < small: continue
            if Bar[0] == 0:
                Bar0Birth.append(Bar[1])
                Bar0Death.append(Bar[2])
            if Bar[0] == 1:
                Bar1Birth.append(Bar[1])
                Bar1Death.append(Bar[2])
            if Bar[0] == 2:
                Bar2Birth.append(Bar[1])
                Bar2Death.append(Bar[2])
        if len(Bar0Death) > 0:
            for i in range(len(Bar0Death)):
                did = int(np.floor(Bar0Death[i]/rs))
                Feature_h_pl[did, 0] += 1.
                Feature_h_pl[:did+1, 1] += 1.
        if len(Bar1Death) > 0:
            for i in range(len(Bar1Death)):
                did = int(np.floor(Bar1Death[i]/rs))
                bid = int(np.floor(Bar1Birth[i]/rs))
                Feature_c_pl[bid, 2] += 1.
                Feature_c_pl[did, 3] += 1.
                Feature_c_pl[bid:did+1, 4] += 1.
        if len(Bar2Death) > 0:
            for i in range(len(Bar2Death)):
                did = int(np.floor(Bar2Death[i]/rs))
                bid = int(np.floor(Bar2Birth[i]/rs))
                Feature_c_pl[bid, 5] += 1.
                Feature_c_pl[did, 6] += 1.
                Feature_c_pl[bid:did+1, 7] += 1.

    Feature = np.concatenate((Feature_c_pl, Feature_c_pl[:,:]-Feature_c_p[:,:], Feature_h_pl, Feature_h_pl[:,:]-Feature_h_p[:,:]), axis=1)
    outfile = open(working_dir+'/'+protein_name+'_feature_complex_alpha_1DCNN.npy', 'w')
    np.save(outfile, Feature)

def GenerateFeature_alpha_2DCNN(protein_name, working_dir):
    thr = 12.0
    rs = 0.1
    lth = int(thr/rs)
    small = 0.01

    OrderedName = np.load('PerformanceOrderAlphaHand.npy')
    X = np.zeros([16, lth, 128], float)

    InFile = open(working_dir+'/'+protein_name+'_alpha.pkl')
    BarCollection = pickle.load(InFile)
    for j in range(len(OrderedName)):
        plname = OrderedName[j]; pname, lname = plname.split('_');
        f_p_i_j = np.zeros([8,lth], float)
        f_pl_i_j = np.zeros([8,lth], float)
        if 'pro_'+pname in BarCollection.keys():
            Bars = BarCollection['pro_'+pname]
            for Bar in Bars:
                if Bar[1] >= thr: continue
                if Bar[2] >= thr: Bar[2] = thr-0.00001
                bid = min(int(np.floor(Bar[1]/rs)), lth-1)
                did = min(int(np.floor(Bar[2]/rs)), lth-1)
                if Bar[0] == 0:
                    f_p_i_j[0,did] += 1.0
                    f_p_i_j[1,bid:did+1] += 1.0
                elif Bar[0] == 1:
                    f_p_i_j[2,bid] += 1.0
                    f_p_i_j[3,did] += 1.0
                    f_p_i_j[4,bid:did+1] += 1.0
                elif Bar[0] == 2:
                    f_p_i_j[5,bid] += 1.0
                    f_p_i_j[6,did] += 1.0
                    f_p_i_j[7,bid:did+1] += 1.0
        if 'com_'+plname in BarCollection.keys():
            Bars = BarCollection['com_'+plname]
            for Bar in Bars:
                if Bar[1] >= thr: continue
                if Bar[2] >= thr: Bar[2] = thr-0.00001
                bid = min(int(np.floor(Bar[1]/rs)), lth-1)
                did = min(int(np.floor(Bar[2]/rs)), lth-1)
                if Bar[0] == 0:
                    f_pl_i_j[0,did] += 1.0
                    f_pl_i_j[1,bid:did+1] += 1.0
                elif Bar[0] == 1:
                    f_pl_i_j[2,bid] += 1.0
                    f_pl_i_j[3,did] += 1.0
                    f_pl_i_j[4,bid:did+1] += 1.0
                elif Bar[0] == 2:
                    f_pl_i_j[5,bid] += 1.0
                    f_pl_i_j[6,did] += 1.0
                    f_pl_i_j[7,bid:did+1] += 1.0
        f_df_i_j = f_pl_i_j[:,:] - f_p_i_j[:,:]
        X[0:8,:,j] = f_pl_i_j[:,:]; X[8:16,:,j] = f_df_i_j[:,:]

    OutFile = open(working_dir+'/'+protein_name+'_feature_complex_alpha_2DCNN.npy', 'w')
    np.save(OutFile, X)
    OutFile.close()
