import numpy as np
import pickle
import os

def GenerateFeature_alpha(ligand_name, working_dir):
    Cut = 12.0

    LIGELE = ['C','N','O','S','CN','CO','CS','NO','NS','OS','CCl','CBr','CP','CF','CNO','CNS','COS','NOS','CNOS','CNOSPFClBrI','H','CH','NH','OH','SH','CNH','COH','CSH','NOH','NSH','OSH','CNOH','CNSH','COSH','NOSH','CNOSH','CNOSPFClBrIH','CClH','CBrH','CPH','CFH']

    Feature_i = []
    pdb = ligand_name
    InFile = open(working_dir+'/'+ligand_name+'_alpha.pkl')
    BarCollection = pickle.load(InFile)
    for el in LIGELE:
        if 'lig_'+el in BarCollection.keys():
            Bars = BarCollection['lig_'+el]
            Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
            for Bar in Bars:
                if Bar[2] < Bar[1]:
                    continue
                if Bar[2] > 12.0 and Bar[0] == 0: continue
                if Bar[2] > 12.0 and Bar[0] > 0: Bar[2] = 12.0
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
        else:
            Feature_i.extend([0.]*42)
    Feature_i = np.asarray(Feature_i, float)

    outfile = open(working_dir+'/'+ligand_name+'_feature_alpha_handcrafted.npy', 'w')
    np.save(outfile, Feature_i)
    outfile.close()

def GenerateFeature_level1(ligand_name, working_dir):
    small = 0.01
    Feature_i = []
    Cut = 12.0

    LIGELE = ['C','N','O','S','CN','CO','CS','NO','NS','OS','CCl','CBr','CP','CF','CNO','CNS','COS','NOS','CNOS','CNOSPFClBrI','H','CH','NH','OH','SH','CNH','COH','CSH','NOH','NSH','OSH','CNOH','CNSH','COSH','NOSH','CNOSH','CNOSPFClBrIH','CClH','CBrH','CPH','CFH']

    pdb = ligand_name

    for el in LIGELE:
        if os.path.exists(working_dir+'/'+ligand_name+'_'+el+'_level1.PH'):
            InFile = open(working_dir+'/'+ligand_name+'_'+el+'_level1.PH')
            lines = InFile.read().splitlines()
            Bars = []
            for line in lines:
                a,b,c = line.split()
                Bars.append([int(a), float(b), float(c)])
            InFile.close()
            Bar0Birth = []; Bar0Death = []; Bar1Birth = []; Bar1Death = []; Bar2Birth = []; Bar2Death = [];
            for Bar in Bars:
                if Bar[2] < Bar[1]:
                    continue
                if Bar[2] > 12.0 and Bar[0] == 0: continue
                if Bar[2] > 12.0 and Bar[0] > 0: Bar[2] = 12.0
                if Bar[0] == 0 and Bar[2]-Bar[1] >= small:
                    Bar0Birth.append(Bar[1])
                    Bar0Death.append(Bar[2])
                if Bar[0] == 1 and Bar[2]-Bar[1] >= small:
                    Bar1Birth.append(Bar[1])
                    Bar1Death.append(Bar[2])
                if Bar[0] == 2 and Bar[2]-Bar[1] >= small:
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
        else:
            Feature_i.extend([0.]*42)
    Feature_i = np.asarray(Feature_i, float)

    outfile = open(working_dir+'/'+ligand_name+'_feature_ligand_level1_handcrafted.npy', 'w')
    np.save(outfile, Feature_i)
    outfile.close()
