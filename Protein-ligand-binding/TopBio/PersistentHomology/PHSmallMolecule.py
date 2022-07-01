import numpy as np
import os
from math import sqrt
import pickle
import sys


def Level1_Rips(s, ligand_name, working_dir):
    LigEleList = []
    LigEleList.extend([['C'],['N'],['O'],['S'],['H']])
    LigEleList.extend([['C','H'],['N','H'],['O','H'],['S','H']])
    LigEleList.extend([['C','N'],['C','O'],['C','S'],['N','O'],['N','S'],['O','S'],['C','Cl'],['C','Br'],['C','P'],['C','F']])
    LigEleList.extend([['C','N','H'],['C','O','H'],['C','S','H'],['N','O','H'],['N','S','H'],['O','S','H'],['C','Cl','H'],['C','Br','H'],['C','P','H'],['C','F','H']])
    LigEleList.extend([['C','N','O'],['C','N','S'],['C','O','S'],['N','O','S'],['C','N','O','S']])
    LigEleList.extend([['C','N','O','H'],['C','N','S','H'],['C','O','S','H'],['N','O','S','H'],['C','N','O','S','H']])
    LigEleList.append(['C','N','O','S','P','F','Cl','Br','I'])
    LigEleList.append(['C','N','O','S','P','F','Cl','Br','I','H']) # 41 combinations in total

    for el in LigEleList:
        elname = ''
        for eel in el:
            elname = elname+eel
        index = []
        for j in range(len(s.pos)):
            if s.atmtyp[j] in el:
                index.append(j)
        if len(index) > 3:
            PtsOutFile = open(working_dir+'/'+ligand_name+'_'+elname+'.pts', 'w')
            AlphaPtsOutFile = open(working_dir+'/'+ligand_name+'_'+elname+'.csv', 'w')
            BondOutFile = open(working_dir+'/'+ligand_name+'_'+elname+'.bds', 'w')
            AlphaPtsOutFile.write('x1,x2,x3\n')
            for j in index:
                PtsOutFile.write(str(s.pos[j][0])+' '+str(s.pos[j][1])+' '+str(s.pos[j][2])+'\n')
                AlphaPtsOutFile.write(str(s.pos[j][0])+','+str(s.pos[j][1])+','+str(s.pos[j][2])+'\n')
            indexarr = np.asarray(index, int)
            for j1 in index:
                for j2 in index:
                    if [j1+1,j2+1] in s.bond or [j2+1,j1+1] in s.bond:
                        id1 = np.where(indexarr==j1)[0][0] + 1
                        id2 = np.where(indexarr==j2)[0][0] + 1
                        BondOutFile.write(str(id1)+' '+str(id2)+'\n')
            PtsOutFile.close()
            AlphaPtsOutFile.close()
            BondOutFile.close()
    line = 'matlab -nodisplay -nodesktop -nosplash -r "DataDir=' + "'"+working_dir+"';pdb="+"'"+ligand_name+"'"+';PH_small_molecule_level1"'
    os.system(line)

def Alpha(s, ligand_name, working_dir):
    LigEleList = []
    LigEleList.extend([['C'],['N'],['O'],['S'],['H']])
    LigEleList.extend([['C','H'],['N','H'],['O','H'],['S','H']])
    LigEleList.extend([['C','N'],['C','O'],['C','S'],['N','O'],['N','S'],['O','S'],['C','Cl'],['C','Br'],['C','P'],['C','F']])
    LigEleList.extend([['C','N','H'],['C','O','H'],['C','S','H'],['N','O','H'],['N','S','H'],['O','S','H'],['C','Cl','H'],['C','Br','H'],['C','P','H'],['C','F','H']])
    LigEleList.extend([['C','N','O'],['C','N','S'],['C','O','S'],['N','O','S'],['C','N','O','S']])
    LigEleList.extend([['C','N','O','H'],['C','N','S','H'],['C','O','S','H'],['N','O','S','H'],['C','N','O','S','H']])
    LigEleList.append(['C','N','O','S','P','F','Cl','Br','I'])
    LigEleList.append(['C','N','O','S','P','F','Cl','Br','I','H'])

    for el in LigEleList:
        elname = ''
        for eel in el:
            elname = elname+eel
        index = []
        for j in range(len(s.pos)):
            if s.atmtyp[j] in el:
                index.append(j)
        if len(index) > 3:
            AlphaPtsOutFile = open(working_dir+'/'+ligand_name+'_'+elname+'.csv', 'w')
            AlphaPtsOutFile.write('x1,x2,x3\n')
            for j in index:
                AlphaPtsOutFile.write(str(s.pos[j][0])+','+str(s.pos[j][1])+','+str(s.pos[j][2])+'\n')
            AlphaPtsOutFile.close()
    LIGELE = ['C','N','O','S','CN','CO','CS','NO','NS','OS','CNO','CNS','COS','NOS','CNOS','CNOSPFClBrI','H','CH','NH','OH','SH','CNH','COH','CSH','NOH','NSH','OSH','CNOH','CNSH','COSH','NOSH','CNOSH','CNOSPFClBrIH','CCl','CClH','CBr','CBrH','CP','CF','CPH','CFH']

    small = 0.01

    BarCollection = {}

    for e in LIGELE:
        if not os.path.exists(working_dir+'/'+ligand_name+'_'+e+'.csv'): continue
        os.system('Rscript PH_Alpha_ligand.R '+working_dir+'/'+ligand_name+'_'+e+'.csv '+working_dir+'/tmp.out')
        tmpinfile = open(working_dir+'/tmp.out')
        lines = tmpinfile.read().splitlines()
        Bars = []
        for line in lines[1:]:
            a,b,c,d = line.split()
            if d!='Inf':
                if float(d)-float(c) >= small:
                    Bars.append([int(b), float(c), float(d)])
        tmpinfile.close()
        BarCollection['lig_'+e] = Bars

    OutFile = open(working_dir+'/'+ligand_name+'_alpha.pkl', 'w')
    pickle.dump(BarCollection, OutFile, protocol=pickle.HIGHEST_PROTOCOL)
    OutFile.close()
