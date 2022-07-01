import numpy as np
import os
from math import sqrt
import pickle
import sys

def Interaction_Rips(cut, protein_name, working_dir):

    print "Computing 0-dimensional persistent homology with interaction distance."

    com = np.load(working_dir+'/'+protein_name+'_'+str(cut)+'.struct')
    PRO = com['PRO']; LIG = com['LIG'];
    for e1 in [['C'],['N'],['O'],['S']]:
        for e2 in [['C'],['N'],['O'],['S'],['P'],['F'],['Cl'],['Br'],['I']]:
            PROpts = []; LIGpts = [];
            for i in range(len(PRO)):
                if PRO[i][0].replace(" ","") in e1:
                    PROpts.append(PRO[i][1])
            for i in range(len(LIG)):
                if LIG[i][0].replace(" ","") in e2:
                    LIGpts.append(LIG[i][1])
            if len(PROpts) > 0 and len(LIGpts) > 0:
                namee1 = ''
                for ie in range(len(e1)):
                    namee1 = namee1+e1[ie]
                namee2 = ''
                for ie in range(len(e2)):
                    namee2 = namee2+e2[ie]
                Name = protein_name+'_'+namee1+'_'+namee2+'_'+str(cut)+'.pts'
                OutFile = open(working_dir+'/'+Name, 'w')
                for i in range(len(PROpts)):
                    OutFile.write('1 '+str(PROpts[i][0])+' '+str(PROpts[i][1])+' '+str(PROpts[i][2])+'\n')
                for i in range(len(LIGpts)):
                    OutFile.write('0 '+str(LIGpts[i][0])+' '+str(LIGpts[i][1])+' '+str(LIGpts[i][2])+'\n')
                OutFile.close()

    line = 'matlab -nodisplay -nodesktop -nosplash -r "DataDir=' + "'"+working_dir+"';pdb="+"'"+protein_name+"'"+';PH_complex_interaction"'
    os.system(line)

def Electrostatics_Rips(cut, protein_name, working_dir):

    print "Computing 0-dimensional persistent homology with electrostatics interaction induced distance."

    LIGELE = ['C','N','O','S','P','F','Cl','Br','I','H']
    PROELE = ['C','N','O','S','H']
    InFile = np.load(working_dir+'/'+protein_name+'_'+str(cut)+'_chg.npz')
    LIG = InFile['LIG']; PRO = InFile['PRO'];
    for e1 in PROELE:
        for e2 in LIGELE:
            OutFile = open(working_dir+'/'+protein_name+'_'+e1+'_'+e2+'_'+str(cut)+'_chg.pts', 'w')
            for j in range(len(LIG)):
                if LIG[j][0] == e2:
                    OutFile.write('0 '+str(LIG[j][1])+' '+str(LIG[j][2])+' '+str(LIG[j][3])+' '+str(LIG[j][4])+'\n')
            for j in range(len(PRO)):
                if PRO[j][0] == e1:
                    OutFile.write('1 '+str(PRO[j][1])+' '+str(PRO[j][2])+' '+str(PRO[j][3])+' '+str(PRO[j][4])+'\n')
            OutFile.close()

    line = 'matlab -nodisplay -nodesktop -nosplash -r "DataDir=' + "'"+working_dir+"';pdb="+"'"+protein_name+"'"+';PH_complex_charge"'
    os.system(line)

def Alpha(cut, protein_name, working_dir):

    print "Computing alpha complex based persistent homology for protein-ligand complex."

    dt = np.dtype([('dim', int), ('birth', float), ('death', float)])
    ProEleCollection = [['C'],['N'],['O'],['S'],['C','N'],['C','O'],['N','O'],['C','N','O'],['C','N','O','S']]
    LigEleCollection = [['C'],['N'],['O'],['S'],['C','N'],['C','O'],['C','S'],['N','O'],['N','S'],['O','S'],['C','N','O'],['C','N','S'],['C','O','S'],['N','O','S'],['C','N','O','S'],['C','N','O','S','P','F','Cl','Br','I']]

    small = 0.01

    com = np.load(working_dir+'/'+protein_name+'_'+str(cut)+'.struct')
    PRO = com['PRO']; LIG = com['LIG'];
    BarCollection = {}
    # alpha computation for complex
    cnt = 0
    for ep in ProEleCollection:
        for el in LigEleCollection:
            propts = []
            for a in range(len(PRO)):
                if PRO[a][0].replace(" ","") in ep:
                    propts.append([ PRO[a][1][0], PRO[a][1][1], PRO[a][1][2] ])
            ligpts = []
            for a in range(len(LIG)):
                if LIG[a][0].replace(" ","") in el:
                    ligpts.append([ LIG[a][1][0], LIG[a][1][1], LIG[a][1][2] ])
            if len(propts) + len(ligpts) > 3:
                pname = ''
                for eep in ep:
                    pname = pname + eep
                lname = ''
                for eel in el:
                    lname = lname + eel
                name = 'com_'+pname+'_'+lname
                pt = propts + ligpts
                Bars = []
                tmpoutfile = open(working_dir+'/pt.csv', 'w')
                tmpoutfile.write('x1,x2,x3\n')
                for pp in pt:
                    tmpoutfile.write(str(pp[0])+','+str(pp[1])+','+str(pp[2])+'\n')
                tmpoutfile.close()
                os.system('Rscript PH_Alpha.R '+working_dir+'/pt.csv '+working_dir+'/tmp.out')
                tmpinfile = open(working_dir+'/tmp.out')
                lines = tmpinfile.read().splitlines()
                for line in lines[1:]:
                    a,b,c,d = line.split()
                    if d!='Inf':
                        if float(d)-float(c) >= small:
                            Bars.append([int(b), float(c), float(d)])
                BarCollection[name] = Bars
                os.system('rm '+working_dir+'/pt.csv '+working_dir+'/tmp.out')

    # alpha computation for protein
    for ep in ProEleCollection:
        propts = []
        for a in range(len(PRO)):
            if PRO[a][0].replace(" ","") in ep:
                propts.append([ PRO[a][1][0], PRO[a][1][1], PRO[a][1][2] ])
        if len(propts) > 3:
            pname = ''
            for eep in ep:
                pname = pname + eep
            name = 'pro_'+pname
            pt = propts
            Bars = []
            tmpoutfile = open(working_dir+'/pt.csv', 'w')
            tmpoutfile.write('x1,x2,x3\n')
            for pp in pt:
                tmpoutfile.write(str(pp[0])+','+str(pp[1])+','+str(pp[2])+'\n')
            tmpoutfile.close()
            os.system('Rscript PH_Alpha.R '+working_dir+'/pt.csv '+working_dir+'/tmp.out')
            tmpinfile = open(working_dir+'/tmp.out')
            lines = tmpinfile.read().splitlines()
            for line in lines[1:]:
                a,b,c,d = line.split()
                if d!='Inf':
                    if float(d)-float(c) >= small:
                        Bars.append([int(b), float(c), float(d)])
            BarCollection[name] = Bars

    OutFile = open(working_dir+'/'+protein_name+'_alpha.pkl', 'w')
    pickle.dump(BarCollection, OutFile, protocol=pickle.HIGHEST_PROTOCOL)
    OutFile.close()
