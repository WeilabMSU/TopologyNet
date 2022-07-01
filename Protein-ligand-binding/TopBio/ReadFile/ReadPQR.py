import numpy as np

def get_pqr_structure(ligand_name, protein_name, working_dir, cut=16.0):

    infile = open(working_dir+'/'+ligand_name+'.mol2')
    lines = infile.read().splitlines()
    LIG = []
    for j in range(len(lines)):
        if len(lines[j]) >= 13:
            if lines[j][0:13] == '@<TRIPOS>ATOM':
                istart = j+1
            if lines[j][0:13] == '@<TRIPOS>BOND':
                iend = j
    for j in range(istart,iend):
        a,b,c,d,e,f,g,h,chg = lines[j].split()
        if '.' in f:
            tmptyp, dummy = f.split('.')
        else:
            tmptyp = f
        if tmptyp == 'CL' or tmptyp == 'Cl':
            tmptyp = 'Cl'
        elif tmptyp == 'BR' or tmptyp == 'Br':
            tmptyp = 'Br'
        LIG.append([tmptyp, float(c), float(d), float(e), float(chg)])
    infile.close()
    PRO = []
    infile = open(working_dir+'/'+protein_name+'.pqr')
    lines = infile.read().splitlines()
    for line in lines:
        if line[0:4] == 'ATOM':
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54]);
            dis = 1000.0
            for j in range(len(LIG)):
                tmpdis = np.linalg.norm(np.array([x,y,z])[:] - np.asarray(LIG[j][1:4])[:])
                if tmpdis < dis: dis = tmpdis
            if dis < cut:
                tmptyp = line[12:17]
                typ = tmptyp.replace(" ","")
                typ = typ[0]
                chg = float(line[54:62])
                PRO.append([typ,x,y,z,chg])
    PRO = np.asarray(PRO); LIG = np.asarray(LIG);
    OutFile = open(working_dir+'/'+protein_name+'_'+str(cut)+'_chg.npz', 'w')
    np.savez(OutFile, PRO=PRO, LIG=LIG)
    OutFile.close()
