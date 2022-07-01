import numpy as np

def get_pdb_structure(s, cut, protein_name, working_dir):

    typatm = np.dtype([('typ', 'S2'), ('pos', float, (3,)), ('rad', float), ('id', int)])
    lig_ele_list = ['C','N','O','S','P','F','Cl','Br','I','H']
    pro_ele_list = ['C','N','O','S','H']
    aa_list = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HSE','HSD','SEC',
               'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','PYL']

    def gettyp( rawtyp ):
        if not rawtyp[1].isalpha():
            typ = ' '+rawtyp[0]
        else:
            typ = rawtyp
        return typ

    LIG = np.zeros([s.natom], dtype = typatm)
    for i in range(s.natom):
        LIG[i]['pos'][:] = s.pos[i,:]
        LIG[i]['typ'] = s.atmtyp[i]
        LIG[i]['id'] = 1

    pronum = 0
    profile = open(working_dir+'/'+protein_name+'.pdb')
    lines = profile.read().splitlines()
    for line in lines:
        if line[0:4] == 'ATOM' and line[17:20] in aa_list:
            typ = line[12:14].replace(" ","")
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            for j in range(0,len(LIG)):
                dis = np.linalg.norm(np.array([x,y,z])-LIG[j]['pos'])
                if dis <= cut and typ in pro_ele_list:
                    pronum += 1
                    break
    PRO = np.zeros([pronum], dtype = typatm)
    j = 0
    for line in lines:
        if line[0:4] == 'ATOM' and line[17:20] in aa_list:
            typ = line[12:14]
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            for k in range(0,len(LIG)):
                dis = np.linalg.norm(np.array([x,y,z])-LIG[k]['pos'])
                if dis <= cut and typ.replace(" ", "") in pro_ele_list:
                    PRO[j]['typ'] = typ; PRO[j]['pos'][:] = np.array([x,y,z]);
                    PRO[j]['id'] = -1; j+=1
                    break
    print 'Number of atoms in ligand/protein with cutoff ', str(cut)
    print len(LIG), '/', len(PRO)
    outname = protein_name+'_'+str(cut)+'.struct';
    outfile = open(working_dir+'/'+outname, 'w')
    np.savez(outfile,PRO=PRO,LIG=LIG)
    outfile.close()
