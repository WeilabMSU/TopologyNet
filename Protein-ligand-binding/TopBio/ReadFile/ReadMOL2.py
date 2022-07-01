import numpy as np
class SmallMolecule:

    def __init__(self, ligand_name, working_dir):
        # Read mol2 file
        mol2file = open(working_dir+'/'+ligand_name+'.mol2')
        lines = mol2file.read().splitlines()
        for i in range(len(lines)):
            if lines[i].replace(" ","") == '@<TRIPOS>MOLECULE':
                self.natom, self.nbond, _, _, _ = lines[i+2].split()
                self.natom = int(self.natom); self.nbond = int(self.nbond);
                break
        self.pos = np.empty([self.natom, 3], float);
        self.atmtyp = [];
        self.chg = np.empty([self.natom], float);
        self.bond = [];
        for i in range(len(lines)):
            if lines[i].replace(" ","") == '@<TRIPOS>ATOM':
                for j in range(self.natom):
                    line = lines[i+j+1]
                    _,_,x,y,z,t,_,_,c = line.split()
                    self.pos[j,:] = np.array([float(x),float(y),float(z)])[:]
                    self.chg[j] = float(c)
                    if '.' in t:
                        tt,_ = t.split('.')
                    else:
                        tt = t
                    self.atmtyp.append(tt)
                break
        for i in range(len(lines)):
            if lines[i].replace(" ","") == '@<TRIPOS>BOND':
                for j in range(self.nbond):
                    line = lines[i+j+1]
                    _,a,b,_ = line.split()
                    self.bond.append([int(a), int(b)]);
                break
