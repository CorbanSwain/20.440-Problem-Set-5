#!\bin\python3

import numpy as np

## Problem 1, Part 1 ##

class Residue:
    AA_MASSES = {'A' : 71.03711,  'R' : 156.10111, 'N' : 114.04293, 'D' : 115.02694, 'C' : 103.00919,
                 'E' : 129.04259, 'Q' : 128.05858, 'G' : 57.02146,  'H' : 137.05891, 'I' : 113.08406,
                 'L' : 113.08406, 'K' : 128.09496, 'M' : 131.04049, 'F' : 147.06841, 'P' : 97.05276,
                 'S' : 87.03203,  'T' : 101.04768, 'W' : 186.07931, 'Y' : 163.06333, 'V' : 99.06841}
    MOD_MASSES = {'p' : 79.97}

    def __init__(self, code):        
        self.code = code

    @property
    def mass(self):
        return self.aa_mass + self.mod_mass

    @property
    def mod_mass(self):
        if self._mod is None:
            return 0
        else:
            return Residue.MOD_MASSES[self._mod]
        
    @property
    def aa_mass(self):
        return Residue.AA_MASSES[self.aa]

    @property
    def code(self):
        return self.aa + self.mod

    @code.setter
    def code(self, value):
        self.aa = value[-1]
        if len(value) == 2:
            self._mod = value[0]
        else:
            self._mod = None

    @property
    def aa(self):
        return self._aa
            
    @aa.setter
    def aa(self, value):
        if value not in Residue.AAS():
            raise ValueError('%s is not a valid amino acid' % value)
        self._aa = value

    @property
    def mod(self):
        return self._mod
            
    @mod.setter
    def mod(self, value):
        if value not in Residue.MODS():
            raise ValueError('%s is not a valid aa modification' % value)
        self._mod = value
        
    @classmethod
    def AAS(cls):
        return cls.AA_MASSES.keys()

    @classmethod
    def MODS(cls):
        return cls.MOD_MASSES.keys()


    
class Protein:
    WATER_MASS = 18.010565
    H_MASS = 1.007276
    residues = []
    
    def __init__(self, residues):
        if type(residues) is list:
            self.residues = residues
        else:
            self.residues = Protein.str2res(residues)

    @staticmethod
    def str2res(seqstr):
        res = []
        ichar = 0
        while ichar < len(seqstr):
            c = seqstr[ichar]
            if c in Residue.MODS():
                next_code = seqstr[ichar:(ichar + 2)] 
                ichar += 2
            else:
                next_code = seqstr[ichar]
                ichar += 1
            res.append(Residue(next_code))
        return res
            
    @property
    def mass(self):
        return sum(r.mass for r in self.residues) + Protein.WATER_MASS

    @property
    def seqstr(self):
        return [r.code for r in self.residues]

    def charged_mass(self, num_protons=1):
        return self.mass + (Protein.H_MASS * num_protons)

    def fragment_mass(self, index):
        y = Protein(self.residues[index:]).charged_mass()
        b = self.charged_mass(2) - y
        return (b, y)
        
    

# protein sequence
seqstr = 'ELFDDPSYVNVQNLDK'
prot = Protein(seqstr)
n_phospho = 1


    
     
