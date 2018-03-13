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
        if self._mod is not None:
            return Residue.MOD_MASSES[self._mod]
        else:
            return 0
        
    @property
    def aa_mass(self):
        return Residue.AA_MASSES[self.aa]

    @property
    def code(self):
        code = self.aa
        if self.mod is not None:
            code = self.mod + code
        return code

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

    @mod.deleter
    def mod(self):
        self._mod = None
        
    @classmethod
    def AAS(cls):
        return cls.AA_MASSES.keys()

    @classmethod
    def MODS(cls):
        return cls.MOD_MASSES.keys()

    def __eq__(self, y):
        return type(y) is Residue and y.code == self.code
        
    
    def __repr__(self):
        return '%s(\'%s\')' % (self.__class__.__name__,
                               self.code)


    
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
        if self:
            return sum(r.mass for r in self.residues) + Protein.WATER_MASS
        else:
            return 0

    @property
    def sequence(self):
        return ''.join(r.code for r in self.residues)

    def charged_mass(self, num_protons=1):
        if self:
            return self.mass + (Protein.H_MASS * num_protons)
        else:
            return 0

    def fragment_mass(self, position):
        y_fragment = Protein(self.residues[position:])
        y = y_fragment.charged_mass()
        if y_fragment and len(y_fragment) < len(self):
            total_charge = 2
        else:
            total_charge = 1
        b = self.charged_mass(total_charge) - y
        return (b, y)
    
    def phosphorylate(self, position):
        self.residues[position].mod = 'p'

    def clear_mods(self, position=None):
        if position is not None:
            del self.residues[position].mod
        else:
            for r in self.residues:
                del r.mod

    def all_ions(self):
        b, y = [], []
        n_frag_sites = len(protein) + 1
        for ifrag in range(n_frag_sites)
            new_b, new_y = protein.fragment_mass(ifrag)
            b.append(new_b)
            y.append(new_y)
        y = list(reversed(y_ion_masses))
        return [(b[i], y[i]) for i in range(n_frag_sites)] 
                
    def __bool__(self):
        return bool(self.residues)

    def __eq__(self, y):
        return type(y) is Protein and y.sequence == self.sequence
                
    def __contains__(self, x):
        return x in self.residues or x in self.sequence
                
    def __len__(self):
        return len(self.residues)

    def __repr__(self):
        return '%s(\'%s\')' % (self.__class__.__name__,
                               self.sequence)

# protein sequence
protein = Protein('ELFDDPSYVNVQNLDK')
p_res = ['R','S','Y']
p_sites = [i for i, r in enumerate(protein.sequence) if r in p_res]

def list_ion_table(ion_masses):
    header = '%10s | %10s | %10s' % ('Num Res.', 'b ion', 'y ion')
    print(header)
    body_fmt = '%10d | %10.2f | %10.2f'
    for i, (b, y) in enumerate(ion_masses):
        print(body_fmt % (i, b, y))

for site in p_sites:
    print('\n Phosphorlating %s-%d' % (protein.sequence[site], site + 1))
    protein.phosphorylate(site)
    list_ion_table(protein.all_ions())
    protein.clear_mods()


    
     
