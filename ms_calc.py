#!\bin\python3

# import numpy as np
import itertools as it

class Residue:
    AA_MASSES = {'A' : 71.03711,  'R' : 156.10111, 'N' : 114.04293, 'D' : 115.02694, 'C' : 103.00919,
                 'E' : 129.04259, 'Q' : 128.05858, 'G' : 57.02146,  'H' : 137.05891, 'I' : 113.08406,
                 'L' : 113.08406, 'K' : 128.09496, 'M' : 131.04049, 'F' : 147.06841, 'P' : 97.05276,
                 'S' : 87.03203,  'T' : 101.04768, 'W' : 186.07931, 'Y' : 163.06333, 'V' : 99.06841}
    MOD_MASSES = {'p' : 79.97}

    def __init__(self, res):
        if res.__class__.__name__ is self.__class__.__name__:
            self.code = res.code
        else:
            self.code = res.strip()

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
    
    def __init__(self, residues=[]):
        if residues.__class__.__name__ is self.__class__.__name__:
            self.residues = residues.residues
        elif residues.__class__.__name__ is 'Residue':
            self.residues = [residues]
        elif type(residues) is list:
            self.residues = [Residue(r) for r in residues]
        else:
            self.residues = Protein.str2res(residues.strip())

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
        for ifrag in range(0, len(self) + 1):
            new_b, new_y = self.fragment_mass(ifrag)
            b.append(new_b)
            y.append(new_y)
        y = list(reversed(y))
        return [(b[i], y[i]) for i in range(0, len(self) + 1)] 
                
    def __bool__(self):
        return bool(self.residues)

    def __eq__(self, y):
        return type(y) is Protein and y.sequence == self.sequence
                
    def __contains__(self, x):
        return x in self.residues or x in self.sequence
                
    def __len__(self):
        return len(self.residues)

    def __add__(self, y):
        x = Protein(self)
        y = Protein(y)
        x.residues = x.residues + y.residues
        return x

    def __radd__(self, x):
        x = Protein(x)
        y = Protein(self)
        x.residues = x.residues + y.residues
        return x
            
    def __repr__(self):
        return '%s(\'%s\')' % (self.__class__.__name__,
                               self.sequence)



# protein sequence
protein = Protein('ELFDDPSYVNVQNLDK')
p_res = ['R','S','Y']
all_res = list(Residue.AAS()) + ['p' + r for r in p_res] + [' ']
p_sites = [i for i, r in enumerate(protein.sequence) if r in p_res]

masses = [120.8, 129.10, 147.11, 185.09, 197.13, 213.09, 215.14, 243.13, 245.11,
          263.10, 280.19, 298.20, 331.07, 344.20, 375.22, 378.13, 390.20, 410.11,
          428.12, 441.21, 459.22, 487.22, 489.27, 499.19, 525.14, 527.19, 543.15,
          569.31, 584.24, 600.30, 617.33, 620.26, 641.23, 651.35, 670.30, 684.33,
          712.30, 716.39, 736.32, 740.30, 742.30, 768.33, 787.23, 813.41, 830.44,
          832.44, 855.34, 868.36, 904.31, 929.51, 931.51, 965.38, 982.41, 1011.35,
          1047.37, 1080.41, 1095.49, 1117.42, 1146.44, 1172.54, 1193.48, 1224.52,
          1242.55, 1259.56, 1276.66, 1321.58, 1339.60, 1356.62, 1392.67, 1454.62,
          1471.65, 1569.64, 1586.68, 1596.65, 1715.73, 1733.73]
check_margin = 0.04
mass_ranges_1 = [(m - check_margin, m + check_margin) for m in masses]
def in_ranges(range_list, query_list):
    out = [False for i in range(len(query_list))]
    for low, hi in range_list:
        for i, q in enumerate(query_list):
            #print('%f --%f-- %f' % (low, q, hi))
            if low <= q <= hi:
                out[i] = True
    return out

def print_ion_table(ion_masses, mass_ranges, silent=False):
    header = '%10s | %12s | %12s' % ('Num Res.', 'b+ ion', 'y+ ion')
    if not silent: print(header)
    body_fmt = '%10d | %2s%10.2f | %2s%10.2f'
    b_count, y_count = 0, 0
    for i, (b, y) in enumerate(ion_masses):
        found_list = in_ranges(mass_ranges, [b, y])
        b_star, y_star = '', ''
        if found_list[0]:
            b_star = '*'
            b_count += 1  
        if found_list[1]:
            y_star = '*'
            y_count += 1    
        # double ions
        b_2plus = (b + Protein.H_MASS) / 2
        y_2plus = (y + Protein.H_MASS) / 2
        found_list = in_ranges(mass_ranges, [b_2plus, y_2plus])
        if found_list[0]:
            b_star += '*'
            b_count += 1  
        if found_list[1]:
            y_star += '*'
            y_count += 1

        if not silent: print(body_fmt % (i, b_star, b, y_star, y))

    if not silent: print('%10s   %12d | %12d   = %d\n\n' % ('N Match. :',
                                                            b_count, y_count,
                                                            b_count + y_count))
    return (b_count, y_count)

for site in p_sites:
    print('\n Phosphorylating %s-%d' % (protein.sequence[site], site + 1))
    protein.phosphorylate(site)
    print_ion_table(protein.all_ions(), mass_ranges_1)
    protein.clear_mods()


## Part 2 ##
    
masses_2 = [147.11, 157.10, 167.08, 175.12, 183.15, 185.09, 197.13, 207.11,
            211.14, 216.04, 236.10, 244.13, 254.11, 266.15, 272.12, 280.17,
            294.18, 298.18, 316.21, 334.22, 357.21, 367.20, 381.21, 385.21,
            403.13, 416.87, 440.16, 452.25, 468.25, 491.68, 505.28, 523.29,
            527.19, 559.15, 569.20, 581.33,
            603.76, 630.18, 638.32, 647.27, # Plot 2
            656.24, 672.78, 684.23, 711.28, 721.30, 730.30, 739.32, 755.27,
            767.36, 786.85, 795.86, 798.30, 826.31, 841.30, 858.29, 893.31,
            895.33, 913.34, 941.35, 955.34, 982.35, 992.38, 1005.98, 1010.39,
            1026.38, 1042.39,
            1068.43, 1097.42, 1109.46, 1127.55, 1139.46, # plot 3
            1166.44, 1188.50, 1206.51, 1275.52, 1293.54, 1300.60, 1344.56,
            1362.56, 1380.57, 1390.55, 1459.61, 1478.62, 1590.713]
check_margin = 0.015
mass_ranges_2 = [(m - check_margin, m + check_margin) for m in masses_2]
target_mz = 795.86
target_mass = target_mz * 2 - (Protein.H_MASS * 2)


def spectra_score(temp_p, mass_ranges, target_mass):
    did_find_candidate = False
    (b, y) = print_ion_table(temp_p.all_ions(), mass_ranges, silent=True)
    if target_mass - temp_p.mass < -0.5:
        score = -(5 * len(temp_p))
    elif len(temp_p) > 12:
        score = (y * 10) + (10 * b) - (2 * 5 * len(temp_p))
        if abs(target_mass - temp_p.mass) < 0.015:
            score += 100
            did_find_candidate = True
    else:
        score = y * 10 - (5 * len(temp_p))
    score = round(score)
    return (score, did_find_candidate)

candidate_seq = ['RDLPVpYEDAASFK', 'LPSSHRDSpYSSFK', 'VGDLPVPSDpYSSFK', 'RDLPVPSDpYSSFK',  'PISSPVPDSpYSSFK',
                 'RPMPVpYEDAASFK', 'RPMPVpYEDAASFK', 'RHSTVPDSpYSSFK']
for c in candidate_seq:
    c_prot = Protein(c)
    c_score, _ = spectra_score(c_prot, mass_ranges_2, target_mass)
    b, y = print_ion_table(c_prot.all_ions(), mass_ranges_2, silent=True)
    print('Sequence: %20s, Score = %4d (%2d, %2d)' % (c, c_score, b, y), end=' ')
    print('Mass Diff: %5.2f' % (target_mass - c_prot.mass))
    # print_ion_table(c_prot.all_ions(), mass_ranges_2)
# PISSPVPDSpYSSFK
#  RDLPVPDSpYSSFK = 290 (6, 13) Score, -0.01 Da mass diff
#  RDLPVPSDpYSSFK = ???         Score, -0.01 Da mass diff
# GVDLPVPSDpYSSFK = 300 (6, 14) Score,  0.00 Da mass diff
# VGDLPVPSDpYSSFK = 300 (6, 14) Score,  0.00 Da mass diff
#  RDLPVpYEDAASFK = 310 (8, 13) Score, -0.01 Da mass diff <--
# LPSSHRDSpYSSFK

input('Press ENTER to continue ...\n')

max_score = None
best_list = []
for i, n_mer in enumerate(it.product(all_res, repeat=4)):
    p = Protein(''.join(list(n_mer)).replace(' ', ''))
    (b, y) = print_ion_table(p.all_ions(), mass_ranges_2, silent=True)
    score, _ = spectra_score(p, mass_ranges_2, target_mass)
    if max_score is None:
        max_score = score
        best_list = [p.sequence]
    if abs(score - max_score) <= 1:
        print('%10d | %20s | %8d' % (i, p.sequence, score))
        best_list.append(p.sequence)
    if score > max_score:
        best_list = [p.sequence]
        max_score = score
        print('%10d | %20s | * %6d' % (i, p.sequence, score))
best_list = list(set([s[1:] for s in best_list]))

candidates = []
for best in best_list:    
    p = Protein(best)
    best_list_old = []
    final_max_score = max_score
    found_something = False
    while not found_something:
        max_score = None
        best_list_2 = []
        i_best = 0
        for i, n_mer in enumerate(it.product(all_res, repeat=4)):
            n_mer_seq = ''.join(list(n_mer)).replace(' ', '') 
            temp_p = 'RH' + n_mer_seq + p
            score, did_find_candidate = spectra_score(temp_p,
                                                      mass_ranges_2,
                                                      target_mass)
            if did_find_candidate:
                found_something = True
            if max_score is None:
                max_score = score
                best_list_2 = [temp_p.sequence[2:]]
            elif abs(score - max_score) < 1:
                print('%10d | %20s | %8d | %10.2f' % (i,temp_p.sequence, score, target_mass - temp_p.mass))
                if did_find_candidate: print_ion_table(temp_p.all_ions(), mass_ranges_2)
                best_list_2.append(temp_p.sequence[2:])
            elif score > max_score:
                best_list_2 = [temp_p.sequence[2:]]
                max_score = score
                print('%10d | %20s | * %6d | %10.2f' % (i,temp_p.sequence, score, target_mass - temp_p.mass))
                if did_find_candidate: print_ion_table(temp_p.all_ions(), mass_ranges_2)
        if max_score > final_max_score:
            best_list_old = list(set([s[1:] for s in best_list_2]))
            candidates = [Protein('RH' + s) for s in best_list_2]
            i_best = 0
            p = Protein(best_list_old[0])
            final_max_score = max_score
        else:
            i_best += 1
            if i_best < len(best_list_old):
                p = Protein(best_list_old[i_best])
            else:
                  break

print('Candidates:')
if candidates:                
    for p in candidates:
        print('  Seq - %s' % p.sequence)
        print('  (M + 2H)++ = %10.2f' % (p.charged_mass(2) / 2))
else:
    print('Failed to find any candidates')

