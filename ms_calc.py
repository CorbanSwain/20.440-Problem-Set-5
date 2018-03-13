#!\bin\python3

import numpy as np

## Problem 1, Part 1 ##

# protein sequence
seq = 'ELFDDPSYVNVQNLDK'
n_aa = len(seq)
n_phospho = 1
# monoisotopic amino acid masses
aa_mass = {'A' : 71.03711,
           'R' : 156.10111,
           'N' : 114.04293,
           'D' : 115.02694,
           'C' : 103.00919,
           'E' : 129.04259,
           'Q' : 128.05858,
           'G' : 57.02146,
           'H' : 137.05891,
           'I' : 113.08406,
           'L' : 113.08406,
           'K' : 128.09496,
           'M' : 131.04049,
           'F' : 147.06841,
           'P' : 97.05276,
           'S' : 87.03203,
           'T' : 101.04768,
           'W' : 186.07931,
           'Y' : 163.06333,
           'V' : 99.06841}

water_mass = 18.010565
h_mass = 1.007276
phospho_mass = 79.97

seq_mass = np.array([aa_mass[aa] for aa in seq])
total_mass = np.sum(seq_mass) + h_mass + water_mass + n_phospho * phospho_mass

print(total_mass)

p_site = -1
y_mass = total_mass
b_mass = 0.0
for i, aa in enumerate(seq):
    delta = aa_mass[aa]
    if b_mass - 0.0 < 0.0001:
        b_mass = h_mass
    y_mass = y_mass - delta
    b_mass = b_mass + delta
    infostr = '%2d | %2s | b_%2d: %8.2f | y_%2d: %8.2f' % \
              (i + +1, aa, i + 1, b_mass, 16 - i, y_mass)
    print(infostr)

    
     
