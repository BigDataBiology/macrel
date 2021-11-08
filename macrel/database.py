GROUPS_SA = ['ALFCGIVW', 'RKQEND', 'MSPTHY'] #solventaccess
GROUPS_HB = ['ILVWAMGT', 'FYSQCN', 'PHKEDR'] # HEIJNE&BLOMBERG1979

CTDD_groups = GROUPS_SA + GROUPS_HB

#' # http://emboss.bioinformatics.nl/cgi-bin/emboss/pepstats
#' # Property      Residues              Number  Mole%
#' # Tiny          (A+C+G+S+T)             4   19.048
#' # Small         (A+B+C+D+G+N+P+S+T+V)   4   19.048
#' # Aliphatic     (A+I+L+V)               5   23.810
#' # Aromatic      (F+H+W+Y)               5   23.810
#' # Non-polar     (A+C+F+G+I+L+M+P+V+W+Y) 11  52.381
#' # Polar         (D+E+H+K+N+Q+R+S+T+Z)   9   42.857
#' # Charged       (B+D+E+H+K+R+Z)         8   38.095
#' # Basic         (H+K+R)                 8   38.095
#' # Acidic        (B+D+E+Z)               0   00.000

_aa_groups = [
    frozenset('ACGST'),          # Tiny
    frozenset('ABCDGNPSTV'),     # Small
    frozenset('AILV'),           # Aliphatic
    frozenset('FHWY'),           # Aromatic
    frozenset('ACFGILMPVWY'),    # Nonpolar
    frozenset('DEHKNQRSTZ'),     # Polar
    frozenset('BDEHKRZ'),        # Charged
    frozenset('HKR'),            # Basic
    frozenset('BDEZ'),           # Acidic
]

pos_pks = {'Nterm': 8.6, 'K': 10.8, 'R': 12.5, 'H': 6.5}

neg_pks = {'Cterm': 3.6, 'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.1}

boman_scale = {'L': -4.92, 'I': -4.92, 'V': -4.04, 'F': -2.98, 'M': -2.35, 'W': -2.33, 'A': -1.81, 'C': -1.28,
               'G': -0.94, 'Y': 0.14, 'T': 2.57, 'S': 3.40, 'H': 4.66, 'Q': 5.54, 'K': 5.55, 'N': 6.64, 'E': 6.81,
               'D': 8.72, 'R': 14.92, 'P': 0., 'X': 0.}

eisenberg = {'I': 1.38, 'F': 1.19, 'V': 1.08, 'L': 1.06, 'W': 0.81, 'M': 0.64, 'A': 0.62,
             'G': 0.48, 'C': 0.29, 'Y': 0.26, 'P': 0.12, 'T': -0.05, 'S': -0.18, 'H': -0.4,
             'E': -0.74, 'N': -0.78, 'Q': -0.85, 'D': -0.9, 'K': -1.5, 'R': -2.53, 'X': 0.}

instability = {
            'A': {'A': 1.0, 'C': 44.94, 'E': 1.0, 'D': -7.49, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': -7.49, 'K': 1.0,
                  'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': 1.0, 'V': 1.0,
                  'Y': 1.0, 'X': 0.},
            'C': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 20.26, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 33.6, 'K': 1.0,
                  'M': 33.6, 'L': 20.26, 'N': 1.0, 'Q': -6.54, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': 33.6, 'W': 24.68,
                  'V': -6.54, 'Y': 1.0, 'X': 0.},
            'E': {'A': 1.0, 'C': 44.94, 'E': 33.6, 'D': 20.26, 'G': 1.0, 'F': 1.0, 'I': 20.26, 'H': -6.54, 'K': 1.0,
                  'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': 1.0, 'T': 1.0, 'W': -14.03,
                  'V': 1.0, 'Y': 1.0, 'X': 0.},
            'D': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0, 'K': -7.49,
                  'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 1.0, 'S': 20.26, 'R': -6.54, 'T': -14.03, 'W': 1.0,
                  'V': 1.0, 'Y': 1.0, 'X': 0.},
            'G': {'A': -7.49, 'C': 1.0, 'E': -6.54, 'D': 1.0, 'G': 13.34, 'F': 1.0, 'I': -7.49, 'H': 1.0, 'K': -7.49,
                  'M': 1.0, 'L': 1.0, 'N': -7.49, 'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0, 'T': -7.49, 'W': 13.34,
                  'V': 1.0, 'Y': -7.49, 'X': 0.},
            'F': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 13.34, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': -14.03,
                  'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': 1.0, 'V': 1.0,
                  'Y': 33.601, 'X': 0.},
            'I': {'A': 1.0, 'C': 1.0, 'E': 44.94, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 13.34, 'K': -7.49,
                  'M': 1.0, 'L': 20.26, 'N': 1.0, 'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': 1.0,
                  'V': -7.49, 'Y': 1.0, 'X': 0.},
            'H': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -9.37, 'F': -9.37, 'I': 44.94, 'H': 1.0, 'K': 24.68,
                  'M': 1.0, 'L': 1.0, 'N': 24.68, 'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0, 'T': -6.54, 'W': -1.88,
                  'V': 1.0, 'Y': 44.94, 'X': 0.},
            'K': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -7.49, 'F': 1.0, 'I': -7.49, 'H': 1.0, 'K': 1.0,
                  'M': 33.6, 'L': -7.49, 'N': 1.0, 'Q': 24.64, 'P': -6.54, 'S': 1.0, 'R': 33.6, 'T': 1.0, 'W': 1.0,
                  'V': -7.49, 'Y': 1.0, 'X': 0.},
            'M': {'A': 13.34, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 58.28, 'K': 1.0,
                  'M': -1.88, 'L': 1.0, 'N': 1.0, 'Q': -6.54, 'P': 44.94, 'S': 44.94, 'R': -6.54, 'T': -1.88, 'W': 1.0,
                  'V': 1.0, 'Y': 24.68, 'X': 0.},
            'L': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': -7.49, 'M': 1.0,
                  'L': 1.0, 'N': 1.0, 'Q': 33.6, 'P': 20.26, 'S': 1.0, 'R': 20.26, 'T': 1.0, 'W': 24.68, 'V': 1.0,
                  'Y': 1.0, 'X': 0.},
            'N': {'A': 1.0, 'C': -1.88, 'E': 1.0, 'D': 1.0, 'G': -14.03, 'F': -14.03, 'I': 44.94, 'H': 1.0, 'K': 24.68,
                  'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': -6.54, 'P': -1.88, 'S': 1.0, 'R': 1.0, 'T': -7.49, 'W': -9.37,
                  'V': 1.0, 'Y': 1.0, 'X': 0.},
            'Q': {'A': 1.0, 'C': -6.54, 'E': 20.26, 'D': 20.26, 'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0, 'K': 1.0,
                  'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 1.0, 'T': 1.0, 'W': 1.0,
                  'V': -6.54, 'Y': -6.54, 'X': 0.},
            'P': {'A': 20.26, 'C': -6.54, 'E': 18.38, 'D': -6.54, 'G': 1.0, 'F': 20.26, 'I': 1.0, 'H': 1.0, 'K': 1.0,
                  'M': -6.54, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': -6.54, 'T': 1.0, 'W': -1.88,
                  'V': 20.26, 'Y': 1.0, 'X': 0.},
            'S': {'A': 1.0, 'C': 33.6, 'E': 20.26, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': 1.0, 'M': 1.0,
                  'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 44.94, 'S': 20.26, 'R': 20.26, 'T': 1.0, 'W': 1.0, 'V': 1.0,
                  'Y': 1.0, 'X': 0.},
            'R': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 20.26, 'K': 1.0,
                  'M': 1.0, 'L': 1.0, 'N': 13.34, 'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 58.28, 'T': 1.0, 'W': 58.28,
                  'V': 1.0, 'Y': -6.54, 'X': 0.},
            'T': {'A': 1.0, 'C': 1.0, 'E': 20.26, 'D': 1.0, 'G': -7.49, 'F': 13.34, 'I': 1.0, 'H': 1.0, 'K': 1.0,
                  'M': 1.0, 'L': 1.0, 'N': -14.03, 'Q': -6.54, 'P': 1.0, 'S': 1.0, 'R': 1.0, 'T': 1.0, 'W': -14.03,
                  'V': 1.0, 'Y': 1.0, 'X': 0.},
            'W': {'A': -14.03, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -9.37, 'F': 1.0, 'I': 1.0, 'H': 24.68, 'K': 1.0,
                  'M': 24.68, 'L': 13.34, 'N': 13.34, 'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0, 'T': -14.03, 'W': 1.0,
                  'V': -7.49, 'Y': 1.0, 'X': 0.},
            'V': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': -14.03, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 1.0, 'K': -1.88,
                  'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0, 'T': -7.49, 'W': 1.0,
                  'V': 1.0, 'Y': -6.54, 'X': 0.},
            'Y': {'A': 24.68, 'C': 1.0, 'E': -6.54, 'D': 24.68, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 13.34, 'K': 1.0,
                  'M': 44.94, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 13.34, 'S': 1.0, 'R': -15.91, 'T': -7.49, 'W': -9.37,
                  'V': 1.0, 'Y': 13.34, 'X': 0.},
            'X': {'A': 0., 'C': 0., 'E': 0., 'D': 0., 'G': 0., 'F': 0., 'I': 0., 'H': 0., 'K': 0.,
                  'M': 0., 'L': 0., 'N': 0., 'Q': 0., 'P': 0., 'S': 0., 'R': 0., 'T': 0., 'W': 0.,
                  'V': 0., 'Y': 0., 'X': 0.}}

# Transform the above into a simpler-to-use form, a single dictionary that maps dimers
instability2 = {}
for k,vs in instability.items():
    for k2,v in vs.items():
        if v:
            instability2[k+k2] = v
